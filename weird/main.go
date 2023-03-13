package main

import (
	"fmt"
	"math"
	"math/rand"

	//"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func main() {

	H := 32768
	LogSlots := 14
	r := rand.New(rand.NewSource(0))

	var Xs distribution.Distribution
	if H != 0 {
		Xs = &distribution.Ternary{H: H}
	} else {
		Xs = nil
	}

	// Simple small sets of parameters, no need to modify it
	var err error
	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     16,
		LogQ:     []int{60, 45},
		LogP:     []int{61},
		LogScale: 45,
		Xs:       Xs,
	}); err != nil {
		panic(err)
	}

	Slots := 1 << LogSlots

	// Kgen (but we generate a zero key)
	kgen := ckks.NewKeyGenerator(params)
	var sk *rlwe.SecretKey
	if H == 0 {
		sk = rlwe.NewSecretKey(params.Parameters)
	} else {
		sk = kgen.GenSecretKeyNew()
	}

	pk := kgen.GenPublicKeyNew(sk)

	// Our evaluator
	eval := ckks.NewEvaluator(params, nil)

	// Our encoder
	ecd := ckks.NewEncoder(params)

	// Our encryptor
	enc := ckks.NewEncryptor(params, pk)

	// Our decryptor
	dec := ckks.NewDecryptor(params, sk)

	// Generates a random plaintext Y = X^{N/(2*slots) = gap} IN THE RING with coefficient uniformely distributed [-1, 1]
	values := make([]complex128, Slots)
	for j := range values {
		values[j] = complex(2*r.Float64()-1, 2*r.Float64()-1)
	}

	// Allocates our plaintext
	pt := ckks.NewPlaintext(params, params.MaxLevel())
	pt.LogSlots = LogSlots // Sets the number of slots since it's not the default LogSlots = LogN-1

	// Encodes on the plaintext
	if err = ecd.Encode(values, pt); err != nil {
		panic(err)
	}

	// Encrypts
	ct := enc.EncryptNew(pt)

	// Instantiates our estimator
	est := estimator.NewEstimator(params.N(), H, params.Q(), params.P())

	stdValues := estimator.GetSTDEncodedVector(ecd, params.N(), LogSlots, values)

	// Sets our base ciphertext
	ctEstimator := estimator.NewCiphertextPK(estimator.NewPlaintext(stdValues[0]*params.DefaultScale().Float64(), stdValues[1]*params.DefaultScale().Float64(), params.MaxLevel()))

	vec := make([]complex128, Slots)
	for j := range vec {
		vec[j] = complex(2*r.Float64()-1, 2*r.Float64()-1)
	}

	// Standard deviation and error IN THE RING of diagonal of the matrix
	stdVec := estimator.GetSTDEncodedVector(ecd, params.N(), LogSlots, vec)

	// estimator plaintext vector
	ptTmp := estimator.NewPlaintext(stdVec[0]*params.DefaultScale().Float64(), stdVec[1]*params.DefaultScale().Float64(), params.MaxLevel())

	// Encrypted
	res := eval.MulNew(ct, vec)

	// Plaintext
	values = mul(values, vec)

	// Estimator
	resEstimator := est.Mul(ctEstimator, ptTmp)

	// Encodes reference values (here we could also skip the previous step and do it by setting pt.EncodingDomain = rlwe.SlotsDomain)
	if err = ecd.Encode(values, pt); err != nil {
		panic(err)
	}

	// Subtracts to result
	eval.Sub(res, pt, res)

	// Decrypt result
	dec.Decrypt(res, pt)

	/// Sets the plaintext scale to 1 (to get the result non-scale)
	pt.Scale = rlwe.NewScale(1)
	pt.EncodingDomain = rlwe.CoefficientsDomain

	// Decrypt IN THE RING
	valuesF64 := make([]float64, params.N())
	if err = ecd.Decode(pt, valuesF64); err != nil {
		panic(err)
	}

	// Logs the estimator (have) and the actual result (want)
	fmt.Println("Have:", est.Std(resEstimator))
	fmt.Println("Want:", math.Log2(estimator.STD(valuesF64)))
}

func mul(a, b []complex128) (res []complex128) {

	res = make([]complex128, len(a))

	for i := 0; i < len(a); i++ {
		res[i] = a[i] * b[i]
	}

	return
}
