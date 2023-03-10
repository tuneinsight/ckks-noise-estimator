package main

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func main() {

	// Simple small sets of parameters, no need to modify it
	var err error
	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     16,
		LogQ:     []int{60, 45},
		LogP:     []int{61},
		LogScale: 45,
	}); err != nil {
		panic(err)
	}

	// This you can modify, it's the log2(slots) parameter
	LogSlots := 13
	Slots := 1 << LogSlots

	// Computes the full factorization of the DFT layers in plaintext
	DFTLayers := GenMatrices(params.LogN(), LogSlots, advanced.Decode)

	// Kgen (but we generate a zero key)
	kgen := ckks.NewKeyGenerator(params)
	sk := rlwe.NewSecretKey(params.Parameters) // sk = 0
	pk := kgen.GenPublicKeyNew(sk)

	// The interface storing the Galois keys
	evk := rlwe.NewEvaluationKeySet()
	for i := range DFTLayers[LogSlots-1] {
		if i != 0 {
			if err = evk.Add(kgen.GenGaloisKeyNew(params.GaloisElementForColumnRotationBy(i), sk)); err != nil {
				panic(err)
			}
		}
	}

	// Our evaluator
	eval := ckks.NewEvaluator(params, evk)

	// Our encoder
	ecd := ckks.NewEncoder(params)

	// Our encryptor
	enc := ckks.NewEncryptor(params, pk)

	// Our decryptor
	dec := ckks.NewDecryptor(params, sk)

	// Generates a random plaintext Y = X^{N/(2*slots) = gap} IN THE RING with coefficient uniformely distributed [-1, 1]
	r := rand.New(rand.NewSource(time.Now().Unix()))
	gap := params.N() / (2 << LogSlots)
	values := make([]float64, params.N())
	for i := 0; i < params.N(); i += gap {
		values[i] = 2 * (r.Float64() - 0.5)
	}

	// Allocates our plaintext
	pt := ckks.NewPlaintext(params, params.MaxLevel())
	pt.LogSlots = LogSlots                      // Sets the number of slots since it's not the default LogSlots = LogN-1
	pt.EncodingDomain = rlwe.CoefficientsDomain // Sets the encoding domain (if you have a better name, like `Ring`, feel free to give feedback)

	// Encodes on the plaintext
	if err = ecd.Encode(values, pt); err != nil {
		panic(err)
	}

	// Computes the standard deviation of the message and its error
	vErr := make([]float64, params.N())
	ecd.Decode(pt, vErr)
	for i := range values {
		vErr[i] -= values[i]
	}

	stdMsg := estimator.STD(values) // Plaintext standard deviation IN THE RING
	stdErr := estimator.STD(vErr)   // Plaintext error standard deviation IN THE RING

	// Encrypts
	ct := enc.EncryptNew(pt)
	ct.EncodingDomain = rlwe.SlotsDomain // We have to do that, else the evaluator will return an error, this will be fixed

	// Allocates the result ciphertext
	res := ckks.NewCiphertext(params, ct.Degree(), ct.Level())

	// To be able to evaluate the plaintext circuit, we need to decode the values that were sampled IN THE RING
	vCmplx := make([]complex128, Slots)

	// R^{N} -> C^{Slots}
	for i, idx, jdx := 0, 0, params.N()>>1; i < Slots; i, idx, jdx = i+1, idx+gap, jdx+gap {
		vCmplx[i] = complex(values[idx], values[jdx])
	}
	// Decoding C^{Slots} -> C^{Slots}
	ecd.FFT(vCmplx, LogSlots)

	// Instantiates our estimator
	est := estimator.NewEstimator(params.N(), 0, params.Q(), params.P())

	// Sets our base ciphertext
	ctEstimator := estimator.NewCiphertextPK(estimator.NewPlaintext(stdMsg*params.DefaultScale().Float64(), stdErr*params.DefaultScale().Float64(), params.MaxLevel()))

	// Sets our receiver ciphertext, with zero error
	resEstimator := estimator.Element{
		Level:   res.Level(),
		Message: estimator.NewFloat(0),
		Noise:   []*big.Float{estimator.NewFloat(0), estimator.NewFloat(0)}, // ModDown by P kills all error
	}

	// Sets our receiver plaintext
	resCmplex := make([]complex128, Slots)

	// Iterates only on the last matrix
	// i: amout to rotate
	// vec: []complex128
	for i, vec := range DFTLayers[LogSlots-1] {

		// Standard deviation and error IN THE RING of diagonal of the matrix

		// >>>>>>>>>>>>>>>>>THIS IS PROBABLY WHERE THE ISSUE IS <<<<<<<<<<<<<<<<<<
		// >>>>>>>>>>>>>>>>>THIS IS PROBABLY WHERE THE ISSUE IS <<<<<<<<<<<<<<<<<<
		// >>>>>>>>>>>>>>>>>THIS IS PROBABLY WHERE THE ISSUE IS <<<<<<<<<<<<<<<<<<
		// >>>>>>>>>>>>>>>>>THIS IS PROBABLY WHERE THE ISSUE IS <<<<<<<<<<<<<<<<<<
		std := estimator.GetSTDEncodedVector(ecd, params.N(), LogSlots, vec)
		// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		// estimator plaintext vector
		ptTmp := estimator.NewPlaintext(std[0]*params.DefaultScale().Float64(), std[1]*params.DefaultScale().Float64(), params.MaxLevel())

		if i == 0 {

			// Encrypted
			tmp := eval.MulNew(ct, vec)
			eval.Add(res, tmp, res)

			// Plaintext
			resCmplex = add(resCmplex, mul(vCmplx, vec))

			// Estimator
			resEstimator = est.Add(resEstimator, est.Mul(ctEstimator, ptTmp))

		} else {

			// Encrypted
			tmp := eval.MulNew(ct, vec)
			eval.Rotate(tmp, i, tmp)
			eval.Add(res, tmp, res)

			// Plaintext
			resCmplex = add(resCmplex, rotate(mul(vCmplx, vec), i))

			// Estimator
			resEstimator = est.Add(resEstimator, est.KeySwitch(est.Mul(ctEstimator, ptTmp)))
		}
	}

	res.LogSlots = LogSlots                      // By default we have Y=X^{N/2} contained in Y=X^{N/4}, ..., contained into Y = X, so we need to overwrite of LogSlots of the reveiver
	res.EncodingDomain = rlwe.CoefficientsDomain // Sets the encoding domain back to what we want
	if err = eval.Rescale(res, params.DefaultScale(), res); err != nil {
		panic(err)
	}

	// We also need to call rescal on the estimated ciphertext
	resEstimator = est.Rescale(resEstimator)

	// Nothing to see here
	vCmplx = resCmplex

	// Re-encode the values to have the result plaintext IN THE RING
	// C^{Slots} -> R^{N}
	ecd.IFFT(vCmplx, LogSlots)
	for i, idx, jdx := 0, 0, params.N()>>1; i < Slots; i, idx, jdx = i+1, idx+gap, jdx+gap {
		values[idx], values[jdx] = real(vCmplx[i]), imag(vCmplx[i])
	}

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

	// Decrypt IN THE RING
	valuesF64 := make([]float64, params.N())
	if err = ecd.Decode(pt, valuesF64); err != nil {
		panic(err)
	}

	// Logs the estimator (have) and the actual result (want)
	fmt.Println("Have:", est.Std(resEstimator))
	fmt.Println("Want:", math.Log2(estimator.STD(valuesF64)))
}

// GenMatrices returns the full factorization IDFT (encoding) or DFT (decoding) matrix.
func GenMatrices(LogN, LogSlots int, Type advanced.DFTType) (DFTLayers []map[int][]complex128) {

	slots := 1 << LogSlots

	roots := ckks.GetRootsFloat64(slots << 2)
	pow5 := make([]int, (slots<<1)+1)
	pow5[0] = 1
	mask := (slots << 2) - 1
	for i := 1; i < (slots<<1)+1; i++ {
		pow5[i] = pow5[i-1] * 5
		pow5[i] &= mask
	}

	var a, b, c [][]complex128
	if Type == advanced.Encode {
		a, b, c = ifftPlainVec(LogSlots, 1<<LogSlots, roots, pow5)
	} else {
		a, b, c = fftPlainVec(LogSlots, 1<<LogSlots, roots, pow5)
	}

	DFTLayers = make([]map[int][]complex128, LogSlots)
	for i := 0; i < LogSlots; i++ {
		DFTLayers[i] = genFFTDiagMatrix(LogSlots, LogSlots-i, a[i], b[i], c[i], Type, false)
	}

	if Type == advanced.Encode {
		for j := range DFTLayers {
			for x := range DFTLayers[j] {
				v := DFTLayers[j][x]
				for i := range v {
					v[i] /= 2.0
				}
			}
		}
	}

	return
}

func genFFTDiagMatrix(logL, fftLevel int, a, b, c []complex128, ltType advanced.DFTType, bitreversed bool) (vectors map[int][]complex128) {

	var rot int

	if ltType == advanced.Encode && !bitreversed || ltType == advanced.Decode && bitreversed {
		rot = 1 << (fftLevel - 1)
	} else {
		rot = 1 << (logL - fftLevel)
	}

	vectors = make(map[int][]complex128)

	if bitreversed {
		ckks.SliceBitReverseInPlaceComplex128(a, 1<<logL)
		ckks.SliceBitReverseInPlaceComplex128(b, 1<<logL)
		ckks.SliceBitReverseInPlaceComplex128(c, 1<<logL)

		if len(a) > 1<<logL {
			ckks.SliceBitReverseInPlaceComplex128(a[1<<logL:], 1<<logL)
			ckks.SliceBitReverseInPlaceComplex128(b[1<<logL:], 1<<logL)
			ckks.SliceBitReverseInPlaceComplex128(c[1<<logL:], 1<<logL)
		}
	}

	addToDiagMatrix(vectors, 0, a)
	addToDiagMatrix(vectors, rot, b)
	addToDiagMatrix(vectors, (1<<logL)-rot, c)

	return
}

func addToDiagMatrix(diagMat map[int][]complex128, index int, vec []complex128) {
	if diagMat[index] == nil {
		diagMat[index] = vec
	} else {
		diagMat[index] = add(diagMat[index], vec)
	}
}

func rotate(x []complex128, n int) (y []complex128) {

	y = make([]complex128, len(x))

	mask := int(len(x) - 1)

	// Rotates to the left
	for i := 0; i < len(x); i++ {
		y[i] = x[(i+n)&mask]
	}

	return
}

func mul(a, b []complex128) (res []complex128) {

	res = make([]complex128, len(a))

	for i := 0; i < len(a); i++ {
		res[i] = a[i] * b[i]
	}

	return
}

func add(a, b []complex128) (res []complex128) {

	res = make([]complex128, len(a))

	for i := 0; i < len(a); i++ {
		res[i] = a[i] + b[i]
	}

	return
}

func fftPlainVec(logN, dslots int, roots []complex128, pow5 []int) (a, b, c [][]complex128) {

	var N, m, index, tt, gap, k, mask, idx1, idx2 int

	N = 1 << logN

	a = make([][]complex128, logN)
	b = make([][]complex128, logN)
	c = make([][]complex128, logN)

	var size int
	if 2*N == dslots {
		size = 2
	} else {
		size = 1
	}

	index = 0
	for m = 2; m <= N; m <<= 1 {

		a[index] = make([]complex128, dslots)
		b[index] = make([]complex128, dslots)
		c[index] = make([]complex128, dslots)

		tt = m >> 1

		for i := 0; i < N; i += m {

			gap = N / m
			mask = (m << 2) - 1

			for j := 0; j < m>>1; j++ {

				k = (pow5[j] & mask) * gap

				idx1 = i + j
				idx2 = i + j + tt

				for u := 0; u < size; u++ {
					a[index][idx1+u*N] = 1
					a[index][idx2+u*N] = -roots[k]
					b[index][idx1+u*N] = roots[k]
					c[index][idx2+u*N] = 1
				}
			}
		}

		index++
	}

	return
}

func ifftPlainVec(logN, dslots int, roots []complex128, pow5 []int) (a, b, c [][]complex128) {

	var N, m, index, tt, gap, k, mask, idx1, idx2 int

	N = 1 << logN

	a = make([][]complex128, logN)
	b = make([][]complex128, logN)
	c = make([][]complex128, logN)

	var size int
	if 2*N == dslots {
		size = 2
	} else {
		size = 1
	}

	index = 0
	for m = N; m >= 2; m >>= 1 {

		a[index] = make([]complex128, dslots)
		b[index] = make([]complex128, dslots)
		c[index] = make([]complex128, dslots)

		tt = m >> 1

		for i := 0; i < N; i += m {

			gap = N / m
			mask = (m << 2) - 1

			for j := 0; j < m>>1; j++ {

				k = ((m << 2) - (pow5[j] & mask)) * gap

				idx1 = i + j
				idx2 = i + j + tt

				for u := 0; u < size; u++ {

					a[index][idx1+u*N] = 1
					a[index][idx2+u*N] = -roots[k]
					b[index][idx1+u*N] = 1
					c[index][idx2+u*N] = roots[k]
				}
			}
		}

		index++
	}

	return
}
