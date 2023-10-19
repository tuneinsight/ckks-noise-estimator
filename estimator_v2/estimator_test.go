package estimator_test

import (
	"fmt"
	"os"
	"runtime/pprof"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	estimator "github.com/tuneinsight/ckks-bootstrapping-precision/estimator_v2"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

type testContext struct {
	params    ckks.Parameters
	estimator estimator.Parameters
	kgen      *rlwe.KeyGenerator
	sk        *rlwe.SecretKey
	pk        *rlwe.PublicKey
	decryptor *rlwe.Decryptor
	encoder   *ckks.Encoder
	evaluator *ckks.Evaluator
}

func newTestContext(params ckks.Parameters) testContext {
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk))

	return testContext{
		params:    params,
		estimator: estimator.NewParameters(params),
		kgen:      kgen,
		sk:        sk,
		pk:        pk,
		decryptor: rlwe.NewDecryptor(params, sk),
		encoder:   ckks.NewEncoder(params, 256),
		evaluator: ckks.NewEvaluator(params, evk),
	}
}

func newTestVector(tc testContext, key rlwe.EncryptionKey, a, b complex128) (values []*bignum.Complex, el estimator.Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prec := tc.encoder.Prec()

	values = make([]*bignum.Complex, tc.params.MaxSlots())

	for i := range values {
		values[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
		}
	}

	el = estimator.NewElement(tc.estimator, values, 1, tc.estimator.DefaultScale())
	el.AddEncodingNoise()

	pt = ckks.NewPlaintext(tc.params, tc.params.MaxLevel())
	tc.encoder.Encode(values, pt)

	switch key := key.(type) {
	case *rlwe.SecretKey:
		ct, _ = rlwe.NewEncryptor(tc.params, key).EncryptNew(pt)
		el.AddEncryptionNoiseSk()
	case *rlwe.PublicKey:
		ct, _ = rlwe.NewEncryptor(tc.params, key).EncryptNew(pt)
		el.AddEncryptionNoisePk()
	}

	return
}

func TestEstimator(t *testing.T) {

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{61, 61, 61},
		LogDefaultScale: 45,
	})

	if err != nil {
		panic(err)
	}

	tc := newTestContext(params)

	testChebyshevPowers(tc, t)

}

func testChebyshevPowers(tc testContext, t *testing.T) {

	t.Run("Chebyshev512Power/sk", func(t *testing.T) {

		values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

		n := 9

		RunProfiled(func() {
			RunTimed("estimator", func() {

				for k := 0; k < n; k++ {
					el.Mul(&el, &el)
					el.Relinearize()
					el.Add(&el, &el)
					el.Add(&el, -1)
					require.NoError(t, el.Rescale())
				}

				el.Normalize()
			})
		})

		RunTimed("encrypted", func() {
			eval := tc.evaluator
			for k := 0; k < n; k++ {
				require.NoError(t, eval.MulRelin(ct, ct, ct))
				require.NoError(t, eval.Mul(ct, 2, ct))
				require.NoError(t, eval.Add(ct, -1, ct))
				require.NoError(t, eval.Rescale(ct, ct))
			}
		})

		mul := bignum.NewComplexMultiplier().Mul

		for k := 0; k < n; k++ {
			for i := range values {

				mul(values[i], values[i], values[i])
				values[i].Add(values[i], values[i])
				values[i].Add(values[i], &bignum.Complex{estimator.NewFloat(-1), estimator.NewFloat(0)})
			}
		}

		ckks.VerifyTestVectors(tc.params, tc.encoder, nil, values, el.Value, 0, nil, true, t)
		ckks.VerifyTestVectors(tc.params, tc.encoder, tc.decryptor, values, ct, 0, nil, true, t)
	})
}

func TestEncryptPK(t *testing.T) {

}

func RunTimed(msg string, f func()) {
	now := time.Now()
	f()
	fmt.Printf("%s: %s\n", msg, time.Since(now))
}

func RunProfiled(f func()) {
	file, err := os.Create("cpu.prof")

	if err != nil {
		panic(err)
	}

	if err := pprof.StartCPUProfile(file); err != nil {
		panic(err)
	}

	f()

	pprof.StopCPUProfile()
	file.Close()
}
