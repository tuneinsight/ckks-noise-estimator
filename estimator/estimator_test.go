package estimator_test

import (
	"fmt"
	"math/big"
	"os"
	"runtime/pprof"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/schemes/ckks"
	"github.com/tuneinsight/lattigo/v4/core/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

type testContext struct {
	params    hefloat.Parameters
	estimator estimator.Parameters
	kgen      *rlwe.KeyGenerator
	sk        *rlwe.SecretKey
	pk        *rlwe.PublicKey
	decryptor *rlwe.Decryptor
	encoder   *hefloat.Encoder
	evaluator *hefloat.Evaluator
}

func newTestContext(params hefloat.Parameters) testContext {
	kgen := hefloat.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk))

	estimator := estimator.NewParameters(params)
	estimator.Heuristic = true

	return testContext{
		params:    params,
		estimator: estimator,
		kgen:      kgen,
		sk:        sk,
		pk:        pk,
		decryptor: hefloat.NewDecryptor(params, sk),
		encoder:   hefloat.NewEncoder(params, 256),
		evaluator: hefloat.NewEvaluator(params, evk),
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

	pt = hefloat.NewPlaintext(tc.params, tc.params.MaxLevel())
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

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            12,
		LogQ:            []int{60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60},
		LogP:            []int{61, 61},
		LogDefaultScale: 60,
	})

	if err != nil {
		panic(err)
	}

	tc := newTestContext(params)

	testChebyshevPowers(tc, t)
	testKeySwitching(tc, t)
}

func newPrecisionStats() ckks.PrecisionStats {
	return ckks.PrecisionStats{
		MaxDelta:        newStats(),
		MinDelta:        newStats(),
		MaxPrecision:    newStats(),
		MinPrecision:    newStats(),
		MeanDelta:       newStats(),
		MeanPrecision:   newStats(),
		MedianDelta:     newStats(),
		MedianPrecision: newStats(),
	}
}

func newStats() ckks.Stats {
	return ckks.Stats{new(big.Float), new(big.Float), new(big.Float)}
}

func testKeySwitching(tc testContext, t *testing.T) {

	t.Run("Rotate/sk", func(t *testing.T) {

		statsHave := newPrecisionStats()
		statsWant := newPrecisionStats()

		d := 8

		gk := tc.kgen.GenGaloisKeyNew(tc.params.GaloisElement(1), tc.sk)

		evk := rlwe.NewMemEvaluationKeySet(nil, gk)

		eval := tc.evaluator.WithKey(evk)

		for w := 0; w < d; w++ {

			values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

			RunTimed("estimator", func() {
				el.Rotate(1)
				el.Decrypt()
				el.Normalize()
			})

			RunTimed("encrypted", func() {
				require.NoError(t, eval.Rotate(ct, 1, ct))
			})

			utils.RotateSliceInPlace(values, 1)

			pWant := hefloat.GetPrecisionStats(tc.params, tc.encoder, nil, values, el.Value[0], 0, false)
			pHave := hefloat.GetPrecisionStats(tc.params, tc.encoder, tc.decryptor, values, ct, 0, false)

			addStats(&statsWant, &pWant, &statsWant)
			addStats(&statsHave, &pHave, &statsHave)
		}

		statsDiv(&statsWant, new(big.Float).SetInt64(int64(d)))
		statsDiv(&statsHave, new(big.Float).SetInt64(int64(d)))

		fmt.Println(statsWant.String())
		fmt.Println(statsHave.String())
	})
}

func testChebyshevPowers(tc testContext, t *testing.T) {

	t.Run("Chebyshev512Power/sk", func(t *testing.T) {

		statsHave := newPrecisionStats()
		statsWant := newPrecisionStats()

		d := 16

		for w := 0; w < d; w++ {

			values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

			n := 11

			RunProfiled(func() {
				RunTimed("estimator", func() {
					pb := estimator.NewPowerBasis(el, bignum.Chebyshev)
					pb.GenPower(1<<n, false)
					el = *pb.Value[1<<n]
					el.Decrypt()
					el.Normalize()
				})
			})

			RunTimed("encrypted", func() {
				eval := tc.evaluator
				pb := hefloat.NewPowerBasis(ct, bignum.Chebyshev)
				if err := pb.GenPower(1<<n, false, eval); err != nil{
					panic(err)
				}
				ct = pb.Value[1<<n]
			})

			mul := bignum.NewComplexMultiplier().Mul

			for k := 0; k < n; k++ {
				for i := range values {
					mul(values[i], values[i], values[i])
					values[i].Add(values[i], values[i])
					values[i].Add(values[i], &bignum.Complex{estimator.NewFloat(-1), estimator.NewFloat(0)})
				}
			}

			pWant := hefloat.GetPrecisionStats(tc.params, tc.encoder, nil, values, el.Value[0], 0, false)
			pHave := hefloat.GetPrecisionStats(tc.params, tc.encoder, tc.decryptor, values, ct, 0, false)

			addStats(&statsWant, &pWant, &statsWant)
			addStats(&statsHave, &pHave, &statsHave)
		}

		statsDiv(&statsWant, new(big.Float).SetInt64(int64(d)))
		statsDiv(&statsHave, new(big.Float).SetInt64(int64(d)))

		fmt.Println(statsWant.String())
		fmt.Println(statsHave.String())
	})
}

func statsDiv(a *ckks.PrecisionStats, b *big.Float) {
	a.MaxDelta.Real.Quo(a.MaxDelta.Real, b)
	a.MaxDelta.Imag.Quo(a.MaxDelta.Imag, b)
	a.MaxDelta.L2.Quo(a.MaxDelta.L2, b)

	a.MinDelta.Real.Quo(a.MinDelta.Real, b)
	a.MinDelta.Imag.Quo(a.MinDelta.Imag, b)
	a.MinDelta.L2.Quo(a.MinDelta.L2, b)

	a.MaxPrecision.Real.Quo(a.MaxPrecision.Real, b)
	a.MaxPrecision.Imag.Quo(a.MaxPrecision.Imag, b)
	a.MaxPrecision.L2.Quo(a.MaxPrecision.L2, b)

	a.MinPrecision.Real.Quo(a.MinPrecision.Real, b)
	a.MinPrecision.Imag.Quo(a.MinPrecision.Imag, b)
	a.MinPrecision.L2.Quo(a.MinPrecision.L2, b)

	a.MeanDelta.Real.Quo(a.MeanDelta.Real, b)
	a.MeanDelta.Imag.Quo(a.MeanDelta.Imag, b)
	a.MeanDelta.L2.Quo(a.MeanDelta.L2, b)

	a.MeanPrecision.Real.Quo(a.MeanPrecision.Real, b)
	a.MeanPrecision.Imag.Quo(a.MeanPrecision.Imag, b)
	a.MeanPrecision.L2.Quo(a.MeanPrecision.L2, b)

	a.MedianDelta.Real.Quo(a.MedianDelta.Real, b)
	a.MedianDelta.Imag.Quo(a.MedianDelta.Imag, b)
	a.MedianDelta.L2.Quo(a.MedianDelta.L2, b)

	a.MedianPrecision.Real.Quo(a.MedianPrecision.Real, b)
	a.MedianPrecision.Imag.Quo(a.MedianPrecision.Imag, b)
	a.MedianPrecision.L2.Quo(a.MedianPrecision.L2, b)
}

func addStats(a, b, c *ckks.PrecisionStats) {
	c.MaxDelta.Real.Add(a.MaxDelta.Real, b.MaxDelta.Real)
	c.MaxDelta.Imag.Add(a.MaxDelta.Imag, b.MaxDelta.Imag)
	c.MaxDelta.L2.Add(a.MaxDelta.L2, b.MaxDelta.L2)

	c.MinDelta.Real.Add(a.MinDelta.Real, b.MinDelta.Real)
	c.MinDelta.Imag.Add(a.MinDelta.Imag, b.MinDelta.Imag)
	c.MinDelta.L2.Add(a.MinDelta.L2, b.MinDelta.L2)

	c.MaxPrecision.Real.Add(a.MaxPrecision.Real, b.MaxPrecision.Real)
	c.MaxPrecision.Imag.Add(a.MaxPrecision.Imag, b.MaxPrecision.Imag)
	c.MaxPrecision.L2.Add(a.MaxPrecision.L2, b.MaxPrecision.L2)

	c.MinPrecision.Real.Add(a.MinPrecision.Real, b.MinPrecision.Real)
	c.MinPrecision.Imag.Add(a.MinPrecision.Imag, b.MinPrecision.Imag)
	c.MinPrecision.L2.Add(a.MinPrecision.L2, b.MinPrecision.L2)

	c.MeanDelta.Real.Add(a.MeanDelta.Real, b.MeanDelta.Real)
	c.MeanDelta.Imag.Add(a.MeanDelta.Imag, b.MeanDelta.Imag)
	c.MeanDelta.L2.Add(a.MeanDelta.L2, b.MeanDelta.L2)

	c.MeanPrecision.Real.Add(a.MeanPrecision.Real, b.MeanPrecision.Real)
	c.MeanPrecision.Imag.Add(a.MeanPrecision.Imag, b.MeanPrecision.Imag)
	c.MeanPrecision.L2.Add(a.MeanPrecision.L2, b.MeanPrecision.L2)

	c.MedianDelta.Real.Add(a.MedianDelta.Real, b.MedianDelta.Real)
	c.MedianDelta.Imag.Add(a.MedianDelta.Imag, b.MedianDelta.Imag)
	c.MedianDelta.L2.Add(a.MedianDelta.L2, b.MedianDelta.L2)

	c.MedianPrecision.Real.Add(a.MedianPrecision.Real, b.MedianPrecision.Real)
	c.MedianPrecision.Imag.Add(a.MedianPrecision.Imag, b.MedianPrecision.Imag)
	c.MedianPrecision.L2.Add(a.MedianPrecision.L2, b.MedianPrecision.L2)
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
