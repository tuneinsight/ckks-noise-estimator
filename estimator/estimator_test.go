package estimator_test

import (
	"fmt"
	"math"
	"math/big"
	"os"
	"runtime/pprof"
	"testing"
	"time"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
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
		encoder:   hefloat.NewEncoder(params, 128),
		evaluator: hefloat.NewEvaluator(params, evk),
	}
}

func newTestVector(tc testContext, key rlwe.EncryptionKey, a, b complex128) (values []*bignum.Complex, el *estimator.Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	prec := tc.encoder.Prec()

	values = make([]*bignum.Complex, tc.params.MaxSlots())

	for i := range values {
		values[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
		}
	}

	values[0][0].SetFloat64(1)
	values[0][1].SetFloat64(0)

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
	default:
		panic("INVALID ENCRYPION KEY")
	}

	return
}

func TestEstimator(t *testing.T) {

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            14,
		LogQ:            []int{55, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	tc := newTestContext(params)

	//testChebyshevPowers(tc, t)
	//testKeySwitching(tc, t)
	//testChebyshevPolynomialEvaluation(tc, t)
	testLinearTransformation(tc, t)
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

func testLinearTransformation(tc testContext, t *testing.T) {

	params := tc.params

	mulCmplx := bignum.NewComplexMultiplier().Mul

	add := func(a, b, c []*bignum.Complex) {
		for i := range c {
			if a[i] != nil && b[i] != nil {
				c[i].Add(a[i], b[i])
			}
		}
	}

	muladd := func(a, b, c []*bignum.Complex) {
		tmp := &bignum.Complex{new(big.Float), new(big.Float)}
		for i := range c {
			if a[i] != nil && b[i] != nil {
				mulCmplx(a[i], b[i], tmp)
				c[i].Add(c[i], tmp)
			}
		}
	}

	prec := tc.encoder.Prec()

	newVec := func(size int) (vec []*bignum.Complex) {
		vec = make([]*bignum.Complex, size)
		for i := range vec {
			vec[i] = &bignum.Complex{new(big.Float).SetPrec(prec), new(big.Float).SetPrec(0)}
		}
		return
	}

	t.Run("LinearTransformation", func(t *testing.T) {

		statsHave := newPrecisionStats()
		statsWant := newPrecisionStats()

		d := 8

		for w := 0; w < d; w++ {
			values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

			slots := ct.Slots()

			nonZeroDiags := []int{-1, 0, 1, 2, 3, 4, 15}

			one := new(big.Float).SetInt64(1)
			zero := new(big.Float)

			diagonals := make(hefloat.Diagonals[*bignum.Complex])
			for _, i := range nonZeroDiags {
				diagonals[i] = make([]*bignum.Complex, slots)
				for j := 0; j < slots; j++ {
					diagonals[i][j] = &bignum.Complex{one, zero}
				}
			}

			for i := range diagonals {
				if i < 0 {
					diagonals[slots+i] = diagonals[i]
					delete(diagonals, i)
				}
			}

			lt := estimator.LinearTransformation{
				Scale:                    params.GetOptimalScalingFactor(ct.Scale, params.DefaultScale(), ct.Level()),
				LogSlots:                 ct.LogDimensions.Cols,
				LogBabyStepGianStepRatio: 1,
				Value:                    diagonals,
			}

			el.EvaluateLinearTransformation(lt)
			el.Decrypt()
			el.Normalize()

			ltparams := hefloat.LinearTransformationParameters{
				DiagonalsIndexList:       diagonals.DiagonalsIndexList(),
				LevelQ:                   ct.Level(),
				LevelP:                   params.MaxLevelP(),
				Scale:                    params.GetOptimalScalingFactor(ct.Scale, params.DefaultScale(), ct.Level()),
				LogDimensions:            ct.LogDimensions,
				LogBabyStepGianStepRatio: 1,
			}

			// Allocate the linear transformation
			linTransf := hefloat.NewLinearTransformation(params, ltparams)

			// Encode on the linear transformation
			require.NoError(t, hefloat.EncodeLinearTransformation[*bignum.Complex](tc.encoder, diagonals, linTransf))

			galEls := hefloat.GaloisElementsForLinearTransformation(params, ltparams)

			evk := rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(galEls, tc.sk)...)

			ltEval := hefloat.NewLinearTransformationEvaluator(tc.evaluator.WithKey(evk))

			require.NoError(t, ltEval.Evaluate(ct, linTransf, ct))

			values = diagonals.Evaluate(values, newVec, add, muladd)

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
					el = pb.Value[1<<n]
					el.Decrypt()
					el.Normalize()
				})
			})

			RunTimed("encrypted", func() {
				eval := tc.evaluator
				pb := hefloat.NewPowerBasis(ct, bignum.Chebyshev)
				if err := pb.GenPower(1<<n, false, eval); err != nil {
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

	t.Run("Chebyshev512Power/pk", func(t *testing.T) {

		statsHave := newPrecisionStats()
		statsWant := newPrecisionStats()

		d := 16

		for w := 0; w < d; w++ {

			values, el, _, ct := newTestVector(tc, tc.pk, -1, 1)

			n := 11

			RunProfiled(func() {
				RunTimed("estimator", func() {
					pb := estimator.NewPowerBasis(el, bignum.Chebyshev)
					pb.GenPower(1<<n, false)
					el = pb.Value[1<<n]
					el.Decrypt()
					el.Normalize()
				})
			})

			RunTimed("encrypted", func() {
				eval := tc.evaluator
				pb := hefloat.NewPowerBasis(ct, bignum.Chebyshev)
				if err := pb.GenPower(1<<n, false, eval); err != nil {
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

func testChebyshevPolynomialEvaluation(tc testContext, t *testing.T) {

	sigmoid := func(x float64) (y float64) {
		return 1 / (math.Exp(-x) + 1)
	}

	poly := hefloat.NewPolynomial(GetChebyshevPoly(1, 15, sigmoid))

	scalar, constant := poly.ChangeOfBasis()

	t.Run("ChebyshevPolynomialEvaluation/sk", func(t *testing.T) {

		statsHave := newPrecisionStats()
		statsWant := newPrecisionStats()

		d := 1

		for w := 0; w < d; w++ {

			values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

			RunProfiled(func() {
				RunTimed("estimator", func() {
					var err error
					var tmp *estimator.Element

					if scalar.Cmp(new(big.Float).SetInt64(1)) != 0 {
						el.Mul(el, scalar)
						el.Add(el, constant)
						el.Rescale()
					}

					if tmp, err = el.EvaluatePolynomial(poly, el.Scale); err != nil {
						t.Fatal(err)
					}
					el = tmp
					el.Decrypt()
					el.Normalize()
				})
			})

			RunTimed("encrypted", func() {

				if scalar.Cmp(new(big.Float).SetInt64(1)) != 0 {
					tc.evaluator.Mul(ct, scalar, ct)
					tc.evaluator.Add(ct, constant, ct)
					tc.evaluator.Rescale(ct, ct)
				}

				polyEval := hefloat.NewPolynomialEvaluator(tc.params, tc.evaluator)
				var err error
				if ct, err = polyEval.Evaluate(ct, poly, ct.Scale); err != nil {
					t.Fatal(err)
				}
			})

			fmt.Println(values[0])

			for i := range values {
				values[i] = poly.Evaluate(values[i])
			}

			fmt.Println(values[0])
			fmt.Println(el.Value[0][0])

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

// GetChebyshevPoly returns the Chebyshev polynomial approximation of f the
// in the interval [-K, K] for the given degree.
func GetChebyshevPoly(K float64, degree int, f64 func(x float64) (y float64)) bignum.Polynomial {

	FBig := func(x *big.Float) (y *big.Float) {
		xF64, _ := x.Float64()
		return new(big.Float).SetPrec(x.Prec()).SetFloat64(f64(xF64))
	}

	var prec uint = 128

	interval := bignum.Interval{
		A:     *bignum.NewFloat(-K, prec),
		B:     *bignum.NewFloat(K, prec),
		Nodes: degree,
	}

	// Returns the polynomial.
	return bignum.ChebyshevApproximation(FBig, interval)
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
