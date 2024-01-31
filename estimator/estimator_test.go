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
		LogQ:            []int{55, 45, 45, 45, 45},
		LogP:            []int{60, 60},
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
	//testLinearTransformation(tc, t)
	testDFT(tc, t)
}



func testDFT(tc testContext, t *testing.T) {
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

	t.Run("CoeffsToSlots", func(t *testing.T) {

		statsHave := estimator.Stats{}
		statsWant := estimator.Stats{}

		DFTMatrixLiteral := hefloat.DFTMatrixLiteral{
			LogSlots:     params.LogMaxSlots(),
			Type:         hefloat.HomomorphicEncode,
			Format:       hefloat.RepackImagAsReal,
			LevelQ:       params.MaxLevelQ(),
			LevelP:       params.MaxLevelP(),
			Levels:       []int{1, 1, 1, 1},
			LogBSGSRatio: 0,
		}

		DFTMatrix := estimator.DFTMatrix{
			DFTMatrixLiteral: DFTMatrixLiteral,
		}

		DFTMatrix.GenMatrices(params.LogN(), 128)

		DFTMatrixHeFloat, err := hefloat.NewDFTMatrixFromLiteral(params, DFTMatrixLiteral, tc.encoder)
		require.NoError(t, err)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(append(DFTMatrixHeFloat.GaloisElements(params), params.GaloisElementForComplexConjugation()), tc.sk)...)

		eval := hefloat.NewEvaluator(params, evk)
		hdftEval := hefloat.NewDFTEvaluator(params, eval)

		d := 1

		for w := 0; w < d; w++ {

			values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

			elReal, elImag := el.CoeffsToSlots(DFTMatrix)

			ctReal, ctImag, err := hdftEval.CoeffsToSlotsNew(ct, DFTMatrixHeFloat)
			require.NoError(t, err)

			for i := range DFTMatrix.Value {
				values = DFTMatrix.Value[i].Evaluate(values, newVec, add, muladd)
			}

			if elImag != nil{

				elReal.Decrypt()
				elReal.Normalize()

				two := new(big.Float).SetInt64(2)
				for i := range values{
					values[i][0].Mul(values[i][0], two)
					values[i][1].Mul(values[i][1], two)
				}

				valuesReal := make([]*bignum.Complex, len(values))
				for i := range valuesReal{
					valuesReal[i] = &bignum.Complex{
						values[i][0],
						new(big.Float),
					}
				}

				pWantReal := hefloat.GetPrecisionStats(tc.params, tc.encoder, nil, valuesReal, elReal.Value[0], 0, false)
				pHaveReal := hefloat.GetPrecisionStats(tc.params, tc.encoder, tc.decryptor, valuesReal, ctReal, 0, false)

				statsWant.Add(pWantReal)
				statsHave.Add(pHaveReal)

				valuesImag := make([]*bignum.Complex, len(values))
				for i := range valuesImag{
					valuesImag[i] = &bignum.Complex{
						values[i][1],
						new(big.Float),
					}
				}

				elImag.Decrypt()
				elImag.Normalize()

				pWantImag := hefloat.GetPrecisionStats(tc.params, tc.encoder, nil, valuesImag, elImag.Value[0], 0, false)
				pHaveImag := hefloat.GetPrecisionStats(tc.params, tc.encoder, tc.decryptor, valuesImag, ctImag, 0, false)
				statsWant.Add(pWantImag)
				statsHave.Add(pHaveImag)

			}else{

				elReal.Decrypt()
				elReal.Normalize()
				
				pWant := hefloat.GetPrecisionStats(tc.params, tc.encoder, nil, values, elReal.Value[0], 0, false)
				pHave := hefloat.GetPrecisionStats(tc.params, tc.encoder, tc.decryptor, values, ctReal, 0, false)
				statsWant.Add(pWant)
				statsHave.Add(pHave)
			}
		}

		statsWant.Finalize()
		statsHave.Finalize()

		fmt.Println(statsWant.String())
		fmt.Println(statsHave.String())
	})
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

		slots := params.MaxSlots()

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
			Scale:                    params.GetOptimalScalingFactor(params.DefaultScale(), params.DefaultScale(), params.MaxLevel()),
			LogSlots:                 params.LogMaxSlots(),
			LogBabyStepGianStepRatio: 1,
			Value:                    diagonals,
		}

		ltparams := hefloat.LinearTransformationParameters{
			DiagonalsIndexList:       diagonals.DiagonalsIndexList(),
			LevelQ:                   params.MaxLevel(),
			LevelP:                   params.MaxLevelP(),
			Scale:                    lt.Scale,
			LogDimensions:            params.LogMaxDimensions(),
			LogBabyStepGianStepRatio: 1,
		}

		// Allocate the linear transformation
		linTransf := hefloat.NewLinearTransformation(params, ltparams)

		// Encode on the linear transformation
		require.NoError(t, hefloat.EncodeLinearTransformation[*bignum.Complex](tc.encoder, diagonals, linTransf))

		galEls := hefloat.GaloisElementsForLinearTransformation(params, ltparams)

		evk := rlwe.NewMemEvaluationKeySet(nil, tc.kgen.GenGaloisKeysNew(galEls, tc.sk)...)

		ltEval := hefloat.NewLinearTransformationEvaluator(tc.evaluator.WithKey(evk))

		statsHave := estimator.Stats{}
		statsWant := estimator.Stats{}

		d := 8

		for w := 0; w < d; w++ {

			values, el, _, ct := newTestVector(tc, tc.sk, -1, 1)

			el.EvaluateLinearTransformation(lt)
			el.Decrypt()
			el.Normalize()

			require.NoError(t, ltEval.Evaluate(ct, linTransf, ct))

			values = diagonals.Evaluate(values, newVec, add, muladd)

			pWant := hefloat.GetPrecisionStats(tc.params, tc.encoder, nil, values, el.Value[0], 0, false)
			pHave := hefloat.GetPrecisionStats(tc.params, tc.encoder, tc.decryptor, values, ct, 0, false)

			statsWant.Add(pWant)
			statsHave.Add(pHave)
		}

		statsWant.Finalize()
		statsHave.Finalize()

		fmt.Println(statsWant.String())
		fmt.Println(statsHave.String())
	})
}

func testKeySwitching(tc testContext, t *testing.T) {

	t.Run("Rotate/sk", func(t *testing.T) {

		statsHave := estimator.Stats{}
		statsWant := estimator.Stats{}

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

			statsWant.Add(pWant)
			statsHave.Add(pHave)
		}

		statsWant.Finalize()
		statsHave.Finalize()

		fmt.Println(statsWant.String())
		fmt.Println(statsHave.String())
	})
}

func testChebyshevPowers(tc testContext, t *testing.T) {

	t.Run("Chebyshev512Power/sk", func(t *testing.T) {

		statsHave := estimator.Stats{}
		statsWant := estimator.Stats{}

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

			statsWant.Add(pWant)
			statsHave.Add(pHave)
		}

		statsWant.Finalize()
		statsHave.Finalize()

		fmt.Println(statsWant.String())
		fmt.Println(statsHave.String())
	})

	t.Run("Chebyshev512Power/pk", func(t *testing.T) {

		statsHave := estimator.Stats{}
		statsWant := estimator.Stats{}

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

			statsWant.Add(pWant)
			statsHave.Add(pHave)
		}

		statsWant.Finalize()
		statsHave.Finalize()

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

		statsHave := estimator.Stats{}
		statsWant := estimator.Stats{}

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

			statsWant.Add(pWant)
			statsHave.Add(pHave)
		}

		statsWant.Finalize()
		statsHave.Finalize()

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
