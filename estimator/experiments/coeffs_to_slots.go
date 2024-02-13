package main

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"

	//"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func main() {

	LogN := 14
	LogScale := 55
	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 55, 55, 55, 55},
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		//Xs: ring.Ternary{H:192},
	})

	mulCmplx := bignum.NewComplexMultiplier().Mul

	var prec uint = 128

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

	newVec := func(size int) (vec []*bignum.Complex) {
		vec = make([]*bignum.Complex, size)
		for i := range vec {
			vec[i] = &bignum.Complex{new(big.Float).SetPrec(prec), new(big.Float).SetPrec(prec)}
		}
		return
	}

	ecd := hefloat.NewEncoder(params)

	kgen := hefloat.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	enc := hefloat.NewEncryptor(params, pk)
	dec := hefloat.NewDecryptor(params, sk)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

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

	DFTMatrix.GenMatrices(params.LogN(), prec)

	DFTMatrixHeFloat, err := hefloat.NewDFTMatrixFromLiteral(params, DFTMatrixLiteral, ecd)

	if err != nil {
		panic(err)
	}

	evk := rlwe.NewMemEvaluationKeySet(nil, kgen.GenGaloisKeysNew(append(DFTMatrixHeFloat.GaloisElements(params), params.GaloisElementForComplexConjugation()), sk)...)

	eval := hefloat.NewEvaluator(params, evk)
	hdftEval := hefloat.NewDFTEvaluator(params, eval)

	for i := 0; i < 1; i++ {

		fmt.Println(i)

		values, el, _, ct := est.NewTestVector(ecd, enc, -1-1i, 1+1i)

		elReal, elImag, err := est.CoeffsToSlotsNew(el, DFTMatrix)

		if err != nil {
			panic(err)
		}

		ctReal, ctImag, err := hdftEval.CoeffsToSlotsNew(ct, DFTMatrixHeFloat)

		if err != nil {
			panic(err)
		}

		for i := range DFTMatrix.Value {
			values = DFTMatrix.Value[i].Evaluate(values, newVec, add, muladd)
		}

		if elImag != nil {
			two := new(big.Float).SetInt64(2)
			for i := range values {
				values[i][0].Mul(values[i][0], two)
				values[i][1].Mul(values[i][1], two)
			}

			valuesReal := make([]*bignum.Complex, len(values))
			for i := range valuesReal {
				valuesReal[i] = &bignum.Complex{
					values[i][0],
					new(big.Float),
				}
			}

			pWantReal := hefloat.GetPrecisionStats(params, ecd, nil, valuesReal, est.Decrypt(elReal), 0, false)
			pHaveReal := hefloat.GetPrecisionStats(params, ecd, dec, valuesReal, ctReal, 0, false)

			statsWant.Add(pWantReal)
			statsHave.Add(pHaveReal)

			valuesImag := make([]*bignum.Complex, len(values))
			for i := range valuesImag {
				valuesImag[i] = &bignum.Complex{
					values[i][1],
					new(big.Float),
				}
			}

			pWantImag := hefloat.GetPrecisionStats(params, ecd, nil, valuesImag, est.Decrypt(elImag), 0, false)
			pHaveImag := hefloat.GetPrecisionStats(params, ecd, dec, valuesImag, ctImag, 0, false)
			statsWant.Add(pWantImag)
			statsHave.Add(pHaveImag)
		} else {

			pWantReal := hefloat.GetPrecisionStats(params, ecd, nil, values, est.Decrypt(elReal), 0, false)
			pHaveReal := hefloat.GetPrecisionStats(params, ecd, dec, values, ctReal, 0, false)

			statsWant.Add(pWantReal)
			statsHave.Add(pHaveReal)
		}
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())
}
