package main

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)


func main(){

	LogN := 10
	LogScale := 45
	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 39, 39, 39}, 
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		Xs:              ring.Ternary{H: 192},
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

	estParams := estimator.NewParameters(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	DFTMatrixLiteral := hefloat.DFTMatrixLiteral{
		LogSlots:     params.LogMaxSlots(),
		Type:         hefloat.HomomorphicDecode,
		Format:       hefloat.RepackImagAsReal,
		LevelQ:       params.MaxLevelQ(),
		LevelP:       params.MaxLevelP(),
		Levels:       []int{1, 1, 1},
		LogBSGSRatio: 0,
	}

	DFTMatrix := estimator.DFTMatrix{
		DFTMatrixLiteral: DFTMatrixLiteral,
	}

	DFTMatrix.GenMatrices(params.LogN(), prec)

	DFTMatrixHeFloat, err := hefloat.NewDFTMatrixFromLiteral(params, DFTMatrixLiteral, ecd)
	
	if err != nil{
		panic(err)
	}

	evk := rlwe.NewMemEvaluationKeySet(nil, kgen.GenGaloisKeysNew(append(DFTMatrixHeFloat.GaloisElements(params)), sk)...)

	eval := hefloat.NewEvaluator(params, evk)
	hdftEval := hefloat.NewDFTEvaluator(params, eval)

	for i := 0; i < 16; i++ {

		fmt.Println(i)

		valuesReal, elReal, _, ctReal := estParams.NewTestVector(ecd, enc, -1, 1)
		valuesImag, elImag, _, ctImag := estParams.NewTestVector(ecd, enc, -1, 1)
		
		elReal = elReal.SlotsToCoeffs(elReal, elImag, DFTMatrix)

		ctReal, err := hdftEval.SlotsToCoeffsNew(ctReal, ctImag, DFTMatrixHeFloat)
		
		if err != nil{
			panic(err)
		}

		for i := range valuesReal{
			valuesReal[i][1].Set(valuesImag[i][0])
		}

		for i := range DFTMatrix.Value {
			valuesReal = DFTMatrix.Value[i].Evaluate(valuesReal, newVec, add, muladd)
		}

		elReal.Decrypt()
		elReal.Normalize()

		pWantReal := hefloat.GetPrecisionStats(params, ecd, nil, valuesReal, elReal.Value[0], 0, false)
		pHaveReal := hefloat.GetPrecisionStats(params, ecd, dec, valuesReal, ctReal, 0, false)

		statsWant.Add(pWantReal)
		statsHave.Add(pHaveReal)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())
}