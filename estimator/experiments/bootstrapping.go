package main

import(
	"fmt"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	bootEst "github.com/tuneinsight/ckks-bootstrapping-precision/estimator/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	//"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
)

func main(){
	LogN := 16
	LogScale := 40
	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 40}, 
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		//Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	btpParametersLit := bootstrapping.ParametersLiteral{
		LogN: utils.Pointy(params.LogN()),
		LogP: []int{61, 61, 61, 61},
		//SlotsToCoeffsFactorizationDepthAndLogScales: [][]int{{45}, {45}, {45}},
		Xs: params.Xs(),
	}

	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	ecd := hefloat.NewEncoder(params)
	enc := rlwe.NewEncryptor(params, pk)
	dec := rlwe.NewDecryptor(params, sk)

	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}

	fmt.Println(len(evk.GaloisKeys))
	
	var eval *bootstrapping.Evaluator
	if eval, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	evalEst := bootEst.NewEvaluator(btpParams)

	estParamsResidual := evalEst.ResidualParameters

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 128; i++{

		fmt.Println(i)

		values, el, _, ct := estParamsResidual.NewTestVector(ecd, enc, -1-1i, 1+1i)
		
		// Encrypted
		ct, _ = eval.Evaluate(ct)

		// Simulated
		evalEst.ScaleDown(el)

		evalEst.ModUp(el)

		elReal, elImag := evalEst.CoeffsToSlots(el)

		dd := elReal.CopyNew()
		dd.Decrypt()
		dd.Normalize()
		fmt.Println(dd.Value[0][0])

		elReal = evalEst.EvalMod(elReal)

		dd = elReal.CopyNew()
		dd.Decrypt()
		dd.Normalize()
		fmt.Println(dd.Value[0][0])

		if elImag != nil{
			elImag = evalEst.EvalMod(elImag)
		}

		el = evalEst.SlotsToCoeffs(elReal, elImag)

		el.Decrypt()
		el.Normalize()

		fmt.Println(el.Value[0][0])
		fmt.Println(values[0])

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values, el.Value[0], 0, false)
		pHave := hefloat.GetPrecisionStats(params, ecd, dec, values, ct, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())

	fmt.Println(estimator.ToLaTeXTable(LogN, LogScale, statsWant, statsHave))

}