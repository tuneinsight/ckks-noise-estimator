package main

import(
	"fmt"
	"math"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	bootEst "github.com/tuneinsight/ckks-bootstrapping-precision/estimator/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func main(){
	LogN := 15
	LogScale := 40
	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 40}, 
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	btpParametersLit := bootstrapping.ParametersLiteral{
		LogN: utils.Pointy(params.LogN()),
		LogP: []int{61, 61, 61, 61},
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
	
	var eval *bootstrapping.Evaluator
	if eval, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	evalEst := bootEst.NewEvaluator(btpParams)

	estParamsResidual := evalEst.ResidualParameters

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 8; i++{

		fmt.Println(i)

		values, el, _, ct := estParamsResidual.NewTestVector(ecd, enc, -1-1i, 1+1i)
		
		// Encrypted
		ct, _, _ = eval.ScaleDown(ct)
		ct, _ = eval.ModUp(ct)
		ctReal, ctImag, _ := eval.CoeffsToSlots(ct)

		ctReal, _ = eval.EvalMod(ctReal)
		ctImag, _ = eval.EvalMod(ctImag)

		ct, _ = eval.SlotsToCoeffs(ctReal, ctImag)
		
		// Simulated
		evalEst.ScaleDown(el)
		evalEst.ModUp(el)
		elReal, elImag := evalEst.CoeffsToSlots(el)
		elReal, elImag = evalEst.EvalMod(elReal), evalEst.EvalMod(elImag)

		el = evalEst.SlotsToCoeffs(elReal, elImag)

		elReal.Decrypt()
		elReal.Normalize()

		el.Decrypt()
		el.Normalize()

		eval.Evaluator.DropLevel(ctReal, ctReal.Level()-params.MaxLevel())

		vv := make([]*bignum.Complex, ctReal.Slots())
		ecd.Decode(dec.DecryptNew(ctReal), vv)

		zel, zct := 0, 0
		for i := range vv{
			if i < 8{
				fmt.Println(i, elReal.Value[0][i], vv[i])
			}
			
			
			if math.Round(math.Abs(real(elReal.Value[0][i].Complex128()))) == 0{
				zel++
			}

			if math.Round(math.Abs(real(vv[i].Complex128()))) == 0{
				zct++
			}
		}

		std0, std1 := estimator.Log2STD(elReal.Value[0])
		std2, std3 := estimator.Log2STD(vv)
		avg0, avg1 := estimator.Log2MAX(elReal.Value[0])
		avg2, avg3 := estimator.Log2MAX(vv)

		fmt.Println("Zero", zel, zct)
		fmt.Println("STD", std0, std2)
		fmt.Println("STD", std1, std3)
		fmt.Println("MAX", avg0, avg2)
		fmt.Println("MAX", avg1, avg3)

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