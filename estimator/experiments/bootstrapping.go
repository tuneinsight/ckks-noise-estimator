package main

import (
	"fmt"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	bootEst "github.com/tuneinsight/ckks-bootstrapping-precision/estimator/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"

	//"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

func main() {
	LogN := 14
	LogScale := 45
	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 45},
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
		Xs:   params.Xs(),
	}

	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}

	fmt.Println(btpParams.BootstrappingParameters.LogQ(), btpParams.BootstrappingParameters.LogP())
	fmt.Println(btpParams.BootstrappingParameters.MaxLevel())

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

	for i := 0; i < 1; i++ {

		fmt.Println(i)

		values, el, _, ct := estParamsResidual.NewTestVector(ecd, enc, -1-1i, 1+1i)

		// Encrypted
		if ct, err = eval.Evaluate(ct); err != nil {
			panic(err)
		}

		// Simulated
		if el, err = evalEst.Bootstrap(el); err != nil {
			panic(err)
		}

		valuesWant := evalEst.ResidualParameters.Decrypt(el)
		fmt.Println(valuesWant[0])
		fmt.Println(values[0])

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values, valuesWant, 0, false)
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
