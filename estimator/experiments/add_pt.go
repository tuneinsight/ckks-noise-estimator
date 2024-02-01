package main

import (
	"fmt"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
)

func main() {

	LogN := 16
	LogScale := 45

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55},
		LogP:            []int{60},
		LogDefaultScale: LogScale,
	})

	if err != nil {
		panic(err)
	}

	ecd := hefloat.NewEncoder(params)

	kgen := hefloat.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	enc := hefloat.NewEncryptor(params, pk)
	dec := hefloat.NewDecryptor(params, sk)
	eval := hefloat.NewEvaluator(params, nil)

	estParams := estimator.NewParameters(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 128; i++ {

		values0, el0, _, ct0 := estParams.NewTestVector(ecd, enc, -1-1i, 1+1i)
		values1, el1, pt1, _ := estParams.NewTestVector(ecd, nil, -1-1i, 1+1i)

		for j := range values0 {
			values0[j].Add(values0[j], values1[j])
		}

		if err := eval.Add(ct0, pt1, ct0); err != nil {
			panic(err)
		}

		el0.Add(el0, el1)

		el0.Decrypt()
		el0.Normalize()

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values0, el0.Value[0], 0, false)
		pHave := hefloat.GetPrecisionStats(params, ecd, dec, values0, ct0, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())

	fmt.Println(estimator.ToLaTeXTable(LogN, LogScale, statsWant, statsHave))
}
