package main

import (
	"fmt"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func main() {

	LogN := 16
	LogScale := 45

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 45},
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

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	mul := bignum.NewComplexMultiplier().Mul

	for i := 0; i < 1; i++ {

		fmt.Println(i)

		values0, el0, _, ct0 := est.NewTestVector(ecd, enc, -1, 1)
		values1, el1, pt1, _ := est.NewTestVector(ecd, nil, -1, 1)

		for j := range values0 {
			mul(values0[j], values1[j], values0[j])
		}

		if err := eval.Mul(ct0, pt1, ct0); err != nil {
			panic(err)
		}

		if err = est.Mul(el0, el1, el0); err != nil {
			panic(err)
		}

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values0, est.Decrypt(el0), 0, false)
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
