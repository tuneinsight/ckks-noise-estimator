package main

import (
	"fmt"

	"estimator"
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
	sk := kgen.GenSecretKeyNew()
	enc := hefloat.NewEncryptor(params, sk)
	dec := hefloat.NewDecryptor(params, sk)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 1; i++ {

		fmt.Println(i)

		values, el, _, ct := est.NewTestVector(ecd, enc, -1-1i, 1+1i)

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values, est.Decrypt(el), 0, false)
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
