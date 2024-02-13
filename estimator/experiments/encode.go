package main

import (
	"fmt"

	"estimator"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
)

func main() {

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{55},
		LogP:            []int{60},
		LogDefaultScale: 45,
	})

	if err != nil {
		panic(err)
	}

	ecd := hefloat.NewEncoder(params)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 1; i++ {

		values, el, pt, _ := est.NewTestVector(ecd, nil, -1-1i, 1+1i)

		pWant := hefloat.GetPrecisionStats(params, ecd, nil, values, est.Decrypt(el), 0, false)
		pHave := hefloat.GetPrecisionStats(params, ecd, nil, values, pt, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())
}
