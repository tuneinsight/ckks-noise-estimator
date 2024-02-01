package main

import (
	"fmt"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
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

	estParams := estimator.NewParameters(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 128; i++ {

		values, el, pt, _ := estParams.NewTestVector(ecd, nil, -1-1i, 1+1i)

		el.Normalize()

		pWant := hefloat.GetPrecisionStats(params, ecd, nil, values, el.Value[0], 0, false)
		pHave := hefloat.GetPrecisionStats(params, ecd, nil, values, pt, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())
}
