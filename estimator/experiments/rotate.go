package main

import (
	"fmt"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils"
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

	k := 1

	evk := rlwe.NewMemEvaluationKeySet(nil, kgen.GenGaloisKeyNew(params.GaloisElement(k), sk))

	eval := hefloat.NewEvaluator(params, evk)

	estParams := estimator.NewParameters(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 128; i++ {

		values, el, _, ct := estParams.NewTestVector(ecd, enc, -1-1i, 1+1i)

		utils.RotateSliceInPlace(values, k)

		if err := eval.Rotate(ct, k, ct); err != nil {
			panic(err)
		}

		el.Rotate(k)
		el.Decrypt()
		el.Normalize()

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
