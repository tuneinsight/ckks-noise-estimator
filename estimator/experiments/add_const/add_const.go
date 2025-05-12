package main

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/ckks-noise-estimator"

	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func main() {

	LogN := 16
	LogScale := 45

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55},
		LogP:            []int{60},
		LogDefaultScale: LogScale,
	})

	if err != nil {
		panic(err)
	}

	ecd := ckks.NewEncoder(params)

	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	dec := ckks.NewDecryptor(params, sk)
	eval := ckks.NewEvaluator(params, nil)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	source := estimator.NewTestRand()

	for i := 0; i < 1; i++ {

		scalar := complex(source.Float64(-1, 1), source.Float64(-1, 1))

		values, el, _, ct := est.NewTestVector(ecd, pk, -1-1i, 1+1i)

		scalarR := new(big.Float).SetFloat64(real(scalar))
		scalarI := new(big.Float).SetFloat64(imag(scalar))

		for j := range values {
			values[j][0].Add(values[j][0], scalarR)
			values[j][1].Add(values[j][1], scalarI)
		}

		if err := eval.Add(ct, scalar, ct); err != nil {
			panic(err)
		}

		if err := est.Add(el, scalar, el); err != nil {
			panic(err)
		}

		pWant := ckks.GetPrecisionStats(params, ecd, nil, values, est.Decrypt(el), 0, false)
		pHave := ckks.GetPrecisionStats(params, ecd, dec, values, ct, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())

	fmt.Println(estimator.ToLaTeXTable(LogN, LogScale, statsWant, statsHave))
}
