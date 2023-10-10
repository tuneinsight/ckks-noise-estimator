package estimator_test

import (
	"time"
	"testing"

	estimator "github.com/tuneinsight/ckks-bootstrapping-precision/estimator_v2"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func TestEstimator(t *testing.T) {

	paramsCKKS, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: 14,
		LogQ: []int{55, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45},
		LogP: []int{61, 61, 61},
		LogDefaultScale: 45,
	})

	if err != nil{
		panic(err)
	}

	p := estimator.NewParameters(paramsCKKS)

	values0 := make([]*bignum.Complex, p.MaxSlots())

	a := complex(-1, 0)
	b := complex(1, 0)

	for i := range values0 {
		values0[i] = &bignum.Complex{
			estimator.NewFloat(sampling.RandFloat64(real(a), real(b))),
			estimator.NewFloat(sampling.RandFloat64(imag(a), imag(b))),
		}
	}


	el0 := estimator.NewElement(p, values0, 1, p.DefaultScale())
	el0.AddEncodingNoise()
	el0.AddEncryptionNoiseSk()

	now := time.Now()

	n := 9

	for k := 0; k < n; k++{
		el0.Mul(&el0, &el0)
		el0.Add(&el0, &el0)
		el0.Add(&el0, -1)
		el0.Relinearize()

		if err := el0.Rescale(); err != nil{
			panic(err)
		}
	}
	
	el0.Normalize()

	t.Log(el0.Value[0])

	t.Logf("Since: %s\n", time.Since(now))

	mul := bignum.NewComplexMultiplier().Mul

	for k := 0; k < n; k++{
		for i := range values0 {
			
			mul(values0[i], values0[i], values0[i])
			values0[i].Add(values0[i], values0[i])
			values0[i].Add(values0[i], &bignum.Complex{estimator.NewFloat(-1), estimator.NewFloat(0)})
		}
	}
	
	ckks.VerifyTestVectors(paramsCKKS, ckks.NewEncoder(paramsCKKS, 256), nil, values0, el0.Value, paramsCKKS.LogDefaultScale(), nil, true, t)
}
