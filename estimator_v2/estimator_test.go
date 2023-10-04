package estimator_test

import (
	"testing"

	estimator "github.com/tuneinsight/ckks-bootstrapping-precision/estimator_v2"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
	"github.com/tuneinsight/lattigo/v4/utils/sampling"
)

func TestEstimator(t *testing.T) {

	LogN := 14
	Sigma := 3.2
	var prec uint = 256
	H := 512
	LogScale := 90

	p := estimator.NewParameters(LogN, Sigma, H, LogScale, prec)

	values0 := make([]*bignum.Complex, p.MaxSlots())
	values1 := make([]*bignum.Complex, p.MaxSlots())

	a := complex(-1, -1)
	b := complex(1, 1)

	for i := range values0 {
		values0[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
		}
		values1[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
		}
	}

	pt0 := estimator.NewPlaintext(p, values0)
	pt0.Add(&pt0, p.RoundingNoise())
	pt0.Add(&pt0, p.EncryptionNoiseSk())

	pt1 := estimator.NewPlaintext(p, values1)
	pt1.Add(&pt1, p.RoundingNoise())
	pt1.Add(&pt1, p.EncryptionNoiseSk())

	pt0.Mul(&pt0, &pt1)

	mul := bignum.NewComplexMultiplier().Mul
	for i := range values0 {
		mul(values0[i], values1[i], values0[i])
	}

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{60},
		LogDefaultScale: LogScale,
	})

	if err != nil {
		t.Fatal(err)
	}

	ckks.VerifyTestVectors(params, ckks.NewEncoder(params, prec), nil, values0, pt0.Value, params.LogDefaultScale(), nil, true, t)
}
