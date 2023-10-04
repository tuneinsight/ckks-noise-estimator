package estimator

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Plaintext struct {
	LogScale int
	Value    []*bignum.Complex
	Noise    []*bignum.Complex
}

func NewPlaintext[T ckks.Float](p Parameters, v []T) Plaintext {

	if len(v) > p.MaxSlots() {
		panic("len(v) > p.MaxSlots()")
	}

	Value := make([]*bignum.Complex, p.MaxSlots())

	prec := p.Prec

	for i := range v {
		Value[i] = bignum.ToComplex(v[i], prec)
	}

	for i := len(v); i < p.MaxSlots(); i++ {
		Value[i] = bignum.ToComplex(0, prec)
	}

	Noise := make([]*bignum.Complex, p.MaxSlots())
	for i := range Noise{
		Noise[i] = bignum.ToComplex(0, prec)
	}

	return Plaintext{
		LogScale: p.LogScale,
		Value:    Value,
	}
}
