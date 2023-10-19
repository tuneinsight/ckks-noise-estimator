package estimator

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Element struct {
	*Parameters
	Degree int
	Level  int
	Scale  big.Float
	Value  []*bignum.Complex
}

func NewElement[T ckks.Float](p Parameters, v []T, Degree int, scale big.Float) Element {

	if len(v) > p.MaxSlots() {
		panic("len(v) > p.MaxSlots()")
	}

	Value := make([]*bignum.Complex, p.MaxSlots())

	for i := range v {
		Value[i] = bignum.ToComplex(v[i], prec)
		Value[i][0].Mul(Value[i][0], &scale)
		Value[i][1].Mul(Value[i][1], &scale)
	}

	for i := len(v); i < p.MaxSlots(); i++ {
		Value[i] = bignum.ToComplex(0, prec)
	}

	return Element{
		Parameters: &p,
		Degree:     Degree,
		Level:      p.MaxLevel(),
		Scale:      scale,
		Value:      Value,
	}
}

func (p *Element) AddEncodingNoise() {
	e := p.RoundingNoise()
	value := p.Value
	for j := range value {
		value[j].Add(value[j], e[j])
	}
}

func (p *Element) AddRoundingNoise() {
	for i := 0; i < p.Degree+1; i++ {

		e := p.RoundingNoise()

		if i > 0 {
			mul := bignum.NewComplexMultiplier().Mul

			sk := p.Sk[i-1]

			for j := range e {
				mul(e[j], sk[j], e[j])
			}
		}

		value := p.Value
		for j := range value {
			value[j].Add(value[j], e[j])
		}
	}
}

func (p *Element) AddEncryptionNoiseSk() {
	e := p.EncryptionNoiseSk()
	value := p.Value
	for i := range value {
		value[i].Add(value[i], e[i])
	}
}

func (p *Element) AddEncryptionNoisePk() {
	e := p.EncryptionNoisePk()
	value := p.Value
	for i := range value {
		value[i].Add(value[i], e[i])
	}
}

func (p *Element) AddKeySwitchingNoise(eCt []*bignum.Complex) {
	eKS := p.KeySwitchingNoise(p.Level, eCt, p.Sk[0])
	value := p.Value
	for i := range value {
		value[i].Add(value[i], eKS[i])
	}
}

func (p *Element) Relinearize() {
	if p.Degree != 2 {
		panic("cannot Relinearize: element.Degree != 2")
	}
	p.AddRoundingNoise()
	p.Degree = 1
}

func (p *Element) Rescale() (err error) {
	if p.Level == 0 {
		return fmt.Errorf("cannot Rescale: element already at level 0")
	}

	Q := p.Q[p.Level]

	p.Scale = *p.Scale.Quo(&p.Scale, Q)

	m := p.Value
	for i := range m {
		m[i][0].Quo(m[i][0], Q)
		m[i][1].Quo(m[i][1], Q)
	}

	p.AddRoundingNoise()

	return
}

func (p *Element) Normalize() {
	Value := p.Value
	scale := &p.Scale
	for i := range Value {
		Value[i][0].Quo(Value[i][0], scale)
		Value[i][1].Quo(Value[i][1], scale)
	}
}
