package estimator

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

type Element struct {
	*Parameters
	Degree int
	Level  int
	Scale  rlwe.Scale
	Value  [3][]*bignum.Complex //(m + e0, e1, e2)
}

func (p Element) CopyNew() *Element {

	Value := [3][]*bignum.Complex{}

	for i := range Value {

		v0 := p.Value[i]
		v1 := make([]*bignum.Complex, len(v0))

		for j := range v0 {
			v1[j] = &bignum.Complex{new(big.Float), new(big.Float)}
			v1[j].Set(v0[j])
		}

		Value[i] = v1
	}

	return &Element{
		Parameters: p.Parameters,
		Degree:     p.Degree,
		Level:      p.Level,
		Scale:      p.Scale,
		Value:      Value,
	}
}

func NewElement[T ckks.Float](p Parameters, v []T, Degree int, scale rlwe.Scale) *Element {

	if len(v) > p.MaxSlots() {
		panic("len(v) > p.MaxSlots()")
	}

	e0 := make([]*bignum.Complex, p.MaxSlots())
	e1 := make([]*bignum.Complex, p.MaxSlots())
	e2 := make([]*bignum.Complex, p.MaxSlots())

	for i := range v {
		e0[i] = bignum.ToComplex(v[i], prec)
		e0[i][0].Mul(e0[i][0], &scale.Value)
		e0[i][1].Mul(e0[i][1], &scale.Value)
	}

	for i := len(v); i < p.MaxSlots(); i++ {
		e0[i] = bignum.ToComplex(0, prec)
	}

	for i := range e1 {
		e1[i] = bignum.ToComplex(0, prec)
		e2[i] = bignum.ToComplex(0, prec)
	}

	return &Element{
		Parameters: &p,
		Degree:     Degree,
		Level:      p.MaxLevel(),
		Scale:      scale,
		Value:      [3][]*bignum.Complex{e0, e1, e2},
	}
}

// Decrypt decrypts the element by evaluating <(el0, el1, el2), (1, sk, sk^2)>.
// Also sets the degree of the elemen to 0.
func (p *Element) Decrypt() {

	v := p.Value[0]

	for i := 1; i < p.Degree+1; i++ {
		mul := bignum.NewComplexMultiplier().Mul
		sk := p.Sk[i-1]
		ei := p.Value[i]
		for j := range ei {
			mul(ei[j], sk[j], ei[j])
		}

		for j := range v {
			v[j].Add(v[j], ei[j])
		}
	}

	p.Degree = 0
}

// AddEncodingNoise adds the encoding noise, which is
// {round(1/2), 0}.
func (p *Element) AddEncodingNoise() {
	e := p.RoundingNoise()
	value := p.Value[0]
	for j := range value {
		value[j].Add(value[j], e[j])
	}
}

// AddRoundingNoise adds the rounding noise,
// which is {round(1/2), round(1/2)}.
func (p *Element) AddRoundingNoise() {
	for i := 0; i < p.Degree+1; i++ {
		r := p.RoundingNoise()
		ei := p.Value[i]
		for j := range ei {
			ei[j].Add(ei[j], r[j])
		}
	}
}

// AddEncryptionNoiseSk adds the encryption noise
// from SK encryption, which is {sigma, 0}.
func (p *Element) AddEncryptionNoiseSk() {
	e := p.NormalNoise(p.Sigma)
	value := p.Value[0]
	for i := range value {
		value[i].Add(value[i], e[i])
	}
}

// AddEncryptionNoisePk adds the encryption noise
// from PK encryption, which is {round(1/2), round(1/2)}
// This assumes that P != 0 and that |N * e * sk| < P.
func (p *Element) AddEncryptionNoisePk() {
	p.Degree = utils.Max(1, p.Degree)
	p.AddRoundingNoise()
}

// AddAutomorphismNoise sets the noise to the key-switching noise, which is
// (p[0], p[1]) = (p[0] + p[1] * sk + round(sum(e_i * qalphai)/P), round(1/2))
func (p *Element) AddAutomorphismNoise() {
	e0, e1 := p.KeySwitchingNoise(p.Level, p.Value[1], p.Sk[0])
	m0, m1 := p.Value[0], p.Value[1]
	for i := range m0 {
		m0[i].Add(m0[i], e0[i])
		m1[i].Set(e1[i])
	}
}

// AddAutomorphismNoiseRaw sets the noise to the key-switching noise scaled by P, which is
// (p[0], p[1]) = (p[0] + p[1] * sk * P + sum(e_i * qalphai), 0)
// This assumes that p[0] is scaled by P.
func (p *Element) AddAutomorphismNoiseRaw() {
	e0 := p.KeySwitchingNoiseRaw(p.Level, p.Value[1], p.Sk[0])
	m0 := p.Value[0]
	m1 := p.Value[1]
	for i := range m0 {
		m0[i].Add(m0[i], e0[i])
		m1[i][0].SetFloat64(0)
		m1[i][1].SetFloat64(0)
	}
}

// AddRelinearizationNoise adds the relinearization noise, which is
// (p[0], p[1]) = (p[0] + p[2]*sk + round(sum(e_i qalphai)/P)), p[1] + round(1/2))
func (p *Element) AddRelinearizationNoise() {
	e0, e1 := p.KeySwitchingNoise(p.Level, p.Value[2], p.Sk[1])
	m0 := p.Value[0]
	m1 := p.Value[1]
	for i := range m0 {
		m0[i].Add(m0[i], e0[i])
		m1[i].Add(m1[i], e1[i])
	}
}

func (p *Element) Normalize() {
	Value := p.Value[0]
	scale := &p.Scale.Value
	for i := range Value {
		Value[i][0].Quo(Value[i][0], scale)
		Value[i][1].Quo(Value[i][1], scale)
	}
}
