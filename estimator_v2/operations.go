package estimator

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Add adds op0 to op1 and writes the result on p, also returns p.
func (p *Element) Add(op0 *Element, op1 rlwe.Operand) *Element {

	m0 := op0.Value
	m2 := p.Value

	switch op1 := op1.(type) {
	case *Element:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		m1 := op1.Value

		for i := range m2 {
			m2[i].Add(m0[i], m1[i])
		}

		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = utils.Max(op0.Degree, op1.Degree)

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		bComplex := bignum.ToComplex(op1, prec)

		bComplex[0].Mul(bComplex[0], &p.Scale)
		bComplex[1].Mul(bComplex[1], &p.Scale)

		for i := range m2 {
			m2[i].Add(m0[i], bComplex)
		}

	default:
		panic(fmt.Errorf("invalid op1.(type): must be *Element, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}

	return p
}

func (p *Element) Mul(op0 *Element, op1 rlwe.Operand) *Element {
	m0 := op0.Value
	m2 := p.Value

	mul := bignum.NewComplexMultiplier().Mul

	switch op1 := op1.(type) {
	case *Element:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic(fmt.Errorf("invalid input dimensions: do not match receiver dimension"))
		}

		if op0.Degree+op1.Degree > 2 {
			panic(fmt.Errorf("invalid input degree: sum cannot exceed 2"))
		}

		// m0 * m1
		m1 := op1.Value
		for i := range m2 {
			mul(m0[i], m1[i], m2[i])
		}

		p.Scale = *p.Scale.Mul(&op0.Scale, &op1.Scale)
		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = op0.Degree + op1.Degree

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		bComplex := bignum.ToComplex(op1, m2[0].Prec())

		for i := range m2 {
			mul(m0[i], bComplex, m2[i])
		}

		p.Scale = *p.Scale.Mul(&op0.Scale, p.Q[p.Level])

	default:
		panic(fmt.Errorf("invalid op1.(type): must be *Element, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}

	return p
}
