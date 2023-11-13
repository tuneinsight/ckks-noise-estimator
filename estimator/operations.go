package estimator

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/core/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Add adds op0 to op1 and writes the result on p, also returns p.
func (p *Element) Add(op0 *Element, op1 rlwe.Operand) *Element {

	switch op1 := op1.(type) {
	case *Element:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = utils.Max(op0.Degree, op1.Degree)

		for i := 0; i < p.Degree+1; i++ {
			m0 := op0.Value[i]
			m1 := op1.Value[i]
			m2 := p.Value[i]
			for j := range m2 {
				m2[j].Add(m0[j], m1[j])
			}
		}

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		m0 := op0.Value[0]
		m2 := p.Value[0]

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
		m00 := op0.Value[0] // (m0 + e00)
		m10 := op1.Value[0] // (m1 + e10)
		m20 := p.Value[0]   // (m2 + e20)

		e01 := op0.Value[1] // e01
		e11 := op1.Value[1] // e11
		e21 := p.Value[1]   // e21

		e22 := p.Value[2] // e22

		tmp := bignum.NewComplex()

		for i := range m00 {

			// e22 = e01 * e11
			mul(e01[i], e11[i], e22[i])

			// (m0 + e00) * e11 + (m1 + e10) * e01
			mul(m00[i], e11[i], tmp)
			mul(m10[i], e01[i], e21[i])
			e21[i].Add(e21[i], tmp)

			// (m0 + e00) * (m1 + e11)
			mul(m00[i], m10[i], m20[i])
		}

		p.Scale = *p.Scale.Mul(&op0.Scale, &op1.Scale)
		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = op0.Degree + op1.Degree

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		bComplex := bignum.ToComplex(op1, p.Value[0][0].Prec())

		if !bComplex.IsInt(){
			bComplex[0].Mul(bComplex[0], p.Q[p.Level])
			bComplex[1].Mul(bComplex[1], p.Q[p.Level])
			p.Scale = *p.Scale.Mul(&op0.Scale, p.Q[p.Level])
		}else{
			p.Scale = op0.Scale
		}

		for i := 0; i < op0.Degree+1; i++ {
			m0 := op0.Value[i]
			m2 := p.Value[i]
			for i := range m2 {
				mul(m0[i], bComplex, m2[i])
			}
		}
		
		p.Degree = op0.Degree

	default:
		panic(fmt.Errorf("invalid op1.(type): must be *Element, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}

	return p
}

func (p *Element) Rotate(k int) *Element {

	if p.Degree != 1 {
		panic("cannot Rotate: Degree != 1")
	}

	utils.RotateSliceInPlace(p.Value[0], k)

	// p.Value[1]: noise of the second component (s term)
	// p.Sk[0]: sk^1
	p.AddKeySwitchingNoise(p.Value[1], p.Sk[0])

	return p
}

func (p *Element) Relinearize() {
	if p.Degree != 2 {
		panic("cannot Relinearize: element.Degree != 2")
	}

	// p.Value[2]: noise of the third component (s^2 term)
	// p.Sk[1]: Sk^2
	p.AddKeySwitchingNoise(p.Value[2], p.Sk[1])
	p.Degree = 1
}

// RotateHoisted applies a rotation without ModDown.
// Returned element is scaled by P.
func (p *Element) RotateHoisted(k int) {

	if p.Degree != 1 {
		panic("cannot Rotate: Degree != 1")
	}

	utils.RotateSliceInPlace(p.Value[0], k)

	// Scales first term by P
	value := p.Value[0]
	P := p.P
	for i := range value {
		value[i][0].Mul(value[i][0], P)
		value[i][1].Mul(value[i][1], P)
	}

	// p.Value[1]: noise of the second component (s term)
	// p.Sk[0]: sk^1
	// Added noise is scaled by P
	p.AddKeySwitchingNoiseRaw(p.Value[1], p.Sk[0])
}

// ModDown divides by P and adds rounding noise.
func (p *Element) ModDown() {
	p.DivideAndAddRoundingNoise(p.P)
}

// Rescale divides by Q[level] and adds rounding noise.
// Returns an error if already at level 0.
func (p *Element) Rescale() {
	if p.Level == 0 {
		panic(fmt.Errorf("cannot Rescale: element already at level 0"))
	}

	Q := p.Q[p.Level]

	p.Scale = *p.Scale.Quo(&p.Scale, Q)
	p.Level--

	p.DivideAndAddRoundingNoise(Q)
}

// DivideAndRound by P and adds rounding noise
func (p *Element) DivideAndAddRoundingNoise(P *big.Float) {

	for i := 0; i < p.Degree+1; i++ {
		m := p.Value[i]
		for j := range m {
			m[j][0].Quo(m[j][0], P)
			m[j][1].Quo(m[j][1], P)
		}
	}

	p.AddRoundingNoise()
}
