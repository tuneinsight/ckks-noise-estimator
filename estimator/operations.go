package estimator

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// Add adds op0 to op1 and writes the result on p, also returns p.
func (p *Element) Add(op0 *Element, op1 rlwe.Operand) *Element {

	switch op1 := op1.(type) {
	case *Element:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		var tmp0, tmp1 *Element

		switch op0.Scale.Cmp(op1.Scale){
		case -1:

			ratio := op1.Scale.Div(op0.Scale)
			// Only scales up if int(ratio) >= 2
			if ratio.Float64() >= 2.0 {
				tmp0 = op0.CopyNew()
				tmp0.Mul(tmp0, ratio.BigInt())
				p.Scale = op1.Scale
			}else{
				tmp0 = op0
			}

			tmp1 = op1

		case 0:
			tmp0 = op0
			tmp1 = op1
			p.Scale = tmp0.Scale
		case 1:
			ratio := op0.Scale.Div(op1.Scale)
			// Only scales up if int(ratio) >= 2
			if ratio.Float64() >= 2.0 {
				tmp1 = op1.CopyNew()
				tmp1.Mul(tmp1, ratio.BigInt())
				p.Scale = op0.Scale
			}else{
				tmp1 = op1
			}

			tmp0 = op0
		}

		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = utils.Max(op0.Degree, op1.Degree)

		for i := 0; i < p.Degree+1; i++ {
			m0 := tmp0.Value[i]
			m1 := tmp1.Value[i]
			m2 := p.Value[i]
			for j := range m2 {
				m2[j].Add(m0[j], m1[j])
			}
		}

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		m0 := op0.Value[0]
		m2 := p.Value[0]

		bComplex := bignum.ToComplex(op1, prec)

		bComplex[0].Mul(bComplex[0], &p.Scale.Value)
		bComplex[1].Mul(bComplex[1], &p.Scale.Value)

		for i := range m2 {
			m2[i].Add(m0[i], bComplex)
		}

	default:
		panic(fmt.Errorf("invalid op1.(type): must be *Element, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}

	return p
}

// Sub subtracts op0 to op1 and writes the result on p, also returns p.
func (p *Element) Sub(op0 *Element, op1 rlwe.Operand) *Element {

	switch op1 := op1.(type) {
	case *Element:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		var tmp0, tmp1 *Element

		switch op0.Scale.Cmp(op1.Scale){
		case -1:

			ratio := op1.Scale.Div(op0.Scale)

			// Only scales up if int(ratio) >= 2
			if ratio.Float64() >= 2.0 {
				tmp0 = op0.CopyNew()
				tmp0.Mul(tmp0, ratio.BigInt())
				tmp0.Scale = op1.Scale
				p.Scale = op1.Scale
			}else{
				tmp0 = op0
			}

			tmp1 = op1

		case 0:
			tmp0 = op0
			tmp1 = op1
			p.Scale = tmp0.Scale
		case 1:
			ratio := op0.Scale.Div(op1.Scale)

			// Only scales up if int(ratio) >= 2
			if ratio.Float64() >= 2.0 {
				tmp1 = op1.CopyNew()
				tmp1.Mul(tmp1, ratio.BigInt())
				tmp1.Scale = op0.Scale
				p.Scale = op0.Scale
			}else{
				tmp1 = op1
			}

			tmp0 = op0
		}

		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = utils.Max(op0.Degree, op1.Degree)

		for i := 0; i < p.Degree+1; i++ {
			m0 := tmp0.Value[i]
			m1 := tmp1.Value[i]
			m2 := p.Value[i]
			for j := range m2 {
				m2[j].Sub(m0[j], m1[j])
			}
		}

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		m0 := op0.Value[0]
		m2 := p.Value[0]

		bComplex := bignum.ToComplex(op1, prec)

		bComplex[0].Mul(bComplex[0], &p.Scale.Value)
		bComplex[1].Mul(bComplex[1], &p.Scale.Value)

		for i := range m2 {
			m2[i].Sub(m0[i], bComplex)
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

		p.Scale = op0.Scale.Mul(op1.Scale)
		p.Level = utils.Min(op0.Level, op1.Level)
		p.Degree = op0.Degree + op1.Degree

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		bComplex := bignum.ToComplex(op1, p.Value[0][0].Prec())

		if !bComplex.IsInt(){
			bComplex[0].Mul(bComplex[0], p.Q[p.Level])
			bComplex[1].Mul(bComplex[1], p.Q[p.Level])
			p.Scale = op0.Scale.Mul(rlwe.NewScale(p.Q[p.Level]))
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

func (p *Element) MulThenAdd(op0 *Element, op1 rlwe.Operand) *Element {

	mul := bignum.NewComplexMultiplier().Mul

	switch op1 := op1.(type) {
	case *Element:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic(fmt.Errorf("invalid input dimensions: do not match receiver dimension"))
		}

		if op0.Degree+op1.Degree > 2 {
			panic(fmt.Errorf("invalid input degree: sum cannot exceed 2"))
		}

		resScale := op0.Scale.Mul(op1.Scale)
		if p.Scale.Cmp(resScale) == -1 {
			ratio := resScale.Div(p.Scale)
			// Only scales up if int(ratio) >= 2
			if ratio.Float64() >= 2.0 {
				p.Mul(p, &ratio.Value)
				p.Scale = resScale
			}
		}

		if op0.Degree == 1 && op1.Degree == 1{
			// m0 * m1
			m00 := op0.Value[0] // (m0 + e00)
			e01 := op0.Value[1] // e01

			m10 := op1.Value[0] // (m1 + e10)
			e11 := op1.Value[1] // e11

			m20 := p.Value[0]   // (m2 + e20)
			e21 := p.Value[1]   // e21
			e22 := p.Value[2] // e22

			tmp := bignum.NewComplex()
			acc := bignum.NewComplex()

			for i := range m00 {

				// e22 = e01 * e11
				mul(e01[i], e11[i], acc)
				e22[i].Add(e22[i], acc)

				// (m0 + e00) * e11 + (m1 + e10) * e01
				mul(m00[i], e11[i], tmp)
				mul(m10[i], e01[i], acc)

				e21[i].Add(e21[i], acc)
				e21[i].Add(e21[i], tmp)

				// (m0 + e00) * (m1 + e11)
				mul(m00[i], m10[i], acc)
				m20[i].Add(m20[i], acc)
			}

			p.Scale = op0.Scale.Mul(op1.Scale)
			p.Level = utils.Min(op0.Level, op1.Level)
			p.Degree = op0.Degree + op1.Degree

		}else{

			acc := bignum.NewComplex()

			p.Degree = utils.Max(op0.Degree, op1.Degree)
			p.Scale = op0.Scale.Mul(op1.Scale)
			p.Level = utils.Min(op0.Level, op1.Level)

			m := op1.Value[0]

			for i := 0; i < p.Degree+1; i++{
				v0 := p.Value[i]
				v1 := op0.Value[i]

				for j := range m{
					mul(v1[j], m[j], acc)
					v0[j].Add(v0[j], acc)
				}
			}
		}

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:

		p.Degree = utils.Max(p.Degree, op0.Degree)
		p.Level = utils.Min(p.Level, op0.Level)

		bComplex := bignum.ToComplex(op1, p.Value[0][0].Prec())

		var scaleRLWE rlwe.Scale

		// If op0 and opOut scales are identical, but the op1 is not a Gaussian integer then multiplies opOut by scaleRLWE.
		// This ensures noiseless addition with opOut = scaleRLWE * opOut + op0 * round(scalar * scaleRLWE).
		if cmp := op0.Scale.Cmp(p.Scale); cmp == 0 {

			if bComplex.IsInt() {
				scaleRLWE = rlwe.NewScale(1)
			} else {
				scaleRLWE = rlwe.NewScale(p.Q[p.Level])

				for i := 1; i < p.Parameters.Parameters.LevelsConsumedPerRescaling(); i++ {
					scaleRLWE = scaleRLWE.Mul(rlwe.NewScale(p.Q[p.Level-i]))
				}

				scaleInt := new(big.Int)
				scaleRLWE.Value.Int(scaleInt)
				p.Mul(p, scaleInt)
				p.Scale = p.Scale.Mul(scaleRLWE)
			}

		} else if cmp == -1 { // opOut.Scale > op0.Scale then the scaling factor for op1 becomes the quotient between the two scales
			scaleRLWE = p.Scale.Div(op0.Scale)
		} else {
			panic(fmt.Errorf("cannot MulThenAdd: op0.Scale > opOut.Scale is not supported"))
		}

		bComplex[0].Mul(bComplex[0], &scaleRLWE.Value)
		bComplex[1].Mul(bComplex[1], &scaleRLWE.Value)

		acc := bignum.NewComplex()

		for i := 0; i < p.Degree+1; i++ {
			m0 := op0.Value[i]
			m2 := p.Value[i]
			for i := range m2 {
				mul(m0[i], bComplex, acc)
				m2[i].Add(m2[i], acc)
			}
		}
		
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

	p.Scale = p.Scale.Div(rlwe.NewScale(Q))
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
