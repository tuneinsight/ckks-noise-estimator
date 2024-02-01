package estimator

import(
	"fmt"
	"math/big"
	"math/cmplx"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// EvaluteMod1 applies a homomorphic mod Q on a vector scaled by Delta, scaled down to mod 1 :
//
//  1. Delta * (Q/Delta * I(X) + m(X)) (Delta = scaling factor, I(X) integer poly, m(X) message)
//  2. Delta * (I(X) + Delta/Q * m(X)) (divide by Q/Delta)
//  3. Delta * (Delta/Q * m(X)) (x mod 1)
//  4. Delta * (m(X)) (multiply back by Q/Delta)
//
// Since Q is not a power of two, but Delta is, then does an approximate division by the closest
// power of two to Q instead. Hence, it assumes that the input plaintext is already scaled by
// the correcting factor Q/2^{round(log(Q))}.
//
// !! Assumes that the input is normalized by 1/K for K the range of the approximation.
//
// Scaling back error correction by 2^{round(log(Q))}/Q afterward is included in the polynomial
func (el *Element) EvaluateMod1(evm hefloat.Mod1Parameters) (*Element, error){
	return el.EvaluateMod1AndScale(evm, 1)
}

func (el *Element) EvaluateMod1AndScale(evm hefloat.Mod1Parameters, scaling complex128) (*Element, error){

	var err error

	if el.Level < evm.LevelQ {
		return nil, fmt.Errorf("cannot Evaluate: ct.Level() < Mod1Parameters.LevelQ")
	}

	el.Level = evm.LevelQ

	// Stores default scales
	prevScaleCt := el.Scale

	// Normalize the modular reduction to mod by 1 (division by Q)
	el.Scale = evm.ScalingFactor()

	// Compute the scales that the ciphertext should have before the double angle
	// formula such that after it it has the scale it had before the polynomial
	// evaluation

	params := el.Parameters.Parameters

	Qi := params.Q()



	targetScale := el.Scale
	for i := 0; i < evm.DoubleAngle; i++ {
		targetScale = targetScale.Mul(rlwe.NewScale(Qi[el.Level-evm.Mod1Poly.Depth()-evm.DoubleAngle+i+1]))
		targetScale.Value.Sqrt(&targetScale.Value)
	}

	// Division by 1/2^r and change of variable for the Chebyshev evaluation
	if evm.Mod1Type == hefloat.CosDiscrete || evm.Mod1Type == hefloat.CosContinuous {
		offset := new(big.Float).Sub(&evm.Mod1Poly.B, &evm.Mod1Poly.A)
		offset.Mul(offset, new(big.Float).SetFloat64(evm.IntervalShrinkFactor()))
		offset.Quo(new(big.Float).SetFloat64(-0.5), offset)
		el.Add(el, offset)
	}

	// Double angle
	sqrt2pi := complex(evm.Sqrt2Pi, 0)

	var mod1Poly bignum.Polynomial
	if evm.Mod1InvPoly == nil {

		scaling := cmplx.Pow(scaling, complex(1/evm.IntervalShrinkFactor(), 0))

		mul := bignum.NewComplexMultiplier().Mul

		mod1Poly = evm.Mod1Poly.Clone()

		scalingPowBig := bignum.NewComplex().SetComplex128(scaling)

		for i := range mod1Poly.Coeffs {
			if mod1Poly.Coeffs[i] != nil {
				mul(mod1Poly.Coeffs[i], scalingPowBig, mod1Poly.Coeffs[i])
			}
		}

		sqrt2pi *= scaling

	} else {
		mod1Poly = evm.Mod1Poly
	}

	// Chebyshev evaluation
	if el, err = el.EvaluatePolynomial(mod1Poly, rlwe.NewScale(targetScale)); err != nil {
		return nil, fmt.Errorf("cannot Evaluate: %w", err)
	}

	for i := 0; i < evm.DoubleAngle; i++ {
		sqrt2pi *= sqrt2pi
		el.Mul(el, el)
		el.Relinearize()
		el.Mul(el, 2)
		el.Add(el, -sqrt2pi)
		el.Rescale()
	}

	// ArcSine
	if evm.Mod1InvPoly != nil {

		mul := bignum.NewComplexMultiplier().Mul

		mod1InvPoly := evm.Mod1InvPoly.Clone()

		scalingBig := bignum.NewComplex().SetComplex128(scaling)

		for i := range mod1InvPoly.Coeffs {
			if mod1InvPoly.Coeffs[i] != nil {
				mul(mod1InvPoly.Coeffs[i], scalingBig, mod1InvPoly.Coeffs[i])
			}
		}

		if el, err = el.EvaluatePolynomial(mod1InvPoly, el.Scale); err != nil {
			return nil, fmt.Errorf("cannot Evaluate: %w", err)
		}
	}

	// Multiplies back by q
	el.Scale = prevScaleCt
	return el, nil
}
