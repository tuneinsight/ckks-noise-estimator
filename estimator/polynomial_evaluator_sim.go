package estimator

import (
	"math/big"
	"math/bits"

	he "github.com/tuneinsight/lattigo/v6/circuits/common/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// simEvaluator is a struct used to pre-computed the scaling
// factors of the polynomial coefficients used by the inlined
// polynomial evaluation by running the polynomial evaluation
// with dummy operands.
// This struct implements the interface he.SimEvaluator.
type SimEvaluator struct {
	params                     ckks.Parameters
	levelsConsumedPerRescaling int
}

// PolynomialDepth returns the depth of the polynomial.
func (d SimEvaluator) PolynomialDepth(degree int) int {
	return d.levelsConsumedPerRescaling * (bits.Len64(uint64(degree)) - 1)
}

// Rescale rescales the target he.SimOperand n times and returns it.
func (d SimEvaluator) Rescale(op0 *he.SimOperand) {
	for i := 0; i < d.levelsConsumedPerRescaling; i++ {
		op0.Scale = op0.Scale.Div(rlwe.NewScale(d.params.Q()[op0.Level]))
		op0.Level--
	}
}

// MulNew multiplies two he.SimOperand, stores the result the target he.SimOperand and returns the result.
func (d SimEvaluator) MulNew(op0, op1 *he.SimOperand) (opOut *he.SimOperand) {
	opOut = new(he.SimOperand)
	opOut.Level = min(op0.Level, op1.Level)
	opOut.Scale = op0.Scale.Mul(op1.Scale)
	return
}

// UpdateLevelAndScaleBabyStep returns the updated level and scale for a baby-step.
func (d SimEvaluator) UpdateLevelAndScaleBabyStep(lead bool, tLevelOld int, tScaleOld rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	tLevelNew = tLevelOld
	tScaleNew = tScaleOld

	if lead {
		for i := 0; i < d.levelsConsumedPerRescaling; i++ {
			tScaleNew = tScaleNew.Mul(rlwe.NewScale(d.params.Q()[tLevelNew-i]))
		}
	}

	return
}

// UpdateLevelAndScaleGiantStep returns the updated level and scale for a giant-step.
func (d SimEvaluator) UpdateLevelAndScaleGiantStep(lead bool, tLevelOld int, tScaleOld, xPowScale rlwe.Scale) (tLevelNew int, tScaleNew rlwe.Scale) {

	Q := d.params.Q()

	var qi *big.Int
	if lead {
		qi = bignum.NewInt(Q[tLevelOld])
		for i := 1; i < d.levelsConsumedPerRescaling; i++ {
			qi.Mul(qi, bignum.NewInt(Q[tLevelOld-i]))
		}
	} else {
		qi = bignum.NewInt(Q[tLevelOld+d.levelsConsumedPerRescaling])
		for i := 1; i < d.levelsConsumedPerRescaling; i++ {
			qi.Mul(qi, bignum.NewInt(Q[tLevelOld+d.levelsConsumedPerRescaling-i]))
		}
	}

	tLevelNew = tLevelOld + d.levelsConsumedPerRescaling
	tScaleNew = tScaleOld.Mul(rlwe.NewScale(qi))
	tScaleNew = tScaleNew.Div(xPowScale)

	return
}
