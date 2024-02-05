package estimator

import (
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

type DFTMatrix struct {
	hefloat.DFTMatrixLiteral
	Value []hefloat.Diagonals[*bignum.Complex]
}

func (d *DFTMatrix) GenMatrices(LogN int, prec uint) {
	d.Value = d.DFTMatrixLiteral.GenMatrices(LogN, prec)
}

func (el *Element) DFT(mat DFTMatrix) *Element {
	for i := range mat.Value {
		scale := el.Parameters.Parameters.GetOptimalScalingFactor(el.Scale, el.DefaultScale(), el.Level)
		el.EvaluateLinearTransformation(LinearTransformation{
			Scale:                    scale,
			LogSlots:                 mat.LogSlots,
			LogBabyStepGianStepRatio: mat.LogBSGSRatio,
			Value:                    mat.Value[i],
		})
		el.Rescale()
	}
	return el
}

func (el *Element) SlotsToCoeffs(elReal, elImag *Element, mat DFTMatrix) (*Element) {

	if elImag != nil {
		el.MulThenAdd(elImag, 1i)
	} else {
		el = elReal
	}

	el = el.DFT(mat)

	return el
}

func (el *Element) CoeffsToSlots(mat DFTMatrix) (elReal, elImag *Element) {

	if mat.Format == hefloat.RepackImagAsReal || mat.Format == hefloat.SplitRealAndImag {

		zV := el.CopyNew()
		zV.DFT(mat)

		elReal = zV.CopyNew().Conjugate()
		elImag = NewElement[*bignum.Complex](*el.Parameters, nil, 1, el.Scale)

		elImag.Sub(zV, elReal)
		elImag.Mul(elImag, -1i)
		elReal.Add(elReal, zV)


		// If repacking, then ct0 and ct1 right n/2 slots are zero.
		if mat.Format == hefloat.RepackImagAsReal && mat.LogSlots < el.LogMaxSlots() {
			elImag.Rotate(1 << mat.LogSlots)
			elReal.Add(elReal, elImag)
			return elReal, nil
		}

		return elReal, elImag
	}

	return el.CopyNew().DFT(mat), nil
}
