package estimator

import (
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
)

func (e *Estimator) HomomorphicEncoding(ct Element, params ckks.Parameters, ecd ckks.Encoder, encodingMatrixLiteral advanced.EncodingMatrixLiteral) Element {

	LTs := GetEncodingMatrixSTD(params, ecd, encodingMatrixLiteral)

	// Real/Imag packing
	if encodingMatrixLiteral.LinearTransformType == advanced.SlotsToCoeffs && params.LogSlots() == params.LogN()-1 {
		ct = e.Add(ct, ct)
	}

	// DFT
	for i := range LTs {
		ct = e.LinearTransform(ct, LTs[i])
		ct = e.Rescale(ct)
	}

	// Real/Imag extraction
	if encodingMatrixLiteral.LinearTransformType == advanced.CoeffsToSlots {
		ct = e.Add(ct, e.KeySwitch(ct)) // Extract real/imag part
		if params.LogSlots() < params.LogN()-1 {
			ct = e.Add(ct, e.KeySwitch(ct)) // Rotates by slots/2 and adds
		}
	}

	return ct
}

func (e *Estimator) DFT(ct Element, LTs []LinearTransform) Element {

	// DFT
	for i := range LTs {
		ct = e.LinearTransform(ct, LTs[i])
		ct = e.Rescale(ct)
	}

	return ct
}

func GetEncodingMatrixSTD(params ckks.Parameters, ecd ckks.Encoder, encodingMatrixLiteral advanced.EncodingMatrixLiteral) (LTs []LinearTransform) {

	encodingMatrices := encodingMatrixLiteral.ComputeDFTMatrices()

	LTs = make([]LinearTransform, len(encodingMatrices))

	logdSlots := params.LogSlots()
	if logdSlots < encodingMatrixLiteral.LogN-1 && encodingMatrixLiteral.RepackImag2Real {
		logdSlots++
	}

	for i, matrix := range encodingMatrices {

		m := make(map[int][2]float64)

		for j, diag := range matrix {

			if encodingMatrixLiteral.LinearTransformType == advanced.CoeffsToSlots {
				if i == len(encodingMatrices)-1 {
					m[j] = GetSTDEncodedVector(ecd, params.N(), logdSlots, diag)
				} else {
					m[j] = GetSTDEncodedVector(ecd, params.N(), params.LogSlots(), diag)
				}
			} else {
				m[j] = GetSTDEncodedVector(ecd, params.N(), logdSlots, diag)
			}
		}

		Level := encodingMatrixLiteral.LevelStart - i

		LTs[i] = NewLinearTransform(m, params.Q()[Level], Level, params.LogSlots(), encodingMatrixLiteral.LogBSGSRatio)
	}

	return
}

func GetSTDEncodedVector(ecd ckks.Encoder, N, LogSlots int, a []complex128) [2]float64 {

	vec := make([]complex128, 1<<LogSlots)

	copy(vec, a)

	ecd.IFFT(vec, LogSlots)

	b := make([]float64, N)

	slots := 1 << LogSlots
	gap := N / (2 * slots)
	for i, j := 0, N>>1; i < slots; i, j = i+gap, j+gap {
		b[i] = real(vec[i])
		b[j] = imag(vec[i])
	}

	params := ecd.Parameters()
	ringQ := params.RingQ().AtLevel(0)
	pt := ckks.NewPlaintext(params, ringQ.Level())
	ecd.EncodeCoeffs(b, pt)
	c := ecd.DecodeCoeffs(pt)

	for i := range c {
		c[i] -= b[i]
	}

	return [2]float64{STD(b), STD(c)}
}
