package estimator

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func (e *Estimator) HomomorphicEncoding(ct Element, params ckks.Parameters, ecd *ckks.Encoder, encodingMatrixLiteral advanced.HomomorphicDFTMatrixLiteral) (Element, error) {

	var err error

	LTs := GetEncodingMatrixSTD(params, ecd, encodingMatrixLiteral)

	// Real/Imag packing
	if encodingMatrixLiteral.Type == advanced.Decode && encodingMatrixLiteral.LogSlots == params.LogN()-1 {
		ct = e.Add(ct, ct)
	}

	// DFT
	for i := range LTs {
		ct = e.LinearTransform(ct, LTs[i])
		if ct, err = e.Rescale(ct); err != nil {
			return ct, err
		}
	}

	// Real/Imag extraction
	if encodingMatrixLiteral.Type == advanced.Decode {
		ct = e.Add(ct, e.KeySwitch(ct)) // Extract real/imag part
		if encodingMatrixLiteral.LogSlots < params.LogN()-1 {
			ct = e.Add(ct, e.KeySwitch(ct)) // Rotates by slots/2 and adds
		}
	}

	return ct, nil
}

func (e *Estimator) DFT(ct Element, LTs []LinearTransform) (Element, error) {

	var err error

	// DFT
	for i := range LTs {
		ct = e.LinearTransform(ct, LTs[i])
		if ct, err = e.Rescale(ct); err != nil {
			return ct, err
		}
	}

	return ct, nil
}

func GetEncodingMatrixSTD(params ckks.Parameters, ecd *ckks.Encoder, encodingMatrixLiteral advanced.HomomorphicDFTMatrixLiteral) (LTs []LinearTransform) {

	encodingMatrices := encodingMatrixLiteral.GenMatrices(params.LogN(), params.DefaultPrecision())

	LTs = make([]LinearTransform, len(encodingMatrices))

	logdSlots := encodingMatrixLiteral.LogSlots
	if logdSlots < params.LogN()-1 && encodingMatrixLiteral.RepackImag2Real {
		logdSlots++
	}

	for i, matrix := range encodingMatrices {

		m := make(map[int][2]float64)

		for j, diag := range matrix {

			if encodingMatrixLiteral.Type == advanced.Encode {
				if i == len(encodingMatrices)-1 {
					m[j] = GetSTDEncodedVector(ecd, params.N(), logdSlots, diag)
				} else {
					m[j] = GetSTDEncodedVector(ecd, params.N(), encodingMatrixLiteral.LogSlots, diag)
				}
			} else {
				m[j] = GetSTDEncodedVector(ecd, params.N(), logdSlots, diag)
			}
		}

		Level := encodingMatrixLiteral.LevelStart - i

		LTs[i] = NewLinearTransform(m, params.Q()[Level], Level, encodingMatrixLiteral.LogSlots, encodingMatrixLiteral.LogBSGSRatio)
	}

	return
}

func GetSTDEncodedVector(ecd *ckks.Encoder, N, LogSlots int, a []*bignum.Complex) [2]float64 {

	prec := a[0].Prec()

	vec := make([]*bignum.Complex, 1<<LogSlots)

	for i := range vec {
		vec[i] = a[i].Copy()
	}

	if err := ecd.IFFT(vec, LogSlots); err != nil {
		panic(err)
	}

	b := make([]*big.Float, N)
	for i := range b {
		b[i] = new(big.Float).SetPrec(prec)
	}

	slots := 1 << LogSlots
	gap := N / (2 * slots)
	for i, j, k := 0, N>>1, 0; k < slots; i, j, k = i+gap, j+gap, k+1 {
		b[i].Set(vec[k][0])
		b[j].Set(vec[k][1])
	}

	params := ecd.Parameters()
	ringQ := params.RingQ().AtLevel(0)
	pt := ckks.NewPlaintext(params, ringQ.Level())
	pt.EncodingDomain = rlwe.CoefficientsDomain

	if err := ecd.Encode(b, pt); err != nil {
		panic(err)
	}

	c := make([]*big.Float, len(b))

	if err := ecd.Decode(pt, c); err != nil {
		panic(err)
	}

	for i := range c {
		c[i].Sub(c[i], b[i])
	}

	return [2]float64{STD(b), STD(c)} // a sqrt(2)^{degree of sparsity} seems to resolve the issue, but where does it come from o7
}
