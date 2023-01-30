package estimator

import(
	"math"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
)


func (e *Estimator) HomomorphicEncoding(ct Element, params ckks.Parameters, ecd ckks.Encoder, encodingMatrixLiteral advanced.EncodingMatrixLiteral) (Element){

	LTs := GetEncodingMatrixSTD(params, ecd, encodingMatrixLiteral)

	// Real/Imag packing
	if encodingMatrixLiteral.LinearTransformType == advanced.SlotsToCoeffs && params.LogSlots() == params.LogN()-1{
		ct = e.Add(ct, ct)
	}

	// DFT
	for i := range LTs[:]{
		ct = e.LinearTransform(ct, LTs[i])
		ct = e.Rescale(ct)
	}

	// Real/Imag extraction
	if encodingMatrixLiteral.LinearTransformType == advanced.CoeffsToSlots && encodingMatrixLiteral.RepackImag2Real{
		ct = e.Add(ct, e.KeySwitch(ct)) // Extract real/imag part
		if params.LogSlots() < params.LogN()-1{
			ct = e.Add(ct, e.KeySwitch(ct)) // Rotates by slots/2 and adds
		}
	}
	
	return ct
}


func GetEncodingMatrixSTD(params ckks.Parameters, ecd ckks.Encoder, encodingMatrixLiteral advanced.EncodingMatrixLiteral) (LTs []LinearTransform){

	encodingMatrices := encodingMatrixLiteral.ComputeDFTMatrices()

	LTs = make([]LinearTransform, len(encodingMatrices))

	logdSlots := params.LogSlots()
	if logdSlots < encodingMatrixLiteral.LogN-1 && encodingMatrixLiteral.RepackImag2Real {
		logdSlots++
	}

	ringQ := params.RingQ().AtLevel(0)
	pt := ckks.NewPlaintext(params, ringQ.Level())

	for i, matrix := range encodingMatrices{

		m := make(map[int][2]float64)

		for j, diag := range matrix{

			if encodingMatrixLiteral.LinearTransformType == advanced.CoeffsToSlots && encodingMatrixLiteral.RepackImag2Real{
				if  i == len(encodingMatrices)-1{
					m[j] = GetSTDEncodedVector(ecd, params.N(), logdSlots, diag, pt)
				}else{
					m[j] = GetSTDEncodedVector(ecd, params.N(), params.LogSlots(), diag, pt)
				}
			}else{
				m[j] = GetSTDEncodedVector(ecd, params.N(), logdSlots, diag, pt)
			}
		}

		Level := encodingMatrixLiteral.LevelStart-i

		LTs[i] = NewLinearTransform(m, params.Q()[Level], Level, params.LogSlots(), encodingMatrixLiteral.LogBSGSRatio)
	}

	return
}

func GetSTDEncodedVector(ecd ckks.Encoder, N, LogSlots int, a []complex128, pt *rlwe.Plaintext) ([2]float64){

	vec := make([]complex128, 1<<LogSlots)

	copy(vec, a)

	ecd.IFFT(vec, LogSlots)

	bb := make([]float64, 2<<LogSlots)

	for i, j := 0, 1<<LogSlots; i < 1<<LogSlots; i, j = i+1, j+1{
		bb[i], bb[j] = real(vec[i]), imag(vec[i])
	}

	b := make([]float64, N)

	slots := 1<<LogSlots
	gap := N/(2*slots)
	for i, j := 0, N>>1; i < slots; i, j = i+gap, j+gap{
		b[i], b[j] = real(vec[i]), imag(vec[i])
	}

	ecd.EncodeCoeffs(b, pt)
	c := ecd.DecodeCoeffs(pt)

	for i := range c{
		c[i] -= b[i]
	}

	// This seems to give much better results than taking the STD of the 
	// encoded plaintext when the Hamming weight of the secret is small. 
	// A possible reason is that the encoded plaintext is a polynomial 
	// in Y = X^{N/n} and terms which are not a multiples of N/n do 
	// not contribute to the plaintext encoding/decoding error.
	correction := math.Pow(math.Sqrt(2), math.Log2(float64(gap))) // = sqrt(2) ** log2(gap)
	return [2]float64{STD(bb)/correction, STD(c)}
}