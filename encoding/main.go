package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var (
	LogN     = 16                     // Log2 ring degree
	LogSlots = 15                     // Log2 #slots
	LogScale = 45                     // Log2 scaling factor
	NbRuns   = 8                      // Number of recorded events
	LtType   = advanced.CoeffsToSlots //advanced.SlotsToCoeffs //
)

func main() {

	var s string
	switch LtType {
	case advanced.CoeffsToSlots:
		s = "c2s"
	case advanced.SlotsToCoeffs:
		s = "s2c"
	}

	f, err := os.Create(fmt.Sprintf("data/experiment_%s_%d_%d_%d.csv", s, LogSlots, NbRuns, time.Now().Unix()))
	if err != nil {
		panic(err)
	}
	defer f.Close()

	w := csv.NewWriter(f)

	// CSV Header
	if err := w.Write(stats.Header); err != nil {
		panic(err)
	}

	w.Flush()

	// H: Hamming weight
	// Depth: matrix decomposition
	for H := 32; H <= 32768; H <<= 1 {
		for Depth := 4; Depth > 2; Depth-- {

			fmt.Printf("H:%d - Depth: %d\n", H, Depth)

			c := NewContext(H, Depth, LogSlots, LtType)

			for i := 0; i < NbRuns; i++ {
				// Generate a new set of keys
				// Evaluates the linear transform
				// Records the precision stats
				c.ComputeStats()
			}

			c.Finalize()

			if err := w.Write(c.ToCSV()); err != nil {
				panic(err)
			}

			w.Flush()
		}
	}
}

type Context struct {
	params         ckks.Parameters
	ecd1N          ckks.Encoder            // Encoder in degree N
	ecd2N          ckks.Encoder            // Encoder in degree 2N
	enc            rlwe.Encryptor          // Encryptor
	dec            rlwe.Decryptor          // Decryptor
	kgen           rlwe.KeyGenerator       // KeyGenerator
	eval           advanced.Evaluator      // Evaluator
	rotations      []int                   // Rotations needed for the linear transform
	encodingMatrix advanced.EncodingMatrix // Encoded linear transform
	stats          *stats.PrecisionStats   // Precision stats
}

func NewContext(H, depth, logSlots int, ltType advanced.LinearTransformType) (c *Context) {

	var err error

	LogQ := make([]int, depth+1)

	LogQ[0] = 60
	for i := 1; i < depth+1; i++ {
		LogQ[i] = LogScale
	}

	LogP := []int{61, 61}

	var params1N ckks.Parameters
	if params1N, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     LogQ,
		LogP:     LogP,
		H:        H,
		LogSlots: logSlots,
		LogScale: LogScale,
	}); err != nil {
		panic(err)
	}

	var params2N ckks.Parameters
	if params2N, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN + 1,
		LogQ:     LogQ,
		LogP:     LogP,
		H:        H,
		LogSlots: logSlots,
		LogScale: LogScale,
	}); err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params1N)
	ecd1N := ckks.NewEncoder(params1N)
	ecd2N := ckks.NewEncoder(params2N)

	Levels := make([]int, params1N.MaxLevel())

	for i := range Levels {
		Levels[i] = 1
	}

	encodingMatrixLiteral := advanced.EncodingMatrixLiteral{
		LinearTransformType: ltType,
		LogN:                params1N.LogN(),
		LogSlots:            params1N.LogSlots(),
		LevelStart:          params1N.MaxLevel(),
		Levels:              Levels,
		RepackImag2Real:     true,
		LogBSGSRatio:        2,
	}

	// Gets the rotations indexes for CoeffsToSlots
	rotations := encodingMatrixLiteral.Rotations()

	// Generates the encoding matrices
	encodingMatrix := advanced.NewHomomorphicEncodingMatrixFromLiteral(encodingMatrixLiteral, ecd1N)

	eval := advanced.NewEvaluator(params1N, rlwe.EvaluationKey{})

	return &Context{
		params:         params1N,
		ecd1N:          ecd1N,
		ecd2N:          ecd2N,
		kgen:           kgen,
		eval:           eval,
		rotations:      rotations,
		encodingMatrix: encodingMatrix,
		stats:          stats.NewPrecisionStats(H),
	}
}

func GetEncodingMatricesSTD(params ckks.Parameters, ecd ckks.Encoder, encodingMatrix advanced.EncodingMatrix){

	ringQ := params.RingQ().AtLevel(1)
	tmp := ringQ.NewPoly()

	variances := [][]float64{}

	for _, matrix := range encodingMatrix.Matrices{

		s := []float64{}

		for _, vec := range matrix.Vec{
			ringQ.INTT(vec.Q, tmp)
			ringQ.IMForm(tmp, tmp)

			values := ecd.DecodeCoeffs(&rlwe.Plaintext{Value:tmp, MetaData:rlwe.MetaData{Scale:rlwe.NewScale(1)}})

			s = append(s, STD(values))
		}

		variances = append(variances, s)
	}

	fmt.Printf("[\n")
	for i := range variances{
		fmt.Printf("[")
		for j := range variances[i]{
			fmt.Printf("%f, ", variances[i][j])
		}
		fmt.Printf("],\n")
	}
	fmt.Println()
}

func STD(slice []float64)(std float64){

	n := float64(len(slice))

	var mean float64
	for _, c := range slice{
		mean += c
	}

	mean /= n

	var tmp float64
	for _, c := range slice{
		tmp = (c - mean)
		std += tmp * tmp
	}

	std /= (n-1)

	return math.Sqrt(std)
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {
	sk := c.kgen.GenSecretKey()
	c.enc = ckks.NewEncryptor(c.params, sk)
	c.dec = ckks.NewDecryptor(c.params, sk)
	rtks := c.kgen.GenRotationKeysForRotations(c.rotations, true, sk)
	c.eval = c.eval.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rtks})
}

// ComputeStats generates a new set of keys, evaluates the linear transform and records the precision stats.
func (c *Context) ComputeStats() {

	c.GenKeys() // Generates a new set of keys for each event

	params := c.params
	ecd := c.ecd1N
	enc := c.enc
	dec := c.dec
	eval := c.eval

	switch c.encodingMatrix.LinearTransformType {
	case advanced.CoeffsToSlots:

		// Generates the vector of random complex values
		values := make([]complex128, params.Slots())
		for i := range values {
			values[i] = complex(utils.RandFloat64(-1, 1), utils.RandFloat64(-1, 1))
		}

		// Splits between real and imaginary
		valuesReal := make([]float64, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = real(values[i])
		}

		valuesImag := make([]float64, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = imag(values[i])
		}

		// Applies bit-reverse on the original complex vector
		ckks.SliceBitReverseInPlaceComplex128(values, params.Slots())

		// Maps to a float vector
		// Add gaps if sparse packing
		valuesFloat := make([]float64, params.N())
		gap := params.N() / (2 * params.Slots())
		for i, jdx, idx := 0, params.N()>>1, 0; i < params.Slots(); i, jdx, idx = i+1, jdx+gap, idx+gap {
			valuesFloat[idx] = real(values[i])
			valuesFloat[jdx] = imag(values[i])
		}

		// Encodes and encrypts in the ring
		plaintext := ckks.NewPlaintext(params, params.MaxLevel())
		ecd.EncodeCoeffs(valuesFloat, plaintext) // Encodes in the ring
		ciphertext := enc.EncryptNew(plaintext)

		// Applies the homomorphic DFT
		ct0, ct1 := CoeffsToSlotsNew(eval, ciphertext, c.encodingMatrix)

		var diff []float64

		// Checks against the original coefficients

		// Sparse packing
		if params.LogSlots() < params.LogN()-1 {

			slots := params.Slots()
			twoslots := slots<<1

			// Plaintext circuit
			vec := make([]complex128, twoslots)

			// Embed real vector into the complex vector (trivial)
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				vec[i] = complex(valuesReal[i], 0)
				vec[j] = complex(valuesImag[i], 0)
			}

			// IFFT
			// Twice the number of slots because we encode
			// a real vector with the imaginary part repacked
			// in the real part.
			ecd.IFFT(vec, params.LogSlots()+1)

			// Extract complex vector into real vector
			vecFloat := make([]float64, params.N())
			for i, idx, jdx := 0, 0, params.N()>>1; i < twoslots; i, idx, jdx = i+1, idx+gap/2, jdx+gap/2 {
				vecFloat[idx] = real(vec[i])
				vecFloat[jdx] = imag(vec[i])
			}

			eval.Sub(ct0, ecd.EncodeCoeffsNew(vecFloat, ct0.Level(), ct0.Scale), ct0)

			ct0.Scale = rlwe.NewScale(1)

			diff = ecd.DecodeCoeffs(dec.DecryptNew(ct0)) // Decodes in the ring

		} else {

			// Embed the reverence vectors into the complex vector (trivial)
			vec0 := make([]complex128, params.Slots())
			vec1 := make([]complex128, params.Slots())
			for i := 0; i < params.Slots(); i++ {
				vec0[i] = complex(valuesReal[i], 0)
				vec1[i] = complex(valuesImag[i], 0)
			}

			// IFFT
			ecd.IFFT(vec0, params.LogSlots())
			ecd.IFFT(vec1, params.LogSlots())

			// Extract complex vectors into real vectors
			vecReal := make([]float64, params.N())
			vecImag := make([]float64, params.N())
			for i, j := 0, params.Slots(); i < params.Slots(); i, j = i+1, j+1 {
				vecReal[i], vecReal[j] = real(vec0[i]), imag(vec0[i])
				vecImag[i], vecImag[j] = real(vec1[i]), imag(vec1[i])
			}

			eval.Sub(ct0, ecd.EncodeCoeffsNew(vecReal, ct0.Level(), ct0.Scale), ct0)
			eval.Sub(ct1, ecd.EncodeCoeffsNew(vecImag, ct1.Level(), ct1.Scale), ct1)

			ct0.Scale = rlwe.NewScale(1)
			ct1.Scale = rlwe.NewScale(1)

			// Full packing = 2 ciphertext, one for the real and one for the imaginary part
			diff = append(ecd.DecodeCoeffs(dec.DecryptNew(ct0)), ecd.DecodeCoeffs(dec.DecryptNew(ct1))...)
		}

		// Compares
		c.stats.Update(diff)

	case advanced.SlotsToCoeffs:

		slots := params.Slots()

		// Generates the n first slots of the test vector (real part to encode)
		valuesReal := make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = complex(utils.RandFloat64(-1, 1), 0)
		}

		// Generates the n first slots of the test vector (imaginary part to encode)
		valuesImag := make([]complex128, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = complex(utils.RandFloat64(-1, 1), 0)
		}

		// If sparse, there there is the space to store both vectors in one
		if params.LogSlots() < params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Encodes with IFFT and encrypts the test vectors
		logSlots := params.LogSlots()
		if params.LogSlots() < params.LogN()-1 {
			logSlots++
		}

		// Apply DFT on the reference vector and encode in the ring = regular encoding
		plaintext := ckks.NewPlaintext(params, params.MaxLevel())
		ecd.Encode(valuesReal, plaintext, logSlots)
		ct0 := enc.EncryptNew(plaintext)
		var ct1 *rlwe.Ciphertext
		if params.LogSlots() == params.LogN()-1 {
			ecd.Encode(valuesImag, plaintext, logSlots)
			ct1 = enc.EncryptNew(plaintext)
		}

		// Applies the homomorphic DFT
		res := SlotsToCoeffsNew(eval, ct0, ct1, c.encodingMatrix)

		// Decrypt and decode in the ring
		

		// The result is always returned as a single complex vector, so if full-packing (2 initial vectors)
		// then repacks both vectors together
		
		if params.LogSlots() == params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
			ckks.SliceBitReverseInPlaceComplex128(valuesReal, slots)
		}else{
			ckks.SliceBitReverseInPlaceComplex128(valuesReal, slots)
			ckks.SliceBitReverseInPlaceComplex128(valuesImag, slots)
		}

		valuesFloat := make([]float64, params.N())
		gap := params.N() / (2*slots)
		for i, idx, jdx := 0, 0, params.N()>>1; i < slots; i, idx, jdx = i+1, idx+gap, jdx+gap {
			valuesFloat[idx] = real(valuesReal[i])
			valuesFloat[jdx] = imag(valuesReal[i])
		}

		eval.Sub(res, ecd.EncodeCoeffsNew(valuesFloat, res.Level(), res.Scale), res)
		res.Scale = rlwe.NewScale(1)

		diff := ecd.DecodeCoeffs(dec.DecryptNew(res))

		c.stats.Update(diff)

	default:
		panic("how did you get there?")

	}
}

// Computes the final precision stats
func (c *Context) Finalize() {
	c.stats.Finalize()
}

// Produces a CSV friendly string.
func (c *Context) ToCSV() []string {
	return c.stats.ToCSV()
}

// CoeffsToSlotsNew applies the homomorphic encoding and returns the result on new ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func CoeffsToSlotsNew(eval advanced.Evaluator, ctIn *rlwe.Ciphertext, ctsMatrices advanced.EncodingMatrix) (ctReal, ctImag *rlwe.Ciphertext) {

	params := eval.Parameters()

	ctReal = ckks.NewCiphertext(params, 1, ctsMatrices.LevelStart)

	if params.LogSlots() == params.LogN()-1 {
		ctImag = ckks.NewCiphertext(params, 1, ctsMatrices.LevelStart)
	}

	CoeffsToSlots(eval, ctIn, ctsMatrices, ctReal, ctImag)
	return
}

// CoeffsToSlots applies the homomorphic encoding and returns the results on the provided ciphertexts.
// Homomorphically encodes a complex vector vReal + i*vImag of size n on a real vector of size 2n.
// If the packing is sparse (n < N/2), then returns ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then returns ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func CoeffsToSlots(eval advanced.Evaluator, ctIn *rlwe.Ciphertext, ctsMatrices advanced.EncodingMatrix, ctReal, ctImag *rlwe.Ciphertext) {

	if ctsMatrices.RepackImag2Real {

		params := eval.Parameters()

		zV := ctIn.CopyNew()

		dft(eval, ctIn, ctsMatrices.Matrices, zV)

		eval.Conjugate(zV, ctReal)

		var tmp *rlwe.Ciphertext
		if ctImag != nil {
			tmp = ctImag
		} else {
			tmp = rlwe.NewCiphertextAtLevelFromPoly(ctReal.Level(), eval.BuffCt().Value[:2])
			tmp.IsNTT = true
		}

		// Imag part
		eval.Sub(zV, ctReal, tmp)
		eval.MultByConst(tmp, complex(0, -1), tmp)

		// Real part
		eval.Add(ctReal, zV, ctReal)

		// If repacking, then ct0 and ct1 right n/2 slots are zero.
		if params.LogSlots() < params.LogN()-1 {
			eval.Rotate(tmp, params.Slots(), tmp)
			eval.Add(ctReal, tmp, ctReal)
		}

		zV = nil

	} else {
		dft(eval, ctIn, ctsMatrices.Matrices, ctReal)
	}
}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on a new ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func SlotsToCoeffsNew(eval advanced.Evaluator, ctReal, ctImag *rlwe.Ciphertext, stcMatrices advanced.EncodingMatrix) (ctOut *rlwe.Ciphertext) {

	params := eval.Parameters()

	if ctReal.Level() < stcMatrices.LevelStart || (ctImag != nil && ctImag.Level() < stcMatrices.LevelStart) {
		panic("ctReal.Level() or ctImag.Level() < EncodingMatrix.LevelStart")
	}

	ctOut = ckks.NewCiphertext(params, 1, stcMatrices.LevelStart)
	SlotsToCoeffs(eval, ctReal, ctImag, stcMatrices, ctOut)
	return

}

// SlotsToCoeffsNew applies the homomorphic decoding and returns the result on the provided ciphertext.
// Homomorphically decodes a real vector of size 2n on a complex vector vReal + i*vImag of size n.
// If the packing is sparse (n < N/2) then ctReal = Ecd(vReal || vImag) and ctImag = nil.
// If the packing is dense (n == N/2), then ctReal = Ecd(vReal) and ctImag = Ecd(vImag).
func SlotsToCoeffs(eval advanced.Evaluator, ctReal, ctImag *rlwe.Ciphertext, stcMatrices advanced.EncodingMatrix, ctOut *rlwe.Ciphertext) {
	// If full packing, the repacking can be done directly using ct0 and ct1.
	if ctImag != nil {
		eval.MultByConst(ctImag, complex(0, 1), ctOut)
		eval.Add(ctOut, ctReal, ctOut)
		dft(eval, ctOut, stcMatrices.Matrices, ctOut)
	} else {
		dft(eval, ctReal, stcMatrices.Matrices, ctOut)
	}
}

func dft(eval advanced.Evaluator, ctIn *rlwe.Ciphertext, plainVectors []ckks.LinearTransform, ctOut *rlwe.Ciphertext) {

	// Sequentially multiplies w with the provided dft matrices.
	scale := ctIn.Scale
	var in, out *rlwe.Ciphertext
	for i, plainVector := range plainVectors {
		in, out = ctOut, ctOut
		if i == 0 {
			in, out = ctIn, ctOut
		}
		eval.LinearTransform(in, plainVector, []*rlwe.Ciphertext{out})

		if err := eval.Rescale(out, scale, out); err != nil {
			panic(err)
		}
	}
}
