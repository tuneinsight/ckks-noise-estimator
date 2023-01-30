package main

import (
	"encoding/csv"
	"os"
	"time"
	"fmt"
	"runtime"
	"math"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var (
	// Scheme Parameters
	LogN     = 16                     // Log2 ring degree
	LogSlots = 14                     // Log2 #slots
	LogScale = 45                     // Log2 scaling factor
	KeyType  = "pk"                   // "pk"

	// Homomorphic Encoding/Decoding Parameters:
	LtType   = advanced.SlotsToCoeffs //advanced.CoeffsToSlots //
	RepackImag2Real = false           // Additional lineart transform [a, b, c, d] -> [(a+ci), (b+di), (a+ci), (b+di)] = [(a+ci), (b+di)]
	LogBSGSRatio = 2                  // Ratio N2/N1

	// Others
	NbRuns   = 1                      // Number of recorded events
	Record   = false                  // Record in CSV
)

func main() {

	var s string
	switch LtType {
	case advanced.CoeffsToSlots:
		s = "c2s"
	case advanced.SlotsToCoeffs:
		s = "s2c"
	}

	var f *os.File
	var w *csv.Writer
	var err error
	
	if Record{
		if f, err = os.Create(fmt.Sprintf("data/experiment_%s_%s_logslots_%d_runs_%d_id_%d.csv", s, KeyType, LogSlots, NbRuns, time.Now().Unix())); err != nil{
			panic(err)
		}

		defer f.Close()

		w = csv.NewWriter(f)

		Header := []string{
			"H",
			"Depth",
			stats.Header[0],
			stats.Header[1],
			stats.Header[2],
			"Estimation",
		}

		// CSV Header
		if err := w.Write(Header); err != nil {
			panic(err)
		}

		w.Flush()
	}
	
	// H: Hamming weight
	// Depth: matrix decomposition
	for H := 32; H <= 32; H <<= 1 {
		for Depth := LogSlots; Depth >= 4; Depth-- {

			fmt.Printf("%s - %s - H:%d - Depth: %d\n", s, KeyType, H, Depth)

			c := NewContext(H, Depth, LogSlots, LtType)

			est := estimator.NewEstimator(c.params.N(), 0, c.params.Q(), c.params.P())

			scale := c.params.DefaultScale().Float64()

			var stdEstimated float64
			for i := 0; i < NbRuns; i++ {
				// Generate a new set of keys
				// Evaluates the linear transform
				// Records the precision stats
				stdMsg, stdErr := c.ComputeStats()

				fmt.Println(stdMsg, stdErr)

				pt := estimator.NewPlaintext(stdMsg * scale, stdErr * scale, c.params.MaxLevel())

				var ct estimator.Element
				switch KeyType{
				case "sk":
					ct = estimator.NewCiphertextSK(pt)
				case "pk":
					ct = estimator.NewCiphertextPK(pt)
				}

				stdEstimated += est.Std(est.HomomorphicEncoding(ct, c.params, c.ecd1N, c.encodingMatrixLiteral))

				runtime.GC()
			}

			stdEstimated /= float64(NbRuns)

			c.Finalize()

			stats := c.ToCSV()

			data := []string{
				fmt.Sprintf("%d", H),
				fmt.Sprintf("%d", Depth),
				stats[0],
				stats[1],
				stats[2],
				fmt.Sprintf("%.5f", stdEstimated),
			}

			fmt.Println(data)

			if Record{
				if err := w.Write(data); err != nil {
					panic(err)
				}

				w.Flush()
			}
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
	encodingMatrixLiteral advanced.EncodingMatrixLiteral
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
		RepackImag2Real:     RepackImag2Real,
		LogBSGSRatio:        LogBSGSRatio,
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
		encodingMatrixLiteral: encodingMatrixLiteral,
		stats:          stats.NewPrecisionStats(),
	}
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {

	sk := c.kgen.GenSecretKey()

	sk = rlwe.NewSecretKey(c.params.Parameters)

	switch KeyType{
	case "sk":
		c.enc = ckks.NewEncryptor(c.params, sk)
	case "pk":
		c.enc = ckks.NewEncryptor(c.params, c.kgen.GenPublicKey(sk))
	}
	
	c.dec = ckks.NewDecryptor(c.params, sk)
	rtks := c.kgen.GenRotationKeysForRotations(c.rotations, true, sk)
	c.eval = c.eval.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rtks})
}

func NewTestVector(params ckks.Parameters, ecd ckks.Encoder) (msgFFT []complex128, msgIFFT []float64, stdMsg, stdErr float64){
	return
}

// ComputeStats generates a new set of keys, evaluates the linear transform and records the precision stats.
func (c *Context) ComputeStats() (stdMsg, stdErr float64){

	c.GenKeys() // Generates a new set of keys for each event

	params := c.params
	ecd := c.ecd1N
	enc := c.enc
	dec := c.dec
	eval := c.eval

	switch c.encodingMatrix.LinearTransformType {
	case advanced.CoeffsToSlots:

		// Circuit:
		//
		// Case 1: Full/Sparse packing IFFT only
		//
		// Input: [(a + ei), (b + fi), (c + gi), (d + hi)]
		// Output: IFFT([(a + ei), (b + fi), (c + gi), (d + hi)])
		//
		// Case 2: Full packing with [(a + bi), 0] -> [a, b] repacking
		//
		// Input: [(a + ei), (b + fi), (c + gi), (d + hi)]
		// Output ct0 = IFFT([a, b, c, d]) ct1 = IFFT([e, f, g, h])
		//
		// Case2: Sparse packing with 
		//
		// Input: [(a + ei), (b + fi), (c + gi), (d + hi)]
		// Output: ct0 = [IFFT([a, b, c, d]), IFFT([e, f, g, h])]

		// Generates the vector of random complex values
		values := make([]complex128, params.Slots())
		for i := range values {
			values[i] = complex(utils.RandFloat64(-1, 1), utils.RandFloat64(-1, 1))
		}

		// Splits between real and imaginary
		valuesReal, valuesImag := make([]float64, params.Slots()), make([]float64, params.Slots())
		for i := range valuesReal {
			valuesReal[i], valuesImag[i] = real(values[i]), imag(values[i])
		}

		stdMsg = estimator.STD(append(valuesReal, valuesImag...))

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

		valuesErr := ecd.DecodeCoeffs(plaintext)
		for i := range valuesErr{
			valuesErr[i] -= valuesFloat[i]
		}
		stdErr = estimator.STD(valuesErr)

		// Applies the homomorphic DFT
		ct0, ct1 := CoeffsToSlotsNew(eval, ciphertext, c.encodingMatrix)

		var diff []float64

		// Plain DFT (sparse & dense packing)
		if !c.encodingMatrix.RepackImag2Real{

			vec := make([]complex128, params.Slots())
			for i := 0; i < params.Slots(); i++ {
				vec[i] = complex(valuesReal[i], valuesImag[i])
			}

			ckks.SliceBitReverseInPlaceComplex128(values, len(values))
			ecd.IFFT(values, params.LogSlots())

			vecReal := make([]float64, params.N())
			gap := params.N() / (2*params.Slots())
			for i, idx, jdx := 0, 0, params.N()>>1; i < params.Slots(); i, idx, jdx = i+1, idx+gap, jdx+gap {
				vecReal[idx], vecReal[jdx] = real(values[i]), imag(values[i])
			}

			eval.Sub(ct0, ecd.EncodeCoeffsNew(vecReal, ct0.Level(), ct0.Scale), ct0)
			ct0.Scale = rlwe.NewScale(1)
			real := ecd.DecodeCoeffs(dec.DecryptNew(ct0))

			// We only consider the plaintext slots
			diff = make([]float64, 2*params.Slots())
			for i, j := 0, 0; i < 2*params.Slots(); i, j = i+1, j+gap{
				diff[i] = real[j]
			}

		// DFT with sparse packing and real | real -> (real, img) repacking
		}else if params.LogSlots() < params.LogN()-1 {

			slots := params.Slots()
			twoslots := slots<<1
	
			// Embed real vector into the complex vector (trivial)
			vec := make([]complex128, twoslots)
			for i, j := 0, slots; i < slots; i, j = i+1, j+1 {
				vec[i], vec[j] = complex(valuesReal[i], 0), complex(valuesImag[i], 0)
			}

			// IFFT
			// Twice the number of slots because we encode
			// a real vector with the imaginary part repacked
			// in the real part.
			ecd.IFFT(vec, params.LogSlots()+1)

			// Extract complex vector into real vector
			vecFloat := make([]float64, params.N())
			for i, idx, jdx := 0, 0, params.N()>>1; i < twoslots; i, idx, jdx = i+1, idx+gap/2, jdx+gap/2 {
				vecFloat[idx], vecFloat[jdx] = real(vec[i]), imag(vec[i])
			}

			eval.Sub(ct0, ecd.EncodeCoeffsNew(vecFloat, ct0.Level(), ct0.Scale), ct0)
			ct0.Scale = rlwe.NewScale(1)

			real := ecd.DecodeCoeffs(dec.DecryptNew(ct0)) // Decodes in the ring

			// We only consider the plaintext slots
			diff = make([]float64, 2*params.Slots())
			for i, j := 0, 0; i < 2*params.Slots(); i, j = i+1, j+gap{
				diff[i] = real[j]
			}

		// DFT with full packing and real | real -> (real, img) repacking
		} else {

			// Embed the reverence vectors into the complex vector (trivial)
			vec0, vec1 := make([]complex128, params.Slots()), make([]complex128, params.Slots())
			for i := 0; i < params.Slots(); i++ {
				vec0[i], vec1[i] = complex(valuesReal[i], 0), complex(valuesImag[i], 0)
			}

			// IFFT
			ecd.IFFT(vec0, params.LogSlots())
			ecd.IFFT(vec1, params.LogSlots())

			// Extract complex vectors into real vectors
			vecReal, vecImag := make([]float64, params.N()), make([]float64, params.N())
			for i, j := 0, params.Slots(); i < params.Slots(); i, j = i+1, j+1 {
				vecReal[i], vecReal[j] = real(vec0[i]), imag(vec0[i])
				vecImag[i], vecImag[j] = real(vec1[i]), imag(vec1[i])
			}

			eval.Sub(ct0, ecd.EncodeCoeffsNew(vecReal, ct0.Level(), ct0.Scale), ct0)
			eval.Sub(ct1, ecd.EncodeCoeffsNew(vecImag, ct1.Level(), ct1.Scale), ct1)

			ct0.Scale = rlwe.NewScale(1)
			ct1.Scale = rlwe.NewScale(1)

			real := ecd.DecodeCoeffs(dec.DecryptNew(ct0))
			imag := ecd.DecodeCoeffs(dec.DecryptNew(ct1))

			// We only consider the plaintext slots
			diff = make([]float64, 4*params.Slots())
			for i, j, k := 0, params.Slots(), 0; i < 2*params.Slots(); i, j, k = i+1,j+1, k+gap{
				diff[i] = real[k]
				diff[j] = imag[k]
			}
		}

		// Compares
		c.stats.Update(diff)

	case advanced.SlotsToCoeffs:

		// Circuit:
		//
		// Case 1: Sparse Plaintext & Repacking Imaginary into Real
		//
		// Input: [a, b, c, d, e, f, g, h] 
		//
		// Step 1: [a, b, c, d, e, f, g, h] -> [(a + ei), (b + fi), (c + gi), (d + hi), (a + ei), (b + fi), (c + gi), (d + hi)] = [(a + ei), (b + fi), (c + gi), (d + hi)]
		// Step 2: FFT([(a + ei), (b + fi), (c + gi), (d + hi)])
		//
		// Output: FFT([(a + ei), (b + fi), (c + gi), (d + hi)])
		//
		// Case 2: Sparse Plaintext Only
		//
		// Input: [(a + ei), (b + fi), (c + gi), (d + hi)]
		// Output: FFT([(a + ei), (b + fi), (c + gi), (d + hi)])
		//
		// Case 3: Dense Plaintext & Repacking Imaginary into Real
		//
		// Input: [a, b, c, d, e, f, g, h] & [i, j, k, l, m, n, o, p]
		//
		// Step1: [a, b, c, d, e, f, g, h] + i* [i, j, k, l, m, n, o, p] = [(a + ii), (b + ij), (c + ik), (d + il), ...]
		// Step2: FFT([(a + ii), (b + ij), (c + ik), (d + il), ...])
		// 
		// Output: FFT([(a + ii), (b + ij), (c + ik), (d + il), ...])
		//
		// Case 4: Dense Plaintext
		//
		// Input: [(a + ii), (b + ij), (c + ik), (d + il), ...]
		// Output: FFT([(a + ii), (b + ij), (c + ik), (d + il), ...])

		var sparse bool

		// Encodes with IFFT and encrypts the test vectors
		logSlots := params.LogSlots()

		// In this case, we are given the vector [(a + ei), (b + fi), (c + gi), (d + hi)] as a real vector [a, b, c, d, e, f, g, h] and
		// it maps it to [(a + ei), (b + fi), (c + gi), (d + hi), (a + ei), (b + fi), (c + gi), (d + hi)], which is gives the encoding
		// of [(a + ei), (b + fi), (c + gi), (d + hi)].
		//
		// So to store [a, b, c, d, e, f, g, h] we need twice the number of slots.
		if params.LogSlots() < params.LogN()-1 && c.encodingMatrix.RepackImag2Real{
			logSlots++
		}

		if params.LogSlots() < params.LogN()-1{
			sparse = true
		}

		gap := params.N() / (2<<logSlots)

		values := make([]complex128, params.Slots())
	

		// Generates a vector vec such that IFFT(vec) has the desired standard deviation
		corrfactor := math.Sqrt(float64(params.N()/2))
		corrfactor = 1

		// We generate a full complex vector if either the plaintext is sparse or if we aren't given the
		// real and imaginary part in separate ciphertexts
		if sparse || !c.encodingMatrix.RepackImag2Real {
			for i := range values{
				values[i] = complex(utils.RandFloat64(-1, 1) * corrfactor , utils.RandFloat64(-1, 1) * corrfactor)
			}
		}else{
			for i := range values{
				values[i] = complex(utils.RandFloat64(-1, 1) * corrfactor , 0)
			}
		}

		plaintext := ckks.NewPlaintext(params, params.MaxLevel())

		valuesCmplx := make([]complex128, 1<<logSlots)
		copy(valuesCmplx, values)
		ecd.IFFT(valuesCmplx, logSlots)

		valuesReal := make([]float64, params.N())
		for i, idx, jdx := 0, 0, params.N()>>1; i < 1<<logSlots; i, idx, jdx = i+1, idx+gap, jdx+gap {
			valuesReal[idx], valuesReal[jdx] = real(valuesCmplx[i]), imag(valuesCmplx[i])
		}

		b := make([]float64, 2<<logSlots)
		for i, j := 0, params.Slots(); i < 1<<logSlots; i, j = i+1, j+1 {
			b[i], b[j] = real(valuesCmplx[i]), imag(valuesCmplx[i])
		}

		stdMsg = estimator.STD(valuesReal)

		fmt.Println(stdMsg)

		ecd.EncodeCoeffs(valuesReal, plaintext)
		ct0 := enc.EncryptNew(plaintext)

		var ct1 *rlwe.Ciphertext
		
		// We generate a second ciphertext storing the imaginary part in the case where we
		// are given the real and imaginary part in separate ciphertexts.
		if params.LogSlots() == params.LogN()-1 && c.encodingMatrix.RepackImag2Real{

			for i := range values{
				values[i] = complex(utils.RandFloat64(-1, 1) * corrfactor , 0)
			}

			buffImag := make([]complex128, 1<<logSlots)
			copy(buffImag, values)
			ecd.IFFT(buffImag, logSlots)

			for i, idx, jdx := 0, 0, params.N()>>1; i < 1<<logSlots; i, idx, jdx = i+1, idx+gap, jdx+gap {
				valuesReal[idx], valuesReal[jdx] = real(buffImag[i]), imag(buffImag[i])
			}

			ecd.EncodeCoeffs(valuesReal, plaintext)
			ct1 = enc.EncryptNew(plaintext)

			for i := range valuesCmplx{
				valuesCmplx[i] += (buffImag[i] * complex(0, 1))
			}
		}

		// Applies the homomorphic DFT
		res := SlotsToCoeffsNew(eval, ct0, ct1, c.encodingMatrix)

		// =========== Plaintext Circuit ===============

		ecd.FFT(valuesCmplx, logSlots)
		ckks.SliceBitReverseInPlaceComplex128(valuesCmplx, len(valuesCmplx))

		valuesFloat := make([]float64, params.N())
		for i, idx, jdx := 0, 0, params.N()>>1; i < 1<<logSlots; i, idx, jdx = i+1, idx+gap, jdx+gap {
			valuesFloat[idx], valuesFloat[jdx] = real(valuesCmplx[i]), imag(valuesCmplx[i])
		}
		// ==============================================

		eval.Sub(res, ecd.EncodeCoeffsNew(valuesFloat, res.Level(), res.Scale), res)
		res.Scale = rlwe.NewScale(1)

		real := ecd.DecodeCoeffs(dec.DecryptNew(res))


		diff := make([]float64, 2<<logSlots)

		gap = params.N() / (2<<logSlots)
		fmt.Println(gap)

		for i, j:= 0, 0; i < 2<<logSlots; i, j = i+1, j+gap{
			diff[i] = real[j]
		}

		d2 := ecd.DecodeCoeffs(dec.DecryptNew(res))

		c.stats.Update(d2)

	default:
		panic("how did you get there?")
	}

	return
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
	ctReal = ckks.NewCiphertext(eval.Parameters(), 1, ctsMatrices.LevelStart)

	if eval.Parameters().LogSlots() == eval.Parameters().LogN()-1 {
		ctImag = ckks.NewCiphertext(eval.Parameters(), 1, ctsMatrices.LevelStart)
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
		if eval.Parameters().LogSlots() < eval.Parameters().LogN()-1 {
			eval.Rotate(tmp, eval.Parameters().Slots(), tmp)
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

	if ctReal.Level() < stcMatrices.LevelStart || (ctImag != nil && ctImag.Level() < stcMatrices.LevelStart) {
		panic("ctReal.Level() or ctImag.Level() < EncodingMatrix.LevelStart")
	}

	ctOut = ckks.NewCiphertext(eval.Parameters(), 1, stcMatrices.LevelStart)
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
	for i, plainVector := range plainVectors[:] {
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