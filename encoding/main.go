package main

import (
	"encoding/csv"
	"os"
	"time"
	"fmt"
	"runtime"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
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
	UseSK    = false                  // True to use secret-key encryption
	NbRuns   = 8                      // Number of recorded events
	LogBSGSRatio = 2                  // Ratio N2/N1
	LtType   = advanced.SlotsToCoeffs //advanced.CoeffsToSlots //
)

func main() {

	var s string
	switch LtType {
	case advanced.CoeffsToSlots:
		s = "c2s"
	case advanced.SlotsToCoeffs:
		s = "s2c"
	}

	var keyType string
	switch UseSK{
	case true:
		keyType = "sk"
	case false:
		keyType = "pk"
	}

	
	f, err := os.Create(fmt.Sprintf("data/experiment_%s_%s_logslots_%d_runs_%d_id_%d.csv", s, keyType, LogSlots, NbRuns, time.Now().Unix()))
	if err != nil {
		panic(err)
	}
	defer f.Close()

	w := csv.NewWriter(f)

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
	

	// H: Hamming weight
	// Depth: matrix decomposition
	for H := 32; H <= 32768; H <<= 1 {
		for Depth := 15; Depth > 1; Depth-- {

			fmt.Printf("%s - %s - H:%d - Depth: %d\n", s, keyType, H, Depth)

			c := NewContext(H, Depth, LogSlots, LtType)

			for i := 0; i < NbRuns; i++ {
				// Generate a new set of keys
				// Evaluates the linear transform
				// Records the precision stats
				c.ComputeStats()

				runtime.GC()
			}

			c.Finalize()

			stats := c.ToCSV()

			stdPredicted := EstimateHomomorphicEncodingNoise(c.params, c.ecd1N, c.encodingMatrixLiteral)

			data := []string{
				fmt.Sprintf("%d", H),
				fmt.Sprintf("%d", Depth),
				stats[0],
				stats[1],
				stats[2],
				fmt.Sprintf("%.5f", stdPredicted),
			}

			fmt.Println(data)

			
			if err := w.Write(data); err != nil {
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
		RepackImag2Real:     true,
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

	if UseSK{
		c.enc = ckks.NewEncryptor(c.params, sk)
	}else{
		c.enc = ckks.NewEncryptor(c.params, c.kgen.GenPublicKey(sk))
	}
	
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
		valuesReal, valuesImag := make([]float64, params.Slots()), make([]float64, params.Slots())
		for i := range valuesReal {
			valuesReal[i], valuesImag[i] = real(values[i]), imag(values[i])
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
		ct0, ct1 := eval.CoeffsToSlotsNew(ciphertext, c.encodingMatrix)

		var diff []float64

		// Checks against the original coefficients

		// Sparse packing
		if params.LogSlots() < params.LogN()-1 {

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

			diff = ecd.DecodeCoeffs(dec.DecryptNew(ct0)) // Decodes in the ring

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

			// Full packing = 2 ciphertext, one for the real and one for the imaginary part
			diff = append(ecd.DecodeCoeffs(dec.DecryptNew(ct0)), ecd.DecodeCoeffs(dec.DecryptNew(ct1))...)
		}

		// Compares
		c.stats.Update(diff)

	case advanced.SlotsToCoeffs:

		slots := params.Slots()

		// Generates the n first slots of the test vector (real and imag part to encode)
		valuesReal, valuesImag := make([]complex128, params.Slots()), make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i], valuesImag[i] = complex(utils.RandFloat64(-1, 1), 0), complex(utils.RandFloat64(-1, 1), 0)
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
		res := eval.SlotsToCoeffsNew(ct0, ct1, c.encodingMatrix)

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
			valuesFloat[idx], valuesFloat[jdx] = real(valuesReal[i]), imag(valuesReal[i])
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


func EstimateHomomorphicEncodingNoise(params ckks.Parameters, ecd ckks.Encoder, encodingMatrixLiteral advanced.EncodingMatrixLiteral) (std float64){

	values := make([]float64, params.N())
	for i := range values{
		values[i] = utils.RandFloat64(-1, 1)
	}

	est := estimator.NewEstimator(params.N(), params.HammingWeight(), params.Q(), params.P())

	LTs := GetEncodingMatrixSTD(params, ecd, encodingMatrixLiteral)

	plaintexterrorstd := 4.0/12.0 //flooring error

	pt := estimator.NewPlaintext(estimator.STD(values) * params.DefaultScale().Float64(), plaintexterrorstd, params.MaxLevel())

	var ct estimator.Element
	if UseSK{
		ct = estimator.NewCiphertextSK(pt)
	}else{
		ct = estimator.NewCiphertextPK(pt)
	}
	
	// Real/Imag packing
	if LtType == advanced.SlotsToCoeffs && params.LogSlots() == params.LogN()-1{
		ct = est.Add(ct, ct)
	}
	
	// DFT
	for i := range LTs{
		ct = est.LinearTransform(ct, LTs[i])
		ct = est.Rescale(ct)
	}

	// Real/Imag extraction
	if LtType == advanced.CoeffsToSlots{
		ct = est.Add(ct, est.KeySwitch(ct)) // Extract real/imag part
		if params.LogSlots() < params.LogN()-1{
			ct = est.Add(ct, est.KeySwitch(ct)) // Rotates by slots/2 and adds
		}
	}
	
	return est.Std(ct)
}

func GetEncodingMatrixSTD(params ckks.Parameters, ecd ckks.Encoder, encodingMatrixLiteral advanced.EncodingMatrixLiteral) (LTs []estimator.LinearTransform){

	encodingMatrices := encodingMatrixLiteral.ComputeDFTMatrices()

	LTs = make([]estimator.LinearTransform, len(encodingMatrices))

	buff := make([]float64, params.N())

	for i, matrix := range encodingMatrices{

		m := make(map[int][2]float64)

		for j, diag := range matrix{
			m[j] = GetSTDEncodedVector(ecd, params.N(), params.LogSlots(), diag, buff)

			if encodingMatrixLiteral.LinearTransformType == advanced.SlotsToCoeffs{
				m[j] = [2]float64{m[j][0], m[j][1]/181.01933598375618}
			}
		}

		Level := encodingMatrixLiteral.LevelStart-i

		LTs[i] = estimator.NewLinearTransform(m, params.Q()[Level], Level, params.LogSlots(), encodingMatrixLiteral.LogBSGSRatio)
	}

	return
}


func GetSTDEncodedVector(ecd ckks.Encoder, N, LogSlots int, a []complex128, b []float64) ([2]float64){

	ecd.IFFT(a, LogSlots)

	slots := 1<<LogSlots
	gap := N/(2*slots)
	for i, j := 0, N>>1; i < slots; i, j = i+gap, j+gap{
		b[i] = real(a[i])
		b[j] = imag(a[i])
	}

	params := ecd.Parameters()
	ringQ := params.RingQ().AtLevel(0)
	pt := ckks.NewPlaintext(params, ringQ.Level())
	ecd.EncodeCoeffs(b, pt)
	c := ecd.DecodeCoeffs(pt)

	for i := range c{
		c[i] -= b[i]
	}

	return [2]float64{estimator.STD(b), estimator.STD(c)}
}