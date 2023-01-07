package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"sort"
	"time"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var (
	LogN     = 16                     // Log2 ring degree
	LogSlots = 14                     // Log2 #slots
	LogScale = 45                     // Log2 scaling factor
	NbRuns   = 32                     // Number of recorded events
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

	f, err := os.Create(fmt.Sprintf("data/experiment_%s_%d_%d_%d.csv", s, LogSlots, NbRuns, time.Now().Unix()))
	if err != nil {
		panic(err)
	}
	defer f.Close()

	w := csv.NewWriter(f)

	// CSV Header
	if err := w.Write([]string{
		"LogN",
		"H",
		"Depth",
		"LogScale",
		"LogSlots",
		"MIN",
		"AVG",
		"MED",
		"STD",
		"MIN",
		"AVG",
		"MED",
		"STD",
	}); err != nil {
		panic(err)
	}

	w.Flush()

	// H: Hamming weight
	// Depth: matrix decomposition
	for H := 32; H <= 32768; H <<= 1 {
		for Depth := 4; Depth > 1; Depth-- {

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
	stats          *PrecisionStats         // Precision stats
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
		stats:          NewPrecisionStats(LogN, H, params1N.MaxLevel(), params1N.LogSlots(), int(math.Round(math.Log2(params1N.DefaultScale().Float64())))),
	}
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
		ct0, ct1 := eval.CoeffsToSlotsNew(ciphertext, c.encodingMatrix)

		var want, have interface{}

		// Checks against the original coefficients

		// Sparse packing
		if params.LogSlots() < params.LogN()-1 {

			coeffsReal := ecd.DecodeCoeffs(dec.DecryptNew(ct0)) // Decodes in the ring

			// Plaintext circuit
			vec := make([]complex128, 2*params.Slots())

			// Embed real vector into the complex vector (trivial)
			for i, j := 0, params.Slots(); i < params.Slots(); i, j = i+1, j+1 {
				vec[i] = complex(valuesReal[i], 0)
				vec[j] = complex(valuesImag[i], 0)
			}

			// IFFT
			// Twice the number of slots because we encode
			// a real vector with the imaginary part repacked
			// in the real part.
			ecd.IFFT(vec, params.LogSlots()+1)

			// Extract complex vector into real vector
			vecReal := make([]float64, params.N())
			for i, idx, jdx := 0, 0, params.N()>>1; i < 2*params.Slots(); i, idx, jdx = i+1, idx+gap/2, jdx+gap/2 {
				vecReal[idx] = real(vec[i])
				vecReal[jdx] = imag(vec[i])
			}

			want = coeffsReal
			have = vecReal

		} else {

			// Full packing = 2 ciphertext, one for the real and one for the imaginary part
			coeffsReal := ecd.DecodeCoeffs(dec.DecryptNew(ct0))
			coeffsImag := ecd.DecodeCoeffs(dec.DecryptNew(ct1))

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

			tmp0 := make([]complex128, params.N())
			tmp1 := make([]complex128, params.N())

			for i := 0; i < params.N(); i++ {
				tmp0[i] = complex(coeffsReal[i], coeffsImag[i])
				tmp1[i] = complex(vecReal[i], vecImag[i])
			}

			have = tmp0
			want = tmp1
		}

		// Compares
		c.stats.Update(have, want)

	case advanced.SlotsToCoeffs:

		// Generates the n first slots of the test vector (real part to encode)
		valuesReal := make([]complex128, params.Slots())
		for i := range valuesReal {
			valuesReal[i] = complex(float64(i+1)/float64(params.Slots()), 0)
		}

		// Generates the n first slots of the test vector (imaginary part to encode)
		valuesImag := make([]complex128, params.Slots())
		for i := range valuesImag {
			valuesImag[i] = complex(float64(i+1)/float64(params.Slots()), 0)
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
		coeffsFloat := ecd.DecodeCoeffs(dec.DecryptNew(res))

		// Extracts the coefficients and construct the complex vector
		// This is simply coefficient ordering
		valuesTest := make([]complex128, params.Slots())
		gap := params.N() / (2 * params.Slots())
		for i, idx := 0, 0; i < params.Slots(); i, idx = i+1, idx+gap {
			valuesTest[i] = complex(coeffsFloat[idx], coeffsFloat[idx+(params.N()>>1)])
		}

		// The result is always returned as a single complex vector, so if full-packing (2 initial vectors)
		// then repacks both vectors together
		if params.LogSlots() == params.LogN()-1 {
			for i := range valuesReal {
				valuesReal[i] += complex(0, real(valuesImag[i]))
			}
		}

		// Result is bit-reversed, so applies the bit-reverse permutation on the reference vector
		ckks.SliceBitReverseInPlaceComplex128(valuesReal, params.Slots())

		c.stats.Update(valuesReal, valuesTest)

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

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	LogN, H, Depth, LogScale, LogSlots int
	MaxDelta                           Stats
	MinPrecision                       Stats
	MeanDelta                          Stats
	MeanPrecision                      Stats
	MedianDelta                        Stats
	MedianPrecision                    Stats
	StdPrecision                       Stats

	/*
		If we want to produce CDF graphs, we can with this
		RealDist, ImagDist, L2Dist []struct {
			Prec  float64
			Count int
		}

		cdfResol int
	*/

	diff     []Stats
	precReal []float64
	precImag []float64
	precL2   []float64
}

// Stats is a struct storing the real, imaginary and L2 norm (modulus)
// about the precision of a complex value.
type Stats struct {
	Real, Imag, L2 float64
}

func NewPrecisionStats(LogN, H, Depth, LogSlots, LogScale int) (prec *PrecisionStats) {
	return &PrecisionStats{
		LogN:      LogN,
		H:         H,
		Depth:     Depth,
		LogScale:  LogScale,
		LogSlots:  LogSlots,
		MaxDelta:  Stats{0, 0, 0},
		MeanDelta: Stats{0, 0, 0},

		diff:     []Stats{},
		precReal: []float64{},
		precImag: []float64{},
		precL2:   []float64{},
	}
}

func (p *PrecisionStats) Update(want, have interface{}) {

	var valuesTest []complex128

	switch have := have.(type) {
	case []complex128:
		valuesTest = have
	case []float64:
		valuesTest = make([]complex128, len(have))
		for i := range have {
			valuesTest[i] = complex(have[i], 0)
		}
	}

	var valuesWant []complex128
	switch want := want.(type) {
	case []complex128:
		valuesWant = want
	case []float64:
		valuesWant = make([]complex128, len(want))
		for i := range want {
			valuesWant[i] = complex(want[i], 0)
		}
	}

	var deltaReal, deltaImag, deltaL2 float64

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))
		deltaL2 = math.Sqrt(deltaReal*deltaReal + deltaImag*deltaImag)

		p.precReal = append(p.precReal, math.Log2(1/deltaReal))
		p.precImag = append(p.precImag, math.Log2(1/deltaImag))
		p.precL2 = append(p.precL2, math.Log2(1/deltaL2))
		p.diff = append(p.diff, Stats{Real: deltaReal, Imag: deltaImag, L2: deltaL2})

		p.MeanDelta.Real += deltaReal
		p.MeanDelta.Imag += deltaImag
		p.MeanDelta.L2 += deltaL2

		if deltaReal > p.MaxDelta.Real {
			p.MaxDelta.Real = deltaReal
		}

		if deltaImag > p.MaxDelta.Imag {
			p.MaxDelta.Imag = deltaImag
		}

		if deltaL2 > p.MaxDelta.L2 {
			p.MaxDelta.L2 = deltaL2
		}
	}
}

func (p *PrecisionStats) Finalize() {

	p.MeanDelta.Real /= float64(len(p.diff))
	p.MeanDelta.Imag /= float64(len(p.diff))
	p.MeanDelta.L2 /= float64(len(p.diff))

	p.MedianDelta = calcmedian(p.diff)

	p.StdPrecision = calcstandarddeviation(p.diff, p.MeanDelta)

	p.MinPrecision = deltaToPrecision(p.MaxDelta)
	p.MeanPrecision = deltaToPrecision(p.MeanDelta)
	p.MedianPrecision = deltaToPrecision(p.MedianDelta)

	/*
		p.cdfResol = 500

		p.RealDist = make([]struct {
			Prec  float64
			Count int
		}, p.cdfResol)
		p.ImagDist = make([]struct {
			Prec  float64
			Count int
		}, p.cdfResol)
		p.L2Dist = make([]struct {
			Prec  float64
			Count int
		}, p.cdfResol)


		p.calcCDF(p.precReal, p.RealDist)
		p.calcCDF(p.precImag, p.ImagDist)
		p.calcCDF(p.precL2, p.L2Dist)
	*/
}

func (prec *PrecisionStats) ToCSV() []string {
	return []string{
		fmt.Sprintf("%d", prec.LogN),
		fmt.Sprintf("%d", prec.H),
		fmt.Sprintf("%d", prec.Depth),
		fmt.Sprintf("%d", prec.LogScale),
		fmt.Sprintf("%d", prec.LogSlots),
		fmt.Sprintf("%.4f", prec.MinPrecision.Real),
		fmt.Sprintf("%.4f", prec.MeanPrecision.Real),
		fmt.Sprintf("%.4f", prec.MedianPrecision.Real),
		fmt.Sprintf("%.4f", prec.StdPrecision.Real),
		fmt.Sprintf("%.4f", prec.MinPrecision.Imag),
		fmt.Sprintf("%.4f", prec.MeanPrecision.Imag),
		fmt.Sprintf("%.4f", prec.MedianPrecision.Imag),
		fmt.Sprintf("%.4f", prec.StdPrecision.Imag),
	}
}

func calcstandarddeviation(diff []Stats, mean Stats) (std Stats) {
	std = Stats{}

	mReal := math.Log2(1 / mean.Real)
	mImag := math.Log2(1 / mean.Imag)
	mL2 := math.Log2(1 / mean.L2)

	var real, imag, l2 float64

	var v float64
	for _, d := range diff {

		v = math.Log2(1/maxFloat64(d.Real, 1e-16)) - mReal
		real += v * v

		v = math.Log2(1/maxFloat64(d.Imag, 1e-16)) - mImag
		imag += v * v

		v = math.Log2(1/maxFloat64(d.L2, 1e-16)) - mL2
		l2 += v * v
	}

	n := float64(len(diff))

	std.Real = math.Sqrt(real / n)
	std.Imag = math.Sqrt(imag / n)
	std.L2 = math.Sqrt(l2 / n)

	return
}

func deltaToPrecision(c Stats) Stats {
	return Stats{math.Log2(1 / maxFloat64(c.Real, 1e-16)), math.Log2(1 / maxFloat64(c.Imag, 1e-16)), math.Log2(1 / maxFloat64(c.L2, 1e-16))}
}

func maxFloat64(a, b float64) (c float64) {
	if a > b {
		return a
	}

	return b
}

/*
func (prec *PrecisionStats) calcCDF(precs []float64, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]
	maxPrec := sortedPrecs[len(sortedPrecs)-1]
	for i := 0; i < prec.cdfResol; i++ {
		curPrec := minPrec + float64(i)*(maxPrec-minPrec)/float64(prec.cdfResol)
		for countSmaller, p := range sortedPrecs {
			if p >= curPrec {
				res[i].Prec = curPrec
				res[i].Count = countSmaller
				break
			}
		}
	}
}
*/

func calcmedian(values []Stats) (median Stats) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = values[i].Real
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].Real = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].Imag
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].Imag = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].L2
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].L2 = tmp[i]
	}

	index := len(values) / 2

	if len(values)&1 == 1 || index+1 == len(values) {
		return Stats{values[index].Real, values[index].Imag, values[index].L2}
	}

	return Stats{(values[index].Real + values[index+1].Real) / 2,
		(values[index].Imag + values[index+1].Imag) / 2,
		(values[index].L2 + values[index+1].L2) / 2}
}
