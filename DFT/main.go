package main

import (
	"encoding/csv"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var (
	// Scheme Parameters
	LogN     = 16   // Log2 ring degree
	LogSlots = 14   // Log2 #slots
	LogScale = 45   // Log2 scaling factor
	H        = 0    // Hamming weight of the secret
	KeyType  = "pk" // "pk"

	// Homomorphic Encoding/Decoding Parameters:
	LtType       = advanced.Decode //advanced.Encode //
	LogBSGSRatio = 2               // Ratio N2/N1
	Depth        = LogSlots

	// Others
	NbRuns = 1     // Number of recorded events
	Record = false // Record in CSV
)

func main() {

	var s string
	switch LtType {
	case advanced.Encode:
		s = "IDFT"
	case advanced.Decode:
		s = "DFT"
	}

	var f *os.File
	var w *csv.Writer
	var err error

	if Record {
		if f, err = os.Create(fmt.Sprintf("data/experiment_%s_%s_H_%d_logslots_%d_runs_%d_id_%d.csv", s, KeyType, H, LogSlots, NbRuns, time.Now().Unix())); err != nil {
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

	fmt.Printf("%s - %s - H:%d - Depth: %d\n", s, KeyType, H, Depth)

	c := NewContext(H, Depth, LogSlots, LtType)

	est := estimator.NewEstimator(c.params.N(), H, c.params.Q(), c.params.P())

	scale := c.params.DefaultScale().Float64()

	nbMatrices := Depth

	var stdEstimated float64
	for i := 0; i < NbRuns; i++ {
		// Generate a new set of keys
		// Evaluates the linear transform
		// Records the precision stats
		stdMsg, stdErr := c.ComputeStats(nbMatrices, LogSlots)

		fmt.Println(stdMsg, stdErr)

		pt := estimator.NewPlaintext(stdMsg*scale, stdErr*scale, c.params.MaxLevel())

		var ct estimator.Element
		switch KeyType {
		case "sk":
			ct = estimator.NewCiphertextSK(pt)
		case "pk":
			ct = estimator.NewCiphertextPK(pt)
		}

		LTs := estimator.GetEncodingMatrixSTD(c.params, c.ecd, c.encodingMatrixLiteral)

		stdEstimated += est.Std(est.DFT(ct, LTs[nbMatrices-1:nbMatrices]))

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

	if Record {
		if err := w.Write(data); err != nil {
			panic(err)
		}

		w.Flush()
	}
}

type Context struct {
	params                ckks.Parameters
	ecd                   *ckks.Encoder                 // Encoder in degree N
	enc                   rlwe.Encryptor                // Encryptor
	dec                   rlwe.Decryptor                // Decryptor
	kgen                  *rlwe.KeyGenerator            // KeyGenerator
	eval                  *ckks.Evaluator               // Evaluator
	galEls                []uint64                      // Rotations needed for the linear transform
	encodingMatrix        advanced.HomomorphicDFTMatrix // Encoded linear transform
	encodingMatrixLiteral advanced.HomomorphicDFTMatrixLiteral
	stats                 *stats.PrecisionStats // Precision stats
}

func NewContext(H, depth, logSlots int, ltType advanced.DFTType) (c *Context) {

	var err error

	LogQ := make([]int, depth+1)

	LogQ[0] = 60
	for i := 1; i < depth+1; i++ {
		LogQ[i] = LogScale
	}

	LogP := []int{61, 61}

	var Xs distribution.Distribution
	if H > 0 {
		Xs = &distribution.Ternary{H: H}
	}

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     LogQ,
		LogP:     LogP,
		Xs:Xs,
		LogScale: LogScale,
	}); err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params)
	ecd := ckks.NewEncoder(params)

	Levels := make([]int, params.MaxLevel())

	for i := range Levels {
		Levels[i] = 1
	}

	encodingMatrixLiteral := advanced.HomomorphicDFTMatrixLiteral{
		Type:         ltType,
		LogSlots:     logSlots,
		LevelStart:   params.MaxLevel(),
		Levels:       Levels,
		LogBSGSRatio: LogBSGSRatio,
	}

	// Gets the rotations indexes for CoeffsToSlots
	galEls := encodingMatrixLiteral.GaloisElements(params)

	// Generates the encoding matrices
	encodingMatrix := advanced.NewHomomorphicDFTMatrixFromLiteral(encodingMatrixLiteral, ecd)

	eval := ckks.NewEvaluator(params, nil)

	return &Context{
		params:                params,
		ecd:                   ecd,
		kgen:                  kgen,
		eval:                  eval,
		galEls:                galEls,
		encodingMatrix:        encodingMatrix,
		encodingMatrixLiteral: encodingMatrixLiteral,
		stats:                 stats.NewPrecisionStats(),
	}
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {

	sk := c.kgen.GenSecretKeyNew()

	if H == 0 {
		sk = rlwe.NewSecretKey(c.params.Parameters)
	}

	switch KeyType {
	case "sk":
		c.enc = ckks.NewEncryptor(c.params, sk)
	case "pk":
		c.enc = ckks.NewEncryptor(c.params, c.kgen.GenPublicKeyNew(sk))
	}

	c.dec = ckks.NewDecryptor(c.params, sk)

	evk := rlwe.NewEvaluationKeySet()

	for _, galEl := range c.galEls {
		if err := evk.Add(c.kgen.GenGaloisKeyNew(galEl, sk)); err != nil {
			panic(err)
		}
	}

	c.eval = c.eval.WithKey(evk)
}

// ComputeStats generates a new set of keys, evaluates the linear transform and records the precision stats.
func (c *Context) ComputeStats(nbMatrices int, LogSlots int) (stdMsg, stdErr float64) {

	c.GenKeys() // Generates a new set of keys for each event

	params := c.params
	ecd := c.ecd
	enc := c.enc
	dec := c.dec
	eval := c.eval

	Slots := 1 << LogSlots

	gap := params.N() / (2 * Slots)

	r := rand.New(rand.NewSource(0))

	values := make([]float64, params.N())
	for i := 0; i < params.N(); i += gap {
		values[i] = 2 * (r.Float64() - 0.5)
	}

	stdMsg = estimator.STD(values)

	pt := ckks.NewPlaintext(params, params.MaxLevel())
	pt.EncodingDomain = rlwe.CoefficientsDomain

	ecd.Encode(values, pt)

	vErr := make([]float64, params.N())
	ecd.Decode(pt, vErr)
	for i := range values {
		vErr[i] -= values[i]
	}
	stdErr = estimator.STD(vErr)

	ct := enc.EncryptNew(pt)

	DFT(eval, ct, c.encodingMatrix.Matrices[nbMatrices-1:nbMatrices])

	// =============== Plaintext circuit ===============

	// R^{2N} -> C^{N}
	vCmplx := make([]complex128, Slots)
	for i, idx, jdx := 0, 0, params.N()>>1; i < Slots; i, idx, jdx = i+1, idx+gap, jdx+gap {
		vCmplx[i] = complex(values[idx], values[jdx])
	}
	ecd.FFT(vCmplx, LogSlots)

	// DFT: C^{N} -> C^{N}
	ptMatrices := c.encodingMatrixLiteral.GenMatrices(params.LogN())
	for _, pt := range ptMatrices[nbMatrices-1:nbMatrices] {
		vCmplx = EvaluateLinearTransform(vCmplx, pt, c.encodingMatrix.LogBSGSRatio)
	}

	// C^{N} -> R^{2N}
	ecd.IFFT(vCmplx, LogSlots)
	for i, idx, jdx := 0, 0, params.N()>>1; i < Slots; i, idx, jdx = i+1, idx+gap, jdx+gap {
		values[idx], values[jdx] = real(vCmplx[i]), imag(vCmplx[i])
	}

	// =================================================

	ecd.Encode(values, pt)
	eval.Sub(ct, pt, ct)
	ct.Scale = rlwe.NewScale(1)
	ecd.Decode(dec.DecryptNew(ct), values)

	c.stats.Update(values)

	return
}

func DFT(eval *ckks.Evaluator, ct *rlwe.Ciphertext, plainVectors []ckks.LinearTransform) {
	scale := ct.Scale
	for _, plainVector := range plainVectors {
		eval.LinearTransform(ct, plainVector, []*rlwe.Ciphertext{ct})
		if err := eval.Rescale(ct, scale, ct); err != nil {
			panic(err)
		}
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

func EvaluateLinearTransform(values []complex128, diags map[int][]complex128, LogBSGSRatio int) (res []complex128) {

	slots := len(values)

	N1 := ckks.FindBestBSGSRatio(diags, len(values), LogBSGSRatio)

	index, _, _ := ckks.BSGSIndex(diags, slots, N1)

	res = make([]complex128, slots)

	for j := range index {

		rot := -j & (slots - 1)

		tmp := make([]complex128, slots)

		for _, i := range index[j] {

			v, ok := diags[j+i]
			if !ok {
				v = diags[j+i-slots]
			}

			a := utils.RotateComplex128Slice(values, i)

			b := utils.RotateComplex128Slice(v, rot)

			for i := 0; i < slots; i++ {
				tmp[i] += a[i] * b[i]
			}
		}

		tmp = utils.RotateComplex128Slice(tmp, j)

		for i := 0; i < slots; i++ {
			res[i] += tmp[i]
		}
	}

	return
}
