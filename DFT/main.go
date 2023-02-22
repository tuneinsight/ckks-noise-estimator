package main

import (
	"encoding/csv"
	"fmt"
	"os"
	"runtime"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

var (
	// Scheme Parameters
	LogN     = 16   // Log2 ring degree
	LogSlots = 15   // Log2 #slots
	LogScale = 45   // Log2 scaling factor
	H        = 0    // Hamming weight of the secret
	KeyType  = "pk" // "pk"

	// Homomorphic Encoding/Decoding Parameters:
	LtType       = advanced.SlotsToCoeffs //advanced.CoeffsToSlots //
	LogBSGSRatio = 2                      // Ratio N2/N1
	Depth        = 10

	// Others
	NbRuns = 1    // Number of recorded events
	Record = true // Record in CSV
)

func main() {

	var s string
	switch LtType {
	case advanced.CoeffsToSlots:
		s = "IDFT"
	case advanced.SlotsToCoeffs:
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

	var stdEstimated float64
	for i := 0; i < NbRuns; i++ {
		// Generate a new set of keys
		// Evaluates the linear transform
		// Records the precision stats
		stdMsg, stdErr := c.ComputeStats()

		fmt.Println(stdMsg, stdErr)

		pt := estimator.NewPlaintext(stdMsg*scale, stdErr*scale, c.params.MaxLevel())

		var ct estimator.Element
		switch KeyType {
		case "sk":
			ct = estimator.NewCiphertextSK(pt)
		case "pk":
			ct = estimator.NewCiphertextPK(pt)
		}

		ec := c.encodingMatrixLiteral
		//ec.LogSlots++

		stdEstimated += est.Std(est.DFT(ct, c.params, c.ecd, ec))

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
	ecd                   ckks.Encoder            // Encoder in degree N
	enc                   rlwe.Encryptor          // Encryptor
	dec                   rlwe.Decryptor          // Decryptor
	kgen                  rlwe.KeyGenerator       // KeyGenerator
	eval                  ckks.Evaluator          // Evaluator
	rotations             []int                   // Rotations needed for the linear transform
	encodingMatrix        advanced.EncodingMatrix // Encoded linear transform
	encodingMatrixLiteral advanced.EncodingMatrixLiteral
	stats                 *stats.PrecisionStats // Precision stats
}

func NewContext(H, depth, logSlots int, ltType advanced.LinearTransformType) (c *Context) {

	var err error

	LogQ := make([]int, depth+1)

	LogQ[0] = 60
	for i := 1; i < depth+1; i++ {
		LogQ[i] = LogScale
	}

	LogP := []int{61, 61}

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     LogQ,
		LogP:     LogP,
		H:        H,
		LogSlots: logSlots,
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

	encodingMatrixLiteral := advanced.EncodingMatrixLiteral{
		LinearTransformType: ltType,
		LogN:                params.LogN(),
		LogSlots:            params.LogSlots(),
		LevelStart:          params.MaxLevel(),
		Levels:              Levels,
		LogBSGSRatio:        LogBSGSRatio,
	}

	// Gets the rotations indexes for CoeffsToSlots
	rotations := encodingMatrixLiteral.Rotations()

	// Generates the encoding matrices
	encodingMatrix := advanced.NewHomomorphicEncodingMatrixFromLiteral(encodingMatrixLiteral, ecd)

	eval := ckks.NewEvaluator(params, rlwe.EvaluationKey{})

	return &Context{
		params:                params,
		ecd:                   ecd,
		kgen:                  kgen,
		eval:                  eval,
		rotations:             rotations,
		encodingMatrix:        encodingMatrix,
		encodingMatrixLiteral: encodingMatrixLiteral,
		stats:                 stats.NewPrecisionStats(),
	}
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {

	sk := c.kgen.GenSecretKey()

	if H == 0 {
		sk = rlwe.NewSecretKey(c.params.Parameters)
	}

	switch KeyType {
	case "sk":
		c.enc = ckks.NewEncryptor(c.params, sk)
	case "pk":
		c.enc = ckks.NewEncryptor(c.params, c.kgen.GenPublicKey(sk))
	}

	c.dec = ckks.NewDecryptor(c.params, sk)
	rtks := c.kgen.GenRotationKeysForRotations(c.rotations, true, sk)
	c.eval = c.eval.WithKey(rlwe.EvaluationKey{Rlk: nil, Rtks: rtks})
}

// ComputeStats generates a new set of keys, evaluates the linear transform and records the precision stats.
func (c *Context) ComputeStats() (stdMsg, stdErr float64) {

	c.GenKeys() // Generates a new set of keys for each event

	params := c.params
	ecd := c.ecd
	enc := c.enc
	dec := c.dec
	eval := c.eval

	gap := params.N() / (2 * params.Slots())

	values := make([]float64, params.N())
	for i := 0; i < params.N(); i += 32 {
		values[i] = utils.RandFloat64(-1, 1)
	}

	stdMsg = estimator.STD(values)

	pt := ckks.NewPlaintext(params, params.MaxLevel())

	ecd.EncodeCoeffs(values, pt)

	vErr := ecd.DecodeCoeffs(pt)
	for i := range values {
		vErr[i] -= values[i]
	}
	stdErr = estimator.STD(vErr)

	ct := enc.EncryptNew(pt)

	nbMatrices := len(c.encodingMatrix.Matrices)

	DFT(eval, ct, c.encodingMatrix.Matrices[:nbMatrices])

	// =============== Plaintext circuit =============== 
	vCmplx := make([]complex128, params.Slots())

	for i, idx, jdx := 0, 0, params.N()>>1; i < params.Slots(); i, idx, jdx = i+1, idx+gap, jdx+gap {
		vCmplx[i] = complex(values[idx], values[jdx])
	}

	ckks.SliceBitReverseInPlaceComplex128(vCmplx, len(vCmplx))

	ptMatrices := c.encodingMatrixLiteral.ComputeDFTMatrices()

	
	for i := range ptMatrices[:nbMatrices]{
		vCmplx = EvaluateLinearTransform(vCmplx, ptMatrices[i], c.encodingMatrix.LogBSGSRatio)
	}

	ckks.SliceBitReverseInPlaceComplex128(vCmplx, len(vCmplx))

	for i, idx, jdx := 0, 0, params.N()>>1; i < params.Slots(); i, idx, jdx = i+1, idx+gap, jdx+gap {
		values[idx], values[jdx] = real(vCmplx[i]), imag(vCmplx[i])
	}

	// ================================================= 

	/*
	have := values
	want := ecd.DecodeCoeffs(dec.DecryptNew(ct))

	havegap := make([]float64, 2*params.Slots())

	for i := 0; i < 2*params.Slots(); i++{
		havegap[i] = have[i*gap] - want[i*gap]
	}

	fmt.Println(havegap)
	*/

	ecd.EncodeCoeffs(values, pt)
	eval.Sub(ct, pt, ct)
	ct.Scale = rlwe.NewScale(1)
	diff := ecd.DecodeCoeffs(dec.DecryptNew(ct))

	/*
	for i := 0; i < 2*params.Slots(); i++{
		diff[i] = diff[i*gap]
	}

	diff = diff[:2*params.Slots()]
	*/

	c.stats.Update(diff)

	return
}


func DFT(eval ckks.Evaluator, ct *rlwe.Ciphertext, plainVectors []ckks.LinearTransform) {
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