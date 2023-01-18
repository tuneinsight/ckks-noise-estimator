package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
	"time"

	"stats"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

var (
	LogN     = 16   // Log2 ring degree
	LogScale = 45   // Log2 scaling factor
	NbRuns   = 1024 // Number of recorded events
)

func main() {

	f, err := os.Create(fmt.Sprintf("data/experiment_rounding_%d_%d_%d.csv", LogN, NbRuns, time.Now().Unix()))
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

		fmt.Printf("H:%d\n", H)

		c := NewContext(H)

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

type Context struct {
	params ckks.Parameters
	ecd    ckks.Encoder          // Encoder
	enc    rlwe.Encryptor        // Encryptor
	dec    rlwe.Decryptor        // Decryptor
	kgen   rlwe.KeyGenerator     // KeyGenerator
	eval   ckks.Evaluator        // Evaluator
	stats  *stats.PrecisionStats // Precision stats
}

func NewContext(H int) (c *Context) {

	var err error

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     []int{60, 45},
		LogP:     []int{61},
		H:        H,
		LogSlots: LogN - 1,
		LogScale: LogScale,
	}); err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params)
	ecd := ckks.NewEncoder(params)
	eval := ckks.NewEvaluator(params, rlwe.EvaluationKey{})

	return &Context{
		params: params,
		ecd:    ecd,
		kgen:   kgen,
		eval:   eval,
		stats:  stats.NewPrecisionStats(LogN, H, params.MaxLevel(), params.LogSlots(), int(math.Round(math.Log2(params.DefaultScale().Float64())))),
	}
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {
	sk := c.kgen.GenSecretKey()
	c.enc = ckks.NewEncryptor(c.params, sk)
	c.dec = ckks.NewDecryptor(c.params, sk)
}

// ComputeStats generates a new set of keys, evaluates the linear transform and records the precision stats.
func (c *Context) ComputeStats() {

	c.GenKeys() // Generates a new set of keys for each event

	params := c.params
	ecd := c.ecd
	enc := c.enc
	dec := c.dec
	eval := c.eval

	// Generates a vector of zero values
	values := make([]float64, params.N())

	// Encodes and encrypts in the ring
	plaintext := ckks.NewPlaintext(params, params.MaxLevel())
	ecd.EncodeCoeffs(values, plaintext)
	ciphertext := enc.EncryptNew(plaintext)

	params.RingQ().DivRoundByLastModulusManyNTT(1, ciphertext.Value[0], eval.BuffQ()[0], ciphertext.Value[0])
	params.RingQ().DivRoundByLastModulusManyNTT(1, ciphertext.Value[1], eval.BuffQ()[0], ciphertext.Value[1])

	ciphertext.Value[0].Coeffs = ciphertext.Value[0].Coeffs[:1]
	ciphertext.Value[1].Coeffs = ciphertext.Value[1].Coeffs[:1]

	ciphertext.Scale = rlwe.NewScale(1)

	c.stats.Update(values, ecd.DecodeCoeffs(dec.DecryptNew(ciphertext)))
}

// Computes the final precision stats
func (c *Context) Finalize() {
	c.stats.Finalize()
}

// Produces a CSV friendly string.
func (c *Context) ToCSV() []string {
	return c.stats.ToCSV()
}
