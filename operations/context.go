package operations

import (
	"fmt"
	"math"
	"math/rand"
	"time"
	"math/big"

	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/ring/distribution"
)

type Context struct {
	params ckks.Parameters
	kgen   *rlwe.KeyGenerator
	enc    rlwe.Encryptor
	dec    rlwe.Decryptor
	ecd    *ckks.Encoder
	eval   *ckks.Evaluator
	std    float64
	sk     *rlwe.SecretKey
	pk     *rlwe.PublicKey
	stats  *stats.PrecisionStats // Precision stats
}

func NewContext(LogN, H, LogScale int) (c *Context) {

	var err error

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     []int{60, LogScale, LogScale, LogScale, LogScale, LogScale, LogScale, LogScale},
		LogP:     []int{61, 61},
		Xs:       &distribution.Ternary{H:H},
		LogScale: LogScale,
	}); err != nil {
		panic(err)
	}

	kgen := ckks.NewKeyGenerator(params)
	ecd := ckks.NewEncoder(params)
	eval := ckks.NewEvaluator(params, nil)

	return &Context{
		params: params,
		ecd:    ecd,
		kgen:   kgen,
		eval:   eval,
		stats:  stats.NewPrecisionStats(),
	}
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {
	c.sk = c.kgen.GenSecretKeyNew()
	c.pk = c.kgen.GenPublicKeyNew(c.sk)
	c.enc = ckks.NewEncryptor(c.params, c.pk)
	c.dec = ckks.NewDecryptor(c.params, c.sk)
}

// Computes the final precision stats
func (c *Context) Finalize() {
	c.stats.Finalize()
}

// Produces a CSV friendly string.
func (c *Context) ToCSV() []string {
	return c.stats.ToCSV()
}

func (c *Context) NewPlaintextVector(std float64, logSlots int, values []float64, pt *rlwe.Plaintext) (stdPt float64) {

	ecd := ckks.NewEncoder(c.params, 128)

	t := time.Now().UnixMilli()
	r := rand.New(rand.NewSource(t))

	gap := c.params.MaxSlots() / (1<<logSlots)

	for i := range values{
		values[i] = 0
	}

	var v float64
	for i := 0; i < len(values); i+= gap{

		v = r.NormFloat64()

		for v > 6.0 {
			v = r.NormFloat64()
		}

		values[i] = v * std
	}

	pt.EncodingDomain = rlwe.CoefficientsDomain
	ecd.Encode(values, pt)

	return standarddeviation(values)
}

func (c *Context) NewPlaintextSlots(min, max float64, logSlots int) (pt *rlwe.Plaintext, values []*big.Float, stdPt float64){

	ecd := c.ecd

	t := time.Now().UnixMilli()
	r := rand.New(rand.NewSource(t))

	slots := 1<<logSlots

	values = make([]*big.Float, slots)

	for i := range values{
		values[i] = new(big.Float).SetPrec(128).SetFloat64((max-min) * r.Float64() + min)
	}

	pt = ckks.NewPlaintext(c.params, c.params.MaxLevel())
	pt.EncodingDomain = rlwe.SlotsDomain
	pt.LogSlots = logSlots
	ecd.Encode(values, pt)

	valuesf64 := make([]float64, c.params.N())

	pt.EncodingDomain = rlwe.CoefficientsDomain
	ecd.Decode(pt, valuesf64)
	pt.EncodingDomain = rlwe.SlotsDomain

	//fmt.Println(valuesf64)

	stdPt = standarddeviation(valuesf64)
	fmt.Println(stdPt)
	return pt, values, stdPt
}


func standarddeviation(values []float64) (std float64){
	var mean float64
	for i := range values{
		mean += values[i]
	}
	mean /= float64(len(values))

	for i := range values{
		tmp := values[i]-mean
		std += tmp*tmp
	}

	std /= float64(len(values)-1)

	return math.Sqrt(std)
}