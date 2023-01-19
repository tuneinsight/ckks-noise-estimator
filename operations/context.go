package operations

import(
	"math/rand"
	"time"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
)

type Context struct{
	params ckks.Parameters
	kgen rlwe.KeyGenerator
	enc rlwe.Encryptor
	dec rlwe.Decryptor
	ecd ckks.Encoder
	eval ckks.Evaluator
	std float64
	sk *rlwe.SecretKey
	pk *rlwe.PublicKey
	stats  *stats.PrecisionStats // Precision stats
}

func NewContext(LogN, H, LogScale int) (c *Context) {

	var err error

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     []int{60, LogScale, LogScale},
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
		stats:  stats.NewPrecisionStats(H),
	}
}

// GenKeys generates a new set of keys.
func (c *Context) GenKeys() {
	c.sk = c.kgen.GenSecretKey()
	c.pk = c.kgen.GenPublicKey(c.sk)
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

func (c *Context) NewPlaintextVector(std float64, values []float64, pt *rlwe.Plaintext) {

	ecd := c.ecd

	t := time.Now().UnixMilli()
	r := rand.New(rand.NewSource(t))

	var v float64
	for i := range values{

		v = r.NormFloat64() 

		for v > 6.0 {
			v = r.NormFloat64() 
		}

		values[i] = v * std
	}

	ecd.EncodeCoeffs(values, pt)

	return
}