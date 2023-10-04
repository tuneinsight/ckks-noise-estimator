package estimator

import (
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Noise samples noise in R[X]/(X^N+1) according to f(), scales it by 2^-logScale and
// returns it in C^N/2.
func (p Parameters) Noise(f func() float64) (e []*bignum.Complex) {
	e = make([]*bignum.Complex, p.MaxSlots())

	scale := new(big.Int).SetInt64(1)
	scale.Lsh(scale, uint(p.LogScale))
	scaleF := new(big.Float).SetPrec(p.Prec).SetInt(scale)

	for i := range e {

		e[i] = &bignum.Complex{}

		ee := e[i]

		ee[0] = new(big.Float).SetPrec(p.Prec).SetFloat64(f())
		ee[1] = new(big.Float).SetPrec(p.Prec).SetFloat64(f())

		ee[0].Quo(ee[0], scaleF)
		ee[1].Quo(ee[1], scaleF)
	}

	// R[X]/(X^N+1) -> C^N/2
	if err := p.Encoder.FFT(e, p.LogMaxSlots()); err != nil {
		panic(err)
	}

	return
}

// RoundingNoise samples a rounding error in the ring and decode it into the canonical embeding
// Standard deviation: sqrt(1/12)
func (p Parameters) RoundingNoise() (e []*bignum.Complex) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	return p.Noise(func() float64 { return r.Float64() - 0.5 })
}

func (p Parameters) EncryptionNoiseSk() (e []*bignum.Complex) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	sigma := p.Sigma

	f := func() (s float64) {
		s = r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return s
	}

	return p.Noise(f)
}
