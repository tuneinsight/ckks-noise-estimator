package estimator

import (
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Noise samples noise in R[X]/(X^N+1) according to f(), scales it by 2^-logScale and
// returns it in C^N/2.
func (p Parameters) Noise(f func() *big.Float) (e []*bignum.Complex) {
	e = make([]*bignum.Complex, p.MaxSlots())
	for i := range e {
		e[i] = &bignum.Complex{f(), f()}
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
	return p.Noise(func() *big.Float { return NewFloat(r.Float64() - 0.5) })
}

func (p Parameters) EncryptionNoiseSk() (e []*bignum.Complex) {

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	sigma := p.Sigma

	f := func() (*big.Float) {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	return p.Noise(f)
}

// Standard Deviation: sqrt(N * (var(noise_base) * var(SK) * P + var(noise_key) * sum(var(q_alpha_i)))
func (p Parameters) KeySwitchingNoise(){
	//[-as + m + e0, a + e1]
}
