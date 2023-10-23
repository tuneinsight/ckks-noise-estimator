package estimator

import (
	"math"
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

// NoiseRingToCanonical samples a noisy vector with standard deviation
// sigma * sqrt(N/2), which emulates the sampling in the ring followed
// by the encoding of the noisy vector with the canonical embeding.
func (p Parameters) NoiseRingToCanonical(sigma float64) (e []*bignum.Complex) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// R[X]/(X^N+1) -> C^N/2 increases variance by sqrt(N/2)
	sigma *= math.Sqrt(float64(p.N() / 2))

	f := func() *big.Float {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	e = make([]*bignum.Complex, p.MaxSlots())
	for i := range e {
		e[i] = &bignum.Complex{f(), f()}
	}

	return
}

func (p Parameters) AddNoiseRingToCanonical(sigma float64, e []*bignum.Complex) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// R[X]/(X^N+1) -> C^N/2 increases variance by sqrt(N/2)
	sigma *= math.Sqrt(float64(p.N() / 2))

	f := func() *big.Float {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	for i := range e {
		e[i].Add(e[i], &bignum.Complex{f(), f()})
	}
}

// RoundingNoise samples a rounding error in the ring and decode it into the canonical embeding
// Standard deviation: sqrt(1/12)
func (p Parameters) RoundingNoise() (e []*bignum.Complex) {

	if p.Heuristic {
		return p.NoiseRingToCanonical(math.Sqrt(1 / 12.0))
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	return p.Noise(func() *big.Float { return NewFloat(r.Float64() - 0.5) })
}

func (p Parameters) NormalNoise(sigma float64) (e []*bignum.Complex) {

	if p.Heuristic {
		return p.NoiseRingToCanonical(sigma)
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	f := func() *big.Float {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	return p.Noise(f)
}

func (p Parameters) KeySwitchingNoise(levelQ int, eCt, sk []*bignum.Complex) (e0, e1 []*bignum.Complex) {

	e := p.KeySwitchingNoiseRaw(levelQ, eCt, sk)

	P := p.P

	for i := range e {
		e[i][0].Quo(e[i][0], P)
		e[i][1].Quo(e[i][1], P)
	}

	r := p.RoundingNoise()

	for i := range e {
		e[i].Add(e[i], r[i])
	}

	return e, p.RoundingNoise()
}

// Standard Deviation: sqrt(N * (var(noise_base) * var(SK) * P + var(noise_key) * sum(var(q_alpha_i)))
// Returns eCt * sk * P + sum(e_i * qalphai)
func (p Parameters) KeySwitchingNoiseRaw(levelQ int, eCt, sk []*bignum.Complex) (e []*bignum.Complex) {

	P := p.P

	e = make([]*bignum.Complex, p.MaxSlots())

	for i := range e {
		e[i] = bignum.NewComplex()
	}

	// eCt * P * sk
	mul := bignum.NewComplexMultiplier().Mul
	tmp := bignum.NewComplex()
	for i := range e {
		mul(eCt[i], sk[i], tmp)
		tmp[0].Mul(tmp[0], P)
		tmp[1].Mul(tmp[1], P)
		e[i].Add(e[i], tmp)
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// var(noise_ct) * var(H) * P + var(ekey) * sum(var(q_alpha_i)))
	decompRNS := DecompRNS(levelQ, p.LevelP)

	if p.Heuristic {

		// sqrt(sum qalphi^2) * sqrt(1/12) * sqrt(N) * eSWK
		sumQAlhai := new(big.Float)

		for i := 0; i < decompRNS; i++ {

			start := i * (p.LevelP + 1)
			end := (i + 1) * (p.LevelP + 1)

			if i == decompRNS-1 {
				end = levelQ + 1
			}

			// prod[qi * ... * ]
			qalphai := NewFloat(1)
			for i := start; i < end; i++ {
				qalphai.Mul(qalphai, p.Q[i])
			}

			// variances (delays)
			qalphai.Mul(qalphai, qalphai)

			// sum of variances
			sumQAlhai.Add(sumQAlhai, qalphai)
		}

		// Uniform distribution
		sumQAlhai.Mul(sumQAlhai, NewFloat(1/12.0))

		// Ring expansion
		sumQAlhai.Mul(sumQAlhai, NewFloat(p.N()))

		// Variance -> Standard Deviation
		sumQAlhai.Sqrt(sumQAlhai)

		f64, _ := sumQAlhai.Float64()

		p.AddNoiseRingToCanonical(f64*p.Sigma, e)

	} else {
		for i := 0; i < decompRNS; i++ {

			start := i * (p.LevelP + 1)
			end := (i + 1) * (p.LevelP + 1)

			if i == decompRNS-1 {
				end = levelQ + 1
			}

			//std(q_alpha_i) * eKey
			qalphai := NewFloat(1)
			for i := start; i < end; i++ {
				qalphai.Mul(qalphai, p.Q[i])
			}

			qalphaiHalf := new(big.Float).Quo(qalphai, new(big.Float).SetInt64(2))

			qalphaInt := new(big.Int)
			qalphai.Int(qalphaInt)

			ei := p.NormalNoise(p.Sigma)

			f := func() (x *big.Float) {
				y := bignum.RandInt(r, qalphaInt)
				x = new(big.Float).SetPrec(prec)
				x.SetInt(y)
				x.Sub(x, qalphaiHalf)
				return
			}

			pi := p.Noise(f)

			for i := range ei {
				mul(ei[i], pi[i], ei[i])
				e[i].Add(e[i], ei[i])
			}
		}
	}

	return
}
