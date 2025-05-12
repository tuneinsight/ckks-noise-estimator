package estimator

import (
	"math"
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Noise samples noise in R[X]/(X^N+1) according to f(), scales it by 2^-logScale and
// returns it in C^N/2.
func (e Estimator) Noise(f func() *big.Float) (noise []*bignum.Complex) {
	noise = make([]*bignum.Complex, e.MaxSlots())
	for i := range noise {
		noise[i] = &bignum.Complex{f(), f()}
	}

	// R[X]/(X^N+1) -> C^N/2
	if err := e.Encoder.FFT(noise, e.LogMaxSlots()); err != nil {
		panic(err)
	}

	return
}

// NoiseRingToCanonical samples a noisy vector with standard deviation
// sigma * sqrt(N/2), which emulates the sampling in the ring followed
// by the encoding of the noisy vector with the canonical embeding.
func (e Estimator) NoiseRingToCanonical(sigma float64) (noise []*bignum.Complex) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// R[X]/(X^N+1) -> C^N/2 increases variance by sqrt(N/2)
	sigma *= math.Sqrt(float64(e.N() / 2))

	f := func() *big.Float {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	noise = make([]*bignum.Complex, e.MaxSlots())
	for i := range noise {
		noise[i] = &bignum.Complex{f(), f()}
	}

	return
}

func (e Estimator) AddNoiseRingToCanonical(sigma float64, noise []*bignum.Complex) {
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// R[X]/(X^N+1) -> C^N/2 increases variance by sqrt(N/2)
	sigma *= math.Sqrt(float64(e.N() / 2))

	f := func() *big.Float {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	for i := range noise {
		noise[i].Add(noise[i], &bignum.Complex{f(), f()})
	}
}

// RoundingNoise samples a rounding error in the ring and decode it into the canonical embeding
// Standard deviation: sqrt(1/12)
func (e Estimator) RoundingNoise() (noise []*bignum.Complex) {

	if e.Heuristic {
		return e.NoiseRingToCanonical(math.Sqrt(1 / 12.0))
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	return e.Noise(func() *big.Float { return NewFloat(r.Float64() - 0.5) })
}

func (e Estimator) NormalNoise(sigma float64) (noise []*bignum.Complex) {

	if e.Heuristic {
		return e.NoiseRingToCanonical(sigma)
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	f := func() *big.Float {

		s := r.NormFloat64() * sigma

		if r.Int()&1 == 0 {
			s *= -1
		}

		return NewFloat(s)
	}

	return e.Noise(f)
}

// {eCt * sk + round(sum(e_i * qalphai)/P), round(1/2)}
func (e Estimator) KeySwitchingNoise(levelQ int, eCt, sk []*bignum.Complex) (e0, e1 []*bignum.Complex) {

	e0 = e.KeySwitchingNoiseRaw(levelQ, eCt, sk)

	P := e.P

	for i := range e0 {
		e0[i][0].Quo(e0[i][0], P)
		e0[i][1].Quo(e0[i][1], P)
	}

	r := e.RoundingNoise()

	for i := range e0 {
		e0[i].Add(e0[i], r[i])
	}

	return e0, e.RoundingNoise()
}

// Standard Deviation: sqrt(N * (var(noise_base) * var(SK) * P + var(noise_key) * sum(var(q_alpha_i)))
// Returns eCt * sk * P + sum(e_i * qalphai)
func (e Estimator) KeySwitchingNoiseRaw(levelQ int, eCt, sk []*bignum.Complex) (noise []*bignum.Complex) {

	P := e.P

	noise = make([]*bignum.Complex, e.MaxSlots())

	for i := range noise {
		noise[i] = bignum.NewComplex()
	}

	// eCt * P * sk
	mul := bignum.NewComplexMultiplier().Mul
	tmp := bignum.NewComplex()
	for i := range noise {
		mul(eCt[i], sk[i], tmp)
		tmp[0].Mul(tmp[0], P)
		tmp[1].Mul(tmp[1], P)
		noise[i].Add(noise[i], tmp)
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// var(noise_ct) * var(H) * P + var(ekey) * sum(var(q_alpha_i)))
	decompRNS := DecompRNS(levelQ, e.LevelP)

	if e.Heuristic {

		// sqrt(sum qalphi^2) * sqrt(1/12) * sqrt(N) * eSWK
		sumQAlhai := new(big.Float)

		for i := 0; i < decompRNS; i++ {

			start := i * (e.LevelP + 1)
			end := (i + 1) * (e.LevelP + 1)

			if i == decompRNS-1 {
				end = levelQ + 1
			}

			// prod[qi * ... * ]
			qalphai := NewFloat(1)
			for i := start; i < end; i++ {
				qalphai.Mul(qalphai, &e.Q[i])
			}

			// variances (delays)
			qalphai.Mul(qalphai, qalphai)

			// sum of variances
			sumQAlhai.Add(sumQAlhai, qalphai)
		}

		// Uniform distribution
		sumQAlhai.Mul(sumQAlhai, NewFloat(1/12.0))

		// Ring expansion
		sumQAlhai.Mul(sumQAlhai, NewFloat(e.N()))

		// Variance -> Standard Deviation
		sumQAlhai.Sqrt(sumQAlhai)

		f64, _ := sumQAlhai.Float64()

		e.AddNoiseRingToCanonical(f64*e.Sigma, noise)

		/*
			stdR, stdI := Log2STD(e)
			maxR, maxI := Log2MAX(e)
			avgR, avgI := Log2AVG(e)

			fmt.Println(stdR, stdI)
			fmt.Println(maxR, maxI)
			fmt.Println(avgR, avgI)
			fmt.Println()
		*/

	} else {
		for i := 0; i < decompRNS; i++ {

			start := i * (e.LevelP + 1)
			end := (i + 1) * (e.LevelP + 1)

			if i == decompRNS-1 {
				end = levelQ + 1
			}

			//std(q_alpha_i) * eKey
			qalphai := NewFloat(1)
			for i := start; i < end; i++ {
				qalphai.Mul(qalphai, &e.Q[i])
			}

			qalphaiHalf := new(big.Float).Quo(qalphai, new(big.Float).SetInt64(2))

			qalphaInt := new(big.Int)
			qalphai.Int(qalphaInt)

			ei := e.NormalNoise(e.Sigma)

			f := func() (x *big.Float) {
				y := bignum.RandInt(r, qalphaInt)
				x = new(big.Float).SetPrec(prec)
				x.SetInt(y)
				x.Sub(x, qalphaiHalf)
				return
			}

			pi := e.Noise(f)

			for i := range ei {
				mul(ei[i], pi[i], ei[i])
				noise[i].Add(noise[i], ei[i])
			}
		}

		/*
			stdR, stdI := Log2STD(noise)
			maxR, maxI := Log2MAX(noise)
			avgR, avgI := Log2AVG(noise)

			fmt.Println(stdR, stdI)
			fmt.Println(maxR, maxI)
			fmt.Println(avgR, avgI)
			fmt.Println()
		*/
	}

	return
}
