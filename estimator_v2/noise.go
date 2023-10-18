package estimator

import (
	"fmt"
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

func (p Parameters) EncryptionNoisePk() (e []*bignum.Complex) {
	e1 := p.RoundingNoise()
	e0 := p.RoundingNoise()
	mul := bignum.NewComplexMultiplier().Mul
	for i := range e0{
		mul(e1[i], p.Sk[0][i], e1[i])
		e0[i].Add(e0[i], e1[i])
	}

	return e0
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
// Returns eCt * sk * P + sum(e_i * qalphai)
func (p Parameters) KeySwitchingNoise(levelQ int, eCt, sk []*bignum.Complex) (e []*bignum.Complex){

	P := p.P

	e = make([]*bignum.Complex, p.MaxSlots())

	for i := range e{
		e[i] = bignum.NewComplex()
	}

	// eCt * P * sk
	mul := bignum.NewComplexMultiplier().Mul
	tmp := bignum.NewComplex()
	for i := range e{
		mul(eCt[i], sk[i], tmp)
		tmp[0].Mul(tmp[0], P)
		tmp[1].Mul(tmp[1], P)
		e[i].Add(e[i], tmp)
	}

	/*
	// var(noise_ct) * var(H) * P + var(ekey) * sum(var(q_alpha_i)))
	decompRNS := DecompRNS(levelQ, p.LevelP)
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

		ei := p.EncryptionNoiseSk()
		for i := range ei{
			ei[i][0].Mul(ei[i][0], qalphai)
			ei[i][1].Mul(ei[i][1], qalphai)

			e[i].Add(e[i], ei[i])
		}
	}
	*/

	for i := range e{
		e[i][0].Quo(e[i][0], P)
		e[i][1].Quo(e[i][1], P)
	}

	r := p.RoundingNoise()

	for i := range e{
		e[i].Add(e[i], r[i])
	}

	for i := 0; i < 8; i++{
		fmt.Println(e[i][0], e[i][1])
	}

	return
}
