package estimator

import (
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v4/he/hefloat"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Parameters struct {
	Encoder
	LogN      int
	Sigma     float64
	Scale     big.Float
	H         int
	Q         []*big.Float
	P         *big.Float
	LevelP    int
	Sk        [][]*bignum.Complex
	Heuristic bool
}

func NewParameters(p hefloat.Parameters) Parameters {

	Q := p.Q()
	P := p.P()

	Qi := make([]*big.Float, len(Q))

	for i := range Q {
		Qi[i] = NewFloat(Q[i])
	}

	Pi := NewFloat(1)
	for i := range P {
		Pi.Mul(Pi, NewFloat(P[i]))
	}

	Encoder := NewEncoder(p.LogN(), prec)

	// Samples a secret-key
	H := p.XsHammingWeight()
	N := p.N()
	skF := make([]*big.Float, N)
	for i := 0; i < H; i++ {
		if i&1 == 0 {
			skF[i] = NewFloat(1)
		} else {
			skF[i] = NewFloat(-1)
		}

	}
	for i := H; i < N; i++ {
		skF[i] = NewFloat(0)
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	r.Shuffle(len(skF), func(i, j int) { skF[i], skF[j] = skF[j], skF[i] })

	sk := make([]*bignum.Complex, p.MaxSlots())
	for i := range sk {
		sk[i] = &bignum.Complex{skF[2*i], skF[2*i+1]}
	}

	// R[X]/(X^N+1) -> C^N/2
	if err := Encoder.FFT(sk, p.LogMaxSlots()); err != nil {
		panic(err)
	}

	mul := bignum.NewComplexMultiplier().Mul
	sk2 := make([]*bignum.Complex, p.MaxSlots())
	for i := range sk2 {
		sk2[i] = &bignum.Complex{NewFloat(0), NewFloat(0)}
		mul(sk[i], sk[i], sk2[i])
	}

	return Parameters{
		LogN:    p.LogN(),
		Sigma:   p.NoiseFreshSK(),
		H:       H,
		Encoder: *Encoder,
		Q:       Qi,
		P:       Pi,
		LevelP:  len(P) - 1,
		Scale:   p.DefaultScale().Value,
		Sk:      [][]*bignum.Complex{sk, sk2},
	}
}

func (p Parameters) N() int {
	return 1 << p.LogN
}

func (p Parameters) MaxSlots() int {
	return 1 << p.LogMaxSlots()
}

func (p Parameters) LogMaxSlots() int {
	return p.LogN - 1
}

func (p Parameters) MaxLevel() int {
	return len(p.Q) - 1
}

func (p Parameters) DefaultScale() big.Float {
	return p.Scale
}
