package estimator

import (
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

type Parameters struct {
	Parameters hefloat.Parameters
	Encoder
	LogN      int
	Sigma     float64
	Scale     rlwe.Scale
	H         int
	Q         []*big.Float
	P         *big.Float
	LevelP    int
	Sk        [][]*bignum.Complex
	Heuristic bool
}

func NewParameters(p hefloat.Parameters) (params Parameters) {

	params = Parameters{}
	params.Parameters = p
	params.LogN = p.LogN()
	params.Sigma = p.NoiseFreshSK()
	params.Scale = p.DefaultScale()
	params.Heuristic = true

	Q := p.Q()
	P := p.P()

	Qi := make([]*big.Float, len(Q))
	for i := range Q {
		Qi[i] = NewFloat(Q[i])
	}
	params.Q = Qi

	Pi := NewFloat(1)
	for i := range P {
		Pi.Mul(Pi, NewFloat(P[i]))
	}
	params.P = Pi
	params.LevelP =  len(P) - 1

	params.Encoder = *NewEncoder(p.LogN(), prec)

	params.H = utils.Min(p.N(), p.XsHammingWeight())

	// Samples a secret-key
	sk := params.SampleSecretKey(params.H)

	mul := bignum.NewComplexMultiplier().Mul
	sk2 := make([]*bignum.Complex, p.MaxSlots())
	for i := range sk2 {
		sk2[i] = &bignum.Complex{NewFloat(0), NewFloat(0)}
		mul(sk[i], sk[i], sk2[i])
	}

	params.Sk = [][]*bignum.Complex{sk, sk2}

	return 
}

func (p Parameters) SampleSecretKey(H int) (sk []*bignum.Complex) {

	N := p.N()

	H = utils.Min(N, H)

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

	sk = make([]*bignum.Complex, p.MaxSlots())
	for i := range sk {
		sk[i] = &bignum.Complex{skF[2*i], skF[2*i+1]}
	}

	// R[X]/(X^N+1) -> C^N/2
	if err := p.Encoder.FFT(sk, p.LogMaxSlots()); err != nil {
		panic(err)
	}

	return
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

func (p Parameters) DefaultScale() rlwe.Scale {
	return p.Scale
}
