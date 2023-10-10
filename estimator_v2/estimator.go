package estimator

import(
	"time"
	"math/rand"
	"math/big"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

type Parameters struct {
	Encoder
	LogN     int
	Sigma    float64
	Scale  big.Float
	H      int
	Q      []*big.Float
	P      *big.Float
	Sk [][]*bignum.Complex
}

func NewParameters(p ckks.Parameters) Parameters {

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
	sk := make([]*bignum.Complex, N)
	for i := 0; i < H; i++{
		if i &1 == 0{
			sk[i] = &bignum.Complex{NewFloat(1), NewFloat(0)}
		}else{
			sk[i] = &bignum.Complex{NewFloat(-1), NewFloat(0)}
		}
		
	}
	for i := H; i < N; i++{
		sk[i] = &bignum.Complex{NewFloat(0), NewFloat(0)}
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	r.Shuffle(N, func(i, j int){sk[i], sk[j] = sk[j], sk[i]})

	// R[X]/(X^N+1) -> C^N/2
	if err := Encoder.FFT(sk, p.LogMaxSlots()); err != nil {
		panic(err)
	}

	mul := bignum.NewComplexMultiplier().Mul
	sk2 := make([]*bignum.Complex, N)
	for i := range sk2{
		sk2[i] = &bignum.Complex{NewFloat(0), NewFloat(0)}
		mul(sk[i], sk[i], sk2[i])
	}

	return Parameters{
		LogN:     p.LogN(),
		Sigma:    p.NoiseFreshSK(),
		H:        H,
		Encoder: *Encoder,
		Q:Qi,
		P:Pi,
		Scale:p.DefaultScale().Value,
		Sk: [][]*bignum.Complex{sk, sk2},
	}
}

func (p Parameters) MaxSlots() int {
	return 1 << p.LogMaxSlots()
}

func (p Parameters) LogMaxSlots() int {
	return p.LogN - 1
}

func (p Parameters) MaxLevel() int {
	return len(p.Q)-1
}

func (p Parameters) DefaultScale() big.Float{
	return p.Scale
}

