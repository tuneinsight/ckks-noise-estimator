package estimator

import (
	"math"
	"math/big"
)

var NoiseFresh = NewFloat(3.2)

type Estimator struct {
	N      *big.Float
	H      *big.Float
	Q      []*big.Float
	P      *big.Float
	LevelP int
}

func NewEstimator(N, H int, Q, P []uint64) *Estimator {

	Qi := make([]*big.Float, len(Q))

	for i := range Qi {
		Qi[i] = NewFloat(Q[i])
	}

	Pi := NewFloat(1)
	for i := range P {
		Pi.Mul(Pi, NewFloat(P[i]))
	}

	return &Estimator{
		N:      NewFloat(N),
		H:      NewFloat(H),
		Q:      Qi,
		P:      Pi,
		LevelP: len(P) - 1,
	}
}

func (e *Estimator) Std(el0 Element) (stdf64 float64) {

	std := new(big.Float).Mul(el0.Noise[0], el0.Noise[0])

	sk := new(big.Float).Set(e.H)

	var tmp = new(big.Float)

	for i := 1; i < el0.Degree()+1; i++ {
		tmp.Mul(el0.Noise[i], el0.Noise[i])
		tmp.Mul(sk, tmp)
		std.Add(std, tmp)
		sk.Mul(sk, e.H)
	}

	std.Sqrt(std)

	stdf64, _ = std.Float64()

	return math.Log2(stdf64)
}

func (e *Estimator) Add(el0, el1 Element) (el2 Element) {

	el2 = Element{}

	el2.Level = MinInt(el0.Level, el1.Level)

	el2.Message = AddSTD(el0.Message, el1.Message)

	if el0.Degree() >= el1.Degree() {

		el2.Noise = make([]*big.Float, el0.Degree()+1)

		for i := 0; i < el1.Degree()+1; i++ {
			el2.Noise[i] = AddSTD(el0.Noise[i], el1.Noise[i])
		}

		for i := el1.Degree() + 1; i < el0.Degree()+1; i++ {
			el2.Noise[i] = new(big.Float).Set(el0.Noise[i])
		}

	} else {

		el2.Noise = make([]*big.Float, el1.Degree()+1)

		for i := 0; i < el0.Degree()+1; i++ {
			el2.Noise[i] = AddSTD(el0.Noise[i], el1.Noise[i])
		}

		for i := el0.Degree() + 1; i < el1.Degree()+1; i++ {
			el2.Noise[i] = new(big.Float).Set(el1.Noise[i])
		}
	}

	return
}

func (e *Estimator) Mul(el0, el1 Element) (el2 Element) {

	el2 = Element{}

	el2.Level = MinInt(el0.Level, el1.Level)

	el2.Message = MulSTD(e.N, el0.Message, el1.Message)

	el2.Noise = make([]*big.Float, el0.Degree()+el1.Degree()+1)

	for i := range el2.Noise {
		el2.Noise[i] = new(big.Float)
	}

	e0 := new(big.Float)
	e1 := new(big.Float)

	for i := 0; i < el0.Degree()+1; i++ {

		e0.Set(el0.Noise[i])

		if i == 0 {
			e0 = AddSTD(e0, el0.Message)
		}

		for j := 0; j < el1.Degree()+1; j++ {

			e1.Set(el1.Noise[i])

			if j == 0 {
				e1 = AddSTD(e1, el1.Message)
			}

			el2.Noise[i+j] = AddSTD(el2.Noise[i+j], MulSTD(e.N, e0, e1))
		}
	}

	el2.Noise[0] = SubSTD(el2.Noise[0], el2.Message)

	return
}

func (e *Estimator) Relinearize(el0 Element) (el1 Element) {

	el1 = Element{}

	el1.Message = new(big.Float).Set(el0.Message)

	el1.Noise = []*big.Float{
		new(big.Float).Set(el0.Noise[0]),
		new(big.Float).Set(el0.Noise[1]),
	}

	for i := 2; i < el0.Degree()+1; i++ {

		tmp := e.ModDown(Element{
			Message: NewFloat(0),
			Noise: []*big.Float{
				e.KeySwitchRawLazy(el0.Level, el0.Noise[i], NoiseFresh, e.H),
				NewFloat(0),
			},
		})

		el1 = e.Add(el1, tmp)
	}

	return
}

// KeySwitchRawLazy: raw output of the key-switching dot-product
// without the division by P. Thus only the error is due to multiplication
// of uniform term of the norm of the RNS basis with the error in the key.
// Additional error due to the multiplication with the encrypted key with the error
// in the term given to the key-switch.
//
// Noise: sqrt(N * (var(noise_base) * var(H) * P + var(noise_key) * sum(var(q_alpha_i)))
//
//
func (e *Estimator) KeySwitchRawLazy(level int, eCt, eKey, H *big.Float) (eKeySwitch *big.Float) {

	// var(noise_ct) * var(H) * P
	eKeySwitch = new(big.Float).Set(e.P)
	eKeySwitch.Mul(eKeySwitch, eCt)
	eKeySwitch.Mul(eKeySwitch, eCt)
	eKeySwitch.Mul(eKeySwitch, H)

	// var(noise_ct) * var(H) * P + var(ekey) * sum(var(q_alpha_i)))
	decompRNS := DecompRNS(level, e.LevelP)
	for i := 0; i < decompRNS; i++ {

		start := i * (e.LevelP + 1)
		end := (i + 1) * (e.LevelP + 1)

		if i == decompRNS-1 {
			end = level + 1
		}

		qalpha := NewFloat(1)

		for i := start; i < end; i++ {
			qalpha.Mul(qalpha, e.Q[i])
		}

		qalpha.Mul(qalpha, qalpha)
		qalpha.Quo(qalpha, NewFloat(12))
		qalpha.Mul(qalpha, eKey)
		qalpha.Mul(qalpha, eKey)

		eKeySwitch.Add(eKeySwitch, qalpha)
	}

	// N * (P * var(eCt) + var(ekey) * sum(var(q_alpha_i)))
	eKeySwitch.Mul(eKeySwitch, e.N)

	// sqrt(N * (P * var(eCt) + var(ekey) * sum(var(q_alpha_i))))
	eKeySwitch.Sqrt(eKeySwitch)

	return
}

func (e *Estimator) ModDown(el0 Element) (el1 Element) {
	return e.DivRound(el0, e.P)
}

func (e *Estimator) Rescale(el0 Element) (el1 Element) {
	el1 = e.DivRound(el0, e.Q[el0.Level])
	el1.Level = el0.Level - 1
	return
}

func (e *Estimator) DivRound(el0 Element, q *big.Float) (el1 Element) {

	el1.Message = new(big.Float).Quo(el0.Message, q)

	el1.Noise = make([]*big.Float, el0.Degree()+1)

	oneTwelve := NewFloat(1)
	oneTwelve.Quo(oneTwelve, NewFloat(12))
	oneTwelve.Sqrt(oneTwelve)

	for i := 0; i < el0.Degree()+1; i++ {
		el1.Noise[i] = new(big.Float).Quo(el0.Noise[i], q)
		el1.Noise[i] = AddSTD(el1.Noise[i], oneTwelve)
	}

	return
}
