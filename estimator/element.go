package estimator

import (
	"math/big"
)

type Element struct {
	Level   int
	Message *big.Float
	Noise   []*big.Float
}

func (e Element) Degree() int {
	return len(e.Noise) - 1
}

func (e Element) CopyNew() Element {
	Noise := make([]*big.Float, len(e.Noise))

	for i := range Noise {
		Noise[i] = new(big.Float).Set(e.Noise[i])
	}

	return Element{
		Level:   e.Level,
		Message: new(big.Float).Set(e.Message),
		Noise:   Noise,
	}
}

func NewPlaintext(msg, noise interface{}, level int) Element {

	var Message *big.Float
	if msg == nil {
		Message = NewFloat(0)
	} else {
		Message = NewFloat(msg)
	}

	var Noise *big.Float
	if noise == nil {
		Noise = NewFloat(0)
	} else {
		Noise = NewFloat(noise)
	}

	return Element{
		Level:   level,
		Message: Message,
		Noise:   []*big.Float{Noise},
	}
}

func NewCiphertextSK(pt Element) Element {
	return Element{
		Level:   pt.Level,
		Message: NewFloat(pt.Message),
		Noise:   []*big.Float{AddSTD(NoiseFresh, pt.Noise[0]), NewFloat(0)},
	}
}

func NewCiphertextPK(pt Element) Element {

	oneTwelve := NewFloat(1)
	oneTwelve.Quo(oneTwelve, NewFloat(12))
	oneTwelve.Sqrt(oneTwelve)

	return Element{
		Level:   pt.Level,
		Message: NewFloat(pt.Message),
		Noise:   []*big.Float{AddSTD(oneTwelve, pt.Noise[0]), oneTwelve}, // ModDown by P kills all error
	}
}
