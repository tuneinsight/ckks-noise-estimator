package estimator

import (
	"math/big"
)

type Element struct {
	Level   int
	Message *big.Float
	Noise   []*big.Float
}

func (e *Element) Degree() int {
	return len(e.Noise) - 1
}
