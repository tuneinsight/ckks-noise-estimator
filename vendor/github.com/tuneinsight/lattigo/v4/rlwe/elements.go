package rlwe

import (
	"math/big"

	"github.com/tuneinsight/lattigo/v4/ring"
)

// AdditiveShare is a type for storing additively shared values in Z_Q[X] (RNS domain).
type AdditiveShare struct {
	Value ring.Poly
}

// AdditiveShareBigint is a type for storing additively shared values
// in Z (positional domain).
type AdditiveShareBigint struct {
	Value []*big.Int
}

// NewAdditiveShare instantiates a new additive share struct for the ring defined
// by the given parameters at maximum level.
func NewAdditiveShare(params Parameters) *AdditiveShare {
	return &AdditiveShare{Value: *ring.NewPoly(params.N(), 0)}
}

// NewAdditiveShareAtLevel instantiates a new additive share struct for the ring defined
// by the given parameters at level `level`.
func NewAdditiveShareAtLevel(params Parameters, level int) *AdditiveShare {
	return &AdditiveShare{Value: *ring.NewPoly(params.N(), level)}
}

// NewAdditiveShareBigint instantiates a new additive share struct composed of "n" big.Int elements.
func NewAdditiveShareBigint(params Parameters, n int) *AdditiveShareBigint {
	v := make([]*big.Int, n)
	for i := range v {
		v[i] = new(big.Int)
	}
	return &AdditiveShareBigint{Value: v}
}
