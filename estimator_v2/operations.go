package estimator

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

// Add adds op0 to op1 and writes the result on p, also returns p.
func (p *Plaintext) Add(op0 *Plaintext, op1 rlwe.Operand) *Plaintext {

	av := op0.Value
	cv := p.Value

	switch op1 := op1.(type) {
	case *Plaintext:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		bv := op1.Value
		for i := range cv {
			cv[i].Add(av[i], bv[i])
		}
	case []*bignum.Complex:

		if len(op0.Value) != len(p.Value) || len(op1) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		for i := range cv {
			cv[i].Add(av[i], op1[i])
		}
	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:
		bComplex := bignum.ToComplex(op1, cv[0].Prec())
		for i := range cv {
			cv[i].Add(av[i], bComplex)
		}

	default:
		panic(fmt.Errorf("invalid op1.(type): must be *Plaintext, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}

	return p
}

func (p *Plaintext) Mul(op0 *Plaintext, op1 rlwe.Operand) *Plaintext {
	av := op0.Value
	cv := p.Value

	mul := bignum.NewComplexMultiplier().Mul

	switch op1 := op1.(type) {
	case *Plaintext:

		if len(op0.Value) != len(p.Value) || len(op1.Value) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		bv := op1.Value
		for i := range cv {
			mul(av[i], bv[i], cv[i])
		}
	case []*bignum.Complex:

		if len(op0.Value) != len(p.Value) || len(op1) != len(p.Value) {
			panic("invalid input dimensions: do not match receiver dimension")
		}

		for i := range cv {
			mul(av[i], op1[i], cv[i])
		}

	case complex128, float64, int, int64, uint, *big.Int, *big.Float, *bignum.Complex:
		bComplex := bignum.ToComplex(op1, cv[0].Prec())
		for i := range cv {
			mul(av[i], bComplex, cv[i])
		}

	default:
		panic(fmt.Errorf("invalid op1.(type): must be *Plaintext, complex128, float64, int, int64, uint, uint64, *big.Int, *big.Float, *bignum.Complex, []complex128, []float64, []*big.Float or []*bignum.Complex, but is %T", op1))
	}

	return p
}
