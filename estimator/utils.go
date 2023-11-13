package estimator

import (
	"fmt"
	"math"
	"math/big"
)

const (
	prec = uint(128)
)

func DecompRNS(levelQ, levelP int) int {
	return (levelQ + levelP + 1) / (levelP + 1)
}

func NewFloat(x interface{}) (s *big.Float) {
	switch x := x.(type) {
	case float64:
		if math.IsNaN(x) || math.IsInf(x, 0) {
			panic(fmt.Errorf("x cannot be negative, NaN or Inf, but is %f", x))
		}
		s = new(big.Float).SetPrec(prec)
		s.SetFloat64(x)
		return
	case *big.Int:
		s = new(big.Float).SetPrec(prec)
		s.SetInt(x)
		return
	case *big.Float:
		s = new(big.Float).SetPrec(prec)
		s.Set(x)
		return
	case int:
		return NewFloat(new(big.Int).SetInt64(int64(x)))
	case int64:
		return NewFloat(new(big.Int).SetInt64(x))
	case uint64:
		return NewFloat(new(big.Int).SetUint64(x))
	default:
		panic(fmt.Errorf("invalid x.(type): must be int, int64, uint64, float64, *big.Int or *big.Float but is %T", x))
	}
}
