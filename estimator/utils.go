package estimator

import (
	"fmt"
	"math"
	"math/big"
)

const (
	prec = uint(256)
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

func AddSTD(a, b *big.Float) (c *big.Float) {
	c = new(big.Float).Mul(a, a)
	tmp := new(big.Float).Mul(b, b)
	c.Add(c, tmp)
	c.Sqrt(c)
	return
}

func SubSTD(a, b *big.Float) (c *big.Float) {
	c = new(big.Float).Mul(a, a)
	tmp := new(big.Float).Mul(b, b)
	c.Sub(c, tmp)
	c.Sqrt(c)
	return
}

func MulSTD(N, a, b *big.Float) (c *big.Float) {
	c = new(big.Float).Mul(a, b)
	tmp := new(big.Float).Quo(N, NewFloat(2))
	tmp.Sqrt(tmp)
	c.Mul(c, tmp)
	return
}

func MinInt(a, b int) (c int) {
	if a > b {
		return b
	}

	return a
}


func STD(slice []float64)(std float64){

	n := float64(len(slice))

	var mean float64
	for _, c := range slice{
		mean += c
	}

	mean /= n

	var tmp float64
	for _, c := range slice{
		tmp = (c - mean)
		std += tmp * tmp
	}

	std /= (n-1)

	return math.Sqrt(std)
}