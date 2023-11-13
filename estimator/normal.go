package estimator

import (
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)

func (p Parameters) SampleNormal(sigma *big.Float) (e []*bignum.Complex) {

	e = make([]*bignum.Complex, p.MaxSlots())
	for i := range e {
		e[i] = &bignum.Complex{new(big.Float), new(big.Float)}
	}

	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	tmp := new(big.Float)

	if sigma.Cmp(new(big.Float).SetFloat64(0x20000000000000)) >= 0 {

		norm := new(big.Float)
		normInt := new(big.Int)

		var negative bool
		for i := range e {

			for j := 0; j < 2; j++ {

				norm.Mul(tmp.SetFloat64(r.NormFloat64()), sigma)

				if norm.Cmp(new(big.Float)) == -1 {
					negative = true
					norm.Neg(norm)
				} else {
					negative = false
				}

				norm.Int(normInt)

				e[i][j].SetInt(bignum.RandInt(r, normInt))

				if negative {
					e[i][j].Neg(e[i][j])
				}
			}
		}
	} else {
		for i := range e {
			for j := 0; j < 2; j++ {
				e[i][j].Mul(tmp.SetFloat64(r.NormFloat64()), sigma)
			}
		}
	}

	return
}
