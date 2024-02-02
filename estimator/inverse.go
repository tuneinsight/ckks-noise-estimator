package estimator

import(
	"math"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func (el *Element) GoldschmidtDivision(log2min float64) {

	// 2^{-(prec - LogN + 1)}
	prec := float64(el.N()/2) / el.DefaultScale().Float64()

	// Estimates the number of iterations required to achieve the desired precision, given the interval [min, 2-min]
	start := 1 - math.Exp2(log2min)
	var iters = 1
	for start >= prec {
		start *= start // Doubles the bit-precision at each iteration
		iters++
	}

	iters = utils.Max(iters, 3)

	a := el
	a.Mul(a, -1)
	b := a.CopyNew()
	a.Add(a, 2)
	b.Add(b, 1)

	tmp := NewElement[*bignum.Complex](*el.Parameters, nil, 1, a.Scale)

	for j := 1; j < iters; j++{

		b.Mul(b, b)
		b.Relinearize()
		b.Rescale()

		tmp.Mul(a, b)
		tmp.Relinearize()
		tmp.Rescale()

		a.SetScale(tmp.Scale)

		a.Add(a, tmp)
	}
}