package main

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/ckks-noise-estimator"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func main() {

	LogN := 16
	LogScale := 55

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{60, 55, 55, 55, 55, 55, 55, 55, 55, 55},
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		//Xs: ring.Ternary{H:192},
	})

	if err != nil {
		panic(err)
	}

	fmt.Println(params.LogQ(), params.LogP())

	ecd := ckks.NewEncoder(params)

	kgen := ckks.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	dec := ckks.NewDecryptor(params, sk)

	evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk))

	eval := ckks.NewEvaluator(params, evk)

	polyEval := polynomial.NewEvaluator(params, eval)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	f := func(x float64) (y float64) {
		return 1 / (math.Exp(-x) + 1)
	}

	deg := 255

	k := 32.0

	poly := polynomial.NewPolynomial(GetChebyshevPoly(k, deg, f))

	scalar, constant := poly.ChangeOfBasis()

	for i := 0; i < 1; i++ {

		fmt.Println(i)

		values, el, _, ct := est.NewTestVector(ecd, pk, complex(-k+1, 0), complex(k-1, 0), 20)

		for i := range values {
			values[i] = poly.Evaluate(values[i])
		}

		if scalar.Cmp(new(big.Float).SetInt64(1)) != 0 {

			if err = est.Mul(el, scalar, el); err != nil {
				panic(err)
			}

			if err = est.Add(el, constant, el); err != nil {
				panic(err)
			}

			if err = est.Rescale(el, el); err != nil {
				panic(err)
			}

			if err = polyEval.Mul(ct, scalar, ct); err != nil {
				panic(err)
			}

			if err = polyEval.Add(ct, constant, ct); err != nil {
				panic(err)
			}

			if err = polyEval.Rescale(ct, ct); err != nil {
				panic(err)
			}
		}

		if el, err = est.EvaluatePolynomialNew(el, poly, el.Scale); err != nil {
			panic(err)
		}

		if ct, err = polyEval.Evaluate(ct, poly, ct.Scale); err != nil {
			panic(err)
		}

		pWant := ckks.GetPrecisionStats(params, ecd, dec, values, est.Decrypt(el), 0, false)
		pHave := ckks.GetPrecisionStats(params, ecd, dec, values, ct, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())

	fmt.Println(estimator.ToLaTeXTable(LogN, LogScale, statsWant, statsHave))
}

// GetChebyshevPoly returns the Chebyshev polynomial approximation of f the
// in the interval [-K, K] for the given degree.
func GetChebyshevPoly(K float64, degree int, f64 func(x float64) (y float64)) bignum.Polynomial {

	FBig := func(x *big.Float) (y *big.Float) {
		xF64, _ := x.Float64()
		return new(big.Float).SetPrec(x.Prec()).SetFloat64(f64(xF64))
	}

	var prec uint = 128

	interval := bignum.Interval{
		A:     *bignum.NewFloat(-K, prec),
		B:     *bignum.NewFloat(K, prec),
		Nodes: degree,
	}

	// Returns the polynomial.
	return bignum.ChebyshevApproximation(FBig, interval)
}

/*
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 25.71 │ 25.79 │ 25.21 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 30.67 │ 30.67 │ 30.17 │
│MED Prec │ 30.37 │ 30.37 │ 29.87 │
│STD Prec │  1.85 │  1.85 │  1.85 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │ 19.29 │ 19.21 │ 19.79 │
│AVG Err  │ 14.33 │ 14.33 │ 14.83 │
│MED Err  │ 14.63 │ 14.63 │ 15.13 │
│STD Err  │  1.85 │  1.85 │  1.85 │
└─────────┴───────┴───────┴───────┘


┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 25.76 │ 25.97 │ 25.26 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 30.67 │ 30.67 │ 30.17 │
│MED Prec │ 30.37 │ 30.37 │ 29.87 │
│STD Prec │  1.85 │  1.86 │  1.85 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │ 19.24 │ 19.03 │ 19.74 │
│AVG Err  │ 14.33 │ 14.33 │ 14.83 │
│MED Err  │ 14.63 │ 14.63 │ 15.13 │
│STD Err  │  1.85 │  1.86 │  1.85 │
└─────────┴───────┴───────┴───────┘


\begin{table}[]
    \centering
    \begin{tabular}{|c||c|c|c||c|c|c|}
    \hline
        & \multicolumn{3}{c||}{Predicted} & \multicolumn{3}{c|}{Actual}  \\
        \hline
        $\log_{2}$ & real & imag & l2 & real & imag & l2\\
        \hline
        Min & 25.71 & 25.79 & 25.21 & 25.76 & 25.97 & 25.26 \\
        Max & 45.00 & 45.00 & 45.00 & 45.00 & 45.00 & 45.00 \\
        AVG & 30.67 & 30.67 & 30.17 & 30.67 & 30.67 & 30.17 \\
        MED & 30.37 & 30.37 & 29.87 & 30.37 & 30.37 & 29.87 \\
        STD &  1.85 &  1.85 &  1.85 &  1.85 &  1.86 &  1.85 \\
        \hline
    \end{tabular}
    \caption{$N=2^{16}$ and $\Delta = 2^{45}$}
    \label{tab:my_label}
\end{table}
*/
