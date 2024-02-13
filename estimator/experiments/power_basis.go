package main

import (
	"fmt"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"

	//"github.com/tuneinsight/lattigo/v5/ring"
	"time"

	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func main() {

	LogN := 16
	LogScale := 55

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{60, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55},
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		//Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	ecd := hefloat.NewEncoder(params)

	fmt.Println(params.LogQ(), params.LogP())

	kgen := hefloat.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	enc := hefloat.NewEncryptor(params, pk)
	dec := hefloat.NewDecryptor(params, sk)

	evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk))

	eval := hefloat.NewEvaluator(params, evk)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	mul := bignum.NewComplexMultiplier().Mul

	n := 12

	for i := 0; i < 1; i++ {

		fmt.Println(i)

		values, el, _, ct := est.NewTestVector(ecd, enc, -1, 1)

		for k := 0; k < n; k++ {
			for l := range values {
				mul(values[l], values[l], values[l])
				values[l].Add(values[l], values[l])
				values[l].Add(values[l], &bignum.Complex{estimator.NewFloat(-1), estimator.NewFloat(0)})
			}
		}

		pbCt := hefloat.NewPowerBasis(ct, bignum.Chebyshev)
		runTimed(func() {
			if err := pbCt.GenPower(1<<n, false, eval); err != nil {
				panic(err)
			}
		})

		ct = pbCt.Value[1<<n]

		pbEl := estimator.NewPowerBasis(el, bignum.Chebyshev)

		runTimed(func() {
			pbEl.GenPower(1<<n, false, est)
		})
		el = pbEl.Value[1<<n]

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values, est.Decrypt(el), 0, false)
		pHave := hefloat.GetPrecisionStats(params, ecd, dec, values, ct, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())

	fmt.Println(estimator.ToLaTeXTable(LogN, LogScale, statsWant, statsHave))
}

func runTimed(f func()) {
	now := time.Now()
	f()
	fmt.Println(time.Since(now))
}

/*
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │  7.50 │  7.55 │  7.00 │
│MAX Prec │ 45.00 │ 42.81 │ 45.00 │
│AVG Prec │ 21.69 │ 21.69 │ 21.19 │
│MED Prec │ 21.42 │ 21.42 │ 20.92 │
│STD Prec │  2.47 │  2.47 │  2.47 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  2.19 │  0.00 │
│MAX Err  │ 37.50 │ 37.45 │ 38.00 │
│AVG Err  │ 23.31 │ 23.31 │ 23.81 │
│MED Err  │ 23.58 │ 23.58 │ 24.08 │
│STD Err  │  2.47 │  2.47 │  2.47 │
└─────────┴───────┴───────┴───────┘


┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │  7.69 │  7.48 │  7.19 │
│MAX Prec │ 44.78 │ 45.00 │ 44.28 │
│AVG Prec │ 21.69 │ 21.69 │ 21.19 │
│MED Prec │ 21.43 │ 21.43 │ 20.93 │
│STD Prec │  2.47 │  2.47 │  2.47 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.22 │  0.00 │  0.72 │
│MAX Err  │ 37.31 │ 37.52 │ 37.81 │
│AVG Err  │ 23.31 │ 23.31 │ 23.81 │
│MED Err  │ 23.57 │ 23.57 │ 24.07 │
│STD Err  │  2.47 │  2.47 │  2.47 │
└─────────┴───────┴───────┴───────┘


\begin{table}[]
    \centering
    \begin{tabular}{|c||c|c|c||c|c|c|}
    \hline
        & \multicolumn{3}{c||}{Predicted} & \multicolumn{3}{c|}{Actual}  \\
        \hline
        $\log_{2}$ & real & imag & l2 & real & imag & l2\\
        \hline
        Min &  7.50 &  7.55 &  7.00 &  7.69 &  7.48 &  7.19 \\
        Max & 45.00 & 42.81 & 45.00 & 44.78 & 45.00 & 44.28 \\
        AVG & 21.69 & 21.69 & 21.19 & 21.69 & 21.69 & 21.19 \\
        MED & 21.42 & 21.42 & 20.92 & 21.43 & 21.43 & 20.93 \\
        STD &  2.47 &  2.47 &  2.47 &  2.47 &  2.47 &  2.47 \\
        \hline
    \end{tabular}
    \caption{$N=2^{16}$ and $\Delta = 2^{45}$}
    \label{tab:my_label}
\end{table}
*/
