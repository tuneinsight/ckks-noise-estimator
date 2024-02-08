package main

import (
	"fmt"
	"math"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/ring"
)

func main() {

	LogN := 16
	LogScale := 55

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{60, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55},
		LogP:            []int{61, 61, 61},
		LogDefaultScale: LogScale,
		Xs: ring.Ternary{H:192},
	})

	if err != nil {
		panic(err)
	}

	ecd := hefloat.NewEncoder(params)

	kgen := hefloat.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	enc := hefloat.NewEncryptor(params, pk)
	dec := hefloat.NewDecryptor(params, sk)

	evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk))

	eval := hefloat.NewEvaluator(params, evk)

	inverseEval := hefloat.NewInverseEvaluator(params, eval, nil)

	estParams := estimator.NewParameters(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	log2min := -4.0

	// 2^{-(prec - LogN + 1)}
	prec := float64(params.N()/2) / params.DefaultScale().Float64()

	// Estimates the number of iterations required to achieve the desired precision, given the interval [min, 2-min]
	start := 1 - math.Exp2(log2min)
	var iters = 1
	for start >= prec {
		start *= start // Doubles the bit-precision at each iteration
		iters++
	}

	// Minimum of 3 iterations
	// This minimum is set in the case where log2min is close to 0.
	iters = utils.Max(iters, 3)

	fmt.Println(iters)

	for i := 0; i < 128; i++ {

		fmt.Println(i)

		values, el, _, ct := NewTestVector(estParams, ecd, enc, 0.1, 1.9)

		GoldschmidtDivision(values, iters)

		el.GoldschmidtDivision(log2min)

		if ct, err = inverseEval.GoldschmidtDivisionNew(ct, log2min); err != nil {
			panic(err)
		}

		el.Decrypt()
		el.Normalize()

		pWant := hefloat.GetPrecisionStats(params, ecd, dec, values, el.Value[0], 0, false)
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


func GoldschmidtDivision(values []*bignum.Complex, iters int){

	a := bignum.NewComplex()
	b := bignum.NewComplex()
	tmp := bignum.NewComplex()

	prec := values[0].Prec()

	one := bignum.NewFloat(1, prec)
	two := bignum.NewFloat(2, prec)

	mul := bignum.NewComplexMultiplier().Mul

	for i := range values{

		a[0].Neg(values[i][0])
		a[1].Neg(values[i][1])
		b[0].Add(a[0], one)
		b[1].Set(a[1])
		a[0].Add(a[0], two)

		for j := 1; j < iters; j++{
			mul(b, b, b)
			mul(a, b, tmp)
			a.Add(a, tmp)
		}

		values[i].Set(a)
	}
}

func NewTestVector(p estimator.Parameters, ecd *hefloat.Encoder, enc *rlwe.Encryptor, a, b float64) (values []*bignum.Complex, el *estimator.Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	params := p.Parameters

	prec := ecd.Prec()

	values = make([]*bignum.Complex, params.MaxSlots())

	for i := range values {
		values[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(a, b), prec),
			bignum.NewFloat(0, prec),
		}
	}

	values[0][0].SetFloat64(1)
	values[0][1].SetFloat64(0)

	el = estimator.NewElement(p, values, 1, params.DefaultScale())
	el.AddEncodingNoise()

	pt = hefloat.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(values, pt); err != nil{
		panic(err)
	}

	if enc != nil{
		ct, _ = enc.EncryptNew(pt)
		switch enc.KeyType().(type) {
		case *rlwe.SecretKey:
			el.AddEncryptionNoiseSk()
		case *rlwe.PublicKey:
			el.AddEncryptionNoisePk()
		default:
			panic("INVALID ENCRYPION KEY")
		}
	}
	
	return
}

/*
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 21.80 │ 22.08 │ 21.30 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 30.04 │ 30.04 │ 29.54 │
│MED Prec │ 29.86 │ 29.86 │ 29.36 │
│STD Prec │  2.08 │  2.08 │  2.08 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │ 23.20 │ 22.92 │ 23.70 │
│AVG Err  │ 14.96 │ 14.96 │ 15.46 │
│MED Err  │ 15.14 │ 15.14 │ 15.64 │
│STD Err  │  2.08 │  2.08 │  2.08 │
└─────────┴───────┴───────┴───────┘


┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 21.49 │ 22.10 │ 20.99 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 30.05 │ 30.04 │ 29.55 │
│MED Prec │ 29.86 │ 29.86 │ 29.36 │
│STD Prec │  2.08 │  2.08 │  2.08 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │ 23.51 │ 22.90 │ 24.01 │
│AVG Err  │ 14.95 │ 14.96 │ 15.45 │
│MED Err  │ 15.14 │ 15.14 │ 15.64 │
│STD Err  │  2.08 │  2.08 │  2.08 │
└─────────┴───────┴───────┴───────┘


\begin{table}[]
    \centering
    \begin{tabular}{|c||c|c|c||c|c|c|}
    \hline
        & \multicolumn{3}{c||}{Predicted} & \multicolumn{3}{c|}{Actual}  \\
        \hline
        $\log_{2}$ & real & imag & l2 & real & imag & l2\\
        \hline
        Min & 21.80 & 22.08 & 21.30 & 21.49 & 22.10 & 20.99 \\
        AVG & 30.04 & 30.04 & 29.54 & 30.05 & 30.04 & 29.55 \\
        STD &  2.08 &  2.08 &  2.08 &  2.08 &  2.08 &  2.08 \\
        \hline
    \end{tabular}
    \caption{$N=2^{16}$ and $\Delta = 2^{45}$}
    \label{tab:my_label}
\end{table}
*/