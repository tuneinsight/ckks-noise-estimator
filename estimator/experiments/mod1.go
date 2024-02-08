package main

import (
	"fmt"
	"math"
	"math/big"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func main() {

	LogN := 16
	LogScale := 45

	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 53},
		LogP:            []int{61, 61, 61, 61},
		LogDefaultScale: LogScale,
		Xs: ring.Ternary{H:192},
	})

	if err != nil {
		panic(err)
	}

	fmt.Println(params.LogQ(), params.LogP())

	ecd := hefloat.NewEncoder(params)

	kgen := hefloat.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	enc := hefloat.NewEncryptor(params, pk)
	dec := hefloat.NewDecryptor(params, sk)
	
	evk := rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk))

	eval := hefloat.NewEvaluator(params, evk)

	mod1Parameters := hefloat.Mod1ParametersLiteral{
		LevelQ:          8,
		Mod1Type:        hefloat.CosDiscrete,
		LogMessageRatio: 8,
		K:               12,
		Mod1Degree:      30,
		DoubleAngle:     3,
		LogScale:        60,
	}

	evm, err := hefloat.NewMod1ParametersFromLiteral(params, mod1Parameters)

	if err != nil{
		panic(err)
	}

	mod1Eval := hefloat.NewMod1Evaluator(eval, hefloat.NewPolynomialEvaluator(params, eval), evm)

	estParams := estimator.NewParameters(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 128; i++ {

		values, el, _, ct := NewTestVectorsMod1(estParams, params, ecd, enc, evm, -1, 1)

		// Scale the message to Delta = Q/MessageRatio
		scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / evm.MessageRatio()))))
		scale = scale.Div(ct.Scale)

		eval.ScaleUp(ct, rlwe.NewScale(math.Round(scale.Float64())), ct)
		el.ScaleUp(rlwe.NewScale(math.Round(scale.Float64())))

		// Scale the message up to Sine/MessageRatio
		scale = evm.ScalingFactor().Div(ct.Scale)
		scale = scale.Div(rlwe.NewScale(evm.MessageRatio()))
		
		eval.ScaleUp(ct, rlwe.NewScale(math.Round(scale.Float64())), ct)
		el.ScaleUp(rlwe.NewScale(math.Round(scale.Float64())))

		// Normalization
		if err = eval.Mul(ct, 1/(evm.K*evm.QDiff), ct); err != nil{
			panic(err)
		}

		el.Mul(el, 1/(evm.K*evm.QDiff))

		if err = eval.Rescale(ct, ct); err != nil{
			panic(err)
		}

		el.Rescale()

		if ct, err = mod1Eval.EvaluateNew(ct); err != nil{
			panic(err)
		}

		if el, err = el.EvaluateMod1(evm); err != nil{
			panic(err)
		}
		
		// ===================================
		ratio := new(big.Float).SetPrec(256).SetFloat64(evm.MessageRatio() * evm.QDiff / 6.28318530717958)
		for j := range values{
			values[j][0].Quo(values[j][0], ratio)
			values[j][0] = bignum.Sin(values[j][0])
			values[j][0].Mul(values[j][0], ratio)
		}

		el.Decrypt()
		el.Normalize()

		fmt.Println(i, values[0][0])

		// ===================================
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


func NewTestVectorsMod1(est estimator.Parameters, params hefloat.Parameters, ecd *hefloat.Encoder, enc *rlwe.Encryptor, evm hefloat.Mod1Parameters, a, b float64) (values []*bignum.Complex, el *estimator.Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	slots := params.MaxSlots()

	values = make([]*bignum.Complex, slots)

	K := evm.K-1
	Q := evm.QDiff * evm.MessageRatio()

	for i := 0; i < slots; i++ {
		values[i] = bignum.ToComplex(math.Round(sampling.RandFloat64(-K, K))*Q + sampling.RandFloat64(a, b), 64)
	}

	values[0] = bignum.ToComplex(K*Q + 0.5, 64)

	el = estimator.NewElement(est, values, 1, params.DefaultScale())
	el.AddEncodingNoise()

	pt = hefloat.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(values, pt); err != nil{
		panic(err)
	}

	if enc != nil {
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
│MIN Prec │ 30.28 │ 31.57 │ 29.78 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 36.59 │ 36.75 │ 36.09 │
│MED Prec │ 36.40 │ 36.46 │ 35.90 │
│STD Prec │  2.00 │  1.82 │  2.00 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │ 14.72 │ 13.43 │ 15.22 │
│AVG Err  │  8.41 │  8.25 │  8.91 │
│MED Err  │  8.60 │  8.54 │  9.10 │
│STD Err  │  2.00 │  1.82 │  2.00 │
└─────────┴───────┴───────┴───────┘


┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 30.17 │ 31.69 │ 29.67 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 36.56 │ 36.75 │ 36.06 │
│MED Prec │ 36.38 │ 36.46 │ 35.88 │
│STD Prec │  1.99 │  1.82 │  1.99 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │ 14.83 │ 13.31 │ 15.33 │
│AVG Err  │  8.44 │  8.25 │  8.94 │
│MED Err  │  8.62 │  8.54 │  9.12 │
│STD Err  │  1.99 │  1.82 │  1.99 │
└─────────┴───────┴───────┴───────┘


\begin{table}[]
    \centering
    \begin{tabular}{|c||c|c|c||c|c|c|}
    \hline
        & \multicolumn{3}{c||}{Predicted} & \multicolumn{3}{c|}{Actual}  \\
        \hline
        $\log_{2}$ & real & imag & l2 & real & imag & l2\\
        \hline
        Min & 30.28 & 31.57 & 29.78 & 30.17 & 31.69 & 29.67 \\
        Max & 45.00 & 45.00 & 45.00 & 45.00 & 45.00 & 45.00 \\
        AVG & 36.59 & 36.75 & 36.09 & 36.56 & 36.75 & 36.06 \\
        MED & 36.40 & 36.46 & 35.90 & 36.38 & 36.46 & 35.88 \\
        STD &  2.00 &  1.82 &  2.00 &  1.99 &  1.82 &  1.99 \\
        \hline
    \end{tabular}
    \caption{$N=2^{16}$ and $\Delta = 2^{45}$}
    \label{tab:my_label}
\end{table}
*/