package main

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/ckks-noise-estimator"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func main() {

	LogN := 16
	LogScale := 45

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		LogQ:            []int{55, 60, 60, 60, 60, 60, 60, 60, 60, 53},
		LogP:            []int{61, 61, 61, 61},
		LogDefaultScale: LogScale,
		Xs:              ring.Ternary{H: 192},
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

	mod1Parameters := mod1.ParametersLiteral{
		LevelQ:          8,
		Mod1Type:        mod1.CosDiscrete,
		LogMessageRatio: 8,
		K:               12,
		Mod1Degree:      30,
		DoubleAngle:     3,
		LogScale:        60,
	}

	evm, err := mod1.NewParametersFromLiteral(params, mod1Parameters)

	if err != nil {
		panic(err)
	}

	mod1Eval := mod1.NewEvaluator(eval, polynomial.NewEvaluator(params, eval), evm)

	est := estimator.NewEstimator(params)

	statsHave := estimator.NewStats()
	statsWant := estimator.NewStats()

	for i := 0; i < 1; i++ {

		values, el, _, ct := NewTestVectorsMod1(est, params, ecd, pk, evm, -1, 1)

		// Scale the message to Delta = Q/MessageRatio
		scale := rlwe.NewScale(math.Exp2(math.Round(math.Log2(float64(params.Q()[0]) / evm.MessageRatio()))))
		scale = scale.Div(ct.Scale)

		if err = eval.ScaleUp(ct, rlwe.NewScale(math.Round(scale.Float64())), ct); err != nil {
			panic(err)
		}

		if err = est.ScaleUp(el, rlwe.NewScale(math.Round(scale.Float64()))); err != nil {
			panic(err)
		}

		// Scale the message up to Sine/MessageRatio
		scale = evm.ScalingFactor().Div(ct.Scale)
		scale = scale.Div(rlwe.NewScale(evm.MessageRatio()))

		if err = eval.ScaleUp(ct, rlwe.NewScale(math.Round(scale.Float64())), ct); err != nil {
			panic(err)
		}

		if err = est.ScaleUp(el, rlwe.NewScale(math.Round(scale.Float64()))); err != nil {
			panic(err)
		}

		// Normalization
		if err = eval.Mul(ct, 1/(evm.K*evm.QDiff), ct); err != nil {
			panic(err)
		}

		if err = est.Mul(el, 1/(evm.K*evm.QDiff), el); err != nil {
			panic(err)
		}

		if err = eval.Rescale(ct, ct); err != nil {
			panic(err)
		}

		if err = est.Rescale(el, el); err != nil {
			panic(err)
		}

		if ct, err = mod1Eval.EvaluateNew(ct); err != nil {
			panic(err)
		}

		if el, err = est.EvaluateMod1New(el, evm); err != nil {
			panic(err)
		}

		// ===================================
		ratio := new(big.Float).SetPrec(256).SetFloat64(evm.MessageRatio() * evm.QDiff / 6.28318530717958)
		for j := range values {
			values[j][0].Quo(values[j][0], ratio)
			values[j][0] = bignum.Sin(values[j][0])
			values[j][0].Mul(values[j][0], ratio)
		}

		// ===================================
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

func NewTestVectorsMod1(est estimator.Estimator, params ckks.Parameters, ecd *ckks.Encoder, key rlwe.EncryptionKey, evm mod1.Parameters, a, b float64) (values []*bignum.Complex, el *estimator.Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	slots := params.MaxSlots()

	values = make([]*bignum.Complex, slots)

	K := evm.K - 1
	Q := evm.QDiff * evm.MessageRatio()

	source := estimator.NewTestRand()

	for i := 0; i < slots; i++ {
		values[i] = bignum.ToComplex(math.Round(source.Float64(-K, K))*Q+source.Float64(a, b), 64)
	}

	values[0] = bignum.ToComplex(K*Q+0.5, 64)

	el = est.NewElement(values, 1, est.MaxLevel(), est.DefaultScale())
	est.AddEncodingNoise(el)

	pt = ckks.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(values, pt); err != nil {
		panic(err)
	}

	enc := rlwe.NewEncryptor(params, key)
	ct, _ = enc.EncryptNew(pt)
	switch key.(type) {
	case *rlwe.SecretKey:
		est.AddEncryptionNoiseSk(el)
	case *rlwe.PublicKey:
		est.AddEncryptionNoisePk(el)
	default:
		panic("INVALID ENCRYPION KEY")
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
