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

	f := func(x float64) (y float64) {
		return 1 / (math.Exp(-x) + 1)
	}

	deg := 255

	k := 32.0

	poly := polynomial.NewPolynomial(GetChebyshevPoly(k, deg, f))

	scale := new(big.Float).SetFloat64(math.Exp2(float64(LogScale)))

	truncVal := new(big.Int).SetUint64(3)

	var acc float64

	wantReal := new(big.Int)
	wantImag := new(big.Int)

	for i := range poly.Coeffs {
		Truncate(poly.Coeffs[i], scale, truncVal, wantReal, wantImag)
	}

	seed := int64(0)

	wantInt, wantCmplx := Eval(truncVal, k, est, ecd, pk, dec, polyEval, &poly, seed)

	for i := 0; i < 1024; i++ {

		fmt.Println(i)

		haveInt, haveCmplx := Eval(truncVal, k, est, ecd, pk, dec, polyEval, &poly, seed)

		stats := ckks.GetPrecisionStats(params, ecd, dec, haveCmplx, wantCmplx, 0, false)

		fmt.Println(stats)

		for i := range wantInt {

			if wantInt[i].Cmp(haveInt[i]) != 0 {
				//fmt.Println(i, want[i], have[i])
				acc++
			}
		}
	}

	fmt.Println(acc)
}

func Eval(trunc *big.Int, k float64, est estimator.Estimator, ecd *ckks.Encoder, pk rlwe.EncryptionKey, dec *rlwe.Decryptor, polyEval *polynomial.Evaluator, poly *polynomial.Polynomial, seed ...int64) ([]*big.Int, []*big.Float) {

	source := estimator.NewTestRand(seed...)

	want, _, _, ct := est.NewTestVectorFromSeed(ecd, pk, complex(-k+1, 0), complex(k-1, 0), source, trunc)

	for i := range want {
		want[i] = poly.Evaluate(want[i])
	}

	var err error

	scalar, constant := poly.ChangeOfBasis()

	if scalar.Cmp(new(big.Float).SetInt64(1)) != 0 {

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

	if ct, err = polyEval.Evaluate(ct, *poly, ct.Scale); err != nil {
		panic(err)
	}

	pt := dec.DecryptNew(ct)
	pt.IsBatched = false

	out := make([]*big.Float, 2*pt.Slots())

	ecd.Decode(pt, out)

	return est.TruncatePlaintext(pt, trunc), out
}

func Truncate(x *bignum.Complex, scale *big.Float, trunc, xReal, xImag *big.Int) {
	x[0].Mul(x[0], scale)
	x[0].Int(xReal)
	xReal.Quo(xReal, trunc)
	DivRound(xReal, trunc, xReal)
	xReal.Mul(xReal, trunc)
	x[0].SetInt(xReal)
	x[0].Quo(x[0], scale)
	x[1].Mul(x[1], scale)
	x[1].Int(xImag)
	DivRound(xImag, trunc, xImag)
	xImag.Mul(xImag, trunc)
	x[1].SetInt(xImag)
	x[1].Quo(x[1], scale)
	return
}

// DivRound sets the target i to round(a/b).
func DivRound(a, b, i *big.Int) {
	_a := new(big.Int).Set(a)
	i.Quo(_a, b)
	r := new(big.Int).Rem(_a, b)
	r2 := new(big.Int).Mul(r, bignum.NewInt(2))
	if r2.CmpAbs(b) != -1.0 {
		if _a.Sign() == b.Sign() {
			i.Add(i, bignum.NewInt(1))
		} else {
			i.Sub(i, bignum.NewInt(1))
		}
	}
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
\begin{figure}
\centering
\begin{tabular}{|c|c|c|}
    \hline
    $\log_{2}(\Delta')$ & $\log_{2}(\#x)$ & $f(x_{i}) \neq \textsf{Dec}(f(\textsf{Enc}(x_{i})))$\\
    \hline
    $17$ & 23 & 8121\\
    $18$ & 23 & 4318\\
    $19$ & 23 & 2267\\
    $20$ & 23 & 791\\
    $21$ & 23 & 686\\
    $22$ & 23 & 222 \\
    $23$ & 23 & 68\\
    $24$ & 28 & 0\\
    \hline
\end{tabular}
\caption{Caption}
\label{fig:enter-label}
\end{figure}
*/
