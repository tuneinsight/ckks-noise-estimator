package operations

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
)

func GetNoiseMulPt(LogN, LogScale int, std float64, nbRuns int) {

	fmt.Printf("Noise Plaintext Multiplication\n")

	for H := 32; H <= 32768; H <<= 1 {

		fmt.Printf("H:%d\n", H)

		c := NewContext(LogN, H, LogScale)

		params := c.params

		values := make([]float64, params.N())

		pt1 := ckks.NewPlaintext(params, params.MaxLevel())
		pt2 := ckks.NewPlaintext(params, params.MaxLevel())

		ct := ckks.NewCiphertext(params, 1, params.MaxLevel())

		est := estimator.NewEstimator(params.N(), params.HammingWeight(), params.Q(), params.P())

		for i := 0; i < nbRuns; i++ {

			c.GenKeys()

			enc := c.enc
			dec := c.dec
			eval := c.eval

			// First plaintext
			pt1.Scale = params.DefaultScale()
			c.NewPlaintextVector(std, values, pt1)

			// Encrypt first plaintext
			enc.Encrypt(pt1, ct)

			// Second plaintext
			pt2.Scale = rlwe.NewScale(params.Q()[pt2.Level()])
			c.NewPlaintextVector(std, values, pt2)

			// Enc(pt1) x pt2
			eval.Mul(ct, pt2, ct)

			ct.Scale = rlwe.NewScale(1)

			// pt1 x pt2
			params.RingQ().MForm(pt1.Value, pt1.Value)
			params.RingQ().MulCoeffsMontgomery(pt1.Value, pt2.Value, pt1.Value)
			pt1.Scale = ct.Scale

			// Dec(Enc(pt1) x pt2)
			dec.Decrypt(ct, pt2)

			params.RingQ().Sub(pt1.Value, pt2.Value, pt1.Value)

			estCT := estimator.NewCiphertextPK(estimator.NewPlaintext(std * params.DefaultScale().Float64(), 1/12.0, params.MaxLevel()))
			estPT := estimator.NewPlaintext(std * rlwe.NewScale(params.Q()[pt2.Level()]).Float64(), 1/12.0, params.MaxLevel())
			estCT = est.Mul(estCT, estPT)

			fmt.Println(math.Log2(STDPoly(params.RingQ(), pt1.Value, 1, false)), est.Std(estCT))
		}
		fmt.Println()
	}
}
