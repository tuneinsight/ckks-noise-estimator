package operations

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func GetNoiseMulCt(LogN, LogSlots, LogScale int, std float64, nbRuns int) {

	fmt.Printf("Noise Ciphertext Multiplication\n")

	for H := 32; H <= 32768; H <<= 1 {

		fmt.Printf("H:%d\n", H)

		c := NewContext(LogN, H, LogScale)

		params := c.params

		values := make([]float64, params.N())

		pt1 := ckks.NewPlaintext(params, params.MaxLevel())
		pt2 := ckks.NewPlaintext(params, params.MaxLevel())

		ct1 := ckks.NewCiphertext(params, 1, params.MaxLevel())
		ct2 := ckks.NewCiphertext(params, 1, params.MaxLevel())
		ct3 := ckks.NewCiphertext(params, 2, params.MaxLevel())

		for i := 0; i < nbRuns; i++ {

			c.GenKeys()

			ecd := c.ecd
			enc := c.enc
			dec := c.dec
			eval := c.eval

			// First plaintext
			pt1.Scale = params.DefaultScale()
			c.NewPlaintextVector(std, LogSlots, values, pt1)

			// Encrypt first plaintext
			enc.Encrypt(pt1, ct1)

			// Second plaintext
			pt2.Scale = rlwe.NewScale(params.Q()[pt2.Level()])
			c.NewPlaintextVector(std, LogSlots, values, pt2)

			enc.Encrypt(pt2, ct2)

			// Enc(pt1) x pt2
			eval.Mul(ct1, ct1, ct3)

			ct1.Scale = rlwe.NewScale(1)

			// pt1 x pt2
			params.RingQ().MForm(pt1.Value, pt1.Value)
			params.RingQ().MulCoeffsMontgomery(pt1.Value, pt1.Value, pt1.Value)
			params.RingQ().IMForm(pt1.Value, pt1.Value)
			pt1.Scale = ct1.Scale

			// Dec(Enc(pt1) x pt2)
			dec.Decrypt(ct3, pt2)

			params.RingQ().Sub(pt1.Value, pt2.Value, pt1.Value)

			// Compares in the ring
			values := make([]float64, params.N())
			ecd.Decode(pt1, values)
			c.stats.Update(values)
		}

		c.Finalize()
		fmt.Println(c.ToCSV())
	}
}
