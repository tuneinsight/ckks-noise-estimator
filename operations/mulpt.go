package operations

import (
	"encoding/csv"
	"fmt"
	"os"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func GetNoiseMulPt(LogN, LogScale int, std float64, nbRuns int) {

	f, err := os.Create(fmt.Sprintf("data/mulpt_%d_%d_%f_%d_%d.csv", LogN, LogScale, std, nbRuns, time.Now().Unix()))
	if err != nil {
		panic(err)
	}
	defer f.Close()

	w := csv.NewWriter(f)

	// CSV Header
	if err := w.Write(stats.Header); err != nil {
		panic(err)
	}

	w.Flush()

	fmt.Printf("Noise Plaintext Multiplication\n")

	for H := 32; H <= 32768; H <<= 1 {

		fmt.Printf("H:%d\n", H)

		c := NewContext(LogN, H, LogScale)

		params := c.params

		values := make([]float64, params.N())

		pt1 := ckks.NewPlaintext(params, params.MaxLevel())
		pt2 := ckks.NewPlaintext(params, params.MaxLevel())

		ct := ckks.NewCiphertext(params, 1, params.MaxLevel())

		for i := 0; i < nbRuns; i++ {

			c.GenKeys()

			ecd := c.ecd
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

			// Compares in the ring
			c.stats.Update(ecd.DecodeCoeffs(pt1))
		}

		c.Finalize()

		if err := w.Write(c.ToCSV()); err != nil {
			panic(err)
		}

		w.Flush()
	}
}
