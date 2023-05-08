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

// GetNoisRescale
//
// sqrt(1/12 + N/12 * var(h)) = sqrt(1/12 + N/12 * h/N) = 1/12 + h/12 = (h+1)/12
//
func GetNoisRescale(LogN, LogScale int, std float64, nbRuns int) {

	f, err := os.Create(fmt.Sprintf("data/rescale_%d_%d_%f_%d_%d.csv", LogN, LogScale, std, nbRuns, time.Now().Unix()))
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

	fmt.Printf("Noise Rescale\n")
	for H := 32; H <= 32768; H <<= 1 {

		fmt.Printf("H:%d\n", H)

		c := NewContext(LogN, H, LogScale)

		params := c.params

		pt := ckks.NewPlaintext(params, params.MaxLevel())
		pt.Scale = rlwe.NewScale(params.Q()[params.MaxLevel()])

		ct := ckks.NewCiphertext(params, 1, params.MaxLevel())

		for i := 0; i < nbRuns; i++ {

			c.GenKeys()

			ecd := c.ecd
			enc := c.enc
			dec := c.dec
			eval := c.eval

			ct := rlwe.NewCiphertextAtLevelFromPoly(params.MaxLevel(), ct.Value)

			pt.Value.Zero()

			enc.Encrypt(pt, ct)

			eval.Rescale(ct, rlwe.NewScale(1), ct)

			ptTmp := rlwe.NewPlaintextAtLevelFromPoly(ct.Level(), pt.Value)
			ptTmp.MetaData = ct.MetaData
			ptTmp.EncodingDomain = rlwe.CoefficientsDomain

			dec.Decrypt(ct, ptTmp)

			values := make([]float64, params.N())
			ecd.Decode(ptTmp, values)

			c.stats.Update(values)
		}

		c.Finalize()

		if err := w.Write(c.ToCSV()); err != nil {
			panic(err)
		}

		w.Flush()
	}
}
