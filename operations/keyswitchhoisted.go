package operations

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/big"
	"os"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
)

func KeySwitchHoisted(LogN, H int, nbRuns int) {

	f, err := os.Create(fmt.Sprintf("data/keyswitchhoisted_%d_%d_%d.csv", LogN, nbRuns, time.Now().Unix()))
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

	fmt.Printf("Noise KeySwitchHoisted\n")

	fmt.Printf("H:%d\n", H)

	c := NewContext(LogN, H, 0)

	params := c.params

	k := 1

	pt := ckks.NewPlaintext(params, params.MaxLevel())
	ct := ckks.NewCiphertext(params, 1, params.MaxLevel())

	coeffsBig := make([]*big.Int, params.N())
	for i := range coeffsBig {
		coeffsBig[i] = new(big.Int)
	}

	for i := 0; i < nbRuns; i++ {

		c.GenKeys()

		rtks := c.kgen.GenRotationKeysForRotations([]int{k}, false, c.sk)

		enc := ckks.NewEncryptor(params, c.pk)
		dec := c.dec
		eval := c.eval

		// Encrypt first plaintext
		enc.EncryptZero(ct)

		// Enc(pt1) x pt2

		eval = eval.WithKey(rlwe.EvaluationKey{Rtks: rtks})

		buffQP := eval.GetRLWEEvaluator().BuffDecompQP
		eval.GetRLWEEvaluator().DecomposeNTT(ct.Level(), params.MaxLevelP(), params.PCount(), ct.Value[1], ct.IsNTT, buffQP)

		ctRots := eval.RotateHoistedLazyNew(ct.Level(), []int{k}, ct.Value[0], buffQP)

		ct2 := &rlwe.Ciphertext{
			Value:    []*ring.Poly{ctRots[k].Value[0].Q, ctRots[k].Value[1].Q},
			MetaData: ct.MetaData,
		}

		// Dec(Enc(pt1) x pt2)
		dec.Decrypt(ct2, pt)

		params.RingQ().INTT(pt.Value, pt.Value)
		params.RingQ().PolyToBigintCentered(pt.Value, 1, coeffsBig)

		values := make([]float64, params.N())
		for i := range values {
			values[i], _ = new(big.Float).SetPrec(128).SetInt(coeffsBig[i]).Float64()
		}

		// Compares in the ring
		c.stats.Update(values)
	}

	c.Finalize()

	if err := w.Write(c.ToCSV()); err != nil {
		panic(err)
	}

	w.Flush()
}

func STD(values []*big.Int) (std float64) {

	mean := new(big.Float).SetPrec(128)

	for i := range values {
		mean.Add(mean, new(big.Float).SetInt(values[i]))
	}

	n := new(big.Float).SetPrec(128).SetFloat64(float64(len(values)))

	mean.Quo(mean, n)

	variance := new(big.Float).SetPrec(128)

	for i := range values {
		tmp := new(big.Float).SetInt(values[i])
		tmp.Sub(tmp, mean)
		tmp.Mul(tmp, tmp)
		variance.Add(variance, tmp)
	}

	variance.Quo(variance, n)
	variance.Sqrt(variance)

	std, _ = variance.Float64()

	return math.Log2(std)
}
