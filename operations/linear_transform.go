package operations

import(
	"os"
	"fmt"
	"encoding/csv"
	"time"
	"math"
	"math/rand"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
)

func GetNoiseLinearTransform(LogN, H, LogScale int, nonZeroDiags []int, std float64, nbRuns int){

	f, err := os.Create(fmt.Sprintf("data/linear_transform_%d_%d_%d_%f_%d_%d.csv", LogN, H, LogScale, std, nbRuns, time.Now().Unix()))
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

	fmt.Printf("Noise LinearTransform\n")

	correction := math.Sqrt(float64(int(1<<(LogN-1))))

	fmt.Printf("H:%d\n", H)

	c := NewContext(LogN, H, LogScale)

	params := c.params

	diags := make(map[int][]complex128)
	for _, i := range nonZeroDiags{
		diags[i] = make([]complex128, params.Slots())
	}

	for i := 0; i < nbRuns; i++ {

		c.GenKeys()

		t := time.Now().UnixMilli()
		r := rand.New(rand.NewSource(t))

		ecd := c.ecd
		enc := ckks.NewEncryptor(params, c.sk)
		dec := c.dec

		for i, diag := range diags{
			for j := range diags[i]{
				diag[j] = NormComplex(r, std * correction)
			}
		}

		LT := ckks.GenLinearTransformBSGS(ecd, diags, params.MaxLevel(), rlwe.NewScale(params.Q()[params.MaxLevel()]), 2, params.LogSlots())

		rots := LT.Rotations()

		rotKey := c.kgen.GenRotationKeysForRotations(rots, false, c.sk)

		eval := c.eval.WithKey(rlwe.EvaluationKey{Rtks: rotKey})
		
		values := make([]complex128, params.Slots())
		for i := range values{
			values[i] = NormComplex(r, std * correction)
		}

		pt := ecd.EncodeNew(values, params.MaxLevel(), params.DefaultScale(), params.LogSlots())

		ct := enc.EncryptNew(pt)

		eval.LinearTransform(ct, LT, []*rlwe.Ciphertext{ct})

		want := EvaluateLinearTransform(values, diags, 2)

		pt.Scale = ct.Scale
		ecd.Encode(want, pt, params.LogSlots())

		eval.Sub(ct, pt, ct)

		ct.Scale = rlwe.NewScale(1)
		dec.Decrypt(ct, pt)

		have := ecd.DecodeCoeffs(pt)

		fmt.Println(have[:4])

		// Compares in the ring
		c.stats.Update(have)
	}

	c.Finalize()

	if err := w.Write(c.ToCSV()); err != nil {
		panic(err)
	}

	w.Flush()
}

func NormComplex(r *rand.Rand, std float64) complex128 {
	real := r.NormFloat64() 
	for real > 6.0 {
		real = r.NormFloat64() 
	}

	imag := r.NormFloat64() 
	for imag > 6.0 {
		imag = r.NormFloat64() 
	}

	return complex(real, imag) * complex(std, 0)
}

func EvaluateLinearTransform(values []complex128, diags map[int][]complex128, LogBSGSRatio int) (res []complex128){

	slots := len(values)

	N1 := ckks.FindBestBSGSRatio(diags, len(values), LogBSGSRatio)

	index, _, _ := ckks.BSGSIndex(diags, slots, N1)

	res = make([]complex128, slots)

	for j := range index {

		rot := -j & (slots - 1)

		tmp := make([]complex128, slots)

		for _, i := range index[j] {

			v, ok := diags[j+i]
			if !ok {
				v = diags[j+i-slots]
			}

			a := utils.RotateComplex128Slice(values, i)

			b := utils.RotateComplex128Slice(v, rot)

			for i := 0; i < slots; i++{
				tmp[i] += a[i] * b[i]
			}
		}

		tmp = utils.RotateComplex128Slice(tmp, j) 

		for i := 0; i < slots; i++{
			res[i] += tmp[i]
		}
	}

	return
}
