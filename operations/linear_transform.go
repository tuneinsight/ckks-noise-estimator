package operations

import(
	"os"
	"fmt"
	"encoding/csv"
	"time"
	"math"
	"math/rand"
	"math/big"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/ckks-bootstrapping-precision/stats"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
)

func GetNoiseLinearTransform(LogN, H, LogScale int, nonZeroDiags map[int]float64, std float64, nbRuns int){

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
	for i := range nonZeroDiags{
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
			std := nonZeroDiags[i] * correction
			for j := range diags[i]{
				diag[j] = NormComplex(r, std)
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

		ct2 := rlwe.NewCiphertextAtLevelFromPoly(ct.Level(), ct.Value)
		ct2.MetaData = ct.MetaData

		//eval.Rescale(ct2, params.DefaultScale(), ct2)

		want := EvaluateLinearTransform(values, diags, 2)

		pt.Scale = ct2.Scale
		ecd.Encode(want, pt, params.LogSlots())
		eval.Sub(ct2, pt, ct2)

		dec.Decrypt(ct2, pt)
		fmt.Println(STDPoly(params.RingQ().AtLevel(pt.Level()), pt.Value, 1, false))
	}
}

func STDPoly(r *ring.Ring, poly *ring.Poly, scale float64, montgomery bool) (std float64){
	coeffsBig := make([]*big.Int, r.N())
	for i := range coeffsBig{
		coeffsBig[i] = new(big.Int)
	}
	tmp := r.NewPoly()
	r.INTT(poly, tmp)
	if montgomery{
		r.IMForm(tmp, tmp)
	}
	r.PolyToBigintCentered(tmp, 1, coeffsBig)
	values := make([]float64, r.N())
	for i := range coeffsBig{
		values[i] = float64(coeffsBig[i].Int64()) / scale
	}
	return math.Log2(estimator.STD(values))
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
