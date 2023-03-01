package operations

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"time"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils"
)

func GetNoiseLinearTransform(LogN, H, LogSlots, LogScale int, nonZeroDiags map[int]float64, std float64, nbRuns int) {

	Log2BSGSRatio := 2

	//correction := math.Sqrt(float64(int(1 << (LogN - 1))))

	fmt.Printf("H:%d\n", H)

	c := NewContext(LogN, H, LogScale)

	params := c.params

	diags := make(map[int][]complex128)
	for i := range nonZeroDiags {
		diags[i] = make([]complex128, 1<<LogSlots)
	}

	est := estimator.NewEstimator(params.N(), params.HammingWeight(), params.Q(), params.P())

	for i := 0; i < nbRuns; i++ {

		c.GenKeys()

		t := time.Now().UnixMilli()
		r := rand.New(rand.NewSource(t))

		ecd := c.ecd
		enc := ckks.NewEncryptor(params, c.sk)
		dec := c.dec

		for i, diag := range diags {
			//std := nonZeroDiags[i] * correction 
			for j := range diags[i] {

				/*
				real := 2*(r.Float64()-0.5)
				imag := math.Sqrt(1-real*real)

				if r.Float64() >= 0.5{
					imag = -imag
				}
				*/

				diag[j] = complex(2*(r.Float64()-0.5), 2*(r.Float64()-0.5)) / 32.0
			}
		}

		LT := ckks.GenLinearTransformBSGS(ecd, diags, params.MaxLevel(), rlwe.NewScale(params.Q()[params.MaxLevel()]), Log2BSGSRatio, LogSlots)

		rots := LT.Rotations()

		rotKey := c.kgen.GenRotationKeysForRotations(rots, false, c.sk)

		eval := c.eval.WithKey(rlwe.EvaluationKey{Rtks: rotKey})

		values := make([]complex128, 1<<LogSlots)
		for i := range values {
			values[i] = NormComplex(r, std)
		}

		pt := ecd.EncodeNew(values, params.MaxLevel(), params.DefaultScale(), LogSlots)

		// Get the standard deviation of the plaintext in the ring
		stdPT := STDPoly(params.RingQ().AtLevel(pt.Level()), pt.Value, 1, false)

		// Encrypt and evaluate the linear transform
		ct := enc.EncryptNew(pt)
		eval.LinearTransform(ct, LT, []*rlwe.Ciphertext{ct})

		// Subtract the linear transform in the clear
		pt.Scale = ct.Scale
		ecd.Encode(EvaluateLinearTransform(values, diags, Log2BSGSRatio), pt, LogSlots)
		eval.Sub(ct, pt, ct)

		// Decrypt and log STD
		dec.Decrypt(ct, pt)

		// Extract each diagonal standard deviation of the LT
		estDiags := make(map[int][2]float64)
		for diag, vec := range LT.Vec{
			estDiags[diag] = [2]float64{
				STDPoly(params.RingQ(), vec.Q, 1, true),
				math.Sqrt(1/12.0),
			}
		}

		// Evaluate
		estLT := estimator.NewLinearTransform(estDiags, 1, params.MaxLevel(), LogSlots, Log2BSGSRatio)
		estCT := estimator.NewCiphertextSK(estimator.NewPlaintext(stdPT, math.Sqrt(1/12.0), params.MaxLevel()))
		estCT = est.LinearTransform(estCT, estLT)
		fmt.Println(math.Log2(STDPoly(params.RingQ().AtLevel(pt.Level()), pt.Value, 1, false)), est.Std(estCT))
	}
}

func STDPoly(r *ring.Ring, poly *ring.Poly, scale float64, montgomery bool) (std float64) {
	coeffsBig := make([]*big.Int, r.N())
	for i := range coeffsBig {
		coeffsBig[i] = new(big.Int)
	}
	tmp := r.NewPoly()
	r.INTT(poly, tmp)
	if montgomery {
		r.IMForm(tmp, tmp)
	}
	r.PolyToBigintCentered(tmp, 1, coeffsBig)
	values := make([]float64, r.N())
	for i := range coeffsBig {
		values[i], _ = new(big.Float).SetInt(coeffsBig[i]).Float64()
		values[i] /= scale
	}
	return estimator.STD(values)
}

func NormFloat(r *rand.Rand, std float64) (f float64) {

	f = r.NormFloat64() * std

	if f > 6*std{
		f = r.NormFloat64() * std
	}

	return 
}

func NormComplex(r *rand.Rand, std float64) complex128 {
	real := r.NormFloat64()
	for real > 6.0*std {
		real = r.NormFloat64()
	}

	imag := r.NormFloat64()
	for imag > 6.0*std {
		imag = r.NormFloat64()
	}

	return complex(real, imag) * complex(std, 0)
}

func EvaluateLinearTransform(values []complex128, diags map[int][]complex128, LogBSGSRatio int) (res []complex128) {

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

			for i := 0; i < slots; i++ {
				tmp[i] += a[i] * b[i]
			}
		}

		tmp = utils.RotateComplex128Slice(tmp, j)

		for i := 0; i < slots; i++ {
			res[i] += tmp[i]
		}
	}

	return
}
