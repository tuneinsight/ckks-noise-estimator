package estimator

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

var delta = 1e-15 //float64 precision

const (
	stdUni = 0.5772058896878792 // Standard deviation of uniform vector [-1, 1]
)

func TestEstimator(t *testing.T) {

	LogN := 16

	N := 1 << LogN
	H := 32768

	sqrtN := math.Sqrt(float64(N))

	ptStd := stdUni / (sqrtN / 1.73248500) * float64(int(1<<45)) // Standard deviation of encoded plaintext with uniform distribution [-1, 1] with 2^{45} scale

	Q := []uint64{1152921504606584833, 35184372744193, 35184373006337, 35184368025601, 35184376545281}
	P := []uint64{2305843009211596801, 2305843009210023937} //,

	level := 4

	stdSk := NewFloat(H)
	stdSk.Quo(stdSk, NewFloat(N))
	stdSk.Sqrt(stdSk)

	est := NewEstimator(N, H, Q, P)

	t.Run("Add/Pt/Pt", func(t *testing.T) {

		pt0 := NewPlaintext(ptStd, math.Sqrt(1/12.0), level)
		pt1 := NewPlaintext(ptStd, math.Sqrt(1/12.0), level)

		pt2 := est.Add(pt0, pt1)

		tmp := AddSTD(pt0.Message, pt1.Message)

		have, _ := pt2.Message.Float64()
		want, _ := tmp.Float64()

		require.InDelta(t, want, have, delta)

		have = est.Std(pt2)
		want, _ = AddSTD(pt0.Noise[0], pt1.Noise[0]).Float64()
		want = math.Log2(want)

		require.InDelta(t, want, have, delta)
	})

	t.Run("Add/Ct/Pt", func(t *testing.T) {

		ct0 := NewCiphertextPK(NewPlaintext(ptStd, math.Sqrt(1/12.0), level))
		pt1 := NewPlaintext(ptStd, math.Sqrt(1/12.0), level)

		ct1 := est.Add(ct0, pt1)

		tmp := AddSTD(ct0.Message, pt1.Message)

		have, _ := ct1.Message.Float64()
		want, _ := tmp.Float64()

		require.InDelta(t, want, have, delta)

		have = est.Std(ct1)

		tmp = MulSTD(est.N, ct1.Noise[1], stdSk)
		tmp = AddSTD(tmp, ct1.Noise[0])

		want, _ = tmp.Float64()
		want = math.Log2(want)

		require.InDelta(t, want, have, delta)
	})

	t.Run("Add/Ct/Ct", func(t *testing.T) {

		ct0 := NewCiphertextPK(NewPlaintext(ptStd, math.Sqrt(1/12.0), level))
		ct1 := NewCiphertextPK(NewPlaintext(ptStd, math.Sqrt(1/12.0), level))

		ct2 := est.Add(ct0, ct1)

		tmp := AddSTD(ct0.Message, ct1.Message)

		have, _ := ct2.Message.Float64()
		want, _ := tmp.Float64()

		require.InDelta(t, want, have, delta)

		have = est.Std(ct2)

		tmp = MulSTD(est.N, ct2.Noise[1], stdSk)
		tmp = AddSTD(tmp, ct2.Noise[0])

		want, _ = tmp.Float64()
		want = math.Log2(want)

		require.InDelta(t, want, have, delta)
	})

	t.Run("Mul/Pt/Pt", func(t *testing.T) {

		pt0 := NewPlaintext(ptStd, math.Sqrt(1/12.0), level)
		pt1 := NewPlaintext(ptStd, math.Sqrt(1/12.0), level)

		pt2 := est.Mul(pt0, pt1)

		tmp := MulSTD(est.N, pt0.Message, pt1.Message)

		have, _ := pt2.Message.Float64()
		want, _ := tmp.Float64()

		tmp = MulSTD(est.N, pt0.Noise[0], pt1.Message)
		tmp = AddSTD(tmp, MulSTD(est.N, pt0.Message, pt1.Noise[0]))

		have = est.Std(pt2)

		want, _ = tmp.Float64()
		want = math.Log2(want)

		require.InDelta(t, want, have, delta)
	})

	t.Run("Mul/Ct/Pt", func(t *testing.T) {

		// EMPIRICALLY VERIFIED

		pt0 := NewPlaintext(ptStd, math.Sqrt(1/12.0), level)
		ct1 := NewCiphertextPK(NewPlaintext(ptStd, math.Sqrt(1/12.0), level))

		ct2 := est.Mul(pt0, ct1)

		tmp := MulSTD(est.N, pt0.Message, ct1.Message)

		have, _ := ct2.Message.Float64()
		want, _ := tmp.Float64()

		tmp0 := MulSTD(est.N, pt0.Noise[0], ct1.Message)
		tmp0 = AddSTD(tmp0, MulSTD(est.N, pt0.Message, ct1.Noise[0]))

		tmp1 := MulSTD(est.N, pt0.Noise[0], ct1.Noise[1])
		tmp1 = AddSTD(tmp1, MulSTD(est.N, pt0.Message, ct1.Noise[1]))

		tmp = MulSTD(est.N, tmp1, stdSk)
		tmp = AddSTD(tmp0, tmp)

		have = est.Std(ct2)

		want, _ = tmp.Float64()
		want = math.Log2(want)

		t.Log(have)
		t.Log(want)

		require.InDelta(t, want, have, delta)
	})

	t.Run("Mul/Ct/Ct", func(t *testing.T) {

		// EMPIRICALLY VERIFIED
		ct := NewCiphertextPK(NewPlaintext(ptStd, 0, level))

		t.Log(ptStd)

		t.Log(est.Std(ct), ct.Message)

		var err error

		for i := 0; i < level; i++ {
			ct = est.Mul(ct, NewPlaintext(ptStd, 0, level))
			ct = est.Relinearize(ct)
			if ct, err = est.Rescale(ct); err != nil {
				t.Fatal(err)
			}
			t.Log(est.Std(ct), ct.Message)
		}

	})

	t.Run("Rescale", func(t *testing.T) {

		// EMPIRICALLY VERIFIED

		level := 1

		var err error

		var msg float64 = 1 << 45
		var e float64 = 1 << 46
		var q = float64(Q[level])

		elem := Element{
			Level:   level,
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(e), NewFloat(e), NewFloat(e)},
		}

		if elem, err = est.Rescale(elem); err != nil {
			t.Fatal(err)
		}

		have := est.Std(elem)
		want := math.Log2(math.Sqrt((e/q*e/q + 1/12.0) + ((e/q*e/q+1/12.0)+(e/q*e/q+1/12.0)*float64(H))*float64(H)))

		//require.InDelta(t, elem.Message, msg/q, delta)
		require.Equal(t, elem.Level, level-1)
		require.InDelta(t, want, have, delta)
	})

	t.Run("LinearTransform", func(t *testing.T) {

		// EMPIRICALLY VERIFIED (UP TO THERE) for CT-SK & CT-PK

		diags := []map[int][2]float64{
			map[int][2]float64{
				0:  {0.00390626732046525, 0},
				1:  {0.00365389161348000, 0},
				2:  {0.00338292536788323, 0},
				3:  {0.00308802841010419, 0},
				4:  {0.00276214866848703, 0},
				5:  {0.00239205905698830, 0},
				6:  {0.00195307507231700, 0},
				7:  {0.00138094268654370, 0},
				-7: {0.00138107715173799, 0},
				-6: {0.00195313989560961, 0},
				-5: {0.00239209373419969, 0},
				-4: {0.00276215569281487, 0},
				-3: {0.00308818469796449, 0},
				-2: {0.00338291628098494, 0},
				-1: {0.00365397976622493, 0},
			},
			map[int][2]float64{
				0:     {0.00024414102316630370, 0},
				2048:  {0.00024414233423388263, 0},
				4096:  {0.00024414206118669165, 0},
				6144:  {0.00024414240022099044, 0},
				8192:  {0.00024414248708016384, 0},
				10240: {0.00024414247911935246, 0},
				12288: {0.00024414180919533906, 0},
				14336: {0.00024414248329247690, 0},
				16384: {0.00024414246898360015, 0},
				18432: {0.00024414239107164240, 0},
				20480: {0.00024414226649732900, 0},
				22528: {0.00024414233600733267, 0},
				24576: {0.00024414232259683394, 0},
				26624: {0.00024414248311350925, 0},
				28672: {0.00024414158296428044, 0},
				30720: {0.00024414248708458494, 0},
			},
		}

		level := 4

		LT := NewLinearTransform(diags[0], Q[level], level, LogN-1, 2)

		ct := NewCiphertextSK(NewPlaintext(ptStd, nil, level))

		t.Log(est.Std(ct))

		ct = est.LinearTransform(ct, LT)

		t.Log(est.Std(ct))

	})

	t.Run("KeySwitchHoisted", func(t *testing.T) {

		// EMPIRICALLY VERIFIED FOR CT-SK & CT-PK

		ct := NewCiphertextPK(NewPlaintext(nil, nil, level))

		ct = est.KeySwitchLazy(ct)

		have := est.Std(ct)

		t.Log(have)

	})

	t.Run("PowerBasis", func(t *testing.T) {

		ct := NewCiphertextPK(NewPlaintext(ptStd, 0, level))

		pb := NewPowerBasis(ct, polynomial.Chebyshev)

		var err error
		for i := 2; i < 32; i <<= 1 {
			if err = pb.GenPower(i, false, est); err != nil {
				t.Fatal(err)
			}
		}

		for i := 3; i < 8; i++ {
			if err = pb.GenPower(i, false, est); err != nil {
				t.Fatal(err)
			}
		}

		for i := 0; i < 64; i++ {
			if p, ok := pb.Value[i]; ok {
				t.Log(i, est.Std(p), p.Level)
			}
		}
	})
}
