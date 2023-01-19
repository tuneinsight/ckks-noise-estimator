package estimator

import (
	"math"
	"math/big"
	"testing"

	"github.com/stretchr/testify/require"
)

var delta = 1e-15 //float64 precision

func TestEstimator(t *testing.T) {

	LogN := 16

	N := 1<<LogN
	H := 32768

	Q := []uint64{1152921504606584833, 35184372744193, 35184373006337, 35184368025601, 35184376545281}
	P := []uint64{2305843009211596801, 2305843009210023937}

	stdSk := NewFloat(H)
	stdSk.Quo(stdSk, NewFloat(N))
	stdSk.Sqrt(stdSk)

	est := NewEstimator(N, H, Q, P)

	t.Run("Add/Pt/Pt", func(t *testing.T) {

		var msg float64 = 1 << 45
		var err float64 = 3.2

		pt0 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err)},
		}

		pt1 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err)},
		}

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
		var msg float64 = 1 << 45
		var err float64 = 3.2

		ct0 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err), NewFloat(err)},
		}

		pt1 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err)},
		}

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

		var msg float64 = 1 << 45
		var err float64 = 3.2

		ct0 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err), NewFloat(err)},
		}

		ct1 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err), NewFloat(err)},
		}

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

		var msg float64 = 1 << 45
		var err float64 = 3.2

		pt0 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err)},
		}

		pt1 := Element{
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err)},
		}

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

		var msg float64 = 1 << 45
		var err float64 = 3.2

		pt0 := NewPlaintext(err*msg, nil, 4)
		ct1 := NewCiphertextPK(NewPlaintext(err*msg, nil, 4))

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

		var msg float64 = 1 << 45
		var err float64 = 3.2

		ct0 := Element{
			Level:   4,
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err), NewFloat(err)},
		}

		ct1 := Element{
			Level:   4,
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err), NewFloat(err)},
		}

		ct2 := est.Mul(ct0, ct1)

		// This shouldn't affect the noise in a noticeable way because
		// it adds a rounding error and the ciphertext
		// message (and noise) is scaled roughtly by
		// the scaling factor from the multiplication.
		ct2 = est.Relinearize(ct2)

		tmp := MulSTD(est.N, ct0.Message, ct1.Message)

		have, _ := ct2.Message.Float64()
		want, _ := tmp.Float64()

		tmp0 := MulSTD(est.N, ct0.Message, ct1.Noise[0])
		tmp0 = AddSTD(tmp0, MulSTD(est.N, ct1.Message, ct0.Noise[0]))
		tmp0 = AddSTD(tmp0, MulSTD(est.N, ct0.Noise[0], ct0.Noise[1]))

		tmp1 := MulSTD(est.N, ct0.Message, ct1.Noise[1])
		tmp1 = AddSTD(tmp1, MulSTD(est.N, ct0.Noise[0], ct1.Noise[1]))
		tmp1 = AddSTD(tmp1, MulSTD(est.N, ct1.Message, ct0.Noise[1]))
		tmp1 = AddSTD(tmp1, MulSTD(est.N, ct1.Noise[0], ct0.Noise[1]))

		tmp2 := MulSTD(est.N, ct0.Noise[1], ct1.Noise[1])

		tmp = MulSTD(est.N, tmp2, stdSk)
		tmp = AddSTD(tmp, tmp1)
		tmp = MulSTD(est.N, tmp, stdSk)
		tmp = AddSTD(tmp, tmp0)

		have = est.Std(ct2)

		want, _ = tmp.Float64()
		want = math.Log2(want)

		require.InDelta(t, want, have, delta)
	})

	t.Run("Rescale", func(t *testing.T) {

		// EMPIRICALLY VERIFIED

		level := 1

		var msg float64 = 1 << 45
		var err float64 = 1 << 46
		var q = float64(Q[level])

		elem := Element{
			Level:   level,
			Message: NewFloat(msg),
			Noise:   []*big.Float{NewFloat(err), NewFloat(err), NewFloat(err)},
		}

		elem = est.Rescale(elem)

		have := est.Std(elem)
		want := math.Log2(math.Sqrt((err/q*err/q + 1/12.0) + ((err/q*err/q+1/12.0)+(err/q*err/q+1/12.0)*float64(H))*float64(H)*0.5))

		//require.InDelta(t, elem.Message, msg/q, delta)
		require.Equal(t, elem.Level, level-1)
		require.InDelta(t, want, have, delta)
	})

	t.Run("LinearTransform", func(t *testing.T){

		diags := map[int]float64{
			-5: 3.2,
			-4: 3.2,
			-3: 3.2,
			-2: 3.2,
			-1: 3.2,
			0: 3.2,
			1: 3.2,
			2: 3.2,
			3: 3.2,
			4: 3.2,
			5: 3.2,
		}

		level := 4

		LT := NewLinearTransform(diags, Q[level], level, LogN-1, 2)

		ct := NewCiphertextPK(NewPlaintext(3.2 * (1<<45), nil, level))

		ct = est.LinearTransform(ct, LT)

		have := est.Std(ct)

		t.Log(have)

	})
}
