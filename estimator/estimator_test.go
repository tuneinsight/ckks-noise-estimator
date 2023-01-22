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
	P := []uint64{2305843009211596801, 2305843009210023937} //, 

	stdSk := NewFloat(H)
	stdSk.Quo(stdSk, NewFloat(N))
	stdSk.Sqrt(stdSk)

	est := NewEstimator(N, H, Q, P)

	t.Run("LinearTransform", func(t *testing.T){

		// EMPIRICALLY VERIFIED (UP TO THERE) for CT-SK & CT-PK

		diags := []map[int]float64{
			map[int]float64{
				 0: 0.003906267320465252,
				 1: 0.0036538916134800097,
				 2: 0.003382925367883238,
				 3: 0.0030880284101041902,
				 4: 0.002762148668487031,
				 5: 0.002392059056988301,
				 6: 0.0019530750723170012,
				 7: 0.0013809426865437015,
				-7: 0.0013810771517379955,
				-6: 0.0019531398956096173,
				-5: 0.0023920937341996994,
				-4: 0.002762155692814875,
				-3: 0.0030881846979644947,
				-2: 0.003382916280984945,
				-1: 0.00365397976622493,
			},
				map[int]float64{
				    0:0.00024414102316630370,
				 2048:0.00024414233423388263,
				 4096:0.00024414206118669165,
				 6144:0.00024414240022099044, 
				 8192:0.00024414248708016384, 
				10240:0.00024414247911935246, 
				12288:0.00024414180919533906, 
				14336:0.00024414248329247690, 
				16384:0.00024414246898360015, 
				18432:0.00024414239107164240, 
				20480:0.00024414226649732900, 
				22528:0.00024414233600733267, 
				24576:0.00024414232259683394, 
				26624:0.00024414248311350925, 
				28672:0.00024414158296428044, 
				30720:0.00024414248708458494,
			},
		}

		level := 4

		LT := NewLinearTransform(diags[0], Q[level], level, LogN-1, 2)

		ct := NewCiphertextSK(NewPlaintext(0.5772058896878792 * (1<<45), nil, level))

		t.Log(est.Std(ct))

		ct = est.LinearTransform(ct, LT)

		t.Log(est.Std(ct))

	})

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
		want := math.Log2(math.Sqrt((err/q*err/q + 1/12.0) + ((err/q*err/q+1/12.0)+(err/q*err/q+1/12.0)*float64(H))*float64(H)))

		//require.InDelta(t, elem.Message, msg/q, delta)
		require.Equal(t, elem.Level, level-1)
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

	t.Run("KeySwitchHoisted", func(t *testing.T){

		// EMPIRICALLY VERIFIED FOR CT-SK & CT-PK

		ct := NewCiphertextPK(NewPlaintext(nil, nil, 4))

		ct = est.KeySwitchLazy(ct)

		have := est.Std(ct)

		t.Log(have)

	})
}
