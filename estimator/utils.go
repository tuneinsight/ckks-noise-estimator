package estimator

import (
	"fmt"
	"math"
	"math/big"
	"math/rand"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func (e Estimator) NewTestVector(ecd *hefloat.Encoder, enc *rlwe.Encryptor, a, b complex128, trunc... uint) (values []*bignum.Complex, el *Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	params := e.Parameters

	prec := ecd.Prec()

	values = make([]*bignum.Complex, params.MaxSlots())

	for i := range values {
		values[i] = &bignum.Complex{
			bignum.NewFloat(sampling.RandFloat64(real(a), real(b)), prec),
			bignum.NewFloat(sampling.RandFloat64(imag(a), imag(b)), prec),
		}
	}

	values[0][0].SetFloat64(1)
	values[0][1].SetFloat64(0)

	el = e.NewElement(values, 1, params.MaxLevel(), params.DefaultScale())
	e.AddEncodingNoise(el)

	pt = hefloat.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(values, pt); err != nil{
		panic(err)
	}

	if len(trunc) != 0{
		truncval := trunc[0]
		coeffs := make([]*big.Int, params.N())
		for i := range coeffs{
			coeffs[i] = new(big.Int)
		}
		params.RingQ().INTT(pt.Value, pt.Value)
		params.RingQ().PolyToBigintCentered(pt.Value, 1, coeffs)
		for i := range coeffs{
			coeffs[i].Rsh(coeffs[i], truncval)
			coeffs[i].Lsh(coeffs[i], truncval)
		}
		params.RingQ().SetCoefficientsBigint(coeffs, pt.Value)
		params.RingQ().NTT(pt.Value, pt.Value)
	}

	if enc != nil{
		ct, _ = enc.EncryptNew(pt)
		switch enc.KeyType().(type) {
		case *rlwe.SecretKey:
			e.AddEncryptionNoiseSk(el)
		case *rlwe.PublicKey:
			e.AddEncryptionNoisePk(el)
		default:
			panic("INVALID ENCRYPION KEY")
		}
	}
	
	return
}

func (e Estimator) TruncatePlaintext(pt *rlwe.Plaintext, trunc *big.Int) (coeffs []*big.Int){
	params := e.Parameters
	coeffs = make([]*big.Int, params.N())
	for i := range coeffs{
		coeffs[i] = new(big.Int)
	}
	rQ := params.RingQ().AtLevel(pt.Level())
	rQ.INTT(pt.Value, pt.Value)
	rQ.PolyToBigintCentered(pt.Value, 1, coeffs)
	for i := range coeffs{
		bignum.DivRound(coeffs[i], trunc, coeffs[i])
		coeffs[i].Mul(coeffs[i], trunc)
	}
	rQ.SetCoefficientsBigint(coeffs, pt.Value)
	rQ.NTT(pt.Value, pt.Value)
	return coeffs
}

func (e Estimator) NewTestVectorFromSeed(ecd *hefloat.Encoder, enc *rlwe.Encryptor, a, b complex128, seed int64, trunc... *big.Int) (values []*bignum.Complex, el *Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	params := e.Parameters

	prec := ecd.Prec()

	values = make([]*bignum.Complex, params.MaxSlots())

	r := rand.New(rand.NewSource(seed))

	Float64 := func(min, max float64) float64 {
		return (max-min) * r.Float64() + min
	}

	for i := range values {
		values[i] = &bignum.Complex{
			bignum.NewFloat(Float64(real(a), real(b)), prec),
			bignum.NewFloat(Float64(imag(a), imag(b)), prec),
		}
	}

	values[0][0].SetFloat64(1)
	values[0][1].SetFloat64(0)

	el = e.NewElement(values, 1, params.MaxLevel(), params.DefaultScale())
	e.AddEncodingNoise(el)

	pt = hefloat.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(values, pt); err != nil{
		panic(err)
	}

	if len(trunc) != 0{
		e.TruncatePlaintext(pt, trunc[0])
	}

	if enc != nil{
		ct, _ = enc.EncryptNew(pt)
		switch enc.KeyType().(type) {
		case *rlwe.SecretKey:
			e.AddEncryptionNoiseSk(el)
		case *rlwe.PublicKey:
			e.AddEncryptionNoisePk(el)
		default:
			panic("INVALID ENCRYPION KEY")
		}
	}
	
	return
}

const (
	prec = uint(128)
)

func Round(x *big.Float){
	x.Add(x, new(big.Float).SetFloat64(0.5*float64(x.Sign())))
	xint, _ := x.Int(nil)
	x.SetInt(xint)
}

func Truncate(x *big.Float, delta *big.Float){
	x.Mul(x, delta)
	Round(x)
	x.Quo(x, delta)
}

func DecompRNS(levelQ, levelP int) int {
	return (levelQ + levelP + 1) / (levelP + 1)
}

func NewFloat(x interface{}) (s *big.Float) {
	switch x := x.(type) {
	case float64:
		if math.IsNaN(x) || math.IsInf(x, 0) {
			panic(fmt.Errorf("x cannot be negative, NaN or Inf, but is %f", x))
		}
		s = new(big.Float).SetPrec(prec)
		s.SetFloat64(x)
		return
	case *big.Int:
		s = new(big.Float).SetPrec(prec)
		s.SetInt(x)
		return
	case *big.Float:
		s = new(big.Float).SetPrec(prec)
		s.Set(x)
		return
	case int:
		return NewFloat(new(big.Int).SetInt64(int64(x)))
	case int64:
		return NewFloat(new(big.Int).SetInt64(x))
	case uint64:
		return NewFloat(new(big.Int).SetUint64(x))
	default:
		panic(fmt.Errorf("invalid x.(type): must be int, int64, uint64, float64, *big.Int or *big.Float but is %T", x))
	}
}

func Log2MIN(x []*bignum.Complex) (real, img float64){

	minR := new(big.Float)
	minI := new(big.Float)
	tmp := new(big.Float)

	for i := range x{

		tmp.Abs(x[i][0])

		if minR.Cmp(tmp) == 1{
			minR.Set(tmp)
		}

		tmp.Abs(x[i][1])

		if minI.Cmp(tmp) == 1{
			minI.Set(tmp)
		}
	}

	minRF64, _ := minR.Float64()
	minIF64, _ := minI.Float64()

	return math.Log2(minRF64), math.Log2(minIF64)
}

func Log2MAX(x []*bignum.Complex) (real, img float64){

	maxR := new(big.Float)
	maxI := new(big.Float)
	tmp := new(big.Float)

	for i := range x{

		tmp.Abs(x[i][0])

		if maxR.Cmp(tmp) == -1{
			maxR.Set(tmp)
		}

		tmp.Abs(x[i][1])

		if maxI.Cmp(tmp) == -1{
			maxI.Set(tmp)
		}
	}

	maxRF64, _ := maxR.Float64()
	maxIF64, _ := maxI.Float64()

	return math.Log2(maxRF64), math.Log2(maxIF64)
}

func Log2AVG(x []*bignum.Complex) (real, img float64){
	n := new(big.Float).SetInt64(int64(len(x)))

	meanR := new(big.Float).SetPrec(x[0].Prec())
	meanI := new(big.Float).SetPrec(x[0].Prec())

	for i := range x{
		meanR.Add(meanR, x[i][0])
		meanI.Add(meanI, x[i][1])
	}

	meanR.Quo(meanR, n)
	meanI.Quo(meanI, n)

	avgRF64, _ := meanR.Float64()
	avgIF64, _ := meanI.Float64()

	return avgRF64, avgIF64
}

func Log2STD(x []*bignum.Complex) (real, img float64){

	n := new(big.Float).SetInt64(int64(len(x)))

	meanR := new(big.Float).SetPrec(x[0].Prec())
	meanI := new(big.Float).SetPrec(x[0].Prec())

	for i := range x{
		meanR.Add(meanR, x[i][0])
		meanI.Add(meanI, x[i][1])
	}

	meanR.Quo(meanR, n)
	meanI.Quo(meanI, n)

	stdR := new(big.Float).SetPrec(x[0].Prec())
	stdI := new(big.Float).SetPrec(x[0].Prec())
	tmp := new(big.Float)

	for i := range x{
		tmp.Sub(x[i][0], meanR)
		tmp.Mul(tmp, tmp)
		stdR.Add(stdR, tmp)

		tmp.Sub(x[i][1], meanI)
		tmp.Mul(tmp, tmp)
		stdI.Add(stdI, tmp)
	}

	stdRF64, _ := stdR.Sqrt(stdR).Float64()
	stdIF64, _ := stdI.Sqrt(stdI).Float64()

	return math.Log2(stdRF64), math.Log2(stdIF64)
}
