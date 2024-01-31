package estimator

import (
	"fmt"
	"math"
	"math/big"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

func (p Parameters) NewTestVector(ecd *hefloat.Encoder, enc *rlwe.Encryptor, a, b complex128) (values []*bignum.Complex, el *Element, pt *rlwe.Plaintext, ct *rlwe.Ciphertext) {

	params := p.Parameters

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

	el = NewElement(p, values, 1, params.DefaultScale())
	el.AddEncodingNoise()

	pt = hefloat.NewPlaintext(params, params.MaxLevel())
	if err := ecd.Encode(values, pt); err != nil{
		panic(err)
	}

	if enc != nil{
		ct, _ = enc.EncryptNew(pt)
		switch enc.KeyType().(type) {
		case *rlwe.SecretKey:
			el.AddEncryptionNoiseSk()
		case *rlwe.PublicKey:
			el.AddEncryptionNoisePk()
		default:
			panic("INVALID ENCRYPION KEY")
		}
	}
	
	return
}

const (
	prec = uint(128)
)

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
