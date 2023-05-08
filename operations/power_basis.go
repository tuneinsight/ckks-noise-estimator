package operations

import (
	"fmt"
	"math"
	"math/big"
	"github.com/fatih/color"

	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/rlwe"
	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
)

func GetNoisePowerBasis(LogN, H, LogSlots, LogScale int, nbRuns int) {

	fmt.Printf("Noise Chebyshev Power Basis\n")

	fmt.Printf("H:%d\n", H)

	c := NewContext(LogN, H, LogScale)

	params := c.params

	ct := ckks.NewCiphertext(params, 1, params.MaxLevel())

	est := estimator.NewEstimator(params.N(), params.XsHammingWeight(), params.Q(), params.P())

	ecd := ckks.NewEncoder(c.params, 128)
	ptRef := ckks.NewPlaintext(c.params, c.params.MaxLevel())
	ptRef.LogSlots = LogSlots

	n := 128

	for i := 0; i < nbRuns; i++ {

		c.GenKeys()

		enc := c.enc
		dec := c.dec
		eval := c.eval

		evk := rlwe.NewEvaluationKeySet()
		evk.RelinearizationKey = c.kgen.GenRelinearizationKeyNew(c.sk)

		eval = eval.WithKey(evk)

		// First plaintext
		pt1, values, stdpt1 := c.NewPlaintextSlots(-1, 1, LogSlots)

		// Encrypt first plaintext
		enc.Encrypt(pt1, ct)

		pb := ckks.NewPowerBasis(ct, polynomial.Chebyshev)

		pbPt := NewPowerBasisPlaintext(values, polynomial.Chebyshev)

		estCT := estimator.NewCiphertextPK(estimator.NewPlaintext(stdpt1 * params.DefaultScale().Float64(), 4/12.0, params.MaxLevel()))
		
		pbEst := estimator.NewPowerBasis(estCT, polynomial.Chebyshev)

		var err error
		for i := 1; i <= n; i<<=1 {

			pbPt.GenPower(i)
			if err = pb.GenPower(i, false, params.DefaultScale(), eval); err != nil {
				panic(err)
			}

			if err = pbEst.GenPower(i, false, est); err != nil {
				panic(err)
			}
		}

		for i := 1; i <= n; i<<=1 {
			if p, ok := pb.Value[i]; ok {

				ptRef.Scale = p.Scale
				ecd.Encode(pbPt.Value[i], ptRef)

				eval.Sub(p, ptRef, p)

				dec.Decrypt(p, pt1)

				toprint := fmt.Sprintf("%3d - %d - %14.12f - %14.12f\n", i, p.Level(), math.Log2(STDPoly(params.RingQ().AtLevel(pt1.Level()), pt1.Value, 1, false)),est.Std(pbEst.Value[i]))
				if i & (i-1) == 0{
					color.Red(toprint)
				}else{
					fmt.Printf("%s", toprint)
				}
			}
		}
	}
	fmt.Println()

}


type PowerBasisPlaintext struct {
	Basis polynomial.Basis
	Value map[int][]*big.Float
}

func NewPowerBasisPlaintext(v []*big.Float, basis polynomial.Basis) (p *PowerBasisPlaintext) {
	p = new(PowerBasisPlaintext)
	p.Value = make(map[int][]*big.Float)
	p.Value[1] = v
	p.Basis = basis
	return
}

func (p *PowerBasisPlaintext) GenPower(n int) {
	if _, ok := p.Value[n]; !ok {
		p.genPower(n)
	}
}

func (p *PowerBasisPlaintext) genPower(n int) {

	if _, ok := p.Value[n]; !ok {

		isPow2 := n&(n-1) == 0

		// Computes the index required to compute the asked ring evaluation
		var a, b, c int
		if isPow2 {
			a, b = n/2, n/2 //Necessary for optimal depth
		} else {
			// [Lee et al. 2020] : High-Precision and Low-Complexity Approximate Homomorphic Encryption by Error Variance Minimization
			// Maximize the number of odd terms of Chebyshev basis
			k := int(math.Ceil(math.Log2(float64(n)))) - 1
			a = (1 << k) - 1
			b = n + 1 - (1 << k)

			if p.Basis == polynomial.Chebyshev {
				c = int(math.Abs(float64(a) - float64(b))) // Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			}
		}


		// Recurses on the given indexes
		p.genPower(a)
		p.genPower(b)

		// Computes C[n] = C[a]*C[b]
		p.Value[n] = mul(p.Value[a], p.Value[b])

		if p.Basis == polynomial.Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			add(p.Value[n], p.Value[n], p.Value[n])

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				addscalar(p.Value[n], new(big.Float).SetPrec(128).SetInt64(-1), p.Value[n])
			} else {
				p.GenPower(c)
				sub(p.Value[n], p.Value[c], p.Value[n])
			}
		}
	}
}

func mul(a, b []*big.Float) (c []*big.Float) {
	c = make([]*big.Float, len(a))

	for i := range a{
		c[i] = new(big.Float).Mul(a[i], b[i])
	}

	return
}

func add(a, b, c []*big.Float){
	for i := range a{
		c[i].Add(a[i], b[i])
	}
}

func sub(a, b, c []*big.Float){
	for i := range a{
		c[i].Sub(a[i], b[i])
	}
}

func addscalar(a  []*big.Float, b *big.Float, c []*big.Float){
	for i := range a{
		c[i].Add(a[i], b)
	}
}