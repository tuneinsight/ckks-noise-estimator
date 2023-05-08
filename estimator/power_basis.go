package estimator

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v4/utils/bignum/polynomial"
)

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	polynomial.Basis
	Value map[int]Element
}

// NewPowerBasis creates a new PowerBasis. It takes as input a ciphertext
// and a basistype. The struct treats the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPowerBasis(ct Element, basis polynomial.Basis) (p *PowerBasis) {
	p = new(PowerBasis)
	p.Value = make(map[int]Element)
	p.Value[1] = ct.CopyNew()
	p.Basis = basis
	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PowerBasis) GenPower(n int, lazy bool, eval *Estimator) (err error) {

	if _, ok := p.Value[n]; !ok {

		var rescale bool
		if rescale, err = p.genPower(n, lazy, eval); err != nil {
			return
		}

		fmt.Println(n, eval.Std(p.Value[n]))

		if rescale {
			if p.Value[n], err = eval.Rescale(p.Value[n]); err != nil {
				return
			}

			fmt.Println(n, eval.Std(p.Value[n]))
		}
	}

	return nil
}

func (p *PowerBasis) genPower(n int, lazy bool, eval *Estimator) (rescale bool, err error) {

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

		var rescaleA, rescaleB bool

		// Recurses on the given indexes
		if rescaleA, err = p.genPower(a, lazy && !isPow2, eval); err != nil {
			return false, err
		}
		if rescaleB, err = p.genPower(b, lazy && !isPow2, eval); err != nil {
			return false, err
		}

		// Computes C[n] = C[a]*C[b]
		if lazy {

			if p.Value[a].Degree() == 2 {
				p.Value[a] = eval.Relinearize(p.Value[a])
			}

			if p.Value[b].Degree() == 2 {
				p.Value[b] = eval.Relinearize(p.Value[b])
			}

			if rescaleA {
				if p.Value[a], err = eval.Rescale(p.Value[a]); err != nil {
					return false, err
				}
			}

			if rescaleB {
				if p.Value[b], err = eval.Rescale(p.Value[b]); err != nil {
					return false, err
				}
			}

			p.Value[n] = eval.Mul(p.Value[a], p.Value[b])

		} else {

			if rescaleA {
				if p.Value[a], err = eval.Rescale(p.Value[a]); err != nil {
					return false, err
				}
			}

			if rescaleB {
				if p.Value[b], err = eval.Rescale(p.Value[b]); err != nil {
					return false, err
				}
			}

			if a == b{
				p.Value[n] = eval.Square(p.Value[a])
			}else{
				p.Value[n] = eval.Mul(p.Value[a], p.Value[b])
			}

			p.Value[n] = eval.Relinearize(p.Value[n])
		}

		

		if p.Basis == polynomial.Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			p.Value[n] = eval.Mul(p.Value[n], 2)
			p.Value[n].Message.Quo(p.Value[n].Message, new(big.Float).SetFloat64(1.63))

			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				//p.Value[n] = eval.Add(p.Value[n], -1) // Constant addition = no noise (well in fact (1/12, 0), so = negligible)

				//

			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				if err = p.GenPower(c, lazy, eval); err != nil {
					return false, err
				}
				p.Value[n] = eval.Add(p.Value[n], p.Value[c])
			}
		}

		return true, nil
	}

	return false, nil
}
