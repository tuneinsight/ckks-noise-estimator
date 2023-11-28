package estimator

import (
	"math"

	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

// PowerBasis is a struct storing powers of a ciphertext.
type PowerBasis struct {
	bignum.Basis
	Value map[int]*Element
}

// NewPowerBasis creates a new PowerBasis. It takes as input a ciphertext
// and a basistype. The struct treats the input ciphertext as a monomial X and
// can be used to generates power of this monomial X^{n} in the given BasisType.
func NewPowerBasis(ct Element, basis bignum.Basis) (p *PowerBasis) {
	p = new(PowerBasis)
	p.Value = make(map[int]*Element)
	p.Value[1] = ct.CopyNew()
	p.Basis = basis
	return
}

// GenPower recursively computes X^{n}.
// If lazy = true, the final X^{n} will not be relinearized.
// Previous non-relinearized X^{n} that are required to compute the target X^{n} are automatically relinearized.
// Scale sets the threshold for rescaling (ciphertext won't be rescaled if the rescaling operation would make the scale go under this threshold).
func (p *PowerBasis) GenPower(n int, lazy bool) (err error) {

	if _, ok := p.Value[n]; !ok {
		if p.genPower(n, lazy) {
			p.Value[n].Rescale()
		}
	}

	return nil
}

func (p *PowerBasis) genPower(n int, lazy bool) (rescale bool) {

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

			if p.Basis == bignum.Chebyshev {
				c = int(math.Abs(float64(a) - float64(b))) // Cn = 2*Ca*Cb - Cc, n = a+b and c = abs(a-b)
			}
		}

		// Recurses on the given indexes
		rescaleA := p.genPower(a, lazy && !isPow2)
		rescaleB := p.genPower(b, lazy && !isPow2)

		// Computes C[n] = C[a]*C[b]
		if lazy {

			if p.Value[a].Degree == 2 {
				p.Value[a].Relinearize()
			}

			if p.Value[b].Degree == 2 {
				p.Value[b].Relinearize()
			}

			if rescaleA {
				p.Value[a].Rescale()
			}

			if rescaleB {
				p.Value[b].Rescale()
			}
			p.Value[n] = p.Value[a].CopyNew()
			p.Value[n].Mul(p.Value[a], p.Value[b])

		} else {

			if rescaleA {
				p.Value[a].Rescale()
			}

			if rescaleB {
				p.Value[b].Rescale()
			}

			p.Value[n] = p.Value[a].CopyNew()
			p.Value[n].Mul(p.Value[a], p.Value[b])
			p.Value[n].Relinearize()
		}

		if p.Basis == bignum.Chebyshev {

			// Computes C[n] = 2*C[a]*C[b]
			p.Value[n].Mul(p.Value[n], 2)
				
			// Computes C[n] = 2*C[a]*C[b] - C[c]
			if c == 0 {
				p.Value[n].Add(p.Value[n], -1)
			} else {
				// Since C[0] is not stored (but rather seen as the constant 1), only recurses on c if c!= 0
				p.GenPower(c, lazy)
				p.Value[n].Sub(p.Value[n], p.Value[c])
			}
		}

		return true
	}

	return false
}
