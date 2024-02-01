package estimator

import (
	"fmt"
	"math/bits"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

func (p *Element) EvaluatePolynomial(poly interface{}, targetScale rlwe.Scale) (opOut *Element, err error) {

	var polyVec he.PolynomialVector
	switch poly := poly.(type) {
	case hefloat.Polynomial:
		polyVec = he.PolynomialVector{Value: []he.Polynomial{he.Polynomial(poly)}}
	case hefloat.PolynomialVector:
		polyVec = he.PolynomialVector(poly)
	case bignum.Polynomial:
		polyVec = he.PolynomialVector{Value: []he.Polynomial{{Polynomial: poly, MaxDeg: poly.Degree(), Lead: true, Lazy: false}}}
	case he.Polynomial:
		polyVec = he.PolynomialVector{Value: []he.Polynomial{poly}}
	case he.PolynomialVector:
		polyVec = poly
	default:
		return nil, fmt.Errorf("cannot Polynomial: invalid polynomial type, must be either bignum.Polynomial, he.Polynomial or he.PolynomialVector, but is %T", poly)
	}

	powerbasis := NewPowerBasis(p, polyVec.Value[0].Basis)

	logDegree := bits.Len64(uint64(polyVec.Value[0].Degree()))
	logSplit := bignum.OptimalSplit(logDegree)

	var odd, even = false, false
	for _, p := range polyVec.Value {
		odd, even = odd || p.IsOdd, even || p.IsEven
	}

	// Computes all the powers of two with relinearization
	// This will recursively compute and store all powers of two up to 2^logDegree
	if err = powerbasis.GenPower(1<<(logDegree-1), false); err != nil {
		return nil, err
	}

	// Computes the intermediate powers, starting from the largest, without relinearization if possible
	for i := (1 << logSplit) - 1; i > 2; i-- {
		if !(even || odd) || (i&1 == 0 && even) || (i&1 == 1 && odd) {
			if err = powerbasis.GenPower(i, polyVec.Value[0].Lazy); err != nil {
				return nil, err
			}
		}
	}

	PS := polyVec.GetPatersonStockmeyerPolynomial(*p.Parameters.Parameters.GetRLWEParameters(), powerbasis.Value[1].Level, powerbasis.Value[1].Scale, targetScale, &SimEvaluator{p.Parameters.Parameters, p.Parameters.Parameters.LevelsConsumedPerRescaling()})

	coeffGetter := he.CoefficientGetter[*bignum.Complex](&hefloat.CoefficientGetter{Values: make([]*bignum.Complex, p.Parameters.Parameters.MaxSlots())})

	if opOut, err = EvaluatePatersonStockmeyerPolynomialVector(PS, coeffGetter, *powerbasis); err != nil {
		return nil, err
	}

	return opOut, err
}

// BabyStep is a struct storing the result of a baby-step
// of the Paterson-Stockmeyer polynomial evaluation algorithm.
type BabyStep struct {
	Degree int
	Value  *Element
}

// EvaluatePatersonStockmeyerPolynomialVector evaluates a pre-decomposed PatersonStockmeyerPolynomialVector on a pre-computed power basis [1, X^{1}, X^{2}, ..., X^{2^{n}}, X^{2^{n+1}}, ..., X^{2^{m}}]
func EvaluatePatersonStockmeyerPolynomialVector[T any](poly he.PatersonStockmeyerPolynomialVector, cg he.CoefficientGetter[T], pb PowerBasis) (res *Element, err error) {

	split := len(poly.Value[0].Value)

	babySteps := make([]*BabyStep, split)

	// Small steps
	for i := range babySteps {

		// eval & cg are not thread-safe
		if babySteps[split-i-1], err = EvaluateBabyStep[T](i, poly, cg, pb); err != nil {
			return nil, fmt.Errorf("cannot EvaluateBabyStep: %w", err)
		}
	}

	// Loops as long as there is more than one sub-polynomial
	for len(babySteps) != 1 {

		// Precomputes the ops to apply in the giant steps loop
		giantsteps := make([]int, len(babySteps))
		for i := 0; i < len(babySteps); i++ {
			if i == len(babySteps)-1 {
				giantsteps[i] = 2
			} else if babySteps[i].Degree == babySteps[i+1].Degree {
				giantsteps[i] = 1
				i++
			}
		}

		for i := 0; i < len(babySteps); i++ {

			// eval is not thread-safe
			if err = EvaluateGianStep(i, giantsteps, babySteps, pb); err != nil {
				return nil, err
			}
		}

		// Discards processed sub-polynomials
		var idx int
		for i := range babySteps {
			if babySteps[i] != nil {
				babySteps[idx] = babySteps[i]
				idx++
			}
		}

		babySteps = babySteps[:idx]
	}

	if babySteps[0].Value.Degree == 2 {
		babySteps[0].Value.Relinearize()
	}

	babySteps[0].Value.Rescale()

	return babySteps[0].Value, nil
}

// EvaluateBabyStep evaluates a baby-step of the PatersonStockmeyer polynomial evaluation algorithm, i.e. the inner-product between the precomputed
// powers [1, T, T^2, ..., T^{n-1}] and the coefficients [ci0, ci1, ci2, ..., ci{n-1}].
func EvaluateBabyStep[T any](i int, poly he.PatersonStockmeyerPolynomialVector, cg he.CoefficientGetter[T], pb PowerBasis) (ct *BabyStep, err error) {

	nbPoly := len(poly.Value)

	polyVec := he.PolynomialVector{
		Value:   make([]he.Polynomial, nbPoly),
		Mapping: poly.Mapping,
	}

	// Transposes the polynomial matrix
	for j := 0; j < nbPoly; j++ {
		polyVec.Value[j] = poly.Value[j].Value[i]
	}

	level := poly.Value[0].Value[i].Level
	scale := poly.Value[0].Value[i].Scale

	ct = new(BabyStep)
	ct.Degree = poly.Value[0].Value[i].Degree()
	if ct.Value, err = EvaluatePolynomialVectorFromPowerBasis(level, polyVec, cg, pb, scale); err != nil {
		return ct, fmt.Errorf("cannot EvaluatePolynomialVectorFromPowerBasis: polynomial[%d]: %w", i, err)
	}

	return ct, nil
}

// EvaluateGianStep evaluates a giant-step of the PatersonStockmeyer polynomial evaluation algorithm, which consists
// in combining the baby-steps <[1, T, T^2, ..., T^{n-1}], [ci0, ci1, ci2, ..., ci{n-1}]> together with powers T^{2^k}.
func EvaluateGianStep(i int, giantSteps []int, babySteps []*BabyStep, pb PowerBasis) (err error) {

	// If we reach the end of the list it means we weren't able to combine
	// the last two sub-polynomials which necessarily implies that that the
	// last one has degree smaller than the previous one and that there is
	// no next polynomial to combine it with.
	// Therefore we update it's degree to the one of the previous one.
	if giantSteps[i] == 2 {
		babySteps[i].Degree = babySteps[i-1].Degree

		// If two consecutive sub-polynomials, from ascending degree order, have the
		// same degree, we combine them.
	} else if giantSteps[i] == 1 {

		even, odd := babySteps[i], babySteps[i+1]

		deg := 1 << bits.Len64(uint64(babySteps[i].Degree))

		if err = EvaluateMonomial(even.Value, odd.Value, pb.Value[deg]); err != nil {
			return
		}

		odd.Degree = 2*deg - 1
		babySteps[i] = nil

		i++
	}

	return
}

// EvaluateMonomial evaluates a monomial of the form a + b * X^{pow} and writes the results in b.
func EvaluateMonomial(a, b, xpow *Element) (err error) {

	if b.Degree == 2 {
		b.Relinearize()
	}

	b.Rescale()

	b.Mul(b, xpow)

	if !a.Scale.InDelta(b.Scale, float64(rlwe.ScalePrecision-12)) {
		return fmt.Errorf("evalMonomial: scale discrepency: (rescale(b) * X^{n}).Scale = %v != a.Scale = %v", &a.Scale.Value, &b.Scale.Value)
	}

	b.Add(b, a)

	return
}

// EvaluatePolynomialVectorFromPowerBasis a method that complies to the interface he.PolynomialVectorEvaluator. This method evaluates P(ct) = sum c_i * ct^{i}.
func EvaluatePolynomialVectorFromPowerBasis[T any](targetLevel int, pol he.PolynomialVector, cg he.CoefficientGetter[T], pb PowerBasis, targetScale rlwe.Scale) (res *Element, err error) {

	// Map[int] of the powers [X^{0}, X^{1}, X^{2}, ...]
	X := pb.Value

	params := *X[1].Parameters
	mapping := pol.Mapping
	even := pol.IsEven()
	odd := pol.IsOdd()

	// Retrieve the degree of the highest degree non-zero coefficient
	// TODO: optimize for nil/zero coefficients
	minimumDegreeNonZeroCoefficient := len(pol.Value[0].Coeffs) - 1
	if even && !odd {
		minimumDegreeNonZeroCoefficient--
	}

	// Gets the maximum degree of the ciphertexts among the power basis
	// TODO: optimize for nil/zero coefficients, odd/even polynomial
	maximumCiphertextDegree := 0
	for i := pol.Value[0].Degree(); i > 0; i-- {
		if x, ok := X[i]; ok {
			maximumCiphertextDegree = utils.Max(maximumCiphertextDegree, x.Degree)
		}
	}

	// If an index slot is given (either multiply polynomials or masking)
	if mapping != nil {

		// If the degree of the poly is zero
		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewElement(params, make([]float64, params.MaxSlots()), 1, targetScale)
			res.Level = targetLevel

			if even {
				res.Add(res, cg.GetVectorCoefficient(pol, 0))
			}

			return
		}

		// Allocates the output ciphertext
		res = NewElement(params, make([]float64, params.MaxSlots()), maximumCiphertextDegree, targetScale)
		res.Level = targetLevel

		if even {
			res.Add(res, cg.GetVectorCoefficient(pol, 0))
		}

		// Loops starting from the highest degree coefficient
		for key := pol.Value[0].Degree(); key > 0; key-- {
			if !(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd) {
				res.MulThenAdd(X[key], cg.GetVectorCoefficient(pol, key))
			}
		}

	} else {

		if minimumDegreeNonZeroCoefficient == 0 {

			res = NewElement(params, make([]float64, params.MaxSlots()), 1, targetScale)
			res.Level = targetLevel

			if even {
				res.Add(res, cg.GetSingleCoefficient(pol.Value[0], 0))
			}

			return
		}

		res = NewElement[*bignum.Complex](params, nil, maximumCiphertextDegree, targetScale)
		res.Level = targetLevel

		if even {
			res.Add(res, cg.GetSingleCoefficient(pol.Value[0], 0))
		}

		for key := pol.Value[0].Degree(); key > 0; key-- {
			if key != 0 && (!(even || odd) || (key&1 == 0 && even) || (key&1 == 1 && odd)) {
				// MulScalarAndAdd automatically scales c to match the scale of res.
				res.MulThenAdd(X[key], cg.GetSingleCoefficient(pol.Value[0], key))
			}
		}


	}

	return
}
