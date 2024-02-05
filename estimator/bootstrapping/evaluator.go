package estimator

import (
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
	"github.com/tuneinsight/lattigo/v5/utils/sampling"
)

type Evaluator struct {
	bootstrapping.Parameters

	ResidualParameters      estimator.Parameters
	BootstrappingParameters estimator.Parameters

	S2CDFTMatrix   estimator.DFTMatrix
	C2SDFTMatrix   estimator.DFTMatrix
	Mod1Parameters hefloat.Mod1Parameters

	EphemeralSecret []*bignum.Complex
}

func NewEvaluator(btpParams bootstrapping.Parameters) Evaluator {

	eval := Evaluator{
		Parameters:              btpParams,
		ResidualParameters:      estimator.NewParameters(btpParams.ResidualParameters),
		BootstrappingParameters: estimator.NewParameters(btpParams.BootstrappingParameters),
	}

	if btpParams.EphemeralSecretWeight != 0{
		eval.EphemeralSecret = eval.BootstrappingParameters.SampleSecretKey(btpParams.EphemeralSecretWeight)
	}

	eval.initialize(btpParams)

	return eval
}

func (eval *Evaluator) initialize(btpParams bootstrapping.Parameters) (err error) {

	params := btpParams.BootstrappingParameters

	if eval.Mod1Parameters, err = hefloat.NewMod1ParametersFromLiteral(params, btpParams.Mod1ParametersLiteral); err != nil {
		return
	}

	// [-K, K]
	K := eval.Mod1Parameters.K

	// Correcting factor for approximate division by Q
	// The second correcting factor for approximate multiplication by Q is included in the coefficients of the EvalMod polynomials
	qDiff := eval.Mod1Parameters.QDiff

	// If the scale used during the EvalMod step is smaller than Q0, then we cannot increase the scale during
	// the EvalMod step to get a free division by MessageRatio, and we need to do this division (totally or partly)
	// during the CoeffstoSlots step
	qDiv := eval.Mod1Parameters.ScalingFactor().Float64() / math.Exp2(math.Round(math.Log2(float64(params.Q()[0]))))

	// Sets qDiv to 1 if there is enough room for the division to happen using scale manipulation.
	if qDiv > 1 {
		qDiv = 1
	}

	// CoeffsToSlots vectors
	// Change of variable for the evaluation of the Chebyshev polynomial + cancelling factor for the DFT and SubSum + eventual scaling factor for the double angle formula

	scale := eval.BootstrappingParameters.DefaultScale().Float64()
	offset := eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()

	C2SScaling := new(big.Float).SetFloat64(qDiv / (K * qDiff))
	StCScaling := new(big.Float).SetFloat64(scale / offset)

	eval.C2SDFTMatrix = estimator.DFTMatrix{DFTMatrixLiteral: btpParams.CoeffsToSlotsParameters}
	if eval.C2SDFTMatrix.Scaling == nil {
		eval.C2SDFTMatrix.Scaling = C2SScaling
	} else {
		eval.C2SDFTMatrix.Scaling = new(big.Float).Mul(btpParams.CoeffsToSlotsParameters.Scaling, C2SScaling)
	}

	eval.C2SDFTMatrix.GenMatrices(btpParams.BootstrappingParameters.LogN(), 128)

	eval.S2CDFTMatrix = estimator.DFTMatrix{DFTMatrixLiteral: btpParams.SlotsToCoeffsParameters}
	if eval.S2CDFTMatrix.Scaling == nil {
		eval.S2CDFTMatrix.Scaling = StCScaling
	} else {
		eval.S2CDFTMatrix.Scaling = new(big.Float).Mul(btpParams.SlotsToCoeffsParameters.Scaling, StCScaling)
	}

	eval.S2CDFTMatrix.GenMatrices(btpParams.BootstrappingParameters.LogN(), 128)

	return
}

// checks if the current message ratio is greater or equal to the last prime times the target message ratio.
func checkMessageRatio(el *estimator.Element, msgRatio float64, r *ring.Ring) bool {
	level := el.Level
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[level])
	currentMessageRatio = currentMessageRatio.Div(el.Scale)
	return currentMessageRatio.Cmp(rlwe.NewScale(r.SubRings[level].Modulus).Mul(rlwe.NewScale(msgRatio))) > -1
}

func (eval Evaluator) Bootstrap(el *estimator.Element) (*estimator.Element) {

	eval.ScaleDown(el)

	eval.ModUp(el)

	elReal, elImag := eval.CoeffsToSlots(el)

	elReal = eval.EvalMod(elReal)

	if elImag != nil{
		elImag = eval.EvalMod(elImag)
	}

	el = eval.SlotsToCoeffs(elReal, elImag)

	return el
}

func (eval Evaluator) ScaleDown(el *estimator.Element) *rlwe.Scale {

	params := &eval.BootstrappingParameters.Parameters

	r := params.RingQ()

	// Removes unecessary primes
	for el.Level != 0 && checkMessageRatio(el, eval.Mod1Parameters.MessageRatio(), r) {
		el.Level--
	}

	// Current Message Ratio
	currentMessageRatio := rlwe.NewScale(r.ModulusAtLevel[el.Level])
	currentMessageRatio = currentMessageRatio.Div(el.Scale)

	// Desired Message Ratio
	targetMessageRatio := rlwe.NewScale(eval.Mod1Parameters.MessageRatio())

	// (Current Message Ratio) / (Desired Message Ratio)
	scaleUp := currentMessageRatio.Div(targetMessageRatio)

	if scaleUp.Cmp(rlwe.NewScale(0.5)) == -1 {
		panic(fmt.Errorf("initial Q/Scale = %f < 0.5*Q[0]/MessageRatio = %f", currentMessageRatio.Float64(), targetMessageRatio.Float64()))
	}

	scaleUpBigint := scaleUp.BigInt()

	el.Mul(el, scaleUpBigint)
	el.Scale = el.Scale.Mul(rlwe.NewScale(scaleUpBigint))

	// errScale = CtIn.Scale/(Q[0]/MessageRatio)
	targetScale := new(big.Float).SetPrec(256).SetInt(r.ModulusAtLevel[0])
	targetScale.Quo(targetScale, new(big.Float).SetFloat64(eval.Mod1Parameters.MessageRatio()))

	if el.Level != 0 {
		el.Rescale()
	}

	// Rescaling error (if any)
	errScale := el.Scale.Div(rlwe.NewScale(targetScale))

	return &errScale
}

func (eval Evaluator) ModUp(el *estimator.Element) {

	var H int
	if eval.EphemeralSecret != nil{
		el.KeySwitch(eval.BootstrappingParameters.Sk[0])
		H = eval.EphemeralSecretWeight
	}else{
		H = eval.BootstrappingParameters.H
	}

	K := new(big.Float).SetPrec(128)

	Q := eval.BootstrappingParameters.Q[0]

	irwinHall := func() *big.Float {
		var d float64
		for i := 0; i < H+1; i++ {
			d += sampling.RandFloat64(0, 1)
		}

		K.SetFloat64(math.Floor(d) - float64(H>>1))

		return new(big.Float).Mul(K, Q)
	}

	values := make([]*bignum.Complex, el.MaxSlots())

	for i := range values {
		values[i] = &bignum.Complex{irwinHall(), irwinHall()}
	}

	eval.BootstrappingParameters.Encoder.FFT(values, el.LogMaxSlots())

	coeffs := el.Value[0]
	for i := range coeffs {
		coeffs[i].Add(coeffs[i], values[i])
	}

	el.Level = eval.BootstrappingParameters.MaxLevel()
	el.Parameters = &eval.BootstrappingParameters

	if eval.EphemeralSecret != nil{
		el.KeySwitch(eval.EphemeralSecret)
	}

	// Scale the message from Q0/|m| to QL/|m|, where QL is the largest modulus used during the bootstrapping.
	if scale := (eval.Mod1Parameters.ScalingFactor().Float64() / eval.Mod1Parameters.MessageRatio()) / el.Scale.Float64(); scale > 1 {
		el.ScaleUp(rlwe.NewScale(scale))
	}
}

func (eval Evaluator) CoeffsToSlots(el *estimator.Element) (elReal, elImag *estimator.Element){
	return el.CoeffsToSlots(eval.C2SDFTMatrix)
}

func (eval Evaluator) SlotsToCoeffs(elReal, elImag *estimator.Element) (el *estimator.Element){
	el = estimator.NewElement[*bignum.Complex](eval.BootstrappingParameters, nil, 1, eval.BootstrappingParameters.DefaultScale())
	
	return el.SlotsToCoeffs(elReal, elImag, eval.S2CDFTMatrix)
}

func (eval Evaluator) EvalMod(el *estimator.Element) (*estimator.Element){
	el, _ = el.EvaluateMod1(eval.Mod1Parameters)
	el.Scale = eval.BootstrappingParameters.DefaultScale()
	return el
}