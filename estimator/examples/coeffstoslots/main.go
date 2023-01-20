package main

import(
	"fmt"
	"github.com/tuneinsight/lattigo/v4/utils"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/ckks/advanced"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
)

var (
	LogN     = 16                     // Log2 ring degree
	LogSlots = 15                     // Log2 #slots
	LogScale = 45                     // Log2 scaling factor
	H        = 32768                  // SecretKey Hamming weight
	Depth    = 4                      // Depth of the homomorphic encoding/decoding
	LtType   = advanced.SlotsToCoeffs //advanced.CoeffsToSlots //
)

func main(){
	var err error

	LogQ := make([]int, Depth+1)

	LogQ[0] = 60
	for i := 1; i < Depth+1; i++ {
		LogQ[i] = LogScale
	}

	LogP := []int{61, 61}

	var params ckks.Parameters
	if params, err = ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:     LogN,
		LogQ:     LogQ,
		LogP:     LogP,
		H:        H,
		LogSlots: LogSlots,
		LogScale: LogScale,
	}); err != nil {
		panic(err)
	}

	values := make([]float64, params.N())

	for i := range values{
		values[i] = utils.RandFloat64(-1, 1)
	}

	est := estimator.NewEstimator(params.N(), params.HammingWeight(), params.Q(), params.P())

	LTs := GetSlotsToCoeffsSTD(params, ckks.NewEncoder(params))

	ct := estimator.NewCiphertextSK(estimator.NewPlaintext(estimator.STD(values) * params.DefaultScale().Float64(), nil, params.MaxLevel()))

	
	switch LtType {
	case advanced.SlotsToCoeffs:
		// Repack imaginary
		if params.LogSlots() == params.LogN()-1{
			ct = est.Add(ct, ct)
		}
	}
	

	for i := range LTs{
		ct = est.LinearTransform(ct, LTs[i])
		ct = est.Rescale(ct)
	}

	switch LtType{
	case advanced.CoeffsToSlots:
		ctConj := est.KeySwitch(ct)

		ct = est.Add(ctConj, ct)

		if params.LogSlots() < params.LogN()-1{
			ct = est.Add(ct, ct)
		}
	}
	
	fmt.Println(ct.Message)
	have := est.Std(ct)

	fmt.Println(have)
}

func p(){
	fmt.Println()
}

func GetSlotsToCoeffsSTD(params ckks.Parameters, ecd ckks.Encoder) (LTs []estimator.LinearTransform){

	Levels := make([]int, params.MaxLevel())
	for i := range Levels {
		Levels[i] = 1
	}

	encodingMatrixLiteral := advanced.EncodingMatrixLiteral{
		LinearTransformType: LtType,
		LogN:                params.LogN(),
		LogSlots:            params.LogSlots(),
		LevelStart:          params.MaxLevel(),
		Levels:              Levels,
		RepackImag2Real:     true,
		LogBSGSRatio:        2,
	}

	encodingMatrices := encodingMatrixLiteral.ComputeDFTMatrices()

	LTs = make([]estimator.LinearTransform, len(encodingMatrices))

	values := make([]float64, params.N())

	for i, matrix := range encodingMatrices{

		m := make(map[int]float64)

		for j, diag := range matrix{

			EncodePlaintext(ecd, params.N(), params.LogSlots(), diag, values)

			m[j] = estimator.STD(values)
		}

		//fmt.Println(m)
		//fmt.Println()

		Level := encodingMatrixLiteral.LevelStart-i

		LTs[i] = estimator.NewLinearTransform(m, params.Q()[Level], Level, params.LogSlots(), encodingMatrixLiteral.LogBSGSRatio)
	}

	return
}

func EncodePlaintext(ecd ckks.Encoder, N, LogSlots int, a []complex128, b []float64){

	ecd.IFFT(a, LogSlots)

	slots := 1<<LogSlots
	gap := N/(2*slots)
	for i, j := 0, N>>1; i < slots; i, j = i+gap, j+gap{
		b[i] = real(a[i])
		b[j] = imag(a[i])
	}
}