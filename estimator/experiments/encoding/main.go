package main


import(
	"fmt"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/ring"
)

/*
Average over 128 runs
LogScale = 45
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 37.20 │ 37.20 │ 36.70 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 40.21 │ 40.21 │ 39.71 │
│MED Prec │ 39.86 │ 39.86 │ 39.36 │
│STD Prec │  1.60 │  1.60 │  1.60 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │  7.80 │  7.80 │  8.30 │
│AVG Err  │  4.79 │  4.79 │  5.29 │
│MED Err  │  5.14 │  5.14 │  5.64 │
│STD Err  │  1.60 │  1.60 │  1.60 │
└─────────┴───────┴───────┴───────┘
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 37.20 │ 37.19 │ 36.70 │
│MAX Prec │ 45.00 │ 45.00 │ 45.00 │
│AVG Prec │ 40.21 │ 40.21 │ 39.71 │
│MED Prec │ 39.86 │ 39.86 │ 39.36 │
│STD Prec │  1.60 │  1.60 │  1.60 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │  7.80 │  7.81 │  8.30 │
│AVG Err  │  4.79 │  4.79 │  5.29 │
│MED Err  │  5.14 │  5.14 │  5.64 │
│STD Err  │  1.60 │  1.60 │  1.60 │
└─────────┴───────┴───────┴───────┘

Average over 128 runs
LogScale = 90
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 82.19 │ 82.18 │ 81.69 │
│MAX Prec │ 90.00 │ 90.00 │ 90.00 │
│AVG Prec │ 85.21 │ 85.21 │ 84.71 │
│MED Prec │ 84.86 │ 84.86 │ 84.36 │
│STD Prec │  1.60 │  1.60 │  1.60 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │  7.81 │  7.82 │  8.31 │
│AVG Err  │  4.79 │  4.79 │  5.29 │
│MED Err  │  5.14 │  5.14 │  5.64 │
│STD Err  │  1.60 │  1.60 │  1.60 │
└─────────┴───────┴───────┴───────┘
┌─────────┬───────┬───────┬───────┐
│    Log2 │ REAL  │ IMAG  │ L2    │
├─────────┼───────┼───────┼───────┤
│MIN Prec │ 82.21 │ 82.20 │ 81.71 │
│MAX Prec │ 90.00 │ 90.00 │ 90.00 │
│AVG Prec │ 85.19 │ 85.19 │ 84.69 │
│MED Prec │ 84.86 │ 84.85 │ 84.36 │
│STD Prec │  1.53 │  1.53 │  1.54 │
├─────────┼───────┼───────┼───────┤
│MIN Err  │  0.00 │  0.00 │  0.00 │
│MAX Err  │  7.79 │  7.80 │  8.29 │
│AVG Err  │  4.81 │  4.81 │  5.31 │
│MED Err  │  5.14 │  5.15 │  5.64 │
│STD Err  │  1.53 │  1.53 │  1.54 │
└─────────┴───────┴───────┴───────┘
*/

func main(){
	params, err := hefloat.NewParametersFromLiteral(hefloat.ParametersLiteral{
		LogN:            16,
		LogQ:            []int{55, 45},
		LogP:            []int{60},
		LogDefaultScale: 45,
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil{
		panic(err)
	}

	ecd := hefloat.NewEncoder(params)

	estParams := estimator.NewParameters(params)
	estParams.Heuristic = true

	statsHave := estimator.Stats{}
	statsWant := estimator.Stats{}

	for i := 0; i < 128; i++{

		values, el, pt, _ := estParams.NewTestVector(ecd, nil, -1-1i, 1+1i)

		el.Normalize()

		pWant := hefloat.GetPrecisionStats(params, ecd, nil, values, el.Value[0], 0, false)
		pHave := hefloat.GetPrecisionStats(params, ecd, nil, values, pt, 0, false)

		statsWant.Add(pWant)
		statsHave.Add(pHave)
	}

	statsWant.Finalize()
	statsHave.Finalize()

	fmt.Println(statsWant.String())
	fmt.Println(statsHave.String())
}

