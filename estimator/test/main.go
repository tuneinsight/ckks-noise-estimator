package main


import(
	"fmt"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
)

func main(){

	LogN := 16
	LogSlots := 15
	LogScale := 45

	minR := 0.5
	maxR := 2.0
	minI := 0.5
	maxI := 2.0

	prec := uint(256)

	pt := estimator.NewPlaintextDistribution(LogN, LogSlots, LogScale, minR, maxR, minI, maxI, prec)

	fmt.Println(pt)
}