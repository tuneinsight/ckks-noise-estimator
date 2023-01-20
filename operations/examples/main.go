package main

import(
	"github.com/tuneinsight/ckks-bootstrapping-precision/operations"
)

var LogN = 16
var H    = 32768
var LogScale = 45
var PlaintextSTD = 1.0
var Runs = 16

func main(){
	//operations.GetNoisRescale(LogN, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulPt(LogN, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulCt(LogN, LogScale, PlaintextSTD, Runs)

	nonZeroDiags := []int{0, -1, -2, -3, -4, -5, -6, -7, 1, 2, 3, 4, 5, 6, 4096}

	operations.GetNoiseLinearTransform(LogN, H, LogScale, nonZeroDiags, PlaintextSTD, Runs)

	//operations.KeySwitchHoisted(LogN, H, Runs)
}