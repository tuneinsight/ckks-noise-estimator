package main

import (
	"github.com/tuneinsight/ckks-bootstrapping-precision/operations"
)

var LogN = 16
var LogSlots = 10
var H = 32768
var LogScale = 45
var PlaintextSTD = 0.5772058896878792
var Runs = 16

func main() {
	//operations.GetNoisRescale(LogN, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulPt(LogN, LogSlots, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulCt(LogN, LogScale, PlaintextSTD, Runs)
	//operations.KeySwitchHoisted(LogN, H, Runs)
	
	ddd := map[int]float64{}

	for i := 0; i < 1024; i++{
		ddd[i] = 0.02209306164843668
	}

	operations.GetNoiseLinearTransform(LogN, H, LogSlots, LogScale, ddd, PlaintextSTD, Runs)
}
