package main

import (
	"math"
	"github.com/tuneinsight/ckks-bootstrapping-precision/operations"
)

var LogN = 16
var LogSlots = 15
var H = 32768
var LogScale = 55
var PlaintextSTD = 1/math.Sqrt(float64(uint64(1<<LogN)))
var Runs = 1

func main() {
	//operations.GetNoisRescale(LogN, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulPt(LogN, LogSlots, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulCt(LogN, LogScale, PlaintextSTD, Runs)
	//operations.KeySwitchHoisted(LogN, H, Runs)
	
	/*
	ddd := map[int]float64{}

	for i := 0; i < 1024; i++{
		ddd[i] = 0.02209306164843668
	}

	operations.GetNoiseLinearTransform(LogN, H, LogSlots, LogScale, ddd, PlaintextSTD, Runs)
	*/

	operations.GetNoisePowerBasis(LogN, H, LogSlots, LogScale, Runs)
}
