package main

import(
	"github.com/tuneinsight/ckks-bootstrapping-precision/operations"
)

var LogN = 16
var LogScale = 45
var PlaintextSTD = 1.0
var Runs = 32

func main(){
	//operations.GetNoisRescale(LogN, LogScale, PlaintextSTD, Runs)
	operations.GetNoiseMulPt(LogN, LogScale, PlaintextSTD, Runs)
}