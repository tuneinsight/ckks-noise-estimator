package main

import(
	"github.com/tuneinsight/ckks-bootstrapping-precision/operations"
)

var LogN = 16
var LogScale = 45
var PlaintextSTD = 3.2
var Runs = 1024

func main(){
	//operations.GetNoisRescale(LogN, LogScale, PlaintextSTD, Runs)
	//operations.GetNoiseMulPt(LogN, LogScale, PlaintextSTD, Runs)
	operations.GetNoiseMulCt(LogN, LogScale, PlaintextSTD, Runs)
}