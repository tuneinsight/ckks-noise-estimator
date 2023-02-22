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
	
	diags := []map[int]float64{
		map[int]float64{
			0:     0.00024414102316630370,
			2048:  0.00024414233423388263,
			4096:  0.00024414206118669165,
			6144:  0.00024414240022099044,
			8192:  0.00024414248708016384,
			10240: 0.00024414247911935246,
			12288: 0.00024414180919533906,
			14336: 0.00024414248329247690,
			16384: 0.00024414246898360015,
			18432: 0.00024414239107164240,
			20480: 0.00024414226649732900,
			22528: 0.00024414233600733267,
			24576: 0.00024414232259683394,
			26624: 0.00024414248311350925,
			28672: 0.00024414158296428044,
			30720: 0.00024414248708458494,
		},
		map[int]float64{
			0:  0.003906267320465252,
			1:  0.0036538916134800097,
			2:  0.003382925367883238,
			3:  0.0030880284101041902,
			4:  0.002762148668487031,
			5:  0.002392059056988301,
			6:  0.0019530750723170012,
			7:  0.0013809426865437015,
			-7: 0.0013810771517379955,
			-6: 0.0019531398956096173,
			-5: 0.0023920937341996994,
			-4: 0.002762155692814875,
			-3: 0.0030881846979644947,
			-2: 0.003382916280984945,
			-1: 0.00365397976622493,
		},
	}
	_=diags

	ddd := map[int]float64{}

	for i := 0; i < 1024; i++{
		ddd[i] = 0.02209306164843668
	}

	operations.GetNoiseLinearTransform(LogN, H, LogSlots, LogScale, ddd, PlaintextSTD, Runs)
}
