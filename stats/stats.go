package stats

import (
	"fmt"
	"math"
)

var Header = []string{
	"H",
	"MIN",
	"AVG",
	"STD",
}

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	LogN, H, Depth, LogScale, LogSlots int
	MaxDelta                           float64
	MinPrecision                       float64
	MeanDelta                          float64
	MeanPrecision                      float64
	MedianDelta                        float64
	MedianPrecision                    float64
	StdPrecision                       float64

	/*
		If we want to produce CDF graphs, we can with this
		RealDist, ImagDist, L2Dist []struct {
			Prec  float64
			Count int
		}

		cdfResol int
	*/

	diff     []float64
	precReal []float64
	precImag []float64
	precL2   []float64
}

func NewPrecisionStats(H int) (prec *PrecisionStats) {
	return &PrecisionStats{
		H:         H,
		MaxDelta:  0,
		MeanDelta: 0,

		diff:     []float64{},
		precReal: []float64{},
		precImag: []float64{},
		precL2:   []float64{},
	}
}

func (p *PrecisionStats) Update(values []float64) {

	var abs float64

	for _, v := range values {

		p.diff = append(p.diff, v)

		abs = math.Abs(v)

		p.precReal = append(p.precReal, abs)

		p.MeanDelta += abs

		if abs > p.MaxDelta {
			p.MaxDelta = abs
		}
	}
}

func (p *PrecisionStats) Finalize() {

	p.StdPrecision = p.calcstandarddeviation(p.diff)
	p.MinPrecision = deltaToPrecision(p.MaxDelta)

	p.MeanDelta /= float64(len(p.diff))
	p.MeanPrecision = deltaToPrecision(p.MeanDelta)
}

func (prec *PrecisionStats) ToCSV() []string {
	return []string{
		fmt.Sprintf("%d", prec.H),
		fmt.Sprintf("%.4f", prec.MinPrecision),
		fmt.Sprintf("%.4f", prec.MeanPrecision),
		fmt.Sprintf("%.4f", prec.StdPrecision),
	}
}

func (p *PrecisionStats) calcstandarddeviation(diff []float64) (std float64) {

	var mean float64
	for _, d := range diff {
		mean += d
	}

	n := float64(len(diff))

	mean /= n

	var v float64
	for _, d := range diff {
		v = (d - mean)
		std += v * v
	}

	std = math.Log2(math.Sqrt(std / (n - 1)))

	fmt.Printf("STD Log2: %7.4f\n", std)

	return
}

func deltaToPrecision(c float64) float64 {
	return math.Log2(1 / maxFloat64(c, 1e-16))
}

func maxFloat64(a, b float64) (c float64) {
	if a > b {
		return a
	}

	return b
}