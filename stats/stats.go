package stats

import (
	"fmt"
	"math"
)

var Header = []string{
	"MIN",
	"AVG",
	"STD",
}

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	MaxDelta      float64
	MinPrecision  float64
	MeanDelta     float64
	MeanPrecision float64
	StdPrecision  float64

	diff []float64
	prec []float64
}

func NewPrecisionStats() (prec *PrecisionStats) {
	return &PrecisionStats{
		MaxDelta:  0,
		MeanDelta: 0,

		diff: []float64{},
		prec: []float64{},
	}
}

func (p *PrecisionStats) Update(values []float64) {

	var abs float64

	for _, v := range values {

		p.diff = append(p.diff, v)

		abs = math.Abs(v)

		p.prec = append(p.prec, abs)

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
		fmt.Sprintf("%.5f", prec.MinPrecision),
		fmt.Sprintf("%.5f", prec.MeanPrecision),
		fmt.Sprintf("%.5f", prec.StdPrecision),
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
