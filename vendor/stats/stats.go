package stats

import (
	"fmt"
	"math"
	"sort"
)

// PrecisionStats is a struct storing statistic about the precision of a CKKS plaintext
type PrecisionStats struct {
	LogN, H, Depth, LogScale, LogSlots int
	MaxDelta                           Stats
	MinPrecision                       Stats
	MeanDelta                          Stats
	MeanPrecision                      Stats
	MedianDelta                        Stats
	MedianPrecision                    Stats
	StdPrecision                       Stats

	/*
		If we want to produce CDF graphs, we can with this
		RealDist, ImagDist, L2Dist []struct {
			Prec  float64
			Count int
		}

		cdfResol int
	*/

	diff     []Stats
	precReal []float64
	precImag []float64
	precL2   []float64
}

// Stats is a struct storing the real, imaginary and L2 norm (modulus)
// about the precision of a complex value.
type Stats struct {
	Real, Imag, L2 float64
}

func NewPrecisionStats(LogN, H, Depth, LogSlots, LogScale int) (prec *PrecisionStats) {
	return &PrecisionStats{
		LogN:      LogN,
		H:         H,
		Depth:     Depth,
		LogScale:  LogScale,
		LogSlots:  LogSlots,
		MaxDelta:  Stats{0, 0, 0},
		MeanDelta: Stats{0, 0, 0},

		diff:     []Stats{},
		precReal: []float64{},
		precImag: []float64{},
		precL2:   []float64{},
	}
}

func (p *PrecisionStats) Update(want, have interface{}) {

	var valuesTest []complex128

	switch have := have.(type) {
	case []complex128:
		valuesTest = have
	case []float64:
		valuesTest = make([]complex128, len(have))
		for i := range have {
			valuesTest[i] = complex(have[i], 0)
		}
	}

	var valuesWant []complex128
	switch want := want.(type) {
	case []complex128:
		valuesWant = want
	case []float64:
		valuesWant = make([]complex128, len(want))
		for i := range want {
			valuesWant[i] = complex(want[i], 0)
		}
	}

	var deltaReal, deltaImag, deltaL2 float64

	for i := range valuesWant {

		deltaReal = math.Abs(real(valuesTest[i]) - real(valuesWant[i]))
		deltaImag = math.Abs(imag(valuesTest[i]) - imag(valuesWant[i]))
		deltaL2 = math.Sqrt(deltaReal*deltaReal + deltaImag*deltaImag)

		p.precReal = append(p.precReal, math.Log2(1/deltaReal))
		p.precImag = append(p.precImag, math.Log2(1/deltaImag))
		p.precL2 = append(p.precL2, math.Log2(1/deltaL2))
		p.diff = append(p.diff, Stats{Real: deltaReal, Imag: deltaImag, L2: deltaL2})

		p.MeanDelta.Real += deltaReal
		p.MeanDelta.Imag += deltaImag
		p.MeanDelta.L2 += deltaL2

		if deltaReal > p.MaxDelta.Real {
			p.MaxDelta.Real = deltaReal
		}

		if deltaImag > p.MaxDelta.Imag {
			p.MaxDelta.Imag = deltaImag
		}

		if deltaL2 > p.MaxDelta.L2 {
			p.MaxDelta.L2 = deltaL2
		}
	}
}

func (p *PrecisionStats) Finalize() {

	p.MeanDelta.Real /= float64(len(p.diff))
	p.MeanDelta.Imag /= float64(len(p.diff))
	p.MeanDelta.L2 /= float64(len(p.diff))

	p.MedianDelta = calcmedian(p.diff)

	p.StdPrecision = calcstandarddeviation(p.diff, p.MeanDelta)

	p.MinPrecision = deltaToPrecision(p.MaxDelta)
	p.MeanPrecision = deltaToPrecision(p.MeanDelta)
	p.MedianPrecision = deltaToPrecision(p.MedianDelta)

	/*
		p.cdfResol = 500

		p.RealDist = make([]struct {
			Prec  float64
			Count int
		}, p.cdfResol)
		p.ImagDist = make([]struct {
			Prec  float64
			Count int
		}, p.cdfResol)
		p.L2Dist = make([]struct {
			Prec  float64
			Count int
		}, p.cdfResol)


		p.calcCDF(p.precReal, p.RealDist)
		p.calcCDF(p.precImag, p.ImagDist)
		p.calcCDF(p.precL2, p.L2Dist)
	*/
}

func (prec *PrecisionStats) ToCSV() []string {
	return []string{
		fmt.Sprintf("%d", prec.LogN),
		fmt.Sprintf("%d", prec.H),
		fmt.Sprintf("%d", prec.Depth),
		fmt.Sprintf("%d", prec.LogScale),
		fmt.Sprintf("%d", prec.LogSlots),
		fmt.Sprintf("%.4f", prec.MinPrecision.Real),
		fmt.Sprintf("%.4f", prec.MeanPrecision.Real),
		fmt.Sprintf("%.4f", prec.MedianPrecision.Real),
		fmt.Sprintf("%.4f", prec.StdPrecision.Real),
		fmt.Sprintf("%.4f", prec.MinPrecision.Imag),
		fmt.Sprintf("%.4f", prec.MeanPrecision.Imag),
		fmt.Sprintf("%.4f", prec.MedianPrecision.Imag),
		fmt.Sprintf("%.4f", prec.StdPrecision.Imag),
	}
}

func calcstandarddeviation(diff []Stats, mean Stats) (std Stats) {
	std = Stats{}

	mReal := math.Log2(1 / mean.Real)
	mImag := math.Log2(1 / mean.Imag)
	mL2 := math.Log2(1 / mean.L2)

	var real, imag, l2 float64

	var v float64
	for _, d := range diff {

		v = math.Log2(1/maxFloat64(d.Real, 1e-16)) - mReal
		real += v * v

		v = math.Log2(1/maxFloat64(d.Imag, 1e-16)) - mImag
		imag += v * v

		v = math.Log2(1/maxFloat64(d.L2, 1e-16)) - mL2
		l2 += v * v
	}

	n := float64(len(diff))

	std.Real = math.Sqrt(real / n)
	std.Imag = math.Sqrt(imag / n)
	std.L2 = math.Sqrt(l2 / n)

	return
}

func deltaToPrecision(c Stats) Stats {
	return Stats{math.Log2(1 / maxFloat64(c.Real, 1e-16)), math.Log2(1 / maxFloat64(c.Imag, 1e-16)), math.Log2(1 / maxFloat64(c.L2, 1e-16))}
}

func maxFloat64(a, b float64) (c float64) {
	if a > b {
		return a
	}

	return b
}

/*
func (prec *PrecisionStats) calcCDF(precs []float64, res []struct {
	Prec  float64
	Count int
}) {
	sortedPrecs := make([]float64, len(precs))
	copy(sortedPrecs, precs)
	sort.Float64s(sortedPrecs)
	minPrec := sortedPrecs[0]
	maxPrec := sortedPrecs[len(sortedPrecs)-1]
	for i := 0; i < prec.cdfResol; i++ {
		curPrec := minPrec + float64(i)*(maxPrec-minPrec)/float64(prec.cdfResol)
		for countSmaller, p := range sortedPrecs {
			if p >= curPrec {
				res[i].Prec = curPrec
				res[i].Count = countSmaller
				break
			}
		}
	}
}
*/

func calcmedian(values []Stats) (median Stats) {

	tmp := make([]float64, len(values))

	for i := range values {
		tmp[i] = values[i].Real
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].Real = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].Imag
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].Imag = tmp[i]
	}

	for i := range values {
		tmp[i] = values[i].L2
	}

	sort.Float64s(tmp)

	for i := range values {
		values[i].L2 = tmp[i]
	}

	index := len(values) / 2

	if len(values)&1 == 1 || index+1 == len(values) {
		return Stats{values[index].Real, values[index].Imag, values[index].L2}
	}

	return Stats{(values[index].Real + values[index+1].Real) / 2,
		(values[index].Imag + values[index+1].Imag) / 2,
		(values[index].L2 + values[index+1].L2) / 2}
}
