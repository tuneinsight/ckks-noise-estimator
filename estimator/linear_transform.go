package estimator

import (
	"math/big"
)

type LinearTransform struct {
	N1         int
	Level      int
	Diagonales map[int]Element
	Index      map[int][]int
}

// NewLinearTransform
//
// diags: map of non-zero diagonales, with [std plaintext, std encoding error]
func NewLinearTransform(diags map[int][2]float64, scale interface{}, Level, LogSlots, LogBSGSRatio int) (LT LinearTransform) {

	slots := 1 << LogSlots

	N1 := FindBestBSGSRatio(diags, slots, LogBSGSRatio)

	Index, _, _ := BSGSIndex(diags, slots, N1)

	Diagonales := make(map[int]Element)
	for j := range Index {
		for _, i := range Index[j] {

			v, ok := diags[j+i]
			if !ok {
				v = diags[j+i-slots]
			}

			Diagonales[j+i] = NewPlaintext(new(big.Float).Mul(NewFloat(v[0]), NewFloat(scale)), new(big.Float).Mul(NewFloat(v[1]), NewFloat(scale)), Level)
		}
	}

	return LinearTransform{N1: N1, Level: Level, Diagonales: Diagonales, Index: Index}
}

func (e *Estimator) LinearTransform(el0 Element, LT LinearTransform) (el1 Element) {

	Level := MinInt(el0.Level, LT.Level)

	// Accumulator outer-loop
	el1 = Element{}
	el1.Level = Level
	el1.Message = NewFloat(0)
	el1.Noise = []*big.Float{
		NewFloat(0),
		NewFloat(0),
	}

	// Hoisted n1 Rotations
	el0RotHoisted := e.KeySwitchLazy(el0)

	// el0 * P
	el0P := e.Mul(el0, e.P)

	index := LT.Index

	for j := range index {

		// Accumulator inner-loop
		tmp := NewPlaintext(nil, nil, Level)

		for _, i := range index[j] {

			if i == 0 {
				tmp = e.Add(tmp, e.Mul(el0P, LT.Diagonales[j])) // el0 * P * diag[j]
			} else {
				tmp = e.Add(tmp, e.Mul(el0RotHoisted, LT.Diagonales[i+j])) // Rotate(el0) * P * diag[i+j]
			}
		}

		if j != 0 {
			// Rotate n2 of sum(P * Rotate(el0) * diags) / P
			el1 = e.Add(el1, e.KeySwitchLazy(e.ModDown(tmp)))
		} else {
			el1 = e.Add(el1, tmp)
		}
	}

	el1 = e.ModDown(el1)

	return
}

// FindBestBSGSRatio finds the best N1*N2 = N for the baby-step giant-step algorithm for matrix multiplication.
func FindBestBSGSRatio(diagMatrix interface{}, maxN int, logMaxRatio int) (minN int) {

	maxRatio := float64(int(1 << logMaxRatio))

	for N1 := 1; N1 < maxN; N1 <<= 1 {

		_, rotN1, rotN2 := BSGSIndex(diagMatrix, maxN, N1)

		nbN1, nbN2 := len(rotN1)-1, len(rotN2)-1

		if float64(nbN2)/float64(nbN1) == maxRatio {
			return N1
		}

		if float64(nbN2)/float64(nbN1) > maxRatio {
			return N1 / 2
		}
	}

	return 1
}

// BSGSIndex returns the index map and needed rotation for the BSGS matrix-vector multiplication algorithm.
func BSGSIndex(el interface{}, slots, N1 int) (index map[int][]int, rotN1, rotN2 []int) {
	index = make(map[int][]int)
	rotN1Map := make(map[int]bool)
	rotN2Map := make(map[int]bool)
	var nonZeroDiags []int
	switch element := el.(type) {
	case map[int]float64:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int][2]float64:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]Element:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	case map[int]bool:
		nonZeroDiags = make([]int, len(element))
		var i int
		for key := range element {
			nonZeroDiags[i] = key
			i++
		}
	}

	for _, rot := range nonZeroDiags {
		rot &= (slots - 1)
		idxN1 := ((rot / N1) * N1) & (slots - 1)
		idxN2 := rot & (N1 - 1)
		if index[idxN1] == nil {
			index[idxN1] = []int{idxN2}
		} else {
			index[idxN1] = append(index[idxN1], idxN2)
		}
		rotN1Map[idxN1] = true
		rotN2Map[idxN2] = true
	}

	rotN1 = []int{}
	for i := range rotN1Map {
		rotN1 = append(rotN1, i)
	}

	rotN2 = []int{}
	for i := range rotN2Map {
		rotN2 = append(rotN2, i)
	}

	return
}
