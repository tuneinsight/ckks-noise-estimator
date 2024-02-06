package estimator

import (
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

type LinearTransformation struct {
	LogSlots                 int
	LogBabyStepGianStepRatio int
	Scale                    rlwe.Scale
	Value                    hefloat.Diagonals[*bignum.Complex]
}

// BSGSIndex returns the BSGSIndex of the target linear transformation.
func (lt LinearTransformation) BSGSIndex() (index map[int][]int, n1, n2 []int) {
	cols := 1 << lt.LogSlots
	N1 := he.FindBestBSGSRatio(lt.Value.DiagonalsIndexList(), cols, lt.LogBabyStepGianStepRatio)
	return he.BSGSIndex(utils.GetKeys(lt.Value), cols, N1)
}

func (p *Element) EvaluateLinearTransformation(lt LinearTransformation) *Element {

	index, _, rotN2 := lt.BSGSIndex()

	ctPreRot := map[int][]*bignum.Complex{}

	for _, k := range rotN2 {
		if k != 0{
			if _, ok := ctPreRot[k]; k != 0 && !ok {
				ctPreRot[k] = p.RotateHoisted(k)
			}
		}
	}

	ctPreRot0 := p.CopyNew().Mul(p, p.P)

	keys := utils.GetSortedKeys(index)

	scale := &lt.Scale.Value
	slots := 1 << lt.LogSlots

	acc := NewElement[*bignum.Complex](*p.Parameters, nil, 1, p.Scale.Mul(lt.Scale))
	acc.Level = p.Level
	res := NewElement[*bignum.Complex](*p.Parameters, nil, 1, p.Scale.Mul(lt.Scale))
	res.Level = p.Level
	tmp := make([]*bignum.Complex, slots)

	bComplex := bignum.NewComplex()
	dComplex := bignum.NewComplex()

	mul := bignum.NewComplexMultiplier().Mul

	for _, j := range keys {

		rot := -j & (slots - 1)

		var cnt int
		for _, i := range index[j] {

			utils.RotateSliceAllocFree(lt.Value[i+j], rot, tmp)

			coeffsAcc0 := acc.Value[0]
			coeffsAcc1 := acc.Value[1]

			if cnt == 0 {
				
				if i == 0{
		
					coeffsRot0 := ctPreRot0.Value[0]
					coeffsRot1 := ctPreRot0.Value[1]

					r := p.RoundingNoise()

					for k := range tmp{
						bComplex[0].Mul(tmp[k][0], scale)
						bComplex[1].Mul(tmp[k][1], scale)
						bComplex.Add(bComplex, r[k])
						mul(coeffsRot0[k], bComplex, coeffsAcc0[k])
						mul(coeffsRot1[k], bComplex, coeffsAcc1[k])
					}

				}else{
					
					coeffsRot0 := ctPreRot[i]

					r := p.RoundingNoise()

					for k := range tmp{
						bComplex[0].Mul(tmp[k][0], scale)
						bComplex[1].Mul(tmp[k][1], scale)
						bComplex.Add(bComplex, r[k])
						mul(coeffsRot0[k], bComplex, coeffsAcc0[k])
					}
				}

			} else {

				if i == 0{
		
					coeffsRot0 := ctPreRot0.Value[0]
					coeffsRot1 := ctPreRot0.Value[1]

					r := p.RoundingNoise()

					for k := range tmp{
						bComplex[0].Mul(tmp[k][0], scale)
						bComplex[1].Mul(tmp[k][1], scale)
						bComplex.Add(bComplex, r[k])
						mul(coeffsRot0[k], bComplex, dComplex)
						coeffsAcc0[k].Add(coeffsAcc0[k], dComplex)
						mul(coeffsRot1[k], bComplex, dComplex)
						coeffsAcc1[k].Add(coeffsAcc1[k], dComplex)
					}

				}else{
					
					coeffsRot0 := ctPreRot[i]

					r := p.RoundingNoise()

					for k := range tmp{
						bComplex[0].Mul(tmp[k][0], scale)
						bComplex[1].Mul(tmp[k][1], scale)
						bComplex.Add(bComplex, r[k])
						mul(coeffsRot0[k], bComplex, dComplex)
						coeffsAcc0[k].Add(coeffsAcc0[k], dComplex)
					}
				}
			}

			cnt++
		}

		if j != 0 {
			
			m0 := acc.Value[0]
			e0 := p.KeySwitchingNoiseRaw(p.Level, p.RoundingNoise(), p.Sk[0])
			for k := range m0 {
				m0[k].Add(m0[k], e0[k])
			}

			utils.RotateSliceInPlace(m0, j)
		}

		for i := 0; i < 2; i++{
			m0, m1 := res.Value[i], acc.Value[i]
			for j := range m0{
				m0[j].Add(m0[j], m1[j])
			}
		}
	}

	res.ModDown()

	*p = *res

	return res
}
