package estimator

import(
	"github.com/tuneinsight/lattigo/v5/he"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/utils"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"
)

type LinearTransformation struct{
	*rlwe.MetaData
	LogBabyStepGianStepRatio int
	Value                    hefloat.Diagonals[*bignum.Complex]
}

// BSGSIndex returns the BSGSIndex of the target linear transformation.
func (lt LinearTransformation) BSGSIndex() (index map[int][]int, n1, n2 []int) {
	cols := 1<<lt.LogDimensions.Cols
	N1 := he.FindBestBSGSRatio(lt.Value.DiagonalsIndexList(), cols, lt.LogBabyStepGianStepRatio)

	return he.BSGSIndex(utils.GetKeys(lt.Value), cols, N1)
}

func (p *Element) EvaluateLinearTransformation(lt LinearTransformation) *Element {

	index, _, rotN2 := lt.BSGSIndex()

	ctPreRot := map[int]*Element{}

	for _, k := range rotN2 {
		if _, ok := ctPreRot[k]; k != 0 && !ok {
			ctPreRot[k] = p.CopyNew()
			ctPreRot[k].RotateHoisted(k)
		}
	}

	ctPreRot[0] = p.CopyNew().Mul(p, p.P)

	keys := utils.GetSortedKeys(index)

	acc := NewElement[*bignum.Complex](*p.Parameters, nil, 1, p.Scale)
	res := NewElement[*bignum.Complex](*p.Parameters, nil, 1, p.Scale)

	for _, j := range keys {

		var cnt int
		for _, i := range index[j] {

			pt := NewElement(*p.Parameters, lt.Value[i+j], 0, lt.Scale)
			pt.AddEncodingNoise()

			if cnt == 0{
				acc.Mul(ctPreRot[i], &pt)
			}else{
				acc.MulThenAdd(ctPreRot[i], &pt)
			}

			cnt++
		}

		if j != 0{
			acc.ModDown()
			acc.RotateHoisted(j)
		}

		res.Add(&res, &acc)
	}

	res.ModDown()

	*p = res

	return &res
}