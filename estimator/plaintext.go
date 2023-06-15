package estimator


import(
	"fmt"
	"time"
	"math/rand"
	"math/big"
	"github.com/tuneinsight/lattigo/v4/ckks"
	"github.com/tuneinsight/lattigo/v4/utils/bignum"
)


type PlaintextDistribution struct{
	MaxReal *big.Float
	MinReal *big.Float
	
	MaxImag *big.Float
	MinImag *big.Float
	
	MaxRing *big.Float
	MinRing *big.Float
}

func NewPlaintextDistribution(LogN, LogSlots, LogScale int, minR, maxR, minI, maxI float64, prec uint) (*PlaintextDistribution){

	params, _ := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN: LogN,
		LogQ: []int{60},
		LogScale: 45,
	})

	ecd := ckks.NewEncoder(params, prec)

	
	support := make([][]*bignum.Complex, 1024)

	r := rand.New(rand.NewSource(time.Now().UnixNano()))
	_=r

	for i := range support{

		vec := make([]*bignum.Complex, 1<<LogSlots)
		support[i] = make([]*bignum.Complex, 1<<LogSlots)

		for j := range vec{
			vec[j] = &bignum.Complex{
				new(big.Float).SetPrec(prec).SetFloat64((maxR-minR)*r.Float64()+minR),
				new(big.Float).SetPrec(prec).SetFloat64((maxI-minI)*r.Float64()+minI),
			}
		}
	
		ecd.IFFT(vec, LogSlots)

		support[i] = vec
	}


	// Transpose

	supportT := make([][]*bignum.Complex, len(support[0]))

	for j := range supportT{
		supportT[j] = make([]*bignum.Complex, len(support))
	}

	for i := range support{
		for j := range support[i]{
			supportT[j][i] = support[i][j]
		}
	}

	fmt.Println(supportT[0])

	for i := range supportT[:8]{
		a, b := Mean(supportT[i])
		fmt.Println(i, a, b)
	}

	/*

	fmt.Println("===========")
	fmt.Println(vec[0])
	fmt.Println(vec[1])
	fmt.Println(vec[2])
	fmt.Println(vec[3])

	minRing := new(big.Float).Set(vec[0][0])
	maxRing := new(big.Float).Set(vec[0][0])

	for i := range vec{
		if vec[i][0].Cmp(minRing) == -1{
			minRing.Set(vec[i][0])
		}

		if vec[i][1].Cmp(minRing) == -1{
			minRing.Set(vec[i][1])
		}

		if vec[i][0].Cmp(maxRing) == 1{
			maxRing.Set(vec[i][0])
		}

		if vec[i][1].Cmp(maxRing) == 1{
			maxRing.Set(vec[i][0])
		}
	}

	return &PlaintextDistribution{
		MinReal: new(big.Float).SetPrec(prec).SetFloat64(minR),
		MaxReal: new(big.Float).SetPrec(prec).SetFloat64(maxR),
		MinImag: new(big.Float).SetPrec(prec).SetFloat64(minI),
		MaxImag: new(big.Float).SetPrec(prec).SetFloat64(maxI),

		MinRing: minRing,
		MaxRing: maxRing,
	}
	*/

	return nil
}

func Mean(vec []*bignum.Complex) (mR, mI *big.Float){
	mR = new(big.Float)
	mI = new(big.Float)

	for _, c := range vec {
		mR.Add(mR, c[0])
		mI.Add(mI, c[1])
	}

	mR.Quo(mR, new(big.Float).SetInt64(int64(len(vec))))
	mI.Quo(mI, new(big.Float).SetInt64(int64(len(vec))))

	return
}

// StandardDeviation computes the scaled standard deviation of the input vector.
func StandardDeviation(vec []*bignum.Complex) (stdR, stdI float64) {

	meanR, meanI := Mean(vec)

	errR := new(big.Float)
	errI := new(big.Float)
	tmp := new(big.Float)
	for _, c := range vec {
		tmp.Sub(c[0], meanR)
		tmp.Mul(tmp, tmp)
		errR.Add(errR, tmp)

		tmp.Sub(c[1], meanI)
		tmp.Mul(tmp, tmp)
		errI.Add(errI, tmp)
	}

	errR.Quo(errR, new(big.Float).SetInt64(int64(len(vec)-1)))
	errR.Sqrt(errR)

	errI.Quo(errI, new(big.Float).SetInt64(int64(len(vec)-1)))
	errI.Sqrt(errI)

	stdR, _ = errR.Float64()
	stdI, _ = errI.Float64()


	return
}