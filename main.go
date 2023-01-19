package main


import(
	"fmt"
	"github.com/tuneinsight/ckks-bootstrapping-precision/estimator"
)

func main(){

	N := 65536
	H := 32768

	Q := []uint64{1152921504606584833, 35184372744193, 35184373006337, 35184368025601, 35184376545281}
	P := []uint64{2305843009211596801, 2305843009210023937}

	e := estimator.NewEstimator(N, H, Q, P)

	msg0 := 3.2 * (1<<45)
	msg1 := 3.2 * (1<<45)

	pt := estimator.NewCiphertextPk(estimator.NewPlaintext(msg0, nil, 4))
	ct := estimator.NewCiphertextPk(estimator.NewPlaintext(msg1, nil, 4))

	ct = e.Mul(ct, pt)

	fmt.Println(e.Std(ct))
}