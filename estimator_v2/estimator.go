package estimator

type Parameters struct {
	Encoder
	LogN     int
	Sigma    float64
	LogScale int
	H        int
	Prec     uint
}

func NewParameters(LogN int, Sigma float64, H, LogScale int, prec uint) Parameters {
	return Parameters{
		LogN:     LogN,
		Sigma:    Sigma,
		H:        H,
		Prec:     prec,
		LogScale: LogScale,
		Encoder:  *NewEncoder(LogN, prec),
	}
}

func (p Parameters) MaxSlots() int {
	return 1 << p.LogMaxSlots()
}

func (p Parameters) LogMaxSlots() int {
	return p.LogN - 1
}
