package estimator

import(
	"math"
	"github.com/tuneinsight/lattigo/v5/schemes/ckks"
)

type Stats struct{
	ckks.PrecisionStats
	N float64
}

func (s *Stats) Add(a ckks.PrecisionStats) {
	s.MINLog2Prec.Real += a.MINLog2Prec.Real
	s.MINLog2Prec.Imag += a.MINLog2Prec.Imag
	s.MINLog2Prec.L2 += a.MINLog2Prec.L2

	s.MAXLog2Prec.Real += a.MAXLog2Prec.Real
	s.MAXLog2Prec.Imag += a.MAXLog2Prec.Imag
	s.MAXLog2Prec.L2 += a.MAXLog2Prec.L2

	s.AVGLog2Prec.Real += a.AVGLog2Prec.Real
	s.AVGLog2Prec.Imag += a.AVGLog2Prec.Imag
	s.AVGLog2Prec.L2 += a.AVGLog2Prec.L2

	s.MEDLog2Prec.Real += a.MEDLog2Prec.Real
	s.MEDLog2Prec.Imag += a.MEDLog2Prec.Imag
	s.MEDLog2Prec.L2 += a.MEDLog2Prec.L2

	s.STDLog2Prec.Real += a.STDLog2Prec.Real
	s.STDLog2Prec.Imag += a.STDLog2Prec.Imag
	s.STDLog2Prec.L2 += a.STDLog2Prec.L2

	s.MINLog2Err.Real += a.MINLog2Err.Real
	s.MINLog2Err.Imag += a.MINLog2Err.Imag
	s.MINLog2Err.L2 += a.MINLog2Err.L2

	s.MAXLog2Err.Real += a.MAXLog2Err.Real
	s.MAXLog2Err.Imag += a.MAXLog2Err.Imag
	s.MAXLog2Err.L2 += a.MAXLog2Err.L2

	s.AVGLog2Err.Real += a.AVGLog2Err.Real
	s.AVGLog2Err.Imag += a.AVGLog2Err.Imag
	s.AVGLog2Err.L2 += a.AVGLog2Err.L2

	s.MEDLog2Err.Real += a.MEDLog2Err.Real
	s.MEDLog2Err.Imag += a.MEDLog2Err.Imag
	s.MEDLog2Err.L2 += a.MEDLog2Err.L2

	s.STDLog2Err.Real += a.STDLog2Err.Real*a.STDLog2Err.Real
	s.STDLog2Err.Imag += a.STDLog2Err.Imag*a.STDLog2Err.Imag
	s.STDLog2Err.L2 += a.STDLog2Err.L2*a.STDLog2Err.L2

	s.N++
}

func (s *Stats) Finalize(){
	s.MINLog2Prec.Real /= s.N
	s.MINLog2Prec.Imag /= s.N
	s.MINLog2Prec.L2 /= s.N

	s.MAXLog2Prec.Real /= s.N
	s.MAXLog2Prec.Imag /= s.N
	s.MAXLog2Prec.L2 /= s.N

	s.AVGLog2Prec.Real /= s.N
	s.AVGLog2Prec.Imag /= s.N
	s.AVGLog2Prec.L2 /= s.N

	s.MEDLog2Prec.Real /= s.N
	s.MEDLog2Prec.Imag /= s.N
	s.MEDLog2Prec.L2 /= s.N

	s.STDLog2Prec.Real /= s.N
	s.STDLog2Prec.Imag /= s.N
	s.STDLog2Prec.L2 /= s.N

	s.MINLog2Err.Real /= s.N
	s.MINLog2Err.Imag /= s.N
	s.MINLog2Err.L2 /= s.N

	s.MAXLog2Err.Real /= s.N
	s.MAXLog2Err.Imag /= s.N
	s.MAXLog2Err.L2 /= s.N

	s.AVGLog2Err.Real /= s.N
	s.AVGLog2Err.Imag /= s.N
	s.AVGLog2Err.L2 /= s.N

	s.MEDLog2Err.Real /= s.N
	s.MEDLog2Err.Imag /= s.N
	s.MEDLog2Err.L2 /= s.N

	s.STDLog2Err.Real = math.Sqrt(s.STDLog2Err.Real/s.N)
	s.STDLog2Err.Imag = math.Sqrt(s.STDLog2Err.Imag/s.N)
	s.STDLog2Err.L2 = math.Sqrt(s.STDLog2Err.L2/s.N)
}