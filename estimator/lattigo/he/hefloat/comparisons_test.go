package hefloat_test

import (
	"math/big"
	"testing"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils/bignum"

	"github.com/stretchr/testify/require"
)

func TestComparisons(t *testing.T) {

	paramsLiteral := testInsecurePrec90

	for _, ringType := range []ring.Type{ring.Standard, ring.ConjugateInvariant} {

		paramsLiteral.RingType = ringType

		if testing.Short() {
			paramsLiteral.LogN = 10
		}

		params, err := hefloat.NewParametersFromLiteral(paramsLiteral)
		require.NoError(t, err)

		var tc *testContext
		if tc, err = genTestParams(params); err != nil {
			t.Fatal(err)
		}

		enc := tc.encryptorSk
		sk := tc.sk
		ecd := tc.encoder
		dec := tc.decryptor
		kgen := tc.kgen

		btp := bootstrapping.NewSecretKeyBootstrapper(params, sk)

		var galKeys []*rlwe.GaloisKey
		if params.RingType() == ring.Standard {
			galKeys = append(galKeys, kgen.GenGaloisKeyNew(params.GaloisElementForComplexConjugation(), sk))
		}

		eval := tc.evaluator.WithKey(rlwe.NewMemEvaluationKeySet(kgen.GenRelinearizationKeyNew(sk), galKeys...))

		polys := hefloat.NewMinimaxCompositePolynomial(hefloat.DefaultMinimaxCompositePolynomialForSign)

		CmpEval := hefloat.NewComparisonEvaluator(params, eval, btp, polys)

		t.Run(GetTestName(params, "Sign"), func(t *testing.T) {

			values, _, ct := newTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			var sign *rlwe.Ciphertext
			sign, err = CmpEval.Sign(ct)
			require.NoError(t, err)

			have := make([]*bignum.Complex, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(sign), have))

			want := make([]*bignum.Complex, params.MaxSlots())

			for i := range have {
				want[i] = polys.Evaluate(values[i])
			}

			hefloat.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Step"), func(t *testing.T) {

			values, _, ct := newTestVectors(tc, enc, complex(-1, 0), complex(1, 0), t)

			var step *rlwe.Ciphertext
			step, err = CmpEval.Step(ct)
			require.NoError(t, err)

			have := make([]*bignum.Complex, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(step), have))

			want := make([]*bignum.Complex, params.MaxSlots())

			half := new(big.Float).SetFloat64(0.5)

			for i := range have {
				want[i] = polys.Evaluate(values[i])
				want[i][0].Mul(want[i][0], half)
				want[i][0].Add(want[i][0], half)
			}

			hefloat.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Max"), func(t *testing.T) {

			values0, _, ct0 := newTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)
			values1, _, ct1 := newTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)

			var max *rlwe.Ciphertext
			max, err = CmpEval.Max(ct0, ct1)
			require.NoError(t, err)

			have := make([]*bignum.Complex, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(max), have))

			want := make([]*bignum.Complex, params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == -1 {
					want[i] = &bignum.Complex{values1[i][0], new(big.Float)}
				} else {
					want[i] = &bignum.Complex{values0[i][0], new(big.Float)}
				}
			}

			hefloat.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})

		t.Run(GetTestName(params, "Min"), func(t *testing.T) {

			values0, _, ct0 := newTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)
			values1, _, ct1 := newTestVectors(tc, enc, complex(-0.5, 0), complex(0.5, 0), t)

			var max *rlwe.Ciphertext
			max, err = CmpEval.Min(ct0, ct1)
			require.NoError(t, err)

			have := make([]*bignum.Complex, params.MaxSlots())

			require.NoError(t, ecd.Decode(dec.DecryptNew(max), have))

			want := make([]*bignum.Complex, params.MaxSlots())

			for i := range have {

				if values0[i][0].Cmp(values1[i][0]) == 1 {
					want[i] = &bignum.Complex{values1[i][0], new(big.Float)}
				} else {
					want[i] = &bignum.Complex{values0[i][0], new(big.Float)}
				}
			}

			hefloat.VerifyTestVectors(params, ecd, nil, want, have, params.LogDefaultScale(), 0, *printPrecisionStats, t)
		})
	}
}
