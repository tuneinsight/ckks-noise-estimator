package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// KeyGenerator is an interface implementing the methods of the KeyGenerator.
type KeyGenerator interface {
	GenSecretKey() (sk *SecretKey)
	GenSecretKeyGaussian() (sk *SecretKey)
	GenSecretKeyWithDistrib(p float64) (sk *SecretKey)
	GenSecretKeyWithHammingWeight(hw int) (sk *SecretKey)
	GenPublicKey(sk *SecretKey) (pk *PublicKey)
	GenKeyPair() (sk *SecretKey, pk *PublicKey)
	GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey)
	GenSwitchingKey(skInput, skOutput *SecretKey) (newevakey *SwitchingKey)
	GenSwitchingKeyForGalois(galEl uint64, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForRotations(ks []int, inclueSwapRows bool, sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey)
	GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet)
	GenSwitchingKeysForRingSwap(skCKKS, skCI *SecretKey) (stdToci, ciToStd *SwitchingKey)
}

// KeyGenerator is a structure that stores the elements required to create new keys,
// as well as a memory buffer for intermediate values.
type keyGenerator struct {
	*skEncryptor
}

// NewKeyGenerator creates a new KeyGenerator, from which the secret and public keys, as well as the evaluation,
// rotation and switching keys can be generated.
func NewKeyGenerator(params Parameters) KeyGenerator {
	return &keyGenerator{
		skEncryptor: newSkEncryptor(params, NewSecretKey(params)),
	}
}

// GenSecretKey generates a new SecretKey with the distribution [1/3, 1/3, 1/3].
func (keygen *keyGenerator) GenSecretKey() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.ternarySampler)
}

// GenSecretKey generates a new SecretKey with the error distribution.
func (keygen *keyGenerator) GenSecretKeyGaussian() (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(keygen.gaussianSampler)
}

// GenSecretKeyWithDistrib generates a new SecretKey with the distribution [(p-1)/2, p, (p-1)/2].
func (keygen *keyGenerator) GenSecretKeyWithDistrib(p float64) (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(ring.NewTernarySampler(keygen.prng, keygen.params.RingQ(), p, false))
}

// GenSecretKeyWithHammingWeight generates a new SecretKey with exactly hw non-zero coefficients.
func (keygen *keyGenerator) GenSecretKeyWithHammingWeight(hw int) (sk *SecretKey) {
	return keygen.genSecretKeyFromSampler(ring.NewTernarySamplerWithHammingWeight(keygen.prng, keygen.params.RingQ(), hw, false))
}

// genSecretKeyFromSampler generates a new SecretKey sampled from the provided Sampler.
func (keygen *keyGenerator) genSecretKeyFromSampler(sampler ring.Sampler) (sk *SecretKey) {
	sk = new(SecretKey)
	ringQP := keygen.params.RingQP()
	sk.Value = ringQP.NewPoly()
	sampler.Read(sk.Value.Q)

	if levelP := sk.LevelP(); levelP > -1 {
		ringQP.ExtendBasisSmallNormAndCenter(sk.Value.Q, levelP, nil, sk.Value.P)
	}

	ringQP.NTT(sk.Value, sk.Value)
	ringQP.MForm(sk.Value, sk.Value)

	return
}

// GenPublicKey generates a new public key from the provided SecretKey.
func (keygen *keyGenerator) GenPublicKey(sk *SecretKey) (pk *PublicKey) {
	pk = NewPublicKey(keygen.params)
	pk.IsNTT = true
	pk.IsMontgomery = true
	keygen.WithKey(sk).EncryptZero(&CiphertextQP{Value: pk.Value, MetaData: pk.MetaData})
	return
}

// GenKeyPair generates a new SecretKey with distribution [1/3, 1/3, 1/3] and a corresponding public key.
func (keygen *keyGenerator) GenKeyPair() (sk *SecretKey, pk *PublicKey) {
	sk = keygen.GenSecretKey()
	return sk, keygen.GenPublicKey(sk)
}

// GenRelinKey generates a new EvaluationKey that will be used to relinearize Ciphertexts during multiplication.
func (keygen *keyGenerator) GenRelinearizationKey(sk *SecretKey, maxDegree int) (evk *RelinearizationKey) {

	levelQ := keygen.params.QCount() - 1
	levelP := keygen.params.PCount() - 1

	evk = new(RelinearizationKey)
	evk.Keys = make([]*SwitchingKey, maxDegree)
	for i := range evk.Keys {
		evk.Keys[i] = NewSwitchingKey(keygen.params, levelQ, levelP)
	}

	keygen.buffQP.Q.CopyValues(sk.Value.Q)
	ringQ := keygen.params.RingQ()
	for i := 0; i < maxDegree; i++ {
		ringQ.MulCoeffsMontgomery(keygen.buffQP.Q, sk.Value.Q, keygen.buffQP.Q)
		keygen.genSwitchingKey(keygen.buffQP.Q, sk, evk.Keys[i])
	}

	return
}

// GenRotationKeys generates a RotationKeySet from a list of galois element corresponding to the desired rotations
// See also GenRotationKeysForRotations.
func (keygen *keyGenerator) GenRotationKeys(galEls []uint64, sk *SecretKey) (rks *RotationKeySet) {
	rks = NewRotationKeySet(keygen.params, galEls)
	for _, galEl := range galEls {
		keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galEl), rks.Keys[galEl])
	}
	return rks
}

func (keygen *keyGenerator) GenSwitchingKeyForRotationBy(k int, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.MaxLevelQ(), keygen.params.MaxLevelP())
	galElInv := keygen.params.GaloisElementForColumnRotationBy(-int(k))
	keygen.genrotKey(sk.Value, galElInv, swk)
	return
}

// GenRotationKeysForRotations generates a RotationKeySet supporting left rotations by k positions for all k in ks.
// Negative k is equivalent to a right rotation by k positions
// If includeConjugate is true, the resulting set contains the conjugation key.
func (keygen *keyGenerator) GenRotationKeysForRotations(ks []int, includeConjugate bool, sk *SecretKey) (rks *RotationKeySet) {
	galEls := make([]uint64, len(ks), len(ks)+1)
	for i, k := range ks {
		galEls[i] = keygen.params.GaloisElementForColumnRotationBy(k)
	}
	if includeConjugate {
		galEls = append(galEls, keygen.params.GaloisElementForRowRotation())
	}
	return keygen.GenRotationKeys(galEls, sk)
}

func (keygen *keyGenerator) GenSwitchingKeyForRowRotation(sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.MaxLevelQ(), keygen.params.MaxLevelP())
	keygen.genrotKey(sk.Value, keygen.params.GaloisElementForRowRotation(), swk)
	return
}

func (keygen *keyGenerator) GenSwitchingKeyForGalois(galoisEl uint64, sk *SecretKey) (swk *SwitchingKey) {
	swk = NewSwitchingKey(keygen.params, keygen.params.MaxLevelQ(), keygen.params.MaxLevelP())
	keygen.genrotKey(sk.Value, keygen.params.InverseGaloisElement(galoisEl), swk)
	return
}

// GenRotationKeysForInnerSum generates a RotationKeySet supporting the InnerSum operation of the Evaluator
func (keygen *keyGenerator) GenRotationKeysForInnerSum(sk *SecretKey) (rks *RotationKeySet) {
	return keygen.GenRotationKeys(keygen.params.GaloisElementsForRowInnerSum(), sk)
}

func (keygen *keyGenerator) genrotKey(sk ringqp.Poly, galEl uint64, swk *SwitchingKey) {

	skIn := sk
	skOut := keygen.buffQP

	ringQ := keygen.params.RingQ()
	ringP := keygen.params.RingP()

	index := ringQ.PermuteNTTIndex(galEl)
	ringQ.PermuteNTTWithIndex(skIn.Q, index, skOut.Q)

	if ringP != nil {
		ringP.PermuteNTTWithIndex(skIn.P, index, skOut.P)
	}

	keygen.genSwitchingKey(skIn.Q, &SecretKey{Value: skOut}, swk)
}

// GenSwitchingKeysForRingSwap generates the necessary switching keys to switch from a standard ring to to a conjugate invariant ring and vice-versa.
func (keygen *keyGenerator) GenSwitchingKeysForRingSwap(skStd, skConjugateInvariant *SecretKey) (stdToci, ciToStd *SwitchingKey) {

	skCIMappedToStandard := &SecretKey{Value: keygen.buffQP}
	keygen.params.RingQ().UnfoldConjugateInvariantToStandard(skConjugateInvariant.Value.Q, skCIMappedToStandard.Value.Q)

	if keygen.params.PCount() != 0 {
		keygen.extendQ2P(keygen.params.MaxLevelP(), skCIMappedToStandard.Value.Q, keygen.buffQ[0], skCIMappedToStandard.Value.P)
	}

	return keygen.GenSwitchingKey(skStd, skCIMappedToStandard), keygen.GenSwitchingKey(skCIMappedToStandard, skStd)
}

// GenSwitchingKey generates a new key-switching key, that will re-encrypt a Ciphertext encrypted under the input key into the output key.
// If the ringDegree(skOutput) > ringDegree(skInput),  generates [-a*SkOut + w*P*skIn_{Y^{N/n}} + e, a] in X^{N}.
// If the ringDegree(skOutput) < ringDegree(skInput),  generates [-a*skOut_{Y^{N/n}} + w*P*skIn + e_{N}, a_{N}] in X^{N}.
// Else generates [-a*skOut + w*P*skIn + e, a] in X^{N}.
// The output switching key is always given in max(N, n) and in the moduli of the output switching key.
// When key-switching a Ciphertext from Y^{N/n} to X^{N}, the Ciphertext must first be mapped to X^{N}
// using SwitchCiphertextRingDegreeNTT(ctSmallDim, nil, ctLargeDim).
// When key-switching a Ciphertext from X^{N} to Y^{N/n}, the output of the key-switch is in still X^{N} and
// must be mapped Y^{N/n} using SwitchCiphertextRingDegreeNTT(ctLargeDim, ringQLargeDim, ctSmallDim).
func (keygen *keyGenerator) GenSwitchingKey(skInput, skOutput *SecretKey) (swk *SwitchingKey) {

	// N -> n (swk is to switch to a smaller dimension).
	if len(skInput.Value.Q.Coeffs[0]) > len(skOutput.Value.Q.Coeffs[0]) {

		levelP := skInput.LevelP()

		// Allocates the switching-key.
		swk = NewSwitchingKey(keygen.params, skOutput.Value.Q.Level(), levelP)

		// Maps the smaller key to the largest with Y = X^{N/n}.
		ring.MapSmallDimensionToLargerDimensionNTT(skOutput.Value.Q, keygen.buffQP.Q)

		// Extends the modulus P of skOutput to the one of skInput
		if levelP != -1 {
			keygen.extendQ2P(levelP, keygen.buffQP.Q, keygen.buffQ[0], keygen.buffQP.P)
		}

		keygen.genSwitchingKey(skInput.Value.Q, &SecretKey{Value: keygen.buffQP}, swk)

	} else { // N -> N or n -> N (swk switch to the same or a larger dimension)

		levelP := utils.MinInt(skOutput.LevelP(), keygen.params.MaxLevelP())

		// Allocates the switching-key.
		swk = NewSwitchingKey(keygen.params, skOutput.Value.Q.Level(), levelP)

		// Maps the smaller key to the largest dimension with Y = X^{N/n}.
		ring.MapSmallDimensionToLargerDimensionNTT(skInput.Value.Q, keygen.buffQ[0])

		// Extends the modulus of the input key to the one of the output key
		// if the former is smaller.
		if skInput.Value.Q.Level() < skOutput.Value.Q.Level() {

			ringQ := keygen.params.RingQ().AtLevel(0)

			// Switches out of the NTT and Montgomery domain.
			ringQ.INTT(keygen.buffQ[0], keygen.buffQP.Q)
			ringQ.IMForm(keygen.buffQP.Q, keygen.buffQP.Q)

			// Extends the RNS basis of the small norm polynomial.
			Qi := ringQ.ModuliChain()
			Q := Qi[0]
			QHalf := Q >> 1

			polQ := keygen.buffQP.Q
			polP := keygen.buffQ[0]
			var sign uint64
			N := ringQ.N()
			for j := 0; j < N; j++ {

				coeff := polQ.Coeffs[0][j]

				sign = 1
				if coeff > QHalf {
					coeff = Q - coeff
					sign = 0
				}

				for i := skInput.LevelQ() + 1; i < skOutput.LevelQ()+1; i++ {
					polP.Coeffs[i][j] = (coeff * sign) | (Qi[i]-coeff)*(sign^1)
				}
			}

			// Switches back to the NTT and Montgomery domain.
			for i := skInput.Value.Q.Level() + 1; i < skOutput.Value.Q.Level()+1; i++ {
				ringQ.SubRings[i].NTT(polP.Coeffs[i], polP.Coeffs[i])
				ringQ.SubRings[i].MForm(polP.Coeffs[i], polP.Coeffs[i])
			}
		}

		keygen.genSwitchingKey(keygen.buffQ[0], skOutput, swk)
	}

	return
}

func (keygen *keyGenerator) extendQ2P(levelP int, polQ, buff, polP *ring.Poly) {
	ringQ := keygen.params.RingQ().AtLevel(0)
	ringP := keygen.params.RingP().AtLevel(levelP)

	// Switches Q[0] out of the NTT and Montgomery domain.
	ringQ.INTT(polQ, buff)
	ringQ.IMForm(buff, buff)

	// Reconstruct P from Q
	Q := ringQ.SubRings[0].Modulus
	QHalf := Q >> 1

	P := ringP.ModuliChain()
	N := ringQ.N()

	var sign uint64
	for j := 0; j < N; j++ {

		coeff := buff.Coeffs[0][j]

		sign = 1
		if coeff > QHalf {
			coeff = Q - coeff
			sign = 0
		}

		for i := 0; i < levelP+1; i++ {
			polP.Coeffs[i][j] = (coeff * sign) | (P[i]-coeff)*(sign^1)
		}
	}

	ringP.NTT(polP, polP)
	ringP.MForm(polP, polP)
}

func (keygen *keyGenerator) genSwitchingKey(skIn *ring.Poly, skOut *SecretKey, swk *SwitchingKey) {

	enc := keygen.WithKey(skOut)
	// Samples an encryption of zero for each element of the switching-key.
	for i := 0; i < len(swk.Value); i++ {
		for j := 0; j < len(swk.Value[0]); j++ {
			enc.EncryptZero(&swk.Value[i][j])
		}
	}

	// Adds the plaintext (input-key) to the switching-key.
	AddPolyTimesGadgetVectorToGadgetCiphertext(skIn, []GadgetCiphertext{swk.GadgetCiphertext}, *keygen.params.RingQP(), keygen.params.Pow2Base(), keygen.buffQ[0])
}
