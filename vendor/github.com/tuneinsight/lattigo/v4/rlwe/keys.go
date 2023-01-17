package rlwe

import (
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
)

// SecretKey is a type for generic RLWE secret keys.
// The Value field stores the polynomial in NTT and Montgomery form.
type SecretKey struct {
	Value ringqp.Poly
}

// PublicKey is a type for generic RLWE public keys.
// The Value field stores the polynomials in NTT and Montgomery form.
type PublicKey struct {
	CiphertextQP
}

// SwitchingKey is a type for generic RLWE public switching keys.
// The Value field stores the polynomials in NTT and Montgomery form.
type SwitchingKey struct {
	GadgetCiphertext
}

// RelinearizationKey is a type for generic RLWE public relinearization keys. It stores a slice with a
// switching key per relinearizable degree. The switching key at index i is used to relinearize a degree
// i+2 ciphertexts back to a degree i + 1 one.
type RelinearizationKey struct {
	Keys []*SwitchingKey
}

// RotationKeySet is a type for storing generic RLWE public rotation keys. It stores a map indexed by the
// galois element defining the automorphism.
type RotationKeySet struct {
	Keys map[uint64]*SwitchingKey
}

// EvaluationKey is a type for storing generic RLWE public evaluation keys. An evaluation key is a union
// of a relinearization key and a set of rotation keys.
type EvaluationKey struct {
	Rlk  *RelinearizationKey
	Rtks *RotationKeySet
}

// NewSecretKey generates a new SecretKey with zero values.
func NewSecretKey(params Parameters) *SecretKey {
	return &SecretKey{Value: params.RingQP().NewPoly()}
}

// LevelQ returns the level of the modulus Q of the target.
func (sk *SecretKey) LevelQ() int {
	return sk.Value.Q.Level()
}

// LevelP returns the level of the modulus P of the target.
// Returns -1 if P is absent.
func (sk *SecretKey) LevelP() int {
	if sk.Value.P != nil {
		return sk.Value.P.Level()
	}

	return -1
}

// NewPublicKey returns a new PublicKey with zero values.
func NewPublicKey(params Parameters) (pk *PublicKey) {
	return &PublicKey{CiphertextQP{Value: [2]ringqp.Poly{params.RingQP().NewPoly(), params.RingQP().NewPoly()}, MetaData: MetaData{IsNTT: true, IsMontgomery: true}}}
}

// LevelQ returns the level of the modulus Q of the target.
func (pk *PublicKey) LevelQ() int {
	return pk.Value[0].Q.Level()
}

// LevelP returns the level of the modulus P of the target.
// Returns -1 if P is absent.
func (pk *PublicKey) LevelP() int {
	if pk.Value[0].P != nil {
		return pk.Value[0].P.Level()
	}

	return -1
}

// Equals checks two PublicKey struct for equality.
func (pk *PublicKey) Equals(other *PublicKey) bool {
	if pk == other {
		return true
	}
	return pk.Value[0].Equals(other.Value[0]) && pk.Value[1].Equals(other.Value[1])
}

// NewRotationKeySet returns a new RotationKeySet with pre-allocated switching keys for each distinct galoisElement value.
func NewRotationKeySet(params Parameters, galoisElement []uint64) (rotKey *RotationKeySet) {
	rotKey = new(RotationKeySet)
	rotKey.Keys = make(map[uint64]*SwitchingKey, len(galoisElement))
	for _, galEl := range galoisElement {
		rotKey.Keys[galEl] = NewSwitchingKey(params, params.MaxLevelQ(), params.MaxLevelP())
	}
	return
}

// GetRotationKey return the rotation key for the given galois element or nil if such key is not in the set. The
// second argument is true  iff the first one is non-nil.
func (rtks *RotationKeySet) GetRotationKey(galoisEl uint64) (*SwitchingKey, bool) {
	if rtks.Keys == nil {
		return nil, false
	}
	rotKey, inSet := rtks.Keys[galoisEl]
	return rotKey, inSet
}

// NewSwitchingKey returns a new public switching key with pre-allocated zero-value
func NewSwitchingKey(params Parameters, levelQ, levelP int) *SwitchingKey {
	return &SwitchingKey{GadgetCiphertext: *NewGadgetCiphertext(
		params,
		levelQ,
		levelP,
		params.DecompRNS(levelQ, levelP),
		params.DecompPw2(levelQ, levelP),
	)}
}

// Equals checks two SwitchingKeys for equality.
func (swk *SwitchingKey) Equals(other *SwitchingKey) bool {
	return swk.GadgetCiphertext.Equals(&other.GadgetCiphertext)
}

// CopyNew creates a deep copy of the target SwitchingKey and returns it.
func (swk *SwitchingKey) CopyNew() *SwitchingKey {
	return &SwitchingKey{GadgetCiphertext: *swk.GadgetCiphertext.CopyNew()}
}

// NewRelinearizationKey creates a new EvaluationKey with zero values.
func NewRelinearizationKey(params Parameters, maxRelinDegree int) (evakey *RelinearizationKey) {
	evakey = new(RelinearizationKey)
	evakey.Keys = make([]*SwitchingKey, maxRelinDegree)
	for d := 0; d < maxRelinDegree; d++ {
		evakey.Keys[d] = NewSwitchingKey(params, params.MaxLevelQ(), params.MaxLevelP())
	}

	return
}

// CopyNew creates a deep copy of the receiver secret key and returns it.
func (sk *SecretKey) CopyNew() *SecretKey {
	if sk == nil {
		return nil
	}
	return &SecretKey{sk.Value.CopyNew()}
}

// CopyNew creates a deep copy of the receiver PublicKey and returns it.
func (pk *PublicKey) CopyNew() *PublicKey {
	if pk == nil {
		return nil
	}
	return &PublicKey{*pk.CiphertextQP.CopyNew()}
}

// Equals checks two RelinearizationKeys for equality.
func (rlk *RelinearizationKey) Equals(other *RelinearizationKey) bool {
	if rlk == other {
		return true
	}
	if (rlk == nil) != (other == nil) {
		return false
	}
	if len(rlk.Keys) != len(other.Keys) {
		return false
	}
	for i := range rlk.Keys {
		if !rlk.Keys[i].Equals(other.Keys[i]) {
			return false
		}
	}
	return true
}

// CopyNew creates a deep copy of the receiver RelinearizationKey and returns it.
func (rlk *RelinearizationKey) CopyNew() *RelinearizationKey {
	if rlk == nil || len(rlk.Keys) == 0 {
		return nil
	}
	rlkb := &RelinearizationKey{Keys: make([]*SwitchingKey, len(rlk.Keys))}
	for i, swk := range rlk.Keys {
		rlkb.Keys[i] = swk.CopyNew()
	}
	return rlkb
}

// Equals checks to RotationKeySets for equality.
func (rtks *RotationKeySet) Equals(other *RotationKeySet) bool {
	if rtks == other {
		return true
	}
	if (rtks == nil) || (other == nil) {
		return false
	}
	if len(rtks.Keys) != len(other.Keys) {
		return false
	}
	for galEl, otherKey := range other.Keys {
		if key, inSet := rtks.Keys[galEl]; !inSet || !otherKey.GadgetCiphertext.Equals(&key.GadgetCiphertext) {
			return false
		}
	}
	return true
}

// Includes checks whether the receiver RotationKeySet includes the given other RotationKeySet.
func (rtks *RotationKeySet) Includes(other *RotationKeySet) bool {
	if (rtks == nil) || (other == nil) {
		return false
	}
	for galEl := range other.Keys {
		if _, inSet := rtks.Keys[galEl]; !inSet {
			return false
		}
	}
	return true
}
