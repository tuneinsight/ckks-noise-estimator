package rlwe

import (
	"fmt"
	"reflect"

	"github.com/tuneinsight/lattigo/v4/ring"
	"github.com/tuneinsight/lattigo/v4/rlwe/ringqp"
	"github.com/tuneinsight/lattigo/v4/utils"
)

// Encryptor a generic RLWE encryption interface.
type Encryptor interface {
	Encrypt(pt *Plaintext, ct interface{})
	EncryptZero(ct interface{})

	EncryptZeroNew(level int) (ct *Ciphertext)
	EncryptNew(pt *Plaintext) (ct *Ciphertext)

	ShallowCopy() Encryptor
	WithKey(key interface{}) Encryptor
}

// PRNGEncryptor is an interface for encrypting RLWE ciphertexts from a secret-key and
// a pre-determined PRNG. An Encryptor constructed from a secret-key complies to this
// interface.
type PRNGEncryptor interface {
	Encryptor
	WithPRNG(prng utils.PRNG) PRNGEncryptor
}

type encryptorBase struct {
	params Parameters
	*encryptorBuffers

	prng            utils.PRNG
	gaussianSampler *ring.GaussianSampler
	ternarySampler  *ring.TernarySampler
	basisextender   *ring.BasisExtender
}

type pkEncryptor struct {
	*encryptorBase
	pk *PublicKey
}

type skEncryptor struct {
	encryptorBase
	sk *SecretKey

	uniformSampler ringqp.UniformSampler
}

// NewEncryptor creates a new Encryptor
// Accepts either a secret-key or a public-key.
func NewEncryptor(params Parameters, key interface{}) Encryptor {
	switch key := key.(type) {
	case *PublicKey, PublicKey:
		return newPkEncryptor(params, key)
	case *SecretKey, SecretKey:
		return newSkEncryptor(params, key)
	default:
		panic("cannot NewEncryptor: key must be either *rlwe.PublicKey or *rlwe.SecretKey")
	}
}

// NewPRNGEncryptor creates a new PRNGEncryptor instance.
func NewPRNGEncryptor(params Parameters, key *SecretKey) PRNGEncryptor {
	return newSkEncryptor(params, key)
}

func newEncryptorBase(params Parameters) *encryptorBase {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}

	var bc *ring.BasisExtender
	if params.PCount() != 0 {
		bc = ring.NewBasisExtender(params.RingQ(), params.RingP())
	}

	return &encryptorBase{
		params:           params,
		prng:             prng,
		gaussianSampler:  ring.NewGaussianSampler(prng, params.RingQ(), params.Sigma(), int(6*params.Sigma())),
		ternarySampler:   ring.NewTernarySamplerWithHammingWeight(prng, params.ringQ, params.h, false),
		encryptorBuffers: newEncryptorBuffers(params),
		basisextender:    bc,
	}
}

func newSkEncryptor(params Parameters, key interface{}) (enc *skEncryptor) {

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(fmt.Errorf("cannot newSkEncryptor: could not create PRNG for symmetric encryptor: %s", err))
	}

	enc = &skEncryptor{*newEncryptorBase(params), nil, ringqp.NewUniformSampler(prng, *params.RingQP())}
	if enc.sk, err = enc.checkSk(key); err != nil {
		panic(err)
	}

	return enc
}

func newPkEncryptor(params Parameters, key interface{}) (enc *pkEncryptor) {
	var err error
	enc = &pkEncryptor{newEncryptorBase(params), nil}
	enc.pk, err = enc.checkPk(key)
	if err != nil {
		panic(err)
	}
	return enc
}

type encryptorBuffers struct {
	buffQ  [2]*ring.Poly
	buffP  [3]*ring.Poly
	buffQP ringqp.Poly
}

func newEncryptorBuffers(params Parameters) *encryptorBuffers {

	ringQ := params.RingQ()
	ringP := params.RingP()

	var buffP [3]*ring.Poly
	if params.PCount() != 0 {
		buffP = [3]*ring.Poly{ringP.NewPoly(), ringP.NewPoly(), ringP.NewPoly()}
	}

	return &encryptorBuffers{
		buffQ:  [2]*ring.Poly{ringQ.NewPoly(), ringQ.NewPoly()},
		buffP:  buffP,
		buffQP: params.RingQP().NewPoly(),
	}
}

// Encrypt encrypts the input plaintext using the stored public-key and writes the result on ct.
// The encryption procedure first samples a new encryption of zero under the public-key and
// then adds the Plaintext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
// If a Plaintext is given, then the output Ciphertext MetaData will match the Plaintext MetaData.
func (enc *pkEncryptor) Encrypt(pt *Plaintext, ct interface{}) {

	if pt == nil {
		enc.EncryptZero(ct)
	} else {
		switch ct := ct.(type) {
		case *Ciphertext:

			ct.MetaData = pt.MetaData

			level := utils.MinInt(pt.Level(), ct.Level())

			ct.Resize(ct.Degree(), level)

			enc.EncryptZero(ct)

			enc.params.RingQ().AtLevel(level).Add(ct.Value[0], pt.Value, ct.Value[0])

		default:
			panic(fmt.Sprintf("cannot Encrypt: input ciphertext type %s is not supported", reflect.TypeOf(ct)))
		}
	}
}

// EncryptNew encrypts the input plaintext using the stored public-key and returns the result on a new Ciphertext.
// The encryption procedure first samples a new encryption of zero under the public-key and
// then adds the Plaintext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// If a Plaintext is given, then the output ciphertext MetaData will match the Plaintext MetaData.
func (enc *pkEncryptor) EncryptNew(pt *Plaintext) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, pt.Level())
	enc.Encrypt(pt, ct)
	return
}

// EncryptZeroNew generates an encryption of zero under the stored public-key and returns it on a new Ciphertext.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc *pkEncryptor) EncryptZeroNew(level int) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, level)
	enc.EncryptZero(ct)
	return
}

// EncryptZero generates an encryption of zero under the stored public-key and writes the result on ct.
// The encryption procedure depends on the parameters: If the auxiliary modulus P is defined, the
// encryption of zero is sampled in QP before being rescaled by P; otherwise, it is directly sampled in Q.
// The method accepts only *rlwe.Ciphertext as input.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc *pkEncryptor) EncryptZero(ct interface{}) {
	switch ct := ct.(type) {
	case *Ciphertext:
		if enc.params.PCount() > 0 {
			enc.encryptZero(ct)
		} else {
			enc.encryptZeroNoP(ct)
		}
	default:
		panic(fmt.Sprintf("cannot Encrypt: input ciphertext type %s is not supported", reflect.TypeOf(ct)))
	}
}

func (enc *pkEncryptor) encryptZero(ct *Ciphertext) {

	levelQ := ct.Level()
	levelP := 0

	ringQP := enc.params.RingQP().AtLevel(levelQ, levelP)

	buffQ0 := enc.buffQ[0]
	buffP0 := enc.buffP[0]
	buffP1 := enc.buffP[1]
	buffP2 := enc.buffP[2]

	u := ringqp.Poly{Q: buffQ0, P: buffP2}

	// We sample a RLWE instance (encryption of zero) over the extended ring (ciphertext ring + special prime)
	enc.ternarySampler.AtLevel(levelQ).Read(u.Q)
	ringQP.ExtendBasisSmallNormAndCenter(u.Q, levelP, nil, u.P)

	// (#Q + #P) NTT
	ringQP.NTT(u, u)

	ct0QP := ringqp.Poly{Q: ct.Value[0], P: buffP0}
	ct1QP := ringqp.Poly{Q: ct.Value[1], P: buffP1}

	// ct0 = u*pk0
	// ct1 = u*pk1
	ringQP.MulCoeffsMontgomery(u, enc.pk.Value[0], ct0QP)
	ringQP.MulCoeffsMontgomery(u, enc.pk.Value[1], ct1QP)

	// 2*(#Q + #P) NTT
	ringQP.INTT(ct0QP, ct0QP)
	ringQP.INTT(ct1QP, ct1QP)

	e := ringqp.Poly{Q: buffQ0, P: buffP2}

	enc.gaussianSampler.AtLevel(levelQ).Read(e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.Add(ct0QP, e, ct0QP)

	enc.gaussianSampler.AtLevel(levelQ).Read(e.Q)
	ringQP.ExtendBasisSmallNormAndCenter(e.Q, levelP, nil, e.P)
	ringQP.Add(ct1QP, e, ct1QP)

	// ct0 = (u*pk0 + e0)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct0QP.Q, ct0QP.P, ct.Value[0])

	// ct1 = (u*pk1 + e1)/P
	enc.basisextender.ModDownQPtoQ(levelQ, levelP, ct1QP.Q, ct1QP.P, ct.Value[1])

	if ct.IsNTT {
		ringQP.RingQ.NTT(ct.Value[0], ct.Value[0])
		ringQP.RingQ.NTT(ct.Value[1], ct.Value[1])
	}
}

func (enc *pkEncryptor) encryptZeroNoP(ct *Ciphertext) {

	levelQ := ct.Level()

	ringQ := enc.params.RingQ().AtLevel(levelQ)

	buffQ0 := enc.buffQ[0]

	enc.ternarySampler.AtLevel(levelQ).Read(buffQ0)
	ringQ.NTT(buffQ0, buffQ0)

	c0, c1 := ct.Value[0], ct.Value[1]

	// ct0 = NTT(u*pk0)
	ringQ.MulCoeffsMontgomery(buffQ0, enc.pk.Value[0].Q, c0)
	// ct1 = NTT(u*pk1)
	ringQ.MulCoeffsMontgomery(buffQ0, enc.pk.Value[1].Q, c1)

	// c0
	if ct.IsNTT {
		enc.gaussianSampler.AtLevel(levelQ).Read(buffQ0)
		ringQ.NTT(buffQ0, buffQ0)
		ringQ.Add(c0, buffQ0, c0)
	} else {
		ringQ.INTT(c0, c0)
		enc.gaussianSampler.AtLevel(levelQ).ReadAndAdd(c0)
	}

	// c1
	if ct.IsNTT {
		enc.gaussianSampler.AtLevel(levelQ).Read(buffQ0)
		ringQ.NTT(buffQ0, buffQ0)
		ringQ.Add(c1, buffQ0, c1)

	} else {
		ringQ.INTT(c1, c1)
		enc.gaussianSampler.AtLevel(levelQ).ReadAndAdd(c1)
	}
}

// Encrypt encrypts the input plaintext using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will panic otherwise.
// If a plaintext is given, the encryptor only accepts *rlwe.Ciphertext, and the generated Ciphertext
// MetaData will match the given Plaintext MetaData.
func (enc *skEncryptor) Encrypt(pt *Plaintext, ct interface{}) {
	if pt == nil {
		enc.EncryptZero(ct)
	} else {
		switch ct := ct.(type) {
		case *Ciphertext:
			ct.MetaData = pt.MetaData
			level := utils.MinInt(pt.Level(), ct.Level())
			ct.Resize(ct.Degree(), level)
			enc.EncryptZero(ct)
			enc.params.RingQ().AtLevel(level).Add(ct.Value[0], pt.Value, ct.Value[0])
		default:
			panic(fmt.Sprintf("cannot Encrypt: input ciphertext type %T is not supported", ct))
		}
	}
}

// Encrypt encrypts the input plaintext using the stored secret-key and returns the result on a new Ciphertext.
// MetaData will match the given Plaintext MetaData.
func (enc *skEncryptor) EncryptNew(pt *Plaintext) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, pt.Level())
	enc.Encrypt(pt, ct)
	return
}

// EncryptZero generates an encryption of zero using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will panic otherwise.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc *skEncryptor) EncryptZero(ct interface{}) {
	switch ct := ct.(type) {
	case *Ciphertext:

		var c1 *ring.Poly
		if ct.Degree() == 1 {
			c1 = ct.Value[1]
		} else {
			c1 = enc.buffQ[1]
		}

		enc.uniformSampler.AtLevel(ct.Level(), -1).Read(ringqp.Poly{Q: c1})
		enc.encryptZero(ct, c1)
	case *CiphertextQP:
		enc.encryptZeroQP(*ct)
	default:
		panic(fmt.Sprintf("cannot EncryptZero: input ciphertext type %T is not supported", ct))
	}
}

// EncryptZeroNew generates an encryption of zero using the stored secret-key and writes the result on ct.
// The method accepts only *rlwe.Ciphertext or *rgsw.Ciphertext as input and will panic otherwise.
// The zero encryption is generated according to the given Ciphertext MetaData.
func (enc *skEncryptor) EncryptZeroNew(level int) (ct *Ciphertext) {
	ct = NewCiphertext(enc.params, 1, level)
	enc.EncryptZero(ct)
	return
}

func (enc *skEncryptor) encryptZero(ct *Ciphertext, c1 *ring.Poly) {

	levelQ := ct.Level()

	ringQ := enc.params.RingQ().AtLevel(levelQ)

	c0 := ct.Value[0]

	ringQ.MulCoeffsMontgomery(c1, enc.sk.Value.Q, c0) // c0 = NTT(sc1)
	ringQ.Neg(c0, c0)                                 // c0 = NTT(-sc1)

	if ct.IsNTT {
		enc.gaussianSampler.AtLevel(levelQ).Read(enc.buffQ[0]) // e
		ringQ.NTT(enc.buffQ[0], enc.buffQ[0])                  // NTT(e)
		ringQ.Add(c0, enc.buffQ[0], c0)                        // c0 = NTT(-sc1 + e)
	} else {
		ringQ.INTT(c0, c0) // c0 = -sc1
		if ct.Degree() == 1 {
			ringQ.INTT(c1, c1) // c1 = c1
		}

		enc.gaussianSampler.AtLevel(levelQ).ReadAndAdd(c0) // c0 = -sc1 + e
	}
}

// EncryptZeroSeeded generates en encryption of zero under sk.
// levelQ : level of the modulus Q
// levelP : level of the modulus P
// sk     : secret key
// sampler: uniform sampler; if `sampler` is nil, then the internal sampler will be used.
// montgomery: returns the result in the Montgomery domain.
func (enc *skEncryptor) encryptZeroQP(ct CiphertextQP) {

	c0, c1 := ct.Value[0], ct.Value[1]

	levelQ, levelP := c0.LevelQ(), c1.LevelP()
	ringQP := enc.params.RingQP().AtLevel(levelQ, levelP)

	// ct = (e, 0)
	enc.gaussianSampler.AtLevel(levelQ).Read(c0.Q)
	if levelP != -1 {
		ringQP.ExtendBasisSmallNormAndCenter(c0.Q, levelP, nil, c0.P)
	}

	ringQP.NTT(c0, c0)
	// ct[1] is assumed to be sampled in of the Montgomery domain,
	// thus -as will also be in the Montgomery domain (s is by default), therefore 'e'
	// must be switched to the Montgomery domain.
	ringQP.MForm(c0, c0)

	// ct = (e, a)
	enc.uniformSampler.AtLevel(levelQ, levelP).Read(c1)

	// (-a*sk + e, a)
	ringQP.MulCoeffsMontgomeryThenSub(c1, enc.sk.Value, c0)

	if !ct.IsNTT {
		ringQP.INTT(c0, c0)
		ringQP.INTT(c1, c1)
	}
}

// ShallowCopy creates a shallow copy of this skEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *pkEncryptor) ShallowCopy() Encryptor {
	return NewEncryptor(enc.params, enc.pk)
}

// ShallowCopy creates a shallow copy of this skEncryptor in which all the read-only data-structures are
// shared with the receiver and the temporary buffers are reallocated. The receiver and the returned
// Encryptors can be used concurrently.
func (enc *skEncryptor) ShallowCopy() Encryptor {
	return NewEncryptor(enc.params, enc.sk)
}

// WithKey returns this encryptor with a new key.
func (enc *skEncryptor) WithKey(key interface{}) Encryptor {
	skPtr, err := enc.checkSk(key)
	if err != nil {
		panic(err)
	}
	return &skEncryptor{enc.encryptorBase, skPtr, enc.uniformSampler}
}

// WithKey returns this encryptor with a new key.
func (enc *pkEncryptor) WithKey(key interface{}) Encryptor {
	pkPtr, err := enc.checkPk(key)
	if err != nil {
		panic(err)
	}
	return &pkEncryptor{enc.encryptorBase, pkPtr}
}

// WithPRNG returns this encryptor with prng as its source of randomness for the uniform
// element c1.
func (enc skEncryptor) WithPRNG(prng utils.PRNG) PRNGEncryptor {
	return &skEncryptor{enc.encryptorBase, enc.sk, enc.uniformSampler.WithPRNG(prng)}
}

// checkPk checks that a given pk is correct for the parameters.
func (enc encryptorBase) checkPk(key interface{}) (pk *PublicKey, err error) {

	switch key := key.(type) {
	case PublicKey:
		pk = &key
	case *PublicKey:
		pk = key
	default:
		return nil, fmt.Errorf("key is not a valid public key type %T", key)
	}

	if pk.Value[0].Q.N() != enc.params.N() || pk.Value[1].Q.N() != enc.params.N() {
		return nil, fmt.Errorf("pk ring degree does not match params ring degree")
	}

	return pk, nil
}

// checkPk checks that a given pk is correct for the parameters.
func (enc encryptorBase) checkSk(key interface{}) (sk *SecretKey, err error) {

	switch key := key.(type) {
	case SecretKey:
		sk = &key
	case *SecretKey:
		sk = key
	default:
		return nil, fmt.Errorf("key is not a valid public key type %T", key)
	}

	if sk.Value.Q.N() != enc.params.N() {
		panic("cannot checkSk: sk ring degree does not match params ring degree")
	}

	return sk, nil
}
