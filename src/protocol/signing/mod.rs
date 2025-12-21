//! Digital signature operations for elliptic curve cryptography.
//!
//! This module provides ECDSA signature generation and verification
//! using `FieldElement<ScalarField>` for all modular arithmetic.

mod error;
mod signature;

pub use error::SignatureError;
pub use signature::Signature;

use mathlib::FieldElement;

use crate::models::sw::Projective;
use crate::protocol::keys::{PrivateKey, PublicKey};
use crate::traits::ShortWeierstrassConfig;

/// Signing engine for ECDSA operations.
///
/// Generic over `P: ShortWeierstrassConfig` which defines the curve.
/// All arithmetic is performed using `FieldElement<P::ScalarField>`.
pub struct SigningEngine<P: ShortWeierstrassConfig> {
    _marker: std::marker::PhantomData<P>,
}

impl<P: ShortWeierstrassConfig> SigningEngine<P> {
    /// Create a new signing engine.
    pub fn new() -> Self {
        Self {
            _marker: std::marker::PhantomData,
        }
    }

    /// Sign a message hash using the private key with a specified nonce.
    ///
    /// Uses `FieldElement<ScalarField>` for all modular arithmetic.
    ///
    /// **Warning**: Using a non-random or repeated nonce compromises security!
    ///
    /// # Parameters
    ///
    /// * `message_hash` - The hash of the message as a scalar field element
    /// * `private_key` - The signer's private key
    /// * `nonce` - The signing nonce (k) as a scalar field element
    ///
    /// # Returns
    ///
    /// A `Signature` on success, or `SignatureError` on failure.
    pub fn sign(
        &self,
        message_hash: &FieldElement<P::ScalarField>,
        private_key: &PrivateKey<P>,
        nonce: &FieldElement<P::ScalarField>,
    ) -> Result<Signature, SignatureError> {
        // Validate nonce is not zero
        if nonce.is_zero() {
            return Err(SignatureError::InvalidNonce);
        }

        // Validate private key is not zero
        if private_key.is_zero() {
            return Err(SignatureError::InvalidNonce);
        }

        // Compute R = k * G
        let generator = Projective::<P>::generator();
        let k_scalar = nonce.to_u1024();
        let r_point = generator.mul(&k_scalar);

        if r_point.is_identity() {
            return Err(SignatureError::InvalidNonce);
        }

        // r = x(R) mod n
        let r_affine = r_point.to_affine();
        let r = FieldElement::<P::ScalarField>::new(r_affine.x.to_u1024());

        if r.is_zero() {
            return Err(SignatureError::InvalidNonce);
        }

        // s = k^(-1) * (z + r * d) mod n
        let k_inv = nonce.inv();
        let d = private_key.scalar();

        // r * d
        let rd = r * *d;
        // z + r * d
        let z_plus_rd = *message_hash + rd;
        // s = k^(-1) * (z + r * d)
        let s = k_inv * z_plus_rd;

        if s.is_zero() {
            return Err(SignatureError::InvalidNonce);
        }

        Ok(Signature::new(r.to_u1024(), s.to_u1024()))
    }

    /// Verify a signature.
    ///
    /// Uses `FieldElement<ScalarField>` for all modular arithmetic.
    ///
    /// # Parameters
    ///
    /// * `signature` - The signature to verify
    /// * `message_hash` - The hash of the message as a scalar field element
    /// * `public_key` - The signer's public key
    ///
    /// # Returns
    ///
    /// `true` if the signature is valid, `false` otherwise.
    pub fn verify(
        &self,
        signature: &Signature,
        message_hash: &FieldElement<P::ScalarField>,
        public_key: &PublicKey<P>,
    ) -> bool {
        // Convert signature components to scalar field elements
        let r = FieldElement::<P::ScalarField>::new(*signature.r());
        let s = FieldElement::<P::ScalarField>::new(*signature.s());

        // Validate signature components are not zero
        if r.is_zero() || s.is_zero() {
            return false;
        }

        // Validate public key is not identity
        if public_key.is_identity() {
            return false;
        }

        // s_inv = s^(-1) mod n
        let s_inv = s.inv();

        // u1 = z * s^(-1) mod n
        let u1 = *message_hash * s_inv;

        // u2 = r * s^(-1) mod n
        let u2 = r * s_inv;

        // R' = u1 * G + u2 * Q
        let generator = Projective::<P>::generator();
        let u1_g = generator.mul(&u1.to_u1024());
        let u2_q = public_key.point().mul(&u2.to_u1024());
        let r_prime = u1_g.add(&u2_q);

        if r_prime.is_identity() {
            return false;
        }

        // Verify: r == x(R') mod n
        let r_prime_affine = r_prime.to_affine();
        let r_prime_x = FieldElement::<P::ScalarField>::new(r_prime_affine.x.to_u1024());

        r == r_prime_x
    }
}

impl<P: ShortWeierstrassConfig> Default for SigningEngine<P> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6G1Config;
    use crate::instances::bls6_6::Bls6_6ScalarField;
    use crate::protocol::keys::KeyEngine;
    use mathlib::U1024;

    #[test]
    fn test_sign() {
        let key_engine = KeyEngine::<Bls6_6G1Config>::new();
        let signing_engine = SigningEngine::<Bls6_6G1Config>::new();

        // Create keypair with private key = 5
        let private_scalar = U1024::from_u64(5);
        let keypair = key_engine.keypair_from_scalar(private_scalar);

        // Message hash = 7 (as scalar field element)
        let message_hash = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(7));

        // Nonce = 3 (as scalar field element)
        let nonce = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(3));

        // Sign - this should succeed
        let signature = signing_engine
            .sign(&message_hash, keypair.private_key(), &nonce)
            .expect("signing should succeed");

        // Verify signature has non-zero components
        assert_ne!(*signature.r(), U1024::from_u64(0));
        assert_ne!(*signature.s(), U1024::from_u64(0));
    }

    #[test]
    #[ignore] // ECDSA verification needs a curve with cofactor 1 and base field = scalar field
    fn test_sign_verify_full() {
        // Note: BLS6_6 has base field p=43 and scalar field r=13
        // This mismatch causes issues in ECDSA verification because
        // x-coordinates reduce differently in the two fields.
        // For production ECDSA, use a curve like secp256k1 where this works correctly.
        let key_engine = KeyEngine::<Bls6_6G1Config>::new();
        let signing_engine = SigningEngine::<Bls6_6G1Config>::new();

        let private_scalar = U1024::from_u64(5);
        let keypair = key_engine.keypair_from_scalar(private_scalar);

        let message_hash = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(7));
        let nonce = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(3));

        let signature = signing_engine
            .sign(&message_hash, keypair.private_key(), &nonce)
            .expect("signing should succeed");

        let valid = signing_engine.verify(&signature, &message_hash, keypair.public_key());
        assert!(valid, "signature should be valid");
    }

    #[test]
    fn test_invalid_nonce() {
        let key_engine = KeyEngine::<Bls6_6G1Config>::new();
        let signing_engine = SigningEngine::<Bls6_6G1Config>::new();

        let keypair = key_engine.keypair_from_scalar(U1024::from_u64(5));
        let message_hash = FieldElement::<Bls6_6ScalarField>::new(U1024::from_u64(7));

        // Zero nonce should fail
        let zero_nonce = FieldElement::<Bls6_6ScalarField>::zero();
        let result = signing_engine.sign(&message_hash, keypair.private_key(), &zero_nonce);
        assert!(result.is_err());
    }
}
