use std::fmt::Debug;

use mathlib::field::element::FieldElement;
use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, U1024};

use crate::models::signature::Signature;

pub trait Curve<'a>: Clone + Debug {
    type Point: ProjectivePoint;

    fn identity(&self) -> Self::Point;
    fn is_on_curve(&self, x: &FieldElement, y: &FieldElement) -> bool;
    fn scalar_params(&self) -> &'a MontgomeryParams;
    fn generator(&self) -> Self::Point;

    /// Generates a keypair (private_key, public_key) from a random scalar.
    ///
    /// The private key is the scalar, and the public key is scalar * generator.
    fn generate_keypair(&self, private_key: &U1024) -> Self::Point {
        self.generator().mul(private_key)
    }

    /// Signs a message hash using ECDSA with the given private key and nonce.
    ///
    /// Returns a Signature containing (r, s) components.
    ///
    /// # Panics
    /// Panics if r or s is zero.
    fn sign(&self, message_hash: &U1024, priv_key: &U1024, k_nonce: &U1024) -> Signature {
        let scalar_params = self.scalar_params();
        let generator = self.generator();

        let r_point = generator.mul(k_nonce);
        let (r_x_elem, _) = r_point.to_affine();
        let r_val = r_x_elem.to_u1024();

        let r_elem = FieldElement::new(r_val, scalar_params);
        let r = r_elem.to_u1024();

        if r == U1024::zero() {
            panic!("r cannot be zero");
        }

        let k_elem = FieldElement::new(*k_nonce, scalar_params);
        let z_elem = FieldElement::new(*message_hash, scalar_params);
        let d_elem = FieldElement::new(*priv_key, scalar_params);

        let k_inv = k_elem.inv();
        let s_elem = k_inv * (z_elem + (r_elem * d_elem));

        let s = s_elem.to_u1024();
        if s == U1024::zero() {
            panic!("s cannot be zero");
        }

        Signature::new(r, s)
    }

    /// Verifies an ECDSA signature against a message hash and public key.
    ///
    /// Returns true if the signature is valid, false otherwise.
    fn verify(&self, signature: &Signature, message_hash: &U1024, pub_key: &Self::Point) -> bool {
        let scalar_params = self.scalar_params();
        let generator = self.generator();
        let n = &scalar_params.modulus;

        if signature.r == U1024::zero() || signature.r >= *n {
            return false;
        }
        if signature.s == U1024::zero() || signature.s >= *n {
            return false;
        }

        let s_elem = FieldElement::new(signature.s, scalar_params);
        let z_elem = FieldElement::new(*message_hash, scalar_params);
        let r_elem = FieldElement::new(signature.r, scalar_params);

        let w = s_elem.inv();

        let u1 = (z_elem * w).to_u1024();
        let u2 = (r_elem * w).to_u1024();

        let p1 = generator.mul(&u1);
        let p2 = pub_key.mul(&u2);
        let p = p1.add(&p2);

        if p.is_identity() {
            return false;
        }

        let (p_x, _) = p.to_affine();
        let p_x_mod_n = FieldElement::new(p_x.to_u1024(), scalar_params);

        p_x_mod_n.to_u1024() == signature.r
    }
}

pub trait ProjectivePoint: Sized + Clone + Debug + PartialEq + Eq {
    fn is_identity(&self) -> bool;

    fn add(&self, rhs: &Self) -> Self;

    fn double(&self) -> Self;

    fn to_affine(&self) -> (FieldElement<'_>, FieldElement<'_>);

    fn mul(&self, scalar: &U1024) -> Self;
}
