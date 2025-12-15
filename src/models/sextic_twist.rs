use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, FieldElement};

use crate::algebra::fields::{Fp, Fp2};
use crate::def_weierstrass_curve;
use crate::traits::{Curve, Field, ProjectivePoint};

def_weierstrass_curve!(SexticTwist, Fp2<'a>);

impl<'a> Curve<'a> for SexticTwist<'a> {
    type Point = STPoint<'a>;

    /// Returns the curve's point at infinity represented in Jacobian coordinates.
    ///
    /// The returned point uses Fp2 coordinates x = (1, 0), y = (1, 0), and z = (0, 0)
    /// and is associated with this curve instance.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// let curve = SexticTwist::default(); // obtain a SexticTwist instance appropriate for your context
    /// let inf = curve.identity();
    /// assert!(inf.is_identity());
    /// ```
    fn identity(&self) -> Self::Point {
        let zero_fp = FieldElement::zero(self.params);
        let one_fp = FieldElement::one(self.params);

        let zero_fp2 = Fp2::new(Fp::from(zero_fp), Fp::from(zero_fp));
        let one_fp2 = Fp2::new(Fp::from(one_fp), Fp::from(zero_fp)); // 1 + 0u

        STPoint {
            x: one_fp2,
            y: one_fp2,
            z: zero_fp2,
            curve: self.clone(),
        }
    }

    /// Checks whether an affine point (x, y) satisfies the curve equation Y^2 = X^3 + aX + b in Fp2.
    ///
    /// # Returns
    ///
    /// `true` if the provided affine coordinates satisfy the curve equation, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// // Obtain the curve and its generator, convert the generator to affine coordinates,
    /// // and verify it lies on the curve.
    /// let curve = SexticTwist::generator().curve.clone();
    /// let (x, y) = SexticTwist::generator().to_affine();
    /// assert!(curve.is_on_curve(&x, &y));
    /// ```
    fn is_on_curve(
        &self,
        x: &<Self::Point as ProjectivePoint<'a>>::Field,
        y: &<Self::Point as ProjectivePoint<'a>>::Field,
    ) -> bool {
        // Check affine: Y^2 = X^3 + aX + b (calculated in Fp2)
        let y2 = *y * *y;
        let x2 = *x * *x;
        let x3 = x2 * *x;
        let ax = self.a * *x;

        let rhs = x3 + ax + self.b;
        y2 == rhs
    }

    /// Access the curve's Montgomery scalar parameters.
    ///
    /// # Returns
    ///
    /// A reference to the `MontgomeryParams` used for scalar-field arithmetic.
    ///
    /// # Examples
    ///
    /// ```
    /// // `curve` is a `SexticTwist` instance
    /// let params = curve.scalar_params();
    /// let _ = params;
    /// ```
    fn scalar_params(&self) -> &'a MontgomeryParams {
        self.scalar_params
    }

    /// Returns the curve's generator point expressed in Jacobian-like coordinates.
    ///
    /// The returned `STPoint` uses the curve's stored generator X and Y values and a Z coordinate
    /// representing the affine point (Z = 1 in Fp2).
    ///
    /// # Examples
    ///
    /// ```
    /// // Given an existing `SexticTwist` instance `curve`, obtain its generator:
    /// let curve = /* obtain or construct a SexticTwist instance */ unimplemented!();
    /// let g = curve.generator();
    /// // `g` is the affine generator embedded in projective coordinates (Z = 1).
    /// assert!(!g.is_identity());
    /// ```
    fn generator(&self) -> Self::Point {
        let x = self.generator_x;
        let y = self.generator_y;
        let _zero_fp = FieldElement::zero(self.params);
        let _zero_fp2 = Fp2::new(Fp::zero(self.params), Fp::zero(self.params));
        let z = Fp2::new(
            Fp::from(FieldElement::one(self.params)),
            Fp::from(FieldElement::zero(self.params)),
        );

        STPoint {
            x,
            y,
            z,
            curve: self.clone(),
        }
    }
}

// --- 2. G2 POINT DEFINITION (Jacobian Coordinates) ---
#[derive(Clone, Debug)]
pub struct STPoint<'a> {
    pub x: Fp2<'a>,
    pub y: Fp2<'a>,
    pub z: Fp2<'a>,
    pub curve: SexticTwist<'a>,
}

impl<'a> PartialEq for STPoint<'a> {
    /// Determines whether two projective points represent the same point on the curve.
    ///
    /// Points with Z == 0 are treated as the identity; otherwise the points are compared by
    /// converting both to affine coordinates and checking that their X and Y components are equal.
    ///
    /// # Returns
    /// `true` if both points are the identity or their affine X and Y coordinates are equal, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given STPoint::new and a curve instance:
    /// // let p = STPoint::new(x1, y1, z1, curve.clone());
    /// // let q = STPoint::new(x2, y2, z2, curve.clone());
    /// // assert_eq!(p == q, p.eq(&q));
    /// ```
    fn eq(&self, other: &Self) -> bool {
        if self.is_identity() {
            return other.is_identity();
        }
        if other.is_identity() {
            return self.is_identity();
        }

        let (x1, y1) = self.to_affine();
        let (x2, y2) = other.to_affine();
        x1 == x2 && y1 == y2
    }
}
impl<'a> Eq for STPoint<'a> {}

impl<'a> STPoint<'a> {
    /// Constructs an `STPoint` from the given projective (Jacobian) coordinates and associates it with `curve`.
    ///
    /// # Examples
    ///
    /// ```
    /// // Create a point in Jacobian form and attach it to a curve:
    /// let p = STPoint::new(x, y, z, curve);
    /// // `p.x`, `p.y`, `p.z` hold the provided coordinates and `p.curve` is `curve`.
    /// ```
    pub fn new(x: Fp2<'a>, y: Fp2<'a>, z: Fp2<'a>, curve: SexticTwist<'a>) -> Self {
        Self { x, y, z, curve }
    }

    /// Computes the additive inverse of this point on the curve.
    ///
    /// The inverse of a point P = (X:Y:Z) in Jacobian coordinates is -P = (X:-Y:Z).
    /// If `self` is the identity (point at infinity), a clone of `self` is returned.
    ///
    /// # Returns
    ///
    /// A new `STPoint` representing the additive inverse of `self`.
    pub fn neg(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }
        // -P = (X, -Y, Z)
        let neg_y = -self.y; // Fp2 implements Neg
        Self {
            x: self.x,
            y: neg_y,
            z: self.z,
            curve: self.curve.clone(),
        }
    }
}

impl<'a> ProjectivePoint<'a> for STPoint<'a> {
    type Field = Fp2<'a>;

    /// Checks whether this point is the curve identity (point at infinity).
    ///
    /// The point is considered the identity when its projective Z coordinate is zero in both Fp2 components.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given an STPoint `p`, returns true when `p` represents the point at infinity.
    /// assert!(p.is_identity());
    /// ```
    fn is_identity(&self) -> bool {
        // Z == 0 in Fp2
        self.z.c0.value == mathlib::U1024::zero() && self.z.c1.value == mathlib::U1024::zero()
    }

    /// Adds two points on the curve and returns their sum in Jacobian coordinates.
    ///
    /// Performs elliptic-curve addition with correct handling of the identity element; if the
    /// inputs are inverses the result is the curve identity. The returned point is associated
    /// with the same curve as `self`.
    ///
    /// # Returns
    /// The curve point representing `self + rhs` in Jacobian coordinates.
    fn add(&self, rhs: &Self) -> Self {
        if self.is_identity() {
            return rhs.clone();
        }
        if rhs.is_identity() {
            return self.clone();
        }

        let z1z1 = self.z.square();
        let z2z2 = rhs.z.square();

        let u1 = self.x * z2z2;
        let u2 = rhs.x * z1z1;

        let s1 = self.y * (rhs.z * z2z2); // Y1 * Z2^3
        let s2 = rhs.y * (self.z * z1z1); // Y2 * Z1^3

        if u1 == u2 {
            return if s1 == s2 {
                self.double()
            } else {
                self.curve.identity()
            };
        }

        let h = u2 - u1;
        let r = s2 - s1;
        let hh = h.square();
        let hhh = hh * h;
        let v = u1 * hh;

        let two_val = mathlib::U1024::from_u64(2);
        let zero_fp = FieldElement::zero(self.curve.params);
        let two_fp = FieldElement::new(two_val, self.curve.params);
        let two_fp2 = Fp2::new(Fp::from(two_fp), Fp::from(zero_fp));

        let x3 = r.square() - hhh - (two_fp2 * v);
        let y3 = (r * (v - x3)) - (s1 * hhh);
        let z3 = self.z * rhs.z * h;

        Self {
            x: x3,
            y: y3,
            z: z3,
            curve: self.curve.clone(),
        }
    }

    /// Doubles this point using Jacobian-coordinate doubling on the sextic-twist curve.
    ///
    /// Uses the Jacobian doubling formulas over Fp2 to compute 2P and returns the resulting point
    /// in the same projective representation.
    ///
    /// # Examples
    ///
    /// ```rust,no_run
    /// // Given an existing STPoint `p` (for example obtained from a curve's generator),
    /// // compute its double.
    /// // let curve = SexticTwist::some_constructor(...);
    /// // let p = curve.generator();
    /// // let r = p.double();
    /// // assert_eq!(r, p.add(&p));
    /// ```
    fn double(&self) -> Self {
        if self.is_identity() {
            return self.clone();
        }

        // Jacobian doubling formulas (same as SWPoint but on Fp2)
        let xx = self.x.square();
        let yy = self.y.square();
        let yyyy = yy.square();
        let zz = self.z.square();

        // 2
        let zero_fp = FieldElement::zero(self.curve.params);
        let two_val = mathlib::U1024::from_u64(2);
        let two_fp = FieldElement::new(two_val, self.curve.params);
        let two_fp2 = Fp2::new(Fp::from(two_fp), Fp::from(zero_fp));

        // S = 2 * ((X * YY) * 2) = 4XY^2
        let s = two_fp2 * ((self.x * yy) * two_fp2);

        // M = 3*X^2 + a*Z^4
        let three_val = mathlib::U1024::from_u64(3);
        let three_fp = FieldElement::new(three_val, self.curve.params);
        let three_fp2 = Fp2::new(Fp::from(three_fp), Fp::from(zero_fp));

        let zzzz = zz.square();
        let m = (three_fp2 * xx) + (self.curve.a * zzzz);

        // X' = M^2 - 2S
        let x_new = m.square() - (s * two_fp2);

        // Z' = (Y + Z)^2 - YY - ZZ = 2YZ
        // Or Z' = 2 * Y * Z
        let z_new = (self.y * self.z) * two_fp2;

        // Y' = M(S - X') - 8YYYY
        let eight_val = mathlib::U1024::from_u64(8);
        let eight_fp = FieldElement::new(eight_val, self.curve.params);
        let eight_fp2 = Fp2::new(Fp::from(eight_fp), Fp::from(zero_fp));

        let t = eight_fp2 * yyyy;
        let y_new = (m * (s - x_new)) - t;

        Self {
            x: x_new,
            y: y_new,
            z: z_new,
            curve: self.curve.clone(),
        }
    }

    /// Converts this Jacobian-coordinate point into affine coordinates over Fp2.
    ///
    /// If the point is the identity (point-at-infinity), returns (0, 0) in Fp2.
    /// Otherwise returns the pair (x_aff, y_aff) computed by multiplying X by Z^{-2}
    /// and Y by Z^{-3}, where Z^{-1} is the multiplicative inverse of Z in Fp2.
    ///
    /// # Examples
    ///
    /// ```
    /// // given a projective point `pt: STPoint`
    /// let (x_aff, y_aff) = pt.to_affine();
    /// ```
    fn to_affine(&self) -> (Self::Field, Self::Field) {
        if self.is_identity() {
            let zero = FieldElement::zero(self.curve.params);
            let zero_fp2 = Fp2::new(Fp::from(zero), Fp::from(zero));
            return (zero_fp2, zero_fp2);
        }

        let z_inv = self.z.inv().unwrap(); // Fp2 invert
        let z2_inv = z_inv.square();
        let z3_inv = z2_inv * z_inv;

        let x_aff = self.x * z2_inv;
        let y_aff = self.y * z3_inv;
        (x_aff, y_aff)
    }

    /// Computes the scalar multiple of this point.
    ///
    /// The result is equivalent to adding this point to itself `scalar` times.
    ///
    /// # Returns
    ///
    /// The point resulting from multiplying this point by `scalar`.
    ///
    /// # Examples
    ///
    /// ```
    /// // `p` is an STPoint and `k` is a mathlib::U1024 scalar.
    /// // Multiply p by k:
    /// let r = p.mul(&k);
    /// ```
    fn mul(&self, scalar: &mathlib::U1024) -> Self {
        let mut res = self.curve.identity();
        let num_bits = scalar.bits();
        for i in (0..num_bits).rev() {
            res = res.double();
            if scalar.bit(i) {
                res = res.add(self);
            }
        }
        res
    }
}