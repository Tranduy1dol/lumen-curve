use std::ops::Mul;

use mathlib::field::montgomery::MontgomeryParams;

use crate::algebra::fields::Fp;
use crate::def_fp2;
use crate::traits::Field;

def_fp2!(Fp2, Fp<'a>);

// Macro implements add, Sub, Neg

impl<'a> Mul for Fp2<'a> {
    type Output = Self;
    /// Multiplies two elements of the quadratic extension field Fp2.
    ///
    /// Uses a Karatsuba-like method to compute (a + b·u)·(c + d·u) with three base-field
    /// multiplications: v0 = a·c, v1 = b·d, v2 = (a+b)·(c+d); the result is
    /// (v0 + β·v1) + (v2 - v0 - v1)·u (here β = -1 is encoded by the implementation).
    ///
    /// # Examples
    ///
    /// ```
    /// // Construct two Fp2 values (field construction shown illustratively;
    /// // replace with the crate's actual constructors as needed).
    /// let x = Fp2 { c0: Fp::one(), c1: Fp::zero() }; // represents 1
    /// let y = Fp2 { c0: Fp::zero(), c1: Fp::one() }; // represents u
    /// let z = x * y; // expected to follow quadratic extension multiplication rules
    /// ```
    fn mul(self, rhs: Self) -> Self {
        // Quadratic extension field multiplication: (a + bu)(c + du)
        // Using Karatsuba: need only 3 multiplications instead of 4
        // v0 = a*c, v1 = b*d, v2 = (a+b)(c+d)
        // Result: (v0 + β*v1) + (v2 - v0 - v1)u
        // where β is the quadratic non-residue

        let v0 = self.c0 * rhs.c0; // a * c
        let v1 = self.c1 * rhs.c1; // b * d  
        let v2 = (self.c0 + self.c1) * (rhs.c0 + rhs.c1); // (a+b)(c+d)

        Self {
            c0: v0 - v1,      // Real part: ac + β·bd where β = -1
            c1: v2 - v0 - v1, // Imaginary part
        }
    }
}

// Neg implemented by macro

// Trait Field
impl<'a> Field<'a> for Fp2<'a> {
    /// Constructs the additive identity (zero) of `Fp2` using the provided Montgomery parameters.
    ///
    /// The returned `Fp2` has both `c0` and `c1` set to the zero value of the base `Fp` field
    /// instantiated with `params`.
    ///
    /// # Examples
    ///
    /// ```
    /// // let params: &MontgomeryParams = /* obtain MontgomeryParams for the base field */ ;
    /// // let z = Fp2::zero(params);
    /// // assert!(z.is_zero());
    /// ```
    fn zero(params: &'a MontgomeryParams) -> Self {
        let z = Fp::zero(params);
        Self { c0: z, c1: z }
    }
    /// Returns whether both components of the quadratic extension element are zero.
    ///
    /// # Returns
    /// `true` if both `c0` and `c1` are zero, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// // assuming `params` is a `&MontgomeryParams` and Fp2::zero is available
    /// let z = Fp2::zero(params);
    /// assert!(z.is_zero());
    /// ```
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Constructs the multiplicative identity element of Fp2 for the provided Montgomery parameters.
    ///
    /// The returned element has `c0 = 1` (the base-field one) and `c1 = 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// // `params` is a reference to MontgomeryParams appropriate for the base field.
    /// let one = Fp2::one(params);
    /// assert_eq!(one.c0, Fp::one(params));
    /// assert!(one.c1.is_zero());
    /// ```
    fn one(params: &'a MontgomeryParams) -> Self {
        let z = Fp::zero(params);
        let o = Fp::one(params);
        Self { c0: o, c1: z }
    }

    /// Computes the multiplicative inverse of this quadratic-extension element when it exists.
    ///
    /// The norm is computed as `a² + b²` for an element `a + b·u`. If the norm is nonzero, returns
    /// `Some(inv)` where `inv = (a, -b) * inv_norm` and `inv_norm` is the inverse of the norm.
    /// Returns `None` when the norm is zero (no multiplicative inverse).
    ///
    /// # Examples
    ///
    /// ```
    /// // `params` must be a valid `&MontgomeryParams` for the base field.
    /// let one = Fp2::one(params);
    /// assert_eq!(one.inv().unwrap(), one);
    ///
    /// let zero = Fp2::zero(params);
    /// assert!(zero.inv().is_none());
    /// ```
    fn inv(&self) -> Option<Self> {
        // For quadratic extension (a + bu), inverse is computed as:
        // (a + bu)^-1 = (a' - bu') / norm
        // where norm must match our multiplication law
        // For β=1: norm = a^2 - b^2

        let a = self.c0;
        let b = self.c1;

        let a_sq = a.square();
        let b_sq = b.square();
        let norm = a_sq + b_sq; // For β = -1: norm = a² - β·b² = a² + b²

        let inv_norm = Field::inv(&norm)?; // Return None if norm = 0

        let z = Fp::zero(b.params);
        let neg_b = z - b;

        Some(Self {
            c0: a * inv_norm,
            c1: neg_b * inv_norm,
        })
    }

    /// Doubles this field element by adding it to itself.
    ///
    /// # Examples
    ///
    /// ```
    /// // `a` is an `Fp2` element previously constructed.
    /// let doubled = a.double();
    /// assert_eq!(doubled, a + a);
    /// ```
    fn double(&self) -> Self {
        *self + *self
    }

    /// Multiplies two Fp2 field elements.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given two Fp2 values `a` and `b` (constructed with appropriate MontgomeryParams):
    /// // let params = /* MontgomeryParams instance */;
    /// // let a = Fp2::one(&params);
    /// // let b = Fp2::one(&params);
    /// // The product can be computed either with `*` or with `mul`.
    /// let c = a.mul(&b);
    /// assert_eq!(c, a * b);
    /// ```
    ///
    /// @returns The field element equal to `self` multiplied by `rhs`.
    fn mul(&self, rhs: &Self) -> Self {
        *self * *rhs
    }

    /// Returns the sum of `self` and `rhs`.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given two `Fp2` values `a` and `b`:
    /// let sum = a.add(&b);
    /// assert_eq!(sum, a + b);
    /// ```
    fn add(&self, rhs: &Self) -> Self {
        *self + *rhs
    }

    /// Squares the field element.
    ///
    /// The square of the element.
    ///
    /// # Examples
    ///
    /// ```
    /// let elem = /* an `Fp2` value */ unimplemented!();
    /// let sq = elem.square();
    /// assert_eq!(sq, elem * elem);
    /// ```
    fn square(&self) -> Self {
        *self * *self
    }
}