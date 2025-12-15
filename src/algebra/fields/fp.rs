use std::ops::{Add, Deref, Mul, Neg, Sub};

use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{FieldElement, U1024};

use crate::traits::{Field, ToU1024};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp<'a>(pub FieldElement<'a>);

impl<'a> Fp<'a> {
    /// Constructs an `Fp` by creating an inner `FieldElement` from the given value and Montgomery parameters.
    ///
    /// # Examples
    ///
    /// ```
    /// // Obtain or construct `MontgomeryParams` appropriate for your modulus.
    /// let params = /* MontgomeryParams instance */ unimplemented!();
    /// let value: U1024 = U1024::from(42u64);
    /// let fp = Fp::new(value, &params);
    /// ```
    pub fn new(value: U1024, params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::new(value, params))
    }
}

impl<'a> From<FieldElement<'a>> for Fp<'a> {
    /// Wraps a `FieldElement` value in the `Fp` newtype.
    ///
    /// # Examples
    ///
    /// ```
    /// // `fe` is a FieldElement obtained from library constructors.
    /// let fe: FieldElement<'_> = /* ... */;
    /// let fp: Fp = fe.into();
    /// ```
    fn from(f: FieldElement<'a>) -> Self {
        Self(f)
    }
}

// Deref to a field element for convenience (accessing .params, .value)
impl<'a> Deref for Fp<'a> {
    type Target = FieldElement<'a>;
    /// Accesses the underlying `FieldElement` by reference.
    ///
    /// Returns a shared reference to the wrapped `FieldElement`.
    ///
    /// # Examples
    ///
    /// ```
    /// // `params` and `U1024` construction omitted for brevity.
    /// let params = /* MontgomeryParams instance */;
    /// let fp = Fp::new(U1024::from(1u64), &params);
    /// let inner = fp.deref();
    /// // `inner` can be used to read fields or call methods on `FieldElement`.
    /// ```
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> Add for Fp<'a> {
    type Output = Self;
    /// Adds two `Fp` elements.
    ///
    /// Returns an `Fp` containing the sum of `self` and `rhs`.
    ///
    /// # Examples
    ///
    /// ```
    /// // given two `Fp` values `a` and `b`
    /// let c = a + b;
    /// ```
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0)
    }
}

impl<'a> Sub for Fp<'a> {
    type Output = Self;
    /// Subtracts one `Fp` value from another and returns the difference.
    ///
    /// # Examples
    ///
    /// ```
    /// // Given two `Fp` values `a` and `b`, compute their difference:
    /// let c = a - b;
    /// ```
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0)
    }
}

impl<'a> Mul for Fp<'a> {
    type Output = Self;
    /// Multiplies two field elements and returns their product wrapped in `Fp`.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// let a = Fp::new(2u32.into(), params);
    /// let b = Fp::new(3u32.into(), params);
    /// let c = a * b;
    /// assert_eq!(c.to_u1024_val(), (2u32 * 3u32).into());
    /// ```
    fn mul(self, rhs: Self) -> Self {
        Self(self.0 * rhs.0)
    }
}

impl<'a> Neg for Fp<'a> {
    type Output = Self;
    /// Computes the additive inverse of the field element.
    ///
    /// # Examples
    ///
    /// ```
    /// // let params = ...; // obtain MontgomeryParams for the field
    /// // let a = Fp::new(3u32.into(), params);
    /// // let neg_a = -a;
    /// // assert_eq!(a + neg_a, Fp::zero(params));
    /// ```
    fn neg(self) -> Self {
        let zero = FieldElement::zero(self.0.params);
        Self(zero - self.0)
    }
}

impl<'a> Field<'a> for Fp<'a> {
    /// Create the additive identity (zero) in the field defined by `params`.
    ///
    /// # Examples
    ///
    /// ```
    /// // let params = MontgomeryParams::new(...);
    /// let z = Fp::zero(&params);
    /// assert!(z.is_zero());
    /// ```
    fn zero(params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::zero(params))
    }

    /// Checks whether the element is the additive zero for its Montgomery parameters.
    ///
    /// # Returns
    ///
    /// `true` if the element equals the field zero for its associated `MontgomeryParams`, `false` otherwise.
    ///
    /// # Examples
    ///
    /// ```
    /// let z = Fp::zero(params);
    /// assert!(z.is_zero());
    ///
    /// let one = Fp::one(params);
    /// assert!(!one.is_zero());
    /// ```
    fn is_zero(&self) -> bool {
        self.0 == FieldElement::zero(self.0.params)
    }

    /// Constructs the multiplicative identity (one) for the given Montgomery parameters.
    ///
    /// # Returns
    ///
    /// The field element representing 1 for `params`.
    ///
    /// # Examples
    ///
    /// ```
    /// let one = Fp::one(params);
    /// // multiplicative identity: one * one == one
    /// assert_eq!(one.mul(&one).to_u1024_val(), one.to_u1024_val());
    /// ```
    fn one(params: &'a MontgomeryParams) -> Self {
        Self(FieldElement::one(params))
    }

    /// Compute the multiplicative inverse of this field element.
    ///
    /// Returns `Some(Self)` containing the inverse when the element is not zero, or `None` when the element is zero.
    ///
    /// # Examples
    ///
    /// ```
    /// // assuming `params` is a valid `&MontgomeryParams` in scope
    /// let a = Fp::new(5u128.into(), params);
    /// let inv = a.inv().expect("nonzero element should have an inverse");
    /// assert_eq!(a.mul(&inv).to_u1024_val(), Fp::one(params).to_u1024_val());
    ///
    /// let zero = Fp::zero(params);
    /// assert!(zero.inv().is_none());
    /// ```
    fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            Some(Self(self.0.inv()))
        }
    }

    /// Computes the element added to itself.
    ///
    /// # Returns
    ///
    /// `Self` equal to this element plus itself.
    ///
    /// # Examples
    ///
    /// ```
    /// // given `fp: Fp`
    /// let doubled = fp.double();
    /// assert_eq!(doubled, fp + fp);
    /// ```
    fn double(&self) -> Self {
        Self(self.0 + self.0)
    }

    /// Multiplies two field elements and returns their product.
    ///
    /// # Returns
    ///
    /// The product of `self` and `rhs`.
    ///
    /// # Examples
    ///
    /// ```
    /// // Construct `params` appropriately for your field implementation.
    /// let params = /* MontgomeryParams */ unimplemented!();
    /// let a = Fp::new(3u128.into(), &params);
    /// let b = Fp::new(4u128.into(), &params);
    /// let c = a.mul(&b);
    /// assert_eq!(c.to_u1024_val(), (3u128 * 4u128).into());
    /// ```
    fn mul(&self, rhs: &Self) -> Self {
        Self(self.0 * rhs.0)
    }

    /// Adds two field elements and returns their sum.
    ///
    /// # Examples
    ///
    /// ```
    /// // `a` and `b` are `Fp` field elements with the same Montgomery parameters.
    /// let sum = a.add(&b);
    /// assert_eq!(sum, a + b);
    /// ```
    fn add(&self, rhs: &Self) -> Self {
        Self(self.0 + rhs.0)
    }

    /// Computes the square of the field element.
    ///
    /// # Returns
    ///
    /// `Self` equal to this element multiplied by itself.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// // `params` should be a valid `&MontgomeryParams` for the field.
    /// let params = /* obtain MontgomeryParams */ todo!();
    /// let a = Fp::one(&params);
    /// let s = a.square();
    /// assert_eq!(s, a.mul(&a));
    /// ```
    fn square(&self) -> Self {
        Self(self.0 * self.0)
    }
}

impl<'a> ToU1024 for Fp<'a> {
    /// Convert this field element into its canonical `U1024` integer representation.
    ///
    /// # Examples
    ///
    /// ```
    /// // Construct an `Fp` value (creation omitted) and obtain its integer value:
    /// // let params = /* MontgomeryParams */ ;
    /// // let fp = Fp::new(some_u1024_value, &params);
    /// let int_val = fp.to_u1024_val();
    /// // `int_val` is a `U1024` representing the field element.
    /// ```
    fn to_u1024_val(&self) -> U1024 {
        self.0.to_u1024()
    }
}