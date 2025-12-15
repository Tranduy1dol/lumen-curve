use std::sync::OnceLock;

use mathlib::field::montgomery::MontgomeryParams;
use mathlib::{BigInt, FieldElement, U1024};

use crate::algebra::fields::Fp;
use crate::models::WeierstrassCurve;

static PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();
static SCALAR_PARAMS: OnceLock<MontgomeryParams> = OnceLock::new();

/// Field parameters for the bls6_6 base field (prime p = 43).
///
/// The returned value is a static reference to the Montgomery parameters used for
/// arithmetic in the curve's base field.
///
/// # Examples
///
/// ```
/// let _params = crate::instances::bls6_6::get_params();
/// ```
pub fn get_params() -> &'static MontgomeryParams {
    PARAMS.get_or_init(|| {
        let p = U1024::from_u64(43);
        MontgomeryParams::new(p, U1024::zero())
    })
}

/// Provides the Montgomery parameters for the scalar field of the bls6_6 curve (modulus = 39).
///
/// # Examples
///
/// ```
/// let _params = get_scalar_params();
/// ```
pub fn get_scalar_params() -> &'static MontgomeryParams {
    SCALAR_PARAMS.get_or_init(|| {
        let order = U1024::from_u64(39);
        MontgomeryParams::new(order, U1024::zero())
    })
}

/// Generator point coordinates for the bls6_6 curve.
///
/// The generator point is returned as an (x, y) pair in the base field, represented as `U1024` values.
///
/// # Examples
///
/// ```
/// let (x, y) = get_generator_coords();
/// assert_eq!(x, U1024::from_u64(13));
/// assert_eq!(y, U1024::from_u64(15));
/// ```
pub fn get_generator_coords() -> (U1024, U1024) {
    (U1024::from_u64(13), U1024::from_u64(15))
}

/// Beta constant for the quadratic extension field.
///
/// Returns the beta value as a `U1024`.
///
/// # Examples
///
/// ```
/// let beta = get_beta();
/// assert_eq!(beta, U1024::from_u64(42));
/// ```
pub fn get_beta() -> U1024 {
    U1024::from_u64(42)
}

/// Constructs the BLS6-6 Weierstrass curve over F_43 with equation y^2 = x^3 + 6.
///
/// The curve is created with curve parameter a = 0 and b = 6, uses the base field
/// modulus p = 43 and the scalar field order = 39, and is returned with its generator point
/// (x, y) = (13, 15).
///
/// # Examples
///
/// ```
/// let _curve = get_curve();
/// // curve is ready for use
/// ```
pub fn get_curve() -> WeierstrassCurve<'static> {
    let params = get_params();
    let scalar_params = get_scalar_params();

    // a = 0, b = 6
    let a = Fp::from(FieldElement::zero(params));
    let b = Fp::new(U1024::from_u64(6), params);

    let (gen_x_val, gen_y_val) = get_generator_coords();
    let gen_x = Fp::new(gen_x_val, params);
    let gen_y = Fp::new(gen_y_val, params);

    WeierstrassCurve::new(a, b, params, scalar_params, gen_x, gen_y)
}