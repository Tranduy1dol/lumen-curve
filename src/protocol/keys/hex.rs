//! Hex encoding/decoding for keys.

use crate::traits::ShortWeierstrassConfig;

use super::keypair::{PrivateKey, PublicKey};

/// Error type for hex parsing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum HexError {
    /// Invalid hex string length.
    InvalidLength { expected: usize, got: usize },
    /// Invalid hex character.
    InvalidCharacter(char),
    /// Invalid format or prefix.
    InvalidFormat(String),
}

impl std::fmt::Display for HexError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HexError::InvalidLength { expected, got } => {
                write!(f, "invalid length: expected {}, got {}", expected, got)
            }
            HexError::InvalidCharacter(c) => {
                write!(f, "invalid hex character: '{}'", c)
            }
            HexError::InvalidFormat(msg) => {
                write!(f, "invalid format: {}", msg)
            }
        }
    }
}

impl std::error::Error for HexError {}

/// Trait for types that can be encoded to hex.
pub trait ToHex {
    /// Encode to hex string.
    fn to_hex(&self, compressed: bool) -> String;
}

/// Trait for types that can be decoded from hex.
pub trait FromHex: Sized {
    /// Decode from hex string.
    fn from_hex(hex: &str) -> Result<Self, HexError>;
}

impl<P: ShortWeierstrassConfig> ToHex for PrivateKey<P> {
    fn to_hex(&self, _compressed: bool) -> String {
        format!("0x{}", self.to_hex_string())
    }
}

impl<P: ShortWeierstrassConfig> FromHex for PrivateKey<P> {
    fn from_hex(hex_str: &str) -> Result<Self, HexError> {
        // Remove optional 0x prefix
        let hex_str = hex_str.strip_prefix("0x").unwrap_or(hex_str);
        let hex_str = hex_str.strip_prefix("0X").unwrap_or(hex_str);

        // Validate hex characters
        for c in hex_str.chars() {
            if !c.is_ascii_hexdigit() {
                return Err(HexError::InvalidCharacter(c));
            }
        }

        PrivateKey::<P>::from_hex_string(hex_str)
            .ok_or_else(|| HexError::InvalidFormat("failed to parse private key".to_string()))
    }
}

impl<P: ShortWeierstrassConfig> ToHex for PublicKey<P> {
    fn to_hex(&self, compressed: bool) -> String {
        let affine = self.point().to_affine();

        if self.point().is_identity() {
            return "0x00".to_string();
        }

        // Get x coordinate as hex
        let x = affine.x.to_u1024();
        let x_hex = format_u1024_hex(&x);

        if compressed {
            // Compressed: prefix (02/03) + x coordinate
            // 02 if y is even, 03 if y is odd
            let y = affine.y.to_u1024();
            let prefix = if y.0[0] & 1 == 0 { "02" } else { "03" };
            format!("0x{}{}", prefix, x_hex)
        } else {
            // Uncompressed: 04 + x + y
            let y = affine.y.to_u1024();
            let y_hex = format_u1024_hex(&y);
            format!("0x04{}{}", x_hex, y_hex)
        }
    }
}

fn format_u1024_hex(val: &mathlib::U1024) -> String {
    let mut hex = String::new();
    let mut started = false;
    for limb in val.0.iter().rev() {
        if *limb != 0 || started {
            if started {
                hex.push_str(&format!("{:016x}", limb));
            } else {
                hex.push_str(&format!("{:x}", limb));
                started = true;
            }
        }
    }
    if hex.is_empty() {
        hex.push('0');
    }
    hex
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::bls6_6::Bls6_6G1Config;
    use mathlib::U1024;

    #[test]
    fn test_private_key_hex_roundtrip() {
        let scalar = U1024::from_u64(5); // Within order 13
        let key = PrivateKey::<Bls6_6G1Config>::new(scalar);
        let hex_str = key.to_hex(false);
        let recovered = PrivateKey::<Bls6_6G1Config>::from_hex(&hex_str).unwrap();
        assert_eq!(key.to_u1024(), recovered.to_u1024());
    }

    #[test]
    fn test_private_key_from_hex_with_prefix() {
        let key = PrivateKey::<Bls6_6G1Config>::from_hex("0x5").unwrap();
        assert_eq!(key.to_u1024().0[0], 5);
    }

    #[test]
    fn test_private_key_from_hex_without_prefix() {
        let key = PrivateKey::<Bls6_6G1Config>::from_hex("7").unwrap();
        assert_eq!(key.to_u1024().0[0], 7);
    }

    #[test]
    fn test_invalid_hex_character() {
        let result = PrivateKey::<Bls6_6G1Config>::from_hex("0xGHIJ");
        assert!(matches!(result, Err(HexError::InvalidCharacter(_))));
    }
}
