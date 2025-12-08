use curvelib::{curves::tiny_jubjub, traits::Curve};
use mathlib::{BigInt, U1024};

#[test]
fn test_ecdsa_signature_verification() {
    let curve = tiny_jubjub::get_curve();

    // 1. Generate Keypair
    // Private key d = 5
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);

    // 2. Prepare Message and Nonce
    // Message hash z = 11
    let message_hash = U1024::from_u64(11);
    // Random nonce k = 7 (must be coprime to order, 7 and 20 are coprime)
    let k_nonce = U1024::from_u64(2);

    // 3. Sign
    let signature = curve.sign(&message_hash, &priv_key, &k_nonce);

    println!("Signature: r={:?}, s={:?}", signature.r, signature.s);

    // 4. Verify
    let valid = curve.verify(&signature, &message_hash, &pub_key);

    assert!(valid, "Signature verification failed");
}

#[test]
fn test_ecdsa_invalid_signature() {
    let curve = tiny_jubjub::get_curve();
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);
    let message_hash = U1024::from_u64(11);
    let k_nonce = U1024::from_u64(2);

    let mut signature = curve.sign(&message_hash, &priv_key, &k_nonce);

    // Tamper with signature
    signature.s = signature.s + U1024::one();

    let valid = curve.verify(&signature, &message_hash, &pub_key);

    assert!(!valid, "Invalid signature should not verify");
}

#[test]
fn test_ecdsa_different_message() {
    let curve = tiny_jubjub::get_curve();
    let priv_key = U1024::from_u64(2);
    let pub_key = curve.generate_keypair(&priv_key);
    let message_hash = U1024::from_u64(11);
    let k_nonce = U1024::from_u64(2);

    let signature = curve.sign(&message_hash, &priv_key, &k_nonce);

    // Verify against a different message
    let wrong_message = U1024::from_u64(12);
    let valid = curve.verify(&signature, &wrong_message, &pub_key);

    assert!(!valid, "Signature should not verify for wrong message");
}
