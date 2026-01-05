# Mathematical Foundations of lumen-curve

This document explains the mathematical theory and cryptographic concepts underlying the `lumen-curve` library.

---

## Table of Contents

1. [Finite Fields](#1-finite-fields)
2. [Field Extensions](#2-field-extensions)
3. [Elliptic Curves](#3-elliptic-curves)
4. [Point Representations](#4-point-representations)
5. [Scalar Multiplication](#5-scalar-multiplication)
6. [Pairing-Based Cryptography](#6-pairing-based-cryptography)
7. [Polynomial Commitments (KZG)](#7-polynomial-commitments-kzg)
8. [Digital Signatures](#8-digital-signatures)
9. [Algorithms](#9-algorithms)

---

## 1. Finite Fields

### Prime Fields ($\mathbb{F}\_p$)

A **prime field** $\mathbb{F}\_p$ is the set of integers $\{0, 1, 2, \ldots, p-1\}$ with arithmetic performed modulo a prime $p$.

**Properties:**
- Addition: $(a + b) \mod p$
- Multiplication: $(a \times b) \mod p$
- Every non-zero element has a multiplicative inverse
- The field has exactly $p$ elements

**In lumen-curve:**
- `FieldElement<C>` represents elements of $\mathbb{F}\_p$ where `C: FieldConfig` defines the modulus
- `Bls6_6BaseField` defines $\mathbb{F}\_{43}$ (the field with 43 elements)
- `Bls6_6ScalarField` defines $\mathbb{F}\_{13}$ (the field with 13 elements)

### Montgomery Form

Field elements are stored in **Montgomery form** for efficient modular multiplication:

$$x\_{mont} = x \cdot R \mod p$$

Where $R = 2^k$ for some $k \geq \log\_2(p)$. Montgomery multiplication avoids expensive division operations.

---

## 2. Field Extensions

### Quadratic Extension ($\mathbb{F}\_{p^2}$)

The field $\mathbb{F}\_{p^2}$ is a degree-2 extension of $\mathbb{F}\_p$, analogous to complex numbers:

$$\mathbb{F}\_{p^2} = \mathbb{F}\_p[u] / (u^2 + 1)$$

Elements have the form $a + b \cdot u$ where $a, b \in \mathbb{F}\_p$ and $u^2 = -1$.

**Arithmetic:**
- Addition: $(a + bu) + (c + du) = (a+c) + (b+d)u$
- Multiplication: $(a + bu)(c + du) = (ac - bd) + (ad + bc)u$

**In lumen-curve:** `Fp2<C>` in `src/algebra/fields/fp2.rs`

### Sextic Extension ($\mathbb{F}\_{p^6}$)

The field $\mathbb{F}\_{p^6}$ is built as a cubic extension over $\mathbb{F}\_{p^2}$:

$$\mathbb{F}\_{p^6} = \mathbb{F}\_{p^2}[v] / (v^3 - \xi)$$

Where $\xi = u$ (the imaginary unit from $\mathbb{F}\_{p^2}$), so $v^3 = u$.

Elements have the form $c\_0 + c\_1 \cdot v + c\_2 \cdot v^2$ where $c\_0, c\_1, c\_2 \in \mathbb{F}\_{p^2}$.

**In lumen-curve:** `Fp6<C>` in `src/algebra/fields/fp6.rs`

### Tower of Extensions

$$\mathbb{F}\_p \xrightarrow{\times 2} \mathbb{F}\_{p^2} \xrightarrow{\times 3} \mathbb{F}\_{p^6}$$

This gives $\mathbb{F}\_{p^6}$ as a degree-6 extension of $\mathbb{F}\_p$, which is the embedding degree needed for BLS6 pairings.

---

## 3. Elliptic Curves

### Short Weierstrass Form

The most common curve form:

$$y^2 = x^3 + ax + b$$

**Curve requirements:**
- The discriminant $\Delta = -16(4a^3 + 27b^2) \neq 0$ (ensures no singular points)

**In lumen-curve:**
- `WeierstrassCurve<C>` and `WeierstrassPoint<C>` (legacy API)
- `Projective<P>` where `P: ShortWeierstrassConfig` (new API)

**BLS6\_6 G1 curve:** $y^2 = x^3 + 6$ over $\mathbb{F}\_{43}$ ($a = 0$, $b = 6$)

### Twisted Edwards Form

An alternative representation with unified addition formulas:

$$ax^2 + y^2 = 1 + dx^2y^2$$

**Advantages:**
- Unified formulas work for all point additions (including doublings)
- Complete formulas exist that work even for edge cases
- Generally faster than Weierstrass for some operations

**In lumen-curve:**
- `EdwardsCurve<C>` and `EdwardsPoint<C>` (legacy API)
- `Projective<P>` where `P: TwistedEdwardsConfig` (new API)

**Tiny Jubjub curve:** $3x^2 + y^2 = 1 + 8x^2y^2$ over $\mathbb{F}\_{13}$ ($a = 3$, $d = 8$)

### Curve Groups

An elliptic curve $E$ over $\mathbb{F}\_p$ forms an **abelian group** where:
- The **identity** is the "point at infinity" (denoted $\mathcal{O}$ or $\infty$)
- **Addition** follows the chord-and-tangent rule
- The **inverse** of $(x, y)$ is $(x, -y)$ for Weierstrass or $(-x, y)$ for Edwards

The number of points on the curve is $|E(\mathbb{F}\_p)|$, which by Hasse's theorem satisfies:

$$\\big| |E(\\mathbb{F}\\_p)| - (p + 1) \\big| \\leq 2\\sqrt{p}$$

---

## 4. Point Representations

### Affine Coordinates

Standard $(x, y)$ representation where $y^2 = x^3 + ax + b$.

**Pros:** Compact storage (2 field elements)  
**Cons:** Point addition requires field inversion (expensive)

### Projective Coordinates

Represent $(x, y)$ as $(X : Y : Z)$ where $x = X/Z$ and $y = Y/Z$.

**The identity point** is represented as $(1 : 1 : 0)$.

**Pros:** Addition without inversions (only multiplications)  
**Cons:** 3 field elements per point

**In lumen-curve:** All arithmetic uses projective coordinates internally.

### Extended Coordinates (Edwards)

For Edwards curves, use $(X : Y : Z : T)$ where:
- $x = X/Z$
- $y = Y/Z$
- $T = XY/Z$

The extra coordinate $T$ speeds up addition by precomputing $xy$.

---

## 5. Scalar Multiplication

Given a point $P$ and scalar $k$, compute $k \cdot P = P + P + \cdots + P$ ($k$ times).

### Double-and-Add Algorithm

The standard binary method (currently implemented):

**Input:** $P$, $k = (k\_{n-1}, \ldots, k\_1, k\_0)\_2$  
**Output:** $k \cdot P$

$$
\begin{aligned}
&R \leftarrow \mathcal{O} \text{ (identity)} \\\\
&\textbf{for } i \text{ from } 0 \text{ to } n-1: \\\\
&\quad \textbf{if } k\_i = 1: R \leftarrow R + P \\\\
&\quad P \leftarrow 2P \\\\
&\textbf{return } R
\end{aligned}
$$

**Complexity:** $O(\log k)$ point operations

**Security Warning:** This algorithm is NOT constant-time and leaks the scalar through timing.

### Montgomery Ladder (Recommended)

A constant-time algorithm:

**Input:** $P$, $k = (k\_{n-1}, \ldots, k\_1, k\_0)\_2$  
**Output:** $k \cdot P$

$$
\begin{aligned}
&R\_0 \leftarrow \mathcal{O}, \quad R\_1 \leftarrow P \\\\
&\textbf{for } i \text{ from } n-1 \text{ down to } 0: \\\\
&\quad \textbf{if } k\_i = 0: R\_1 \leftarrow R\_0 + R\_1, \quad R\_0 \leftarrow 2 \cdot R\_0 \\\\
&\quad \textbf{else: } R\_0 \leftarrow R\_0 + R\_1, \quad R\_1 \leftarrow 2 \cdot R\_1 \\\\
&\textbf{return } R\_0
\end{aligned}
$$

Always performs the same operations regardless of bit values.

---

## 6. Pairing-Based Cryptography

### What is a Pairing?

A **bilinear pairing** is a map:

$$e: \mathbb{G}\_1 \times \mathbb{G}\_2 \rightarrow \mathbb{G}\_T$$

Where $\mathbb{G}\_1$, $\mathbb{G}\_2$ are elliptic curve groups and $\mathbb{G}\_T$ is a multiplicative group (usually in a field extension).

**Properties:**
1. **Bilinearity:** $e(aP, bQ) = e(P, Q)^{ab}$
2. **Non-degeneracy:** $e(P, Q) \neq 1$ for generators $P$, $Q$
3. **Computability:** Can be computed efficiently

### BLS Curves

**BLS (Boneh-Lynn-Shacham)** curves are pairing-friendly curves where:
- $\mathbb{G}\_1$ is on $E(\mathbb{F}\_p)$
- $\mathbb{G}\_2$ is on a twisted curve $E'(\mathbb{F}\_{p^{k/d}})$
- $\mathbb{G}\_T$ is a subgroup of $\mathbb{F}\_{p^k}^*$

The **embedding degree** $k$ determines the security and efficiency tradeoff.

**BLS6\_6 in lumen-curve:** $k = 6$, using a sextic twist

### Miller's Algorithm

Computes the Tate pairing through the **Miller loop**:

**Input:** $P \in \mathbb{G}\_1$, $Q \in \mathbb{G}\_2$, $r$ = group order  
**Output:** $f\_{r,P}(Q)$ before final exponentiation

$$
\begin{aligned}
&f \leftarrow 1, \quad T \leftarrow P \\\\
&\textbf{for } i \text{ from } n-2 \text{ down to } 0: \\\\
&\quad f \leftarrow f^2 \cdot \ell\_{T,T}(Q) \quad \text{// line through } T,T \text{ evaluated at } Q \\\\
&\quad T \leftarrow 2T \\\\
&\quad \textbf{if } r\_i = 1: f \leftarrow f \cdot \ell\_{T,P}(Q), \quad T \leftarrow T + P \\\\
&\textbf{return } f
\end{aligned}
$$

**In lumen-curve:** `src/protocol/pairing/miller.rs`

### Final Exponentiation

The Miller loop output must be raised to power $(p^k - 1)/r$ to get a unique pairing result:

$$e(P, Q) = f\_{r,P}(Q)^{(p^k - 1)/r}$$

**In lumen-curve:** `src/protocol/pairing/final_exp.rs`

### Sextic Twist

A **twist** maps points from $E'(\mathbb{F}\_{p^2})$ to $E(\mathbb{F}\_{p^6})$ for efficient pairing computation:

$$\psi: E'(\mathbb{F}\_{p^2}) \rightarrow E(\mathbb{F}\_{p^6})$$
$$\psi(x, y) = (x \cdot \omega^2, y \cdot \omega^3)$$

Where $\omega$ is a primitive 6th root of unity.

**In lumen-curve:** `SexticTwist<C>` and `TwistPoint<C>` in `src/models/sextic_twist.rs`

---

## 7. Polynomial Commitments (KZG)

The **Kate-Zaverucha-Goldberg (KZG)** scheme allows committing to polynomials with constant-size commitments and proofs.

### Trusted Setup

Generate **Structured Reference String (SRS)**:

Given secret $\tau$ (must be destroyed after setup):
- G1 powers: $[G, \tau G, \tau^2 G, \ldots, \tau^n G]$
- G2 element: $\tau H$ (where $H$ is G2 generator)

### Commitment

For polynomial $P(x) = \sum\_i c\_i x^i$:

$$C = \sum\_i c\_i \cdot [\tau^i]G = [P(\tau)]G$$

This is a single G1 point regardless of polynomial degree.

### Opening

To prove $P(z) = y$:

1. Compute quotient $Q(x) = \frac{P(x) - y}{x - z}$
2. Proof: $\pi = [Q(\tau)]G$

### Verification

Using pairings, check:

$$e(\pi, [\tau]H - [z]H) = e(C - [y]G, H)$$

Expands to: $e(\pi, \tau H - zH) = e(C - yG, H)$

**In lumen-curve:** `src/protocol/commitment/kzg/mod.rs`

---

## 8. Digital Signatures

### ECDSA (Conceptual)

The library provides key pair structures for elliptic curve signatures:

**Key Generation:**
1. Choose random scalar $d \in [1, n-1]$
2. Compute $Q = d \cdot G$ (public key)

**Signing message $m$:**
1. Compute $e = \text{Hash}(m)$
2. Choose random $k \in [1, n-1]$
3. Compute $(x\_1, y\_1) = k \cdot G$
4. $r = x\_1 \mod n$ (if $r = 0$, go to step 2)
5. $s = k^{-1}(e + rd) \mod n$ (if $s = 0$, go to step 2)
6. Signature is $(r, s)$

**Verification:**
1. Compute $e = \text{Hash}(m)$
2. $w = s^{-1} \mod n$
3. $u\_1 = ew \mod n$, $u\_2 = rw \mod n$
4. $(x\_1, y\_1) = u\_1 \cdot G + u\_2 \cdot Q$
5. Valid if $r \equiv x\_1 \pmod{n}$

**In lumen-curve:**
- `PrivateKey<C>`: Scalar in scalar field
- `PublicKey<C>`: Point on curve
- `Signature`: $(r, s)$ pair in `src/protocol/signing/signature.rs`

---

## 9. Algorithms

### Tonelli-Shanks (Square Roots)

Computes $\sqrt{n}$ in $\mathbb{F}\_p$ when $n$ is a quadratic residue.

**Algorithm:**
1. Check $n^{(p-1)/2} = 1$ (Euler's criterion)
2. Write $p - 1 = Q \cdot 2^S$ where $Q$ is odd
3. Find quadratic non-residue $z$
4. Initialize: $c = z^Q$, $r = n^{(Q+1)/2}$, $t = n^Q$, $m = S$
5. Loop until $t = 1$:
   - Find least $i$ such that $t^{2^i} = 1$
   - $b = c^{2^{(m-i-1)}}$
   - $r = r \cdot b$, $c = b^2$, $t = t \cdot c$, $m = i$

**Special case:** If $p \equiv 3 \pmod{4}$, then $\sqrt{n} = n^{(p+1)/4}$

**In lumen-curve:** `src/algebra/sqrt_mod.rs`

### Karatsuba Multiplication ($\mathbb{F}\_{p^2}$)

For multiplying $(a + bu)(c + du)$:

Standard: 4 multiplications ($ac$, $ad$, $bc$, $bd$)

**Karatsuba:** 3 multiplications
1. $v\_0 = ac$
2. $v\_1 = bd$
3. $v\_2 = (a+b)(c+d)$
4. Result: $(v\_0 - v\_1) + (v\_2 - v\_0 - v\_1)u$

**In lumen-curve:** Used in `Fp2::mul()` in `src/algebra/fields/fp2.rs`

---

## References

### Foundational Papers

1. **Miller (1986)** - "Use of Elliptic Curves in Cryptography"
2. **Boneh, Lynn, Shacham (2001)** - "Short Signatures from the Weil Pairing" (BLS signatures)
3. **Kate, Zaverucha, Goldberg (2010)** - "Constant-Size Commitments to Polynomials" (KZG)

### Curve Standards

4. **IEEE P1363** - Standard for Public-Key Cryptography
5. **SEC 2** - Recommended Elliptic Curve Domain Parameters

### Implementation References

6. **Handbook of Elliptic and Hyperelliptic Curve Cryptography** (2006)
7. **Guide to Pairing-Based Cryptography** (2017)

### Online Resources

8. [Explicit Formulas Database](https://hyperelliptic.org/EFD/) - Optimized curve formulas
9. [Pairing-Friendly Curves](https://www.ietf.org/archive/id/draft-irtf-cfrg-pairing-friendly-curves-11.html) - IETF draft

---

## Glossary

| Term | Definition |
|------|------------|
| **Affine** | Standard $(x, y)$ point representation |
| **Cofactor** | $h = |E(\mathbb{F}\_p)|/r$ where $r$ is the prime subgroup order |
| **Embedding Degree** | Smallest $k$ where $r \mid p^k - 1$ |
| **Generator** | A point $G$ that generates the curve group |
| **Pairing** | Bilinear map from two curve groups to a target group |
| **Projective** | $(X : Y : Z)$ representation avoiding inversions |
| **Scalar Field** | $\mathbb{F}\_r$ where $r$ is the curve group order |
| **SRS** | Structured Reference String (trusted setup output) |
| **Twist** | Isomorphic curve over an extension field |
