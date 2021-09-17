use crate::calc_sha256_scalar;
use crate::gen_rand_scalar;
use crate::Signature;
use crate::GPK;
use crate::USK;
use alloc::vec::Vec;
use bls12_381::{pairing, G1Projective, Scalar};
use group::{Curve, GroupEncoding};
use rand::RngCore;

pub struct SignCredBeforeHashing {
    y: Scalar,
    r: Scalar,
    r_delta: Scalar,
    t: G1Projective,
    delta: Scalar,
    r_first: G1Projective,
    r_second: G1Projective,
}

pub struct SignCredAfterHashing {
    s: Scalar,
    s_delta: Scalar,
}

pub fn process_sign_before_hash(
    x: &Scalar,
    rx: &Scalar,
    pk: &G1Projective,
    rng: &mut impl RngCore,
) -> SignCredBeforeHashing {
    let y = gen_rand_scalar(rng);
    let r = gen_rand_scalar(rng);
    let r_delta = gen_rand_scalar(rng);
    let t = pk * y;
    let delta = y * x;
    let r_first = pk * r;
    let r_second = t * rx + pk * -r_delta;

    SignCredBeforeHashing {
        y,
        r,
        r_delta,
        t,
        delta,
        r_first,
        r_second,
    }
}

pub fn process_sign_after_hash(y: &SignCredBeforeHashing, c: &Scalar) -> SignCredAfterHashing {
    let s = y.r + c * y.y;
    let s_delta = y.r_delta + c * y.delta;

    SignCredAfterHashing { s, s_delta }
}

pub fn sign(usk: &USK, gpk: &GPK, rng: &mut impl RngCore) -> Signature {
    let USK { a: _, x } = usk;
    let GPK {
        h,
        u,
        v,
        w,
        g1: _,
        g2,
    } = gpk;

    let rx = gen_rand_scalar(rng);
    let a = process_sign_before_hash(x, &rx, u, rng);
    let b = process_sign_before_hash(x, &rx, v, rng);

    let t3 = usk.a + h * (a.y + b.y);

    let a1 = pairing(&t3.to_affine(), &g2.to_affine());
    let a2 = pairing(&h.to_affine(), &w.to_affine());
    let a3 = pairing(&h.to_affine(), &g2.to_affine());

    let r3 = a1 * rx + a2 * (-a.r - b.r) + a3 * (-a.r_delta - b.r_delta);

    let mut c: Vec<u8> = vec![];
    c.append(&mut a.t.to_bytes().as_ref().to_vec());
    c.append(&mut b.t.to_bytes().as_ref().to_vec());
    c.append(&mut t3.to_bytes().as_ref().to_vec());
    c.append(&mut a.r_first.to_bytes().as_ref().to_vec());
    c.append(&mut b.r_first.to_bytes().as_ref().to_vec());
    c.append(&mut r3.to_bytes().as_ref().to_vec());
    c.append(&mut a.r_second.to_bytes().as_ref().to_vec());
    c.append(&mut b.r_second.to_bytes().as_ref().to_vec());

    let c = calc_sha256_scalar(&c);

    let sx = rx + c * x;

    let s_a = process_sign_after_hash(&a, &c);
    let s_b = process_sign_after_hash(&b, &c);

    Signature {
        t1: a.t,
        t2: b.t,
        t3,
        c,
        sa: s_a.s,
        sb: s_b.s,
        sx,
        s_delta1: s_a.s_delta,
        s_delta2: s_b.s_delta,
    }
}
