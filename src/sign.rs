use crate::calc_sha256_scalar;
use crate::gen_rand_scalar;
use crate::Signature;
use crate::GPK;
use crate::USK;
use alloc::vec::Vec;
use bls12_381::{pairing, G1Projective, G2Projective, Scalar};
use group::{Curve, Group, GroupEncoding};
use rand::RngCore;

pub struct SignCredBeforeHashing {
    a: Scalar,
    ra: Scalar,
    r_delta1: Scalar,
    t1: G1Projective,
    delta1: Scalar,
    r1: G1Projective,
    r4: G1Projective,
}

pub struct SignCredAfterHashing {
    sa: Scalar,
    s_delta1: Scalar,
}

pub fn process_sign_before_hash(
    x: &Scalar,
    rx: &Scalar,
    u: &G1Projective,
    rng: &mut impl RngCore,
) -> SignCredBeforeHashing {
    let a = gen_rand_scalar(rng);
    let ra = gen_rand_scalar(rng);
    let r_delta1 = gen_rand_scalar(rng);
    let t1 = u * a;
    let delta1 = a * x;
    let r1 = u * ra;
    let r4 = t1 * rx + u * -r_delta1;

    SignCredBeforeHashing {
        a,
        ra,
        r_delta1,
        t1,
        delta1,
        r1,
        r4,
    }
}

pub fn process_sign_after_hash(a: &SignCredBeforeHashing, c: &Scalar) -> SignCredAfterHashing {
    let sa = a.ra + c * a.a;
    let s_delta1 = a.r_delta1 + c * a.delta1;

    SignCredAfterHashing { sa, s_delta1 }
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

    let t3 = usk.a + h * (a.a + b.a);

    let a1 = pairing(&t3.to_affine(), &g2.to_affine());
    let a2 = pairing(&h.to_affine(), &w.to_affine());
    let a3 = pairing(&h.to_affine(), &g2.to_affine());

    let r3 = a1 * rx + a2 * (-a.ra - b.ra) + a3 * (-a.r_delta1 - b.r_delta1);

    let mut c: Vec<u8> = vec![];
    c.append(&mut a.t1.to_bytes().as_ref().to_vec());
    c.append(&mut b.t1.to_bytes().as_ref().to_vec());
    c.append(&mut t3.to_bytes().as_ref().to_vec());
    c.append(&mut a.r1.to_bytes().as_ref().to_vec());
    c.append(&mut b.r1.to_bytes().as_ref().to_vec());
    c.append(&mut r3.to_bytes().as_ref().to_vec());
    c.append(&mut a.r4.to_bytes().as_ref().to_vec());
    c.append(&mut b.r4.to_bytes().as_ref().to_vec());

    let c = calc_sha256_scalar(&c);

    let sx = rx + c * x;

    let s_a = process_sign_after_hash(&a, &c);
    let s_b = process_sign_after_hash(&b, &c);

    Signature {
        t1: a.t1,
        t2: b.t1,
        t3,
        c,
        sa: s_a.sa,
        sb: s_b.sa,
        sx,
        s_delta1: s_a.s_delta1,
        s_delta2: s_b.s_delta1,
    }
}
