#![no_std]
#[macro_use]
extern crate slice_as_array;
#[macro_use]
extern crate alloc;

pub mod sign;

use alloc::vec::Vec;
use bls12_381::{pairing, G1Projective, G2Projective, Scalar};
use byteorder::{BigEndian, ByteOrder};
use ff::Field;
use group::{Curve, Group, GroupEncoding};
use rand::RngCore;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

#[derive(Serialize, Deserialize)]
pub struct GSK {
    pub xi_1: Scalar,
    pub xi_2: Scalar,
    pub gamma: Scalar,
}

#[derive(Serialize, Deserialize)]
pub struct GPK {
    pub h: G1Projective,
    pub u: G1Projective,
    pub v: G1Projective,
    pub w: G2Projective,
    pub g1: G1Projective,
    pub g2: G2Projective,
}

pub struct SetUpResult {
    pub gpk: GPK,
    pub gsk: GSK,
}

#[derive(Serialize, Deserialize)]
pub struct USK {
    pub x: Scalar,
    pub a: G1Projective,
}

#[derive(Serialize, Deserialize)]
pub struct Signature {
    pub t1: G1Projective,
    pub t2: G1Projective,
    pub t3: G1Projective,
    pub c: Scalar,
    pub sa: Scalar,
    pub sb: Scalar,
    pub sx: Scalar,
    pub s_delta1: Scalar,
    pub s_delta2: Scalar,
}

pub fn setup(rng: &mut impl RngCore) -> SetUpResult {
    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    let xi_1 = gen_rand_scalar(rng);
    let xi_2 = gen_rand_scalar(rng);

    let tmp = gen_rand_g1(rng);

    let u = tmp * xi_1;
    let v = tmp * xi_2;

    let h = u * xi_2;

    let gamma = gen_rand_scalar(rng);
    let w = g2 * gamma;

    let gsk = GSK { xi_1, xi_2, gamma };
    let gpk = GPK { h, u, v, w, g1, g2 };

    SetUpResult { gsk, gpk }
}

pub fn issue(gsk: &GSK, gpk: &GPK, rng: &mut impl RngCore) -> USK {
    let x = gen_rand_scalar(rng);
    let tmp = (gsk.gamma + x).invert().unwrap();
    let a = gpk.g1 * tmp;

    USK { a, x }
}

pub fn verify(signature: &Signature, gpk: &GPK) -> Result<(), ()> {
    let Signature {
        t1,
        t2,
        t3,
        c,
        sa,
        sb,
        sx,
        s_delta1,
        s_delta2,
    } = signature;

    let GPK { h, u, v, w, g1, g2 } = gpk;

    let r1_v = u * sa + t1 * -c;
    let r2_v = v * sb + t2 * -c;

    let a1_v = pairing(&t3.to_affine(), &g2.to_affine());
    let a2_v = pairing(&h.to_affine(), &w.to_affine());
    let a3_v = pairing(&h.to_affine(), &g2.to_affine());
    let a4_v = pairing(&t3.to_affine(), &w.to_affine());
    let a5_v = pairing(&g1.to_affine(), &g2.to_affine());

    let r3_v = a1_v * sx + a2_v * (-sa - sb) + a3_v * (-s_delta1 - s_delta2) + (a4_v - a5_v) * c;

    let r4_v = t1 * sx + u * -s_delta1;
    let r5_v = t2 * sx + v * -s_delta2;

    let mut c_v: Vec<u8> = vec![];
    c_v.append(&mut t1.to_bytes().as_ref().to_vec());
    c_v.append(&mut t2.to_bytes().as_ref().to_vec());
    c_v.append(&mut t3.to_bytes().as_ref().to_vec());
    c_v.append(&mut r1_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut r2_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut r3_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut r4_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut r5_v.to_bytes().as_ref().to_vec());

    let c_v = calc_sha256_scalar(&c_v);

    if c_v == *c {
        Ok(())
    } else {
        Err(())
    }
}

pub fn is_signed_member(usk: &USK, signature: &Signature, gsk: &GSK) -> bool {
    let a_v = signature.t3 - (signature.t1 * gsk.xi_2 + signature.t2 * gsk.xi_1);

    usk.a == a_v
}

fn calc_sha256_scalar(vec: &[u8]) -> Scalar {
    let mut hasher = Sha256::new();
    hasher.update(vec);
    let hashed = hasher.finalize().to_vec();

    let mut schalar: Vec<u64> = vec![0; hashed.len() / 8];
    BigEndian::read_u64_into(&hashed, &mut schalar);
    let schalar = slice_as_array!(&schalar, [u64; 4]).unwrap();

    Scalar::from_raw(*schalar)
}

fn gen_rand_scalar(rng: &mut impl RngCore) -> Scalar {
    Scalar::random(rng)
}

fn gen_rand_g1(rng: &mut impl RngCore) -> G1Projective {
    G1Projective::random(rng)
}

#[test]
fn test_all() {
    use rand::thread_rng;
    let mut rng = thread_rng();
    let SetUpResult { gsk, gpk } = setup(&mut rng);
    let usk = issue(&gsk, &gpk, &mut rng);
    let sig = sign::sign(&usk, &gpk, &mut rng);
    verify(&sig, &gpk).unwrap();

    assert!(is_signed_member(&usk, &sig, &gsk));
}
