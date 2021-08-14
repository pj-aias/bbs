#[macro_use]
extern crate slice_as_array;

use bls12_381::{pairing, G1Projective, G2Projective, Scalar};
use byteorder::{BigEndian, ByteOrder};
use ff::Field;
use group::{Curve, Group, GroupEncoding};
use rand::RngCore;
use sha2::{Digest, Sha256};

fn main() {
    use rand::thread_rng;
    let mut rng = thread_rng();

    let g1 = G1Projective::generator();
    let g2 = G2Projective::generator();

    ////// Set up
    let xi_1 = Scalar::random(&mut rng);
    let xi_2 = Scalar::random(&mut rng);

    let tmp = G1Projective::random(&mut rng);

    let u = tmp * xi_1;
    let v = tmp * xi_2;

    let h = u * xi_2;

    let gamma = Scalar::random(&mut rng);
    let w = G2Projective::generator() * gamma;

    ////// Issue
    let x = Scalar::random(&mut rng);
    let tmp = (gamma + x).invert().unwrap();

    let Ai = g1 * tmp;

    ////// Sign
    let a = Scalar::random(&mut rng);
    let b = Scalar::random(&mut rng);

    let ra = Scalar::random(&mut rng);
    let rb = Scalar::random(&mut rng);
    let rx = Scalar::random(&mut rng);
    let r_delta1 = Scalar::random(&mut rng);
    let r_delta2 = Scalar::random(&mut rng);

    let T1 = u * a;
    let T2 = v * b;
    let T3 = Ai + h * (a + b);

    let delta1 = a * x;
    let delta2 = b * x;

    let R1 = u * ra;
    let R2 = v * rb;

    let a1 = pairing(&T3.to_affine(), &g2.to_affine());
    let a2 = pairing(&h.to_affine(), &w.to_affine());
    let a3 = pairing(&h.to_affine(), &g2.to_affine());

    let R3 = a1 * rx + a2 * (-ra - rb) + a3 * (-r_delta1 - r_delta2);

    let R4 = T1 * rx + u * -r_delta1;
    let R5 = T2 * rx + v * -r_delta2;

    let mut c: Vec<u8> = vec![];
    c.append(&mut T1.to_bytes().as_ref().to_vec());
    c.append(&mut T2.to_bytes().as_ref().to_vec());
    c.append(&mut T3.to_bytes().as_ref().to_vec());
    c.append(&mut R1.to_bytes().as_ref().to_vec());
    c.append(&mut R2.to_bytes().as_ref().to_vec());
    c.append(&mut R3.to_bytes().as_ref().to_vec());
    c.append(&mut R4.to_bytes().as_ref().to_vec());
    c.append(&mut R5.to_bytes().as_ref().to_vec());

    let c = calc_sha256_scalar(&c);

    let sa = ra + c * a;
    let sb = rb + c * b;
    let sx = rx + c * x;
    let s_delta1 = r_delta1 + c * delta1;
    let s_delta2 = r_delta2 + c * delta2;

    ////// Verify
    let R1_v = u * sa + T1 * -c;
    let R2_v = v * sb + T2 * -c;

    let a1_v = pairing(&T3.to_affine(), &g2.to_affine());
    let a2_v = pairing(&h.to_affine(), &w.to_affine());
    let a3_v = pairing(&h.to_affine(), &g2.to_affine());
    let a4_v = pairing(&T3.to_affine(), &w.to_affine());
    let a5_v = pairing(&g1.to_affine(), &g2.to_affine());

    let R3_v = a1_v * sx + a2_v * (-sa - sb) + a3_v * (-s_delta1 - s_delta2) + (a4_v - a5_v) * c;

    let R4_v = T1 * sx + u * -s_delta1;
    let R5_v = T2 * sx + v * -s_delta2;

    let mut c_v: Vec<u8> = vec![];
    c_v.append(&mut T1.to_bytes().as_ref().to_vec());
    c_v.append(&mut T2.to_bytes().as_ref().to_vec());
    c_v.append(&mut T3.to_bytes().as_ref().to_vec());
    c_v.append(&mut R1_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut R2_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut R3_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut R4_v.to_bytes().as_ref().to_vec());
    c_v.append(&mut R5_v.to_bytes().as_ref().to_vec());

    let c_v = calc_sha256_scalar(&c_v);

    assert_eq!(c, c_v);

    ////// Open
    let A_v = T3 - (T1 * xi_2 + T2 * xi_1);
    assert_eq!(Ai, A_v);
}

pub fn calc_sha256_scalar(vec: &[u8]) -> Scalar {
    let mut hasher = Sha256::new();
    hasher.update(vec);
    let hashed = hasher.finalize().to_vec();

    let mut schalar: Vec<u64> = vec![0; hashed.len() / 8];
    BigEndian::read_u64_into(&hashed, &mut schalar);
    let schalar = slice_as_array!(&schalar, [u64; 4]).unwrap();

    Scalar::from_raw(*schalar)
}
