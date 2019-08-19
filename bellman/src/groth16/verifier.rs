use ff::PrimeField;
use group::{CurveAffine, CurveProjective};
use pairing::{Engine, PairingCurveAffine};

use super::{PreparedVerifyingKey, Proof, VerifyingKey};

use crate::SynthesisError;

pub fn prepare_verifying_key<E: Engine>(vk: &VerifyingKey<E>) -> PreparedVerifyingKey<E> {
    let mut gamma = vk.gamma_g2;
    gamma.negate();
    let mut delta = vk.delta_g2;
    delta.negate();
    let mut theta = vk.theta_g2;
    theta.negate();

    PreparedVerifyingKey {
        alpha_g1_beta_g2: E::pairing(vk.alpha_g1, vk.beta_g2),
        neg_gamma_g2: gamma.prepare(),
        neg_delta_g2: delta.prepare(),
        neg_theta_g2: theta.prepare(),
        ic: vk.ic.clone(),
    }
}

pub fn verify_proof<'a, E: Engine>(
    pvk: &'a PreparedVerifyingKey<E>,
    proof: &Proof<E>,
    public_inputs: &[E::Fr],
) -> Result<bool, SynthesisError> {
    if (public_inputs.len() + 1) != pvk.ic.len() {
        return Err(SynthesisError::MalformedVerifyingKey);
    }

    let mut acc = pvk.ic[0].into_projective();

    for (i, b) in public_inputs.iter().zip(pvk.ic.iter().skip(1)) {
        acc.add_assign(&b.mul(i.into_repr()));
    }

    // The original verification equation is:
    // A * B = alpha * beta + inputs * gamma + C * delta + D * theta
    // ... however, we rearrange it so that it is:
    // A * B - inputs * gamma - C * delta - D * theta = alpha * beta
    // or equivalently:
    // A * B + inputs * (-gamma) + C * (-delta) + D * (-theta) = alpha * beta
    // which allows us to do a single final exponentiation.

    let acc_prepared = acc.into_affine().prepare();
    let a_prepared = proof.a.prepare();
    let b_prepared = proof.b.prepare();
    let c_prepared = proof.c.prepare();
    let d_prepared = proof.d.map(|d| d.prepare());

    let mut pairings = Vec::with_capacity(4);
    pairings.push((&a_prepared, &b_prepared));
    pairings.push((&acc_prepared, &pvk.neg_gamma_g2));
    pairings.push((&c_prepared, &pvk.neg_delta_g2));
    if let Some(ref d_prepared) = d_prepared {
        pairings.push((d_prepared, &pvk.neg_theta_g2));
    }

    Ok(E::final_exponentiation(&E::miller_loop(pairings.iter())).unwrap() == pvk.alpha_g1_beta_g2)
}
