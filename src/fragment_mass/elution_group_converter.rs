use super::fragment_mass_builder::FragmentMassBuilder;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::isotopes::peptide_isotopes;
use log::{
    error,
    warn,
};
use rayon::prelude::*;
use rustyms::error::{
    Context,
    CustomError,
};
use rustyms::{
    LinearPeptide,
    MolecularCharge,
    MolecularFormula,
    MultiChemical,
};
use std::collections::HashMap;
use std::ops::RangeInclusive;
use timsquery::models::elution_group::ElutionGroup;

/// Super simple 1/k0 prediction.
///
/// This is a simple prediction of the retention time based on the m/z and charge.
/// On my data it gets MAPE 1.82802 so, this prediction + 10% error is a pretty solid way
/// to set an extraction window for mobility if you dont know anything for the peptide.
///
/// Example:
/// ```
/// use timsseek::fragment_mass::elution_group_converter::supersimpleprediction;
/// let mass = 1810.917339999999;
/// let charge = 2;
/// let out = supersimpleprediction(mass / charge as f64, charge);
/// assert!((out - 1.105151).abs() < 0.001 );
/// ```
pub fn supersimpleprediction(mz: f64, charge: i32) -> f64 {
    let intercept_ = -1.660e+00;
    let log1p_mz = (mz + 1.).ln();
    let sq_mz_over_charge = mz.powi(2) / charge as f64;
    let log1p_sq_mz_over_charge = (sq_mz_over_charge + 1.).ln();

    intercept_
        + (-3.798e-01 * log1p_mz)
        + (-2.389e-04 * mz)
        + (3.957e-01 * log1p_sq_mz_over_charge)
        + (4.157e-07 * sq_mz_over_charge)
        + (1.417e-01 * charge as f64)
}

#[derive(Debug)]
pub struct SequenceToElutionGroupConverter {
    pub precursor_charge_range: RangeInclusive<u8>,
    pub fragment_buildder: FragmentMassBuilder,
    pub max_precursor_mz: f64,
    pub min_precursor_mz: f64,
    pub max_fragment_mz: f64,
    pub min_fragment_mz: f64,
}

impl Default for SequenceToElutionGroupConverter {
    fn default() -> Self {
        Self {
            precursor_charge_range: 2..=3,
            fragment_buildder: FragmentMassBuilder::default(),
            max_precursor_mz: 1000.,
            min_precursor_mz: 400.,
            max_fragment_mz: 2000.,
            min_fragment_mz: 200.,
        }
    }
}

const PROTON_MASS: f64 = 1.007276466;

// TODO: Find right way ...
const NEUTRON_MASS: f64 = 1.00;

fn count_carbon_sulphur(form: &MolecularFormula) -> (u16, u16) {
    let mut ncarbon = 0;
    let mut nsulphur = 0;

    for (elem, count, _extras) in form.elements() {
        match (elem, count) {
            (&rustyms::Element::C, Some(cnt)) => {
                ncarbon += *cnt;
            }
            (&rustyms::Element::S, Some(cnt)) => {
                nsulphur += *cnt;
            }
            _ => {}
        }
    }

    (ncarbon, nsulphur)
}

impl SequenceToElutionGroupConverter {
    pub fn convert_sequence(
        &self,
        sequence: &str,
        id: u64,
    ) -> Result<(Vec<ElutionGroup<SafePosition>>, Vec<u8>), CustomError> {
        let mut peptide = LinearPeptide::pro_forma(sequence)?;
        let pep_formulas = peptide.formulas();
        let (pep_mono_mass, pep_formula) = if pep_formulas.len() > 1 {
            return Err(CustomError::error(
                "Peptide contains more than one formula.",
                "",
                Context::none(),
            ));
        } else {
            let form = pep_formulas[0].clone();
            let mono_mass = pep_formulas[0].mass(rustyms::MassMode::Monoisotopic);
            (mono_mass.value, form)
        };
        let (ncarbon, nsulphur) = count_carbon_sulphur(&pep_formula);
        let pep_isotope = peptide_isotopes(ncarbon, nsulphur);
        let mut expected_prec_inten = vec![1e-3f32; 4];

        for (ii, isot) in pep_isotope.iter().enumerate() {
            expected_prec_inten[1 + ii] = *isot
        }

        let mut out = Vec::new();
        let mut out_charges = Vec::new();

        for charge in self.precursor_charge_range.clone() {
            // Q: Why am I adding the charge here manually instead of using the calculator in the
            // Formula?
            let precursor_mz = (pep_mono_mass + (charge as f64 * PROTON_MASS)) / charge as f64;
            let nmf = NEUTRON_MASS / (charge as f64);

            if precursor_mz < self.min_precursor_mz || precursor_mz > self.max_precursor_mz {
                continue;
            }

            peptide = peptide.charge_carriers(Some(MolecularCharge::proton(charge.into())));

            let mut fragment_mzs = self
                .fragment_buildder
                .fragment_mzs_from_linear_peptide(&peptide)?;
            fragment_mzs
                .retain(|(_pos, mz, _)| *mz > self.min_fragment_mz && *mz < self.max_fragment_mz);

            let mobility = supersimpleprediction(precursor_mz, charge as i32);
            let mut precursor_mzs = vec![precursor_mz; 4];
            precursor_mzs[0] -= nmf;
            precursor_mzs[2] += nmf;
            precursor_mzs[3] += 2. * nmf;

            let fragment_expect_inten =
                HashMap::from_iter(fragment_mzs.iter().map(|(k, _, v)| (*k, *v)));
            let fragment_mzs = HashMap::from_iter(fragment_mzs.iter().map(|(k, v, _)| (*k, *v)));

            out.push(ElutionGroup {
                id,
                precursor_mzs,
                mobility: mobility as f32,
                rt_seconds: 0.0f32,
                // precursor_charge: charge,
                fragment_mzs,
                expected_fragment_intensity: Some(fragment_expect_inten),
                expected_precursor_intensity: Some(expected_prec_inten.clone()),
            });
            out_charges.push(charge);
        }

        Ok((out, out_charges))
    }

    pub fn convert_sequences<'a, T: AsRef<str> + Send + Sync + std::fmt::Debug>(
        &self,
        sequences: &'a [T],
    ) -> Result<(Vec<&'a T>, Vec<ElutionGroup<SafePosition>>, Vec<u8>), CustomError> {
        let (seqs, (eg, crg)) = sequences
            .par_iter()
            .enumerate()
            .flat_map(|(id, sequence)| {
                let tmp = self.convert_sequence(sequence.as_ref(), id as u64);
                match tmp {
                    Ok(x) => {
                        let expanded_sequence: Vec<&T> =
                            (0..(x.0.len())).map(|_x| sequence).collect();
                        Some((expanded_sequence, (x.0, x.1)))
                    }
                    Err(e) => {
                        warn!("Error converting sequence {:?}, err: {:?}", sequence, e);
                        None
                    }
                }
            })
            .flatten()
            .collect();
        Ok((seqs, eg, crg))
    }

    pub fn convert_enumerated_sequences<'a, T: AsRef<str> + Send + Sync + std::fmt::Debug>(
        &self,
        enum_sequences: &'a [(usize, T)],
    ) -> Result<(Vec<&'a T>, Vec<ElutionGroup<SafePosition>>, Vec<u8>), CustomError> {
        let (seqs, (eg, crg)) = enum_sequences
            .par_iter()
            .flat_map(|(i, s)| {
                let tmp = self.convert_sequence(s.as_ref(), *i as u64);
                match tmp {
                    Ok(x) => {
                        let expanded_sequence: Vec<&T> = (0..(x.0.len())).map(|_x| s).collect();
                        Some((expanded_sequence, (x.0, x.1)))
                    }
                    Err(e) => {
                        error!("Error converting sequence {:?}, err: {:?}", s, e);
                        None
                    }
                }
            })
            .flatten()
            .collect();
        Ok((seqs, eg, crg))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustyms::model::{
        Location,
        Model,
    };
    use rustyms::system::f64::MassOverCharge;
    use rustyms::system::mass_over_charge::mz;
    use rustyms::system::{
        e,
        Charge,
    };

    #[test]
    fn test_converter() {
        let seq = "PEPTIDEPINK/2";
        let fragment_mass_builder = FragmentMassBuilder::default();
        let converter = SequenceToElutionGroupConverter {
            precursor_charge_range: 2..=3,
            fragment_buildder: FragmentMassBuilder {
                model: Model {
                    a: (Location::None, Vec::new()),
                    b: (Location::SkipNC(2, 2), vec![]),
                    c: (Location::None, Vec::new()),
                    d: (Location::None, Vec::new()),
                    v: (Location::None, Vec::new()),
                    w: (Location::None, Vec::new()),
                    x: (Location::None, Vec::new()),
                    y: (Location::SkipNC(2, 2), vec![]),
                    z: (Location::None, Vec::new()),
                    precursor: vec![],
                    ppm: MassOverCharge::new::<mz>(20.0),
                    glycan_fragmentation: None,
                },
                max_charge: Charge::new::<e>(2.0),
            },
            max_precursor_mz: 1000.,
            min_precursor_mz: 400.,
            max_fragment_mz: 2000.,
            min_fragment_mz: 200.,
        };
        let seq_slc = vec![seq];
        let out = converter.convert_sequences(&seq_slc).unwrap();
        assert_eq!(out.0.len(), 2);
    }
}
