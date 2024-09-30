use super::fragment_mass_builder::FragmentMassBuilder;
use log::warn;
use rustyms::error::Context;
use rustyms::error::CustomError;
use rustyms::fragment::Position;
use rustyms::LinearPeptide;
use rustyms::MolecularCharge;
use std::collections::HashMap;
use std::ops::{RangeBounds, RangeInclusive};

use rustyms::MultiChemical;
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

impl SequenceToElutionGroupConverter {
    pub fn convert_sequence(
        &self,
        sequence: &str,
        id: u64,
    ) -> Result<Vec<ElutionGroup<Option<Position>>>, CustomError> {
        let mut peptide = LinearPeptide::pro_forma(sequence)?;
        let pep_formulas = peptide.formulas();
        let pep_mono_mass = if pep_formulas.len() > 1 {
            return Err(CustomError::error(
                &"Peptide contains more than one formula.",
                &"",
                Context::none(),
            ));
        } else {
            pep_formulas[0].mass(rustyms::MassMode::Monoisotopic)
        }
        .value;

        let mut out = Vec::new();

        for charge in self.precursor_charge_range.clone() {
            let precursor_mz = (pep_mono_mass + (charge as f64 * PROTON_MASS)) / charge as f64;

            if precursor_mz < self.min_precursor_mz || precursor_mz > self.max_precursor_mz {
                continue;
            }

            peptide = peptide.charge_carriers(Some(MolecularCharge::proton(charge.into())));

            let mut fragment_mzs = self
                .fragment_buildder
                .fragment_mzs_from_linear_peptide(&peptide)?;
            fragment_mzs
                .retain(|(_pos, mz)| *mz > self.min_fragment_mz && *mz < self.max_fragment_mz);

            let mobility = supersimpleprediction(precursor_mz, charge as i32);

            out.push(ElutionGroup {
                id,
                precursor_mz,
                mobility: mobility as f32,
                rt_seconds: 0.0f32,
                precursor_charge: charge,
                fragment_mzs: HashMap::from_iter(fragment_mzs.into_iter()),
            })
        }

        Ok(out)
    }

    pub fn convert_sequences(&self, sequences: &[&str]) -> Vec<ElutionGroup<Option<Position>>> {
        sequences
            .iter()
            .enumerate()
            .flat_map(|(id, sequence)| {
                let tmp = self.convert_sequence(sequence, id as u64);
                match tmp {
                    Ok(x) => Some(x),
                    Err(e) => {
                        warn!("Error converting sequence {}, err: {:?}", sequence, e);
                        None
                    }
                }
            })
            .flatten()
            .collect()
    }
}

// ElutionGroup {
//
// }
// peptide_seq: &str,
// let peptide = LinearPeptide::pro_forma(peptide_seq)?;
//

#[cfg(test)]
mod tests {
    use super::*;
    use rustyms::model::Location;
    use rustyms::model::Model;
    use rustyms::system::{e, f64::MassOverCharge, mass_over_charge::mz, Charge};

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
        let out = converter.convert_sequences(&[seq]);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].precursor_charge, 2);
        assert_eq!(out[1].precursor_charge, 3);
    }
}
