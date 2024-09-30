use rustyms::error::CustomError;
use rustyms::fragment::Position;
use rustyms::model::Location;
use rustyms::spectrum::MassMode;
use rustyms::system::{e, f64::MassOverCharge, mass_over_charge::mz, Charge};
use rustyms::LinearPeptide;
use rustyms::Model;
use rustyms::*;

#[derive(Debug)]
pub struct FragmentMassBuilder {
    pub model: Model,
    pub max_charge: Charge,
}

impl Default for FragmentMassBuilder {
    fn default() -> Self {
        // let b_ions = Model {
        //     a: (Location::None, Vec::new()),
        //     b: (
        //         Location::SkipNC(2, 2),
        //         vec![],
        //         // Do I want to remove
        //         // neutral losses?
        //     ),
        //     c: (Location::None, Vec::new()),
        //     d: (Location::None, Vec::new()),
        //     v: (Location::None, Vec::new()),
        //     w: (Location::None, Vec::new()),
        //     x: (Location::None, Vec::new()),
        //     y: (Location::None, Vec::new()),
        //     z: (Location::None, Vec::new()),
        //     precursor: vec![],
        //     ppm: MassOverCharge::new::<mz>(20.0),
        //     glycan_fragmentation: None,
        // };

        // let y_ions = Model {
        //     a: (Location::None, Vec::new()),
        //     b: (Location::None, Vec::new()),
        //     c: (Location::None, Vec::new()),
        //     d: (Location::None, Vec::new()),
        //     v: (Location::None, Vec::new()),
        //     w: (Location::None, Vec::new()),
        //     x: (Location::None, Vec::new()),
        //     y: (Location::SkipNC(2, 2), vec![]),
        //     z: (Location::None, Vec::new()),
        //     precursor: vec![],
        //     ppm: MassOverCharge::new::<mz>(20.0),
        //     glycan_fragmentation: None,
        // };
        let by_ions = Model {
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
        };
        let max_charge: Charge = Charge::new::<e>(2.0);
        Self {
            model: by_ions,
            max_charge,
        }
    }
}

impl FragmentMassBuilder {
    pub fn fragment_mzs_from_linear_peptide(
        &self,
        peptide: &LinearPeptide,
        // peptide_seq: &str,
    ) -> Result<Vec<(Option<Position>, f64)>, CustomError> {
        // let peptide = LinearPeptide::pro_forma(peptide_seq)?;
        let ions = peptide.generate_theoretical_fragments(self.max_charge, &self.model);
        // Does this generate ions above the charge of the precursor?
        Ok(ions
            .into_iter()
            .map(|x| {
                (
                    x.ion.position().copied(),
                    x.mz(MassMode::Monoisotopic).value,
                )
            })
            .collect::<Vec<(Option<Position>, f64)>>())
    }
}
