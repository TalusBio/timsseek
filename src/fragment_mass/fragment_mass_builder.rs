use rustyms::error::Context;
use rustyms::error::CustomError;
use rustyms::fragment::FragmentType;
use rustyms::model::Location;
use rustyms::spectrum::MassMode;
use rustyms::system::{e, f64::MassOverCharge, mass_over_charge::mz, Charge};
use rustyms::LinearPeptide;
use rustyms::Model;
use serde::Serialize;
use std::fmt::Display;

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SafePosition {
    pub series_id: u8,
    pub series_number: u32,
}

impl Serialize for SafePosition {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&format!("{}", self))
    }
}

impl SafePosition {
    fn new(x: FragmentType) -> Result<Self, CustomError> {
        let (series_id, series_number) = match x {
            FragmentType::a(position) => ('a' as u8, position.series_number as u32),
            FragmentType::b(position) => ('b' as u8, position.series_number as u32),
            FragmentType::c(position) => ('c' as u8, position.series_number as u32),
            FragmentType::d(position) => ('d' as u8, position.series_number as u32),
            FragmentType::x(position) => ('x' as u8, position.series_number as u32),
            FragmentType::y(position) => ('y' as u8, position.series_number as u32),
            FragmentType::z(position) => ('z' as u8, position.series_number as u32),
            FragmentType::precursor => (0, 0),
            _ => {
                return Err(CustomError::error(
                    "Invalid fragment type",
                    x.to_string(),
                    Context::none(),
                ))
            }
        };

        Ok(Self {
            series_id,
            series_number,
        })
    }
}

impl Display for SafePosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.{}", self.series_id as char, self.series_number)
    }
}

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
    ) -> Result<Vec<(SafePosition, f64)>, CustomError> {
        let ions = peptide.generate_theoretical_fragments(self.max_charge, &self.model);
        // Does this generate ions above the charge of the precursor?
        ions.into_iter()
            .map(|x| {
                Ok((
                    // x.ion.position().copied(),
                    SafePosition::new(x.ion.clone())?,
                    x.mz(MassMode::Monoisotopic).value,
                ))
            })
            .collect()
    }
}
