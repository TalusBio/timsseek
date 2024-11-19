use rustyms::error::{
    Context,
    CustomError,
};
use rustyms::fragment::FragmentType;
use rustyms::model::Location;
use rustyms::spectrum::MassMode;
use rustyms::system::f64::MassOverCharge;
use rustyms::system::mass_over_charge::mz;
use rustyms::system::{
    e,
    Charge,
};
use rustyms::{
    Fragment,
    LinearPeptide,
    Model,
};
use serde::Serialize;
use std::fmt::Display;

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct SafePosition {
    pub series_id: u8,
    pub series_number: u32,
    pub charge: u8,
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
    fn new(x: FragmentType, charge: u8) -> Result<Self, CustomError> {
        let (series_id, series_number) = match x {
            FragmentType::a(position) => (b'a', position.series_number as u32),
            FragmentType::b(position) => (b'b', position.series_number as u32),
            FragmentType::c(position) => (b'c', position.series_number as u32),
            FragmentType::d(position) => (b'd', position.series_number as u32),
            FragmentType::x(position) => (b'x', position.series_number as u32),
            FragmentType::y(position) => (b'y', position.series_number as u32),
            FragmentType::z(position) => (b'z', position.series_number as u32),
            FragmentType::precursor => (0, 0),
            _ => {
                return Err(CustomError::error(
                    "Invalid fragment type",
                    x.to_string(),
                    Context::none(),
                ));
            }
        };

        Ok(Self {
            series_id,
            series_number,
            charge,
        })
    }
}

impl Display for SafePosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}.{}^{}",
            self.series_id as char, self.series_number, self.charge
        )
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
    ) -> Result<Vec<(SafePosition, f64, f32)>, CustomError> {
        // NOTE: I have to add this retain bc it generates precursor ions even if they are not
        // defined.
        let ions: Vec<Fragment> = peptide
            .generate_theoretical_fragments(self.max_charge, &self.model)
            .into_iter()
            .filter(|x| match x.ion {
                FragmentType::precursor => false,
                _ => true,
            })
            .collect();

        // Does this generate ions above the charge of the precursor?
        ions.into_iter()
            .map(|x| {
                let intensity = match x.ion {
                    FragmentType::Y(_) => 1.0,
                    FragmentType::B(_) => 0.5,
                    _ => 0.01,
                };
                Ok((
                    SafePosition::new(x.ion.clone(), x.charge.abs().value as u8)?,
                    x.mz(MassMode::Monoisotopic).value,
                    intensity,
                ))
            })
            .collect()
    }
}
