use crate::digest;
use crate::errors::TimsSeekError;
use crate::fragment_mass::fragment_mass_builder::SafePosition;
use crate::models::{
    DecoyMarking,
    DigestSlice,
    NamedQueryChunk,
};
use rayon::prelude::*;
use serde::{
    Deserialize,
    Serialize,
};
use std::path;
use std::sync::Arc;
use timsquery::models::elution_group::ElutionGroup;
use timsrust::TimsRustError;

#[derive(Debug, Clone)]
pub struct Speclib {
    digests: Vec<DigestSlice>,
    charges: Vec<u8>,
    queries: Vec<ElutionGroup<SafePosition>>,
}

pub struct SpeclibIterator {
    speclib: Speclib,
    chunk_size: usize,
    max_iterations: usize,
    iteration_index: usize,
}

impl SpeclibIterator {
    pub fn new(speclib: Speclib, chunk_size: usize) -> Self {
        let max_iters = speclib.digests.len() / chunk_size;
        Self {
            speclib,
            chunk_size,
            max_iterations: max_iters,
            iteration_index: 0,
        }
    }
}

impl Iterator for SpeclibIterator {
    type Item = NamedQueryChunk;

    fn next(&mut self) -> Option<Self::Item> {
        // No need to make decoys when we have a speclib!!
        let out = self
            .speclib
            .get_chunk(self.iteration_index, self.chunk_size);
        self.iteration_index += 1;
        out
    }
}

impl ExactSizeIterator for SpeclibIterator {
    fn len(&self) -> usize {
        self.max_iterations
    }
}

impl Speclib {
    pub fn from_json(json: &str) -> Self {
        let speclib: Vec<SpeclibElement> = serde_json::from_str(json).unwrap();

        let (queries, (charges, digests)): (
            Vec<ElutionGroup<SafePosition>>,
            (Vec<u8>, Vec<DigestSlice>),
        ) = speclib
            .into_par_iter()
            .map(|x| {
                let charge = x.precursor.charge;
                let elution_group = x.elution_group;
                let digest = x.precursor.into();
                (elution_group, (charge, digest))
            })
            .unzip();

        Self {
            digests,
            charges,
            queries,
        }
    }

    pub fn from_ndjson(json: &str) -> Self {
        // Split on newlines and parse each ...
        let lines: Vec<&str> = json.split('\n').collect();
        let mut digests = Vec::new();
        let mut charges = Vec::new();
        let mut queries = Vec::new();

        let mut num_show = 10;
        for line in lines {
            // Continue if the line is empty.
            if line.is_empty() {
                continue;
            }
            let elem: SpeclibElement = match serde_json::from_str(line) {
                Ok(x) => x,
                Err(e) => {
                    panic!("Error parsing line: {:?}", line);
                    // return Err(TimsSeekError::TimsRust(TimsRustError::Serde(e)));
                }
            };

            if num_show > 0 {
                num_show -= 1;
                println!("{:?}", elem);
            }
            charges.push(elem.precursor.charge);
            digests.push(elem.precursor.into());
            queries.push(elem.elution_group);
        }

        if digests.is_empty() {
            panic!("No digests found in speclib file");
        }

        Self {
            digests,
            charges,
            queries,
        }
    }

    pub fn from_ndjson_file(path: &path::Path) -> Result<Self, TimsSeekError> {
        let json = std::fs::read_to_string(path)?;
        Ok(Self::from_ndjson(&json))
    }

    fn get_chunk(&self, chunk_index: usize, chunk_size: usize) -> Option<NamedQueryChunk> {
        let start = chunk_index * chunk_size;
        if start >= self.digests.len() {
            return None;
        }
        let end = start + chunk_size;
        let end = if end > self.digests.len() {
            self.digests.len()
        } else {
            end
        };
        let digests = &self.digests[start..end];
        let charges = &self.charges[start..end];
        let queries = &self.queries[start..end];
        Some(NamedQueryChunk::new(
            digests.to_vec(),
            charges.to_vec(),
            queries.to_vec(),
        ))
    }

    pub fn as_iterator(self, chunk_size: usize) -> SpeclibIterator {
        SpeclibIterator::new(self, chunk_size)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SpeclibElement {
    precursor: PrecursorEntry,
    elution_group: ElutionGroup<SafePosition>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PrecursorEntry {
    sequence: String,
    charge: u8,
    decoy: bool,
}

impl From<PrecursorEntry> for DigestSlice {
    fn from(x: PrecursorEntry) -> Self {
        let decoy = if x.decoy {
            DecoyMarking::ReversedDecoy
        } else {
            DecoyMarking::Target
        };
        let seq: Arc<str> = x.sequence.clone().into();
        let range = 0..seq.as_ref().len();
        DigestSlice::new(seq, range, decoy)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_speclib() {
        let json = r#"[
            {
                "precursor": {
                    "sequence": "PEPTIDEPINK",
                    "charge": 2,
                    "decoy": false
                },
                "elution_group": {
                    "id": 0,
                    "precursor_mzs": [
                        1810.917339999999,
                        1810.917339999999
                    ],
                    "fragment_mzs": {
                        "a1": 123.0,
                        "b1": 123.0,
                        "c1^2": 123.0
                    },
                    "precursor_charge": 2,
                    "mobility": 0.8,
                    "rt_seconds": 0.0,
                    "decoy": false,
                    "expected_precursor_intensity": [
                        1.0,
                        1.0
                    ],
                    "expected_fragment_intensity": {
                        "a1": 1.0,
                        "b1": 1.0,
                        "c1^2": 1.0
                    }
                }
            }
        ]"#;
        let speclib = Speclib::from_json(json);
        assert_eq!(speclib.digests.len(), 1);
        assert_eq!(speclib.charges.len(), 1);
        assert_eq!(speclib.queries.len(), 1);
        println!("{:?}", speclib);

        assert_eq!(speclib.digests[0].decoy, DecoyMarking::Target);
        assert_eq!(speclib.digests[0].len(), 11);
        assert_eq!(speclib.queries[0].fragment_mzs.len(), 3);
    }
}
