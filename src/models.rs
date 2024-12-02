use crate::fragment_mass::fragment_mass_builder::SafePosition;
use rayon::iter::Zip as RayonZip;
use rayon::prelude::*;
use rayon::vec::IntoIter as RayonVecIntoIter;
use serde::{
    Deserialize,
    Serialize,
};
use std::collections::HashSet;
use std::ops::Range;
use std::sync::Arc;
use timsquery::models::elution_group::ElutionGroup;

/// The different labels that denote if a sequence is a decoy or not.
///
/// NOTE: The main difference between the decoy and reversed decoy is that the reversed decoy
/// has already been reversed, thus converting it to a string can be done as-is.
#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, std::hash::Hash, PartialOrd, Ord)]
pub enum DecoyMarking {
    Target,
    Decoy,
    ReversedDecoy,
}
impl DecoyMarking {
    pub fn as_str(&self) -> &'static str {
        match self {
            DecoyMarking::Target => "Target",
            DecoyMarking::Decoy => "Decoy",
            DecoyMarking::ReversedDecoy => "Decoy",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DigestSlice {
    ref_seq: Arc<str>,
    range: Range<usize>,
    pub decoy: DecoyMarking,
}

impl Serialize for DigestSlice {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let local_str = Into::<String>::into(self.clone());
        serializer.serialize_str(local_str.as_str())
    }
}

impl DigestSlice {
    pub fn new(ref_seq: Arc<str>, range: Range<usize>, decoy: DecoyMarking) -> Self {
        Self {
            ref_seq,
            range,
            decoy,
        }
    }

    pub fn as_decoy(&self) -> DigestSlice {
        DigestSlice {
            ref_seq: self.ref_seq.clone(),
            range: self.range.clone(),
            decoy: DecoyMarking::Decoy,
        }
    }

    pub fn as_decoy_string(&self) -> String {
        as_decoy_string(&self.ref_seq.as_ref()[self.range.clone()])
    }

    pub fn len(&self) -> usize {
        self.range.len()
    }

    pub fn is_empty(&self) -> bool {
        self.range.is_empty()
    }
}

pub fn deduplicate_digests(mut digest_slices: Vec<DigestSlice>) -> Vec<DigestSlice> {
    let mut seen = HashSet::new();
    digest_slices.retain(|x| {
        let local_str: String = x.clone().into();
        let is_first = !seen.contains(&local_str);
        seen.insert(local_str);
        is_first
    });
    digest_slices
}

impl From<DigestSlice> for String {
    fn from(x: DigestSlice) -> Self {
        let tmp = &x.ref_seq.as_ref()[x.range.clone()];

        match x.decoy {
            DecoyMarking::Target => tmp.to_string(),
            DecoyMarking::ReversedDecoy => tmp.to_string(),
            DecoyMarking::Decoy => as_decoy_string(tmp),
        }
    }
}

fn as_decoy_string(sequence: &str) -> String {
    let mut sequence = sequence.to_string();
    let inner_rev = sequence[1..(sequence.len() - 1)]
        .chars()
        .rev()
        .collect::<String>();
    sequence.replace_range(1..(sequence.len() - 1), &inner_rev);

    sequence
}

#[derive(Debug, Clone)]
pub struct NamedQueryChunk {
    digests: Vec<DigestSlice>,
    charges: Vec<u8>,
    pub queries: Vec<ElutionGroup<SafePosition>>,
}

impl NamedQueryChunk {
    pub fn new(
        digests: Vec<DigestSlice>,
        charges: Vec<u8>,
        queries: Vec<ElutionGroup<SafePosition>>,
    ) -> Self {
        assert_eq!(digests.len(), charges.len());
        assert_eq!(digests.len(), queries.len());
        Self {
            digests,
            charges,
            queries,
        }
    }

    pub fn into_zip_par_iter(
        self,
    ) -> RayonZip<
        RayonVecIntoIter<ElutionGroup<SafePosition>>,
        RayonZip<RayonVecIntoIter<DigestSlice>, RayonVecIntoIter<u8>>,
    > {
        // IN THEORY I should implement IntoIter for this struct
        // but I failed at it (skill issues?) so this will do for now.
        // JSPP - 2024-11-21
        self.queries.into_par_iter().zip(
            self.digests
                .into_par_iter()
                .zip(self.charges.into_par_iter()),
        )
    }

    pub fn len(&self) -> usize {
        self.queries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.queries.is_empty()
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decoy() {
        let seq: Arc<str> = "PEPTIDEPINK".into();
        let my_digest = DigestSlice {
            ref_seq: seq.clone(),
            range: 0..seq.as_ref().len(),
            decoy: DecoyMarking::Target,
        };
        let decoy = my_digest.as_decoy_string();
        assert_eq!(Into::<String>::into(my_digest.clone()), "PEPTIDEPINK");
        assert_eq!(Into::<String>::into(decoy.clone()), "PNIPEDITPEK");
    }

    #[test]
    fn test_deduplicate_digests() {
        let seq: Arc<str> = "PEPTIDEPINKTOMATOTOMATO".into();
        let seq2: Arc<str> = "PEPTIDEPINKTOMATO".into();
        let digests: Vec<DigestSlice> = vec![
            DigestSlice {
                ref_seq: seq.clone(),
                range: 0..seq.as_ref().len(),
                decoy: DecoyMarking::Target,
            },
            DigestSlice {
                ref_seq: seq.clone(),
                range: 0..seq2.as_ref().len(), // Note the short length
                decoy: DecoyMarking::Target,
            },
            DigestSlice {
                ref_seq: seq2.clone(),
                range: 0..seq2.as_ref().len(),
                decoy: DecoyMarking::Target,
            },
        ];
        let deduped = deduplicate_digests(digests);
        assert_eq!(deduped.len(), 2);
        assert_eq!(deduped[0].len(), seq.as_ref().len());
        assert_eq!(deduped[1].len(), seq2.as_ref().len());
    }
}
