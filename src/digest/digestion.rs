use crate::models::DecoyMarking;
use regex::Regex;
use std::collections::HashSet;
use std::ops::Range;
use std::sync::Arc;

// TODO: Reimplement this using the type-state pattern.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DigestSlice {
    ref_seq: Arc<str>,
    range: Range<usize>,
    decoy: DecoyMarking,
}

impl DigestSlice {
    pub fn as_decoy_string(&self) -> String {
        as_decoy_string(&self.ref_seq.as_ref()[self.range.clone()])
    }

    pub fn len(&self) -> usize {
        self.range.len()
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

impl Into<String> for DigestSlice {
    fn into(self) -> String {
        let tmp = &self.ref_seq.as_ref()[self.range.clone()];

        match self.decoy {
            DecoyMarking::Target => tmp.to_string(),
            DecoyMarking::Decoy => as_decoy_string(tmp),
        }
    }
}

#[derive(Debug, Clone)]
pub enum DigestionEnd {
    CTerm,
    NTerm,
}

#[derive(Debug, Clone)]
pub struct DigestionPattern {
    pub regex: Regex,
    pub skip_suffix: Option<char>,
    pub skip_prefix: Option<char>,
}

impl DigestionPattern {
    pub fn trypsin() -> Self {
        DigestionPattern {
            regex: Regex::new("([KR])").unwrap(),
            skip_suffix: Some('P'),
            skip_prefix: None,
        }
    }

    pub fn trypsin_norestriction() -> Self {
        DigestionPattern {
            regex: Regex::new("([KR])").unwrap(),
            skip_suffix: None,
            skip_prefix: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DigestionParameters {
    pub min_length: usize,
    pub max_length: usize,
    pub pattern: DigestionPattern,
    pub digestion_end: DigestionEnd,
    pub max_missed_cleavages: usize,
}

impl DigestionParameters {
    // This section is NEARLY copy-pasted from the Sage implementation.
    // Mike, you rock! sorry about that.
    fn cleavage_sites(&self, sequence: &str) -> Vec<Range<usize>> {
        let mut sites = Vec::new();
        let mut left = 0;
        for mat in self.pattern.regex.find_iter(sequence) {
            let right = match self.digestion_end {
                DigestionEnd::CTerm => mat.end(),
                DigestionEnd::NTerm => mat.start(),
            };

            // Is this needed? Shouldnt I just use the regex?
            // Fun fact ... lookbehinds are not supported so I do need it ...
            if let Some(skip) = self.pattern.skip_suffix {
                if right < sequence.len() && sequence[right..].starts_with(skip) {
                    continue;
                }
            }

            if let Some(skip) = self.pattern.skip_prefix {
                if left > 0 && sequence[left - 1..].ends_with(skip) {
                    continue;
                }
            }
            sites.push(left..right);
            left = right;
        }
        if left < sequence.len() {
            sites.push(left..sequence.len());
        }
        sites
    }

    pub fn digest(&self, sequence: Arc<str>) -> Vec<DigestSlice> {
        let sites = self.cleavage_sites(sequence.as_ref());
        let num_sites = sites.len();
        (0..sites.len())
            .flat_map(|i| {
                let start = sites[i].start;
                let local_out: Vec<DigestSlice> = (0..(self.max_missed_cleavages + 1))
                    .filter_map(|j| {
                        if i + j > num_sites - 1 {
                            return None;
                        }
                        let end = sites[i + j].end;
                        let span = end - start;

                        if span < self.min_length || span > self.max_length {
                            return None;
                        }
                        Some(DigestSlice {
                            ref_seq: sequence.clone(),
                            range: start..end,
                            decoy: DecoyMarking::Target,
                        })
                    })
                    .collect();
                local_out
            })
            .collect()
    }

    pub fn digest_multiple(&self, sequences: &[Arc<str>]) -> Vec<DigestSlice> {
        sequences
            .iter()
            .flat_map(|seq| self.digest(seq.clone()))
            .collect()
    }
}

pub fn as_decoy_string(sequence: &str) -> String {
    let mut sequence = sequence.to_string();
    let inner_rev = sequence[1..(sequence.len() - 1)]
        .chars()
        .rev()
        .collect::<String>();
    sequence.replace_range(1..(sequence.len() - 1), &inner_rev);

    sequence
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cleavage_sites() {
        let params = DigestionParameters {
            min_length: 3,
            max_length: 7,
            pattern: DigestionPattern::trypsin(),
            digestion_end: DigestionEnd::CTerm,
            max_missed_cleavages: 1,
        };
        let seq = "PEPTIKDEPINK";
        let sites = params.cleavage_sites(seq);
        assert_eq!(
            sites.len(),
            2,
            "There should be 2 cleavage sites got: {:?}",
            sites
        );
        assert_eq!(sites[0].start, 0);
        assert_eq!(sites[0].end, 6);
    }

    #[test]
    fn test_digest() {
        let params = DigestionParameters {
            min_length: 3,
            max_length: 7,
            pattern: DigestionPattern::trypsin(),
            digestion_end: DigestionEnd::CTerm,
            max_missed_cleavages: 0,
        };
        let seq: Arc<str> = "PEPTIKDEPINK".into();
        let digests = params.digest(seq);
        assert_eq!(digests.len(), 2);
        assert_eq!(digests[0].len(), 6);
        assert_eq!(Into::<String>::into(digests[0].clone()), "PEPTIK");
        assert_eq!(Into::<String>::into(digests[1].clone()), "DEPINK");
    }

    #[test]
    fn test_digest_nterm() {
        let params = DigestionParameters {
            min_length: 3,
            max_length: 7,
            pattern: DigestionPattern::trypsin(),
            digestion_end: DigestionEnd::NTerm,
            max_missed_cleavages: 1,
        };
        let seq: Arc<str> = "PEPTIKDEPINK".into();
        let digests = params.digest(seq);
        assert_eq!(digests.len(), 3, "Expected 3 digests, got: {:?}", digests);
        assert_eq!(digests[0].len(), 5);
        assert_eq!(Into::<String>::into(digests[0].clone()), "PEPTI");
        assert_eq!(Into::<String>::into(digests[1].clone()), "KDEPIN");
        assert_eq!(Into::<String>::into(digests[2].clone()), "KDEPINK");
    }

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
