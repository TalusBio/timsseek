use regex::Regex;
use std::ops::Range;

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

    pub fn digest<'a>(&self, sequence: &'a str) -> Vec<Digest<'a>> {
        let sites = self.cleavage_sites(sequence);
        let num_sites = sites.len();
        (0..sites.len())
            .flat_map(|i| {
                let start = sites[i].start;
                let local_out: Vec<Digest> = (0..(self.max_missed_cleavages + 1))
                    .filter_map(|j| {
                        if i + j > num_sites - 1 {
                            return None;
                        }
                        let end = sites[i + j].end;
                        let span = end - start;

                        if span < self.min_length || span > self.max_length {
                            return None;
                        }
                        Some(Digest {
                            sequence: &sequence[start..end],
                        })
                    })
                    .collect();
                local_out
            })
            .collect()
    }

    pub fn digest_multiple<'a>(&self, sequences: &[&'a str]) -> Vec<Digest<'a>> {
        sequences.iter().flat_map(|seq| self.digest(seq)).collect()
    }
}

#[derive(Debug)]
pub struct Digest<'a> {
    pub sequence: &'a str,
}

impl<'a> Digest<'a> {
    pub fn as_decoy_string(&self) -> String {
        as_decoy_string(self.sequence)
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
        let seq = "PEPTIKDEPINK";
        let digests = params.digest(seq);
        assert_eq!(digests.len(), 2);
        assert_eq!(digests[0].sequence.len(), 6);
        assert_eq!(digests[0].sequence, "PEPTIK");
        assert_eq!(digests[1].sequence, "DEPINK");
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
        let seq = "PEPTIKDEPINK";
        let digests = params.digest(seq);
        assert_eq!(digests.len(), 3, "Expected 3 digests, got: {:?}", digests);
        assert_eq!(digests[0].sequence.len(), 5);
        assert_eq!(digests[0].sequence, "PEPTI");
        assert_eq!(digests[1].sequence, "KDEPIN");
        assert_eq!(digests[2].sequence, "KDEPINK");
    }

    #[test]
    fn test_decoy() {
        let my_digest = Digest {
            sequence: "PEPTIDEPINK",
        };
        let decoy = my_digest.as_decoy_string();
        assert_eq!(my_digest.sequence, "PEPTIDEPINK");
        assert_eq!(decoy, "PNIPEDITPEK");
    }
}
