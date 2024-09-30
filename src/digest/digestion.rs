use regex::Regex;
use std::ops::Range;

pub enum DigestionEnd {
    CTerm,
    NTerm,
}

struct DigestionPattern {
    pub regex: Regex,
    pub skip_suffix: Option<char>,
    pub skip_prefix: Option<char>,
}

impl DigestionPattern {
    fn trypsin() -> Self {
        DigestionPattern {
            regex: Regex::new("([KR])").unwrap(),
            skip_suffix: Some('P'),
            skip_prefix: None,
        }
    }

    fn trypsin_norestriction() -> Self {
        DigestionPattern {
            regex: Regex::new("([KR])").unwrap(),
            skip_suffix: None,
            skip_prefix: None,
        }
    }
}

struct DigestionParameters {
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

                        println!("start: {:?}, end: {:?}", start, end);
                        if span < self.min_length || span > self.max_length {
                            println!("span: {:?}", span);
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
}

#[derive(Debug)]
struct Digest<'a> {
    pub sequence: &'a str,
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
}
