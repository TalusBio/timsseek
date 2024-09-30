use super::models::{ProteinSequence, ProteinSequenceBuilder};
use log::*;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

#[derive(Debug)]
pub struct ProteinSequenceCollection {
    pub sequences: Vec<ProteinSequence>,
}

#[derive(Debug)]
pub struct ProteinSequenceNmerIndex {
    nmer_size: usize,
    // Q: Does the hashmap store the string or just the hash?
    index: HashMap<Arc<[u8]>, Vec<usize>>,
    sequences: Vec<ProteinSequence>,
}

impl ProteinSequenceNmerIndex {
    pub fn new(nmer_size: usize, sequences: Vec<ProteinSequence>) -> Self {
        let st = Instant::now();
        let mut index = HashMap::new();
        for (curr_id, sequence) in sequences.iter().enumerate() {
            // let curr_id = sequence.id;
            let sequence = sequence.sequence.as_bytes();

            sequence.windows(nmer_size).for_each(|window| {
                let key = Arc::from(window);
                index
                    .entry(key)
                    .and_modify(|e: &mut Vec<usize>| {
                        e.push(curr_id);
                    })
                    .or_insert(vec![curr_id]);
            });
        }
        let elapsed = st.elapsed();
        info!("Indexing took {:#?}", elapsed);

        Self {
            nmer_size,
            index,
            sequences,
        }
    }

    pub fn from_collection(collection: ProteinSequenceCollection, nmer_size: usize) -> Self {
        Self::new(nmer_size, collection.sequences)
    }

    pub fn query_sequences(&self, query: &[u8]) -> Option<Vec<usize>> {
        let first_window = query.get(0..self.nmer_size)?;
        let key = Arc::from(first_window);
        let mut options = self.index.get(&key)?.to_vec();
        for window in query.windows(self.nmer_size) {
            if options.is_empty() {
                return None;
            }
            let key = Arc::from(window);
            let local_options = self.index.get(&key);

            match local_options {
                Some(local_options) => {
                    options.retain(|&id| local_options.contains(&id));
                }
                None => {
                    return None;
                }
            }
        }
        // Finally filter for the full sequence being contained.
        // For instance if the nmer is 2 and the query seq is "FOOPP", it will
        // match "FOP" (wrong) and "FOOOP" (correct)
        // And we want to preseve only the later.
        options.retain(|&id| {
            self.sequences[id]
                .sequence
                .as_bytes()
                .windows(query.len())
                .any(|w| w == query)
        });

        if options.is_empty() {
            None
        } else {
            Some(options)
        }
    }

    fn get_sequence(&self, id: usize) -> Option<&ProteinSequence> {
        self.sequences.get(id)
    }

    fn len(&self) -> usize {
        self.sequences.len()
    }
}

impl ProteinSequenceCollection {
    pub fn from_fasta(fasta: &str) -> ProteinSequenceCollection {
        let mut sequences = vec![];
        let mut num = 0;
        let mut current_sequence = ProteinSequenceBuilder::new(num);
        for line in fasta.lines() {
            if line.starts_with(">") {
                if !current_sequence.is_empty() {
                    sequences.push(current_sequence.build());
                }
                current_sequence = ProteinSequenceBuilder::new(num);
                num += 1;
                let description = line.trim_start_matches('>').trim();
                current_sequence = current_sequence.with_description(description);
            } else {
                current_sequence = current_sequence.append_sequence(line.trim());
            }
        }
        sequences.push(current_sequence.build());
        ProteinSequenceCollection { sequences }
    }

    pub fn from_fasta_file<P: AsRef<Path>>(
        file: P,
    ) -> Result<ProteinSequenceCollection, std::io::Error> {
        let fasta = std::fs::read_to_string(file);
        match fasta {
            Ok(fasta) => Ok(Self::from_fasta(&fasta)),
            Err(e) => Err(e),
        }
    }
}

type ProteinPeptideIdPair = (u32, u32);

pub struct ProteinPeptideGraph {
    pub edges: Vec<ProteinPeptideIdPair>,
}

// Tests ...
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_parsing() {
        let dummy_fasta_string = r#">mysupercoolprotein
PEPTIDEPINK
PEPTIDEPINKPEPTIDEPINK
PEPTIDEPINK

> mysupercoolprotein2
PEPTIDEPLNK
PEPTIDEPLNK

"#;
        let fasta = ProteinSequenceCollection::from_fasta(dummy_fasta_string);
        println!("{:?}", fasta);
        assert_eq!(fasta.sequences.len(), 2);
        assert_eq!(
            fasta.sequences[0].sequence,
            "PEPTIDEPINKPEPTIDEPINKPEPTIDEPINKPEPTIDEPINK"
        );
        assert_eq!(fasta.sequences[1].sequence, "PEPTIDEPLNKPEPTIDEPLNK");
        assert_eq!(fasta.sequences[0].description, "mysupercoolprotein");
        assert_eq!(fasta.sequences[1].description, "mysupercoolprotein2");
    }
}
