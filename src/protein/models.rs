use std::ops::Range;
use std::sync::Arc;

#[derive(Debug)]
pub struct ProteinSequence {
    pub id: u32, // Self incremental identifier within the fasta file.
    pub description: String,
    pub sequence: Arc<str>,
}


#[derive(Debug)]
pub struct ProteinSequenceBuilder {
    pub id: u32,
    pub description: Option<String>,
    pub sequence: String,
}

impl ProteinSequenceBuilder {
    pub fn new(id: u32) -> Self {
        Self {
            id,
            description: None,
            sequence: String::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    pub fn with_description(mut self, description: &str) -> Self {
        self.description = Some(description.to_string());
        self
    }

    pub fn append_sequence(mut self, sequence: &str) -> Self {
        self.sequence.push_str(sequence);
        self
    }

    pub fn build(self) -> ProteinSequence {
        debug_assert!(self.description.is_some());
        ProteinSequence {
            id: self.id,
            description: self.description.unwrap(),
            sequence: self.sequence.into(),
        }
    }
}
