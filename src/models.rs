use serde::Serialize;

#[derive(Debug, Clone, Copy, Serialize, PartialEq, Eq, std::hash::Hash, PartialOrd, Ord)]
pub enum DecoyMarking {
    Target,
    Decoy,
}
impl DecoyMarking {
    pub fn as_str(&self) -> &'static str {
        match self {
            DecoyMarking::Target => "Target",
            DecoyMarking::Decoy => "Decoy",
        }
    }
}
