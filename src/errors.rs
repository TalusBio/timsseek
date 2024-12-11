use serde_json;
use timsquery::TimsqueryError;
use timsrust::TimsRustError;

#[derive(Debug)]
pub enum TimsSeekError {
    TimsRust(TimsRustError),
    Timsquery(TimsqueryError),
    Io(std::io::Error),
    ParseError { msg: String },
}

impl std::fmt::Display for TimsSeekError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

type Result<T> = std::result::Result<T, TimsSeekError>;

impl From<std::io::Error> for TimsSeekError {
    fn from(x: std::io::Error) -> Self {
        Self::Io(x)
    }
}

impl From<TimsRustError> for TimsSeekError {
    fn from(x: TimsRustError) -> Self {
        Self::TimsRust(x)
    }
}

impl From<TimsqueryError> for TimsSeekError {
    fn from(x: TimsqueryError) -> Self {
        Self::Timsquery(x)
    }
}

impl From<std::num::ParseIntError> for TimsSeekError {
    fn from(x: std::num::ParseIntError) -> Self {
        Self::ParseError { msg: x.to_string() }
    }
}

impl Into<TimsSeekError> for serde_json::Error {
    fn into(self) -> TimsSeekError {
        TimsSeekError::ParseError {
            msg: self.to_string(),
        }
    }
}
