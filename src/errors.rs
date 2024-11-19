use timsquery::TimsqueryError;
use timsrust::TimsRustError;

#[derive(Debug)]
pub enum TimsSeekError {
    TimsRust(TimsRustError),
    Timsquery(TimsqueryError),
    Io(std::io::Error),
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
