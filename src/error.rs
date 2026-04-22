use std::path::PathBuf;
use thiserror::Error;

/// Project-level typed error surface.
///
/// Every variant carries enough context for a human-readable message without
/// exposing raw data values (reads, quality strings, etc.).
#[derive(Debug, Error)]
pub enum RvScreenError {
    /// A filesystem or I/O operation failed.
    #[error("I/O error on `{path}`: {source}")]
    IoError {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// A file or record could not be parsed (FASTQ, BAM, TOML, JSON, …).
    #[error("parse error in `{path}` at line {line}: {reason}")]
    ParseError {
        path: PathBuf,
        line: u64,
        reason: String,
    },

    /// Alignment step failed or produced an unacceptable result.
    #[error("alignment error for contig `{contig}`: {reason}")]
    AlignmentError { contig: String, reason: String },

    /// A value violated a domain constraint (e.g. threshold out of range).
    #[error("validation error in `{field}`: {reason}")]
    ValidationError { field: String, reason: String },

    /// A reference-bundle operation failed (build, load, integrity check).
    #[error("bundle error (version `{version}`): {reason}")]
    BundleError { version: String, reason: String },

    /// Configuration file could not be loaded or is invalid.
    #[error("config error in `{path}`: {reason}")]
    ConfigError { path: PathBuf, reason: String },
}

impl RvScreenError {
    /// Convenience constructor for [`RvScreenError::IoError`].
    pub fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        Self::IoError {
            path: path.into(),
            source,
        }
    }

    /// Convenience constructor for [`RvScreenError::ParseError`].
    pub fn parse(path: impl Into<PathBuf>, line: u64, reason: impl Into<String>) -> Self {
        Self::ParseError {
            path: path.into(),
            line,
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RvScreenError::AlignmentError`].
    pub fn alignment(contig: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::AlignmentError {
            contig: contig.into(),
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RvScreenError::ValidationError`].
    pub fn validation(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::ValidationError {
            field: field.into(),
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RvScreenError::BundleError`].
    pub fn bundle(version: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::BundleError {
            version: version.into(),
            reason: reason.into(),
        }
    }

    /// Convenience constructor for [`RvScreenError::ConfigError`].
    pub fn config(path: impl Into<PathBuf>, reason: impl Into<String>) -> Self {
        Self::ConfigError {
            path: path.into(),
            reason: reason.into(),
        }
    }
}

/// Project-level `Result` alias for library code.
pub type Result<T> = std::result::Result<T, RvScreenError>;

#[cfg(test)]
mod tests {
    use super::*;

    fn display(e: &RvScreenError) -> String {
        e.to_string()
    }

    #[test]
    fn io_error_display_contains_path() {
        let err = RvScreenError::io(
            "/data/sample.fastq.gz",
            std::io::Error::new(std::io::ErrorKind::NotFound, "file not found"),
        );
        let msg = display(&err);
        assert!(
            msg.contains("/data/sample.fastq.gz"),
            "expected path in message, got: {msg}"
        );
        assert!(
            msg.contains("file not found"),
            "expected source reason in message, got: {msg}"
        );
    }

    #[test]
    fn parse_error_display_contains_line() {
        let err = RvScreenError::parse("/data/ref.fa", 42, "unexpected EOF");
        let msg = display(&err);
        assert!(
            msg.contains("42"),
            "expected line number in message, got: {msg}"
        );
        assert!(
            msg.contains("unexpected EOF"),
            "expected reason in message, got: {msg}"
        );
    }

    #[test]
    fn alignment_error_display_contains_contig() {
        let err = RvScreenError::alignment("chr1", "no seed hits");
        let msg = display(&err);
        assert!(
            msg.contains("chr1"),
            "expected contig name in message, got: {msg}"
        );
    }

    #[test]
    fn validation_error_display_contains_field() {
        let err = RvScreenError::validation("min_mapq", "value -1 is below the allowed minimum 0");
        let msg = display(&err);
        assert!(
            msg.contains("min_mapq"),
            "expected field name in message, got: {msg}"
        );
    }

    #[test]
    fn bundle_error_display_contains_version() {
        let err = RvScreenError::bundle("v2.3.1", "SHA-256 mismatch for layer hpv18");
        let msg = display(&err);
        assert!(
            msg.contains("v2.3.1"),
            "expected version in message, got: {msg}"
        );
    }

    #[test]
    fn config_error_display_contains_path() {
        let err =
            RvScreenError::config("/etc/rvscreen.toml", "missing required key `reference_dir`");
        let msg = display(&err);
        assert!(
            msg.contains("/etc/rvscreen.toml"),
            "expected config path in message, got: {msg}"
        );
        assert!(
            msg.contains("reference_dir"),
            "expected key name in message, got: {msg}"
        );
    }

    #[test]
    fn test_error_display() {
        let cases: &[(&str, RvScreenError)] = &[
            (
                "/data/sample.fastq.gz",
                RvScreenError::io(
                    "/data/sample.fastq.gz",
                    std::io::Error::new(std::io::ErrorKind::NotFound, "file not found"),
                ),
            ),
            (
                "unexpected EOF",
                RvScreenError::parse("/data/ref.fa", 42, "unexpected EOF"),
            ),
            ("chr1", RvScreenError::alignment("chr1", "no seed hits")),
            (
                "min_mapq",
                RvScreenError::validation("min_mapq", "value -1 is below the allowed minimum 0"),
            ),
            (
                "v2.3.1",
                RvScreenError::bundle("v2.3.1", "SHA-256 mismatch for layer hpv18"),
            ),
            (
                "/etc/rvscreen.toml",
                RvScreenError::config("/etc/rvscreen.toml", "missing required key `reference_dir`"),
            ),
        ];
        for (needle, err) in cases {
            let msg = err.to_string();
            assert!(
                msg.contains(needle),
                "display for {:?} missing `{needle}`, got: {msg}",
                std::mem::discriminant(err),
            );
        }
    }

    #[test]
    fn result_alias_compiles() {
        let ok: Result<u32> = Ok(1);
        assert_eq!(ok.unwrap(), 1);

        let err: Result<u32> = Err(RvScreenError::validation("x", "too large"));
        assert!(err.is_err());
    }
}
