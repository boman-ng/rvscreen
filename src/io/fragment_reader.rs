use crate::error::{Result, RvScreenError};
use noodles::sam;
use std::path::Path;

/// A paired fragment presented to the pipeline.
///
/// All input formats normalize into this same structure so downstream screening can stay
/// format-agnostic and pull-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FragmentRecord {
    pub fragment_key: String,
    pub r1_seq: Vec<u8>,
    pub r1_qual: Vec<u8>,
    pub r2_seq: Vec<u8>,
    pub r2_qual: Vec<u8>,
}

impl FragmentRecord {
    pub fn new(
        fragment_key: String,
        r1_seq: Vec<u8>,
        r1_qual: Vec<u8>,
        r2_seq: Vec<u8>,
        r2_qual: Vec<u8>,
    ) -> Self {
        Self {
            fragment_key,
            r1_seq,
            r1_qual,
            r2_seq,
            r2_qual,
        }
    }
}

/// A boxed iterator over fragment records with error handling.
///
/// This is the shared streaming contract consumed by the pipeline. Readers are opened once,
/// iterated once, and the caller controls pace through `next()`.
pub type FragmentStream = Box<dyn Iterator<Item = Result<FragmentRecord>>>;

/// Factory trait for opening unified fragment readers.
pub trait FragmentReaderFactory {
    fn open_reader(&self) -> Result<FragmentStream>;

    fn kind_label(&self) -> &'static str;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum MateSlot {
    Read1,
    Read2,
}

#[derive(Debug, Clone)]
pub(crate) struct FragmentMateRecord {
    pub record_number: u64,
    pub qname: String,
    pub fragment_key: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub slot: Option<MateSlot>,
}

impl FragmentMateRecord {
    pub(crate) fn new(
        record_number: u64,
        qname: String,
        fragment_key: String,
        seq: Vec<u8>,
        qual: Vec<u8>,
        slot: Option<MateSlot>,
    ) -> Self {
        Self {
            record_number,
            qname,
            fragment_key,
            seq,
            qual,
            slot,
        }
    }
}

pub(crate) fn mate_slot_from_flags(flags: sam::alignment::record::Flags) -> Option<MateSlot> {
    match (flags.is_first_segment(), flags.is_last_segment()) {
        (true, false) => Some(MateSlot::Read1),
        (false, true) => Some(MateSlot::Read2),
        _ => None,
    }
}

pub(crate) fn build_name_sorted_fragment_record(
    path: &Path,
    format_label: &str,
    first: &FragmentMateRecord,
    second: &FragmentMateRecord,
) -> Result<FragmentRecord> {
    if first.fragment_key != second.fragment_key {
        return Err(RvScreenError::parse(
            path,
            first.record_number,
            format!(
                "adjacent {format_label} records are not a name-sorted pair: record {} `{}` -> `{}`, record {} `{}` -> `{}`",
                first.record_number,
                first.qname,
                first.fragment_key,
                second.record_number,
                second.qname,
                second.fragment_key,
            ),
        ));
    }

    let (r1, r2) = ordered_mates(first, second);

    Ok(FragmentRecord::new(
        first.fragment_key.clone(),
        r1.seq.clone(),
        r1.qual.clone(),
        r2.seq.clone(),
        r2.qual.clone(),
    ))
}

pub(crate) fn orphaned_name_sorted_record_error(
    path: &Path,
    format_label: &str,
    orphan: &FragmentMateRecord,
    guidance: &str,
) -> RvScreenError {
    RvScreenError::parse(
        path,
        orphan.record_number,
        format!(
            "orphaned {format_label} record {} `{}`: adjacent mate is missing; {guidance}",
            orphan.record_number, orphan.qname,
        ),
    )
}

fn ordered_mates<'a>(
    first: &'a FragmentMateRecord,
    second: &'a FragmentMateRecord,
) -> (&'a FragmentMateRecord, &'a FragmentMateRecord) {
    match (first.slot, second.slot) {
        (Some(MateSlot::Read1), Some(MateSlot::Read2)) => (first, second),
        (Some(MateSlot::Read2), Some(MateSlot::Read1)) => (second, first),
        _ => (first, second),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flagged_name_sorted_mates_are_reordered_into_r1_r2_slots() {
        let first = FragmentMateRecord::new(
            1,
            "read/2".to_string(),
            "read".to_string(),
            b"TGCA".to_vec(),
            vec![31, 31, 31, 31],
            Some(MateSlot::Read2),
        );
        let second = FragmentMateRecord::new(
            2,
            "read/1".to_string(),
            "read".to_string(),
            b"ACGT".to_vec(),
            vec![30, 30, 30, 30],
            Some(MateSlot::Read1),
        );

        let fragment = build_name_sorted_fragment_record(
            Path::new("sample.bam"),
            "BAM",
            &first,
            &second,
        )
        .expect("flagged mates should build one fragment");

        assert_eq!(fragment.fragment_key, "read");
        assert_eq!(fragment.r1_seq, b"ACGT");
        assert_eq!(fragment.r2_seq, b"TGCA");
    }

    #[test]
    fn unflagged_name_sorted_mates_keep_adjacency_order() {
        let first = FragmentMateRecord::new(
            1,
            "read".to_string(),
            "read".to_string(),
            b"AAAA".to_vec(),
            vec![10, 10, 10, 10],
            None,
        );
        let second = FragmentMateRecord::new(
            2,
            "read".to_string(),
            "read".to_string(),
            b"CCCC".to_vec(),
            vec![20, 20, 20, 20],
            None,
        );

        let fragment = build_name_sorted_fragment_record(
            Path::new("sample.cram"),
            "CRAM",
            &first,
            &second,
        )
        .expect("unflagged mates should still build one fragment");

        assert_eq!(fragment.r1_seq, b"AAAA");
        assert_eq!(fragment.r2_seq, b"CCCC");
    }
}
