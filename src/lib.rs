#![allow(
    clippy::bool_assert_comparison,
    clippy::duplicate_mod,
    clippy::items_after_test_module,
    clippy::redundant_closure,
    clippy::sliced_string_as_bytes,
    clippy::too_many_arguments,
    clippy::unnecessary_literal_unwrap
)]

#[cfg(test)]
extern crate self as rvscreen;

pub mod adjudicate;
pub mod aggregate;
pub mod align;
pub mod architecture;
pub mod audit;
pub mod calibration;
pub mod cli;
pub mod decision;
pub mod error;
pub mod io;
pub mod pipeline;
pub mod qc;
pub mod reference;
pub mod report;
pub mod sampling;
pub mod types;
