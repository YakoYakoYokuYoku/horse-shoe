//! A Rust Language infinitesimal calculus crate
//!
//! This library has got only numerical differentiation and
//! integration, more sources coming soon.

mod functions;

mod numder;
mod numint;

pub use numder::Diff;
pub use numint::Intg;
