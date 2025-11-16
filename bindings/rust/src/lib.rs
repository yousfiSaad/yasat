//! # YASAT - Yet Another SAT Solver
//!
//! Rust bindings for YASAT, a CDCL (Conflict-Driven Clause Learning) SAT solver.
//!
//! ## Features
//!
//! - CDCL algorithm implementation
//! - Simple, safe Rust API
//! - Support for loading CNF files
//! - Programmatic clause building
//! - Zero-cost abstractions over C API
//!
//! ## Example
//!
//! ```rust,no_run
//! use yasat::{Solver, SatResult};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut solver = Solver::new()?;
//!
//!     // Add clauses: (x1 ∨ x2) ∧ (¬x1 ∨ x3)
//!     solver.add_clause(&[1, 2])?;
//!     solver.add_clause(&[-1, 3])?;
//!
//!     match solver.solve()? {
//!         SatResult::Sat => {
//!             println!("SAT");
//!             let assignments = solver.get_assignments()?;
//!             println!("Solution: {:?}", assignments);
//!         }
//!         SatResult::Unsat => {
//!             println!("UNSAT");
//!         }
//!     }
//!
//!     Ok(())
//! }
//! ```

pub mod sys;

use std::ffi::{CStr, CString};
use std::path::Path;
use std::ptr;
use thiserror::Error;

/// Error type for YASAT operations
#[derive(Debug, Error)]
pub enum Error {
    #[error("Failed to create solver")]
    CreationFailed,

    #[error("Invalid argument: {0}")]
    InvalidArgument(String),

    #[error("Out of memory")]
    OutOfMemory,

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("I/O error: {0}")]
    IoError(String),

    #[error("Solver not initialized")]
    NotInitialized,

    #[error("Solver already solved")]
    AlreadySolved,

    #[error("Invalid variable: {0}")]
    InvalidVariable(i32),

    #[error("UTF-8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),

    #[error("Null error: {0}")]
    NulError(#[from] std::ffi::NulError),

    #[error("Unknown error")]
    Unknown,
}

/// Result of a SAT solving operation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SatResult {
    /// Formula is satisfiable
    Sat,
    /// Formula is unsatisfiable
    Unsat,
}

/// A SAT solver instance
///
/// The solver maintains state including:
/// - Added clauses
/// - Variable assignments (after solving)
/// - Learned clauses (CDCL algorithm)
///
/// ## Example
///
/// ```rust,no_run
/// # use yasat::Solver;
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let mut solver = Solver::new()?;
/// solver.add_clause(&[1, 2])?;
/// solver.add_clause(&[-1, 3])?;
/// let result = solver.solve()?;
/// # Ok(())
/// # }
/// ```
pub struct Solver {
    inner: ptr::NonNull<sys::yasat_solver>,
}

impl Solver {
    /// Create a new solver instance
    ///
    /// # Errors
    ///
    /// Returns `Error::CreationFailed` if the solver cannot be created
    /// (usually due to memory allocation failure).
    pub fn new() -> Result<Self, Error> {
        let ptr = unsafe { sys::yasat_solver_create() };

        ptr::NonNull::new(ptr)
            .map(|inner| Solver { inner })
            .ok_or(Error::CreationFailed)
    }

    /// Create a solver from a CNF file
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the CNF file in DIMACS format
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use yasat::Solver;
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut solver = Solver::from_file("problem.cnf")?;
    /// let result = solver.solve()?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Error> {
        let mut solver = Self::new()?;
        solver.load_cnf_file(path)?;
        Ok(solver)
    }

    /// Create a solver from a CNF string
    ///
    /// # Arguments
    ///
    /// * `cnf` - CNF formula in DIMACS format
    pub fn from_string(cnf: &str) -> Result<Self, Error> {
        let mut solver = Self::new()?;
        solver.load_cnf_string(cnf)?;
        Ok(solver)
    }

    /// Add a clause to the solver
    ///
    /// # Arguments
    ///
    /// * `literals` - Array of literals where positive numbers represent
    ///   variables and negative numbers represent negated variables
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use yasat::Solver;
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut solver = Solver::new()?;
    ///
    /// // Add clause: (x1 ∨ x2)
    /// solver.add_clause(&[1, 2])?;
    ///
    /// // Add clause: (¬x1 ∨ x3)
    /// solver.add_clause(&[-1, 3])?;
    /// # Ok(())
    /// # }
    /// ```
    pub fn add_clause(&mut self, literals: &[i32]) -> Result<(), Error> {
        if literals.is_empty() {
            return Err(Error::InvalidArgument("Empty clause".to_string()));
        }

        let result = unsafe {
            sys::yasat_add_clause(
                self.inner.as_ptr(),
                literals.as_ptr(),
                literals.len(),
            )
        };

        Self::check_error(result)
    }

    /// Load a CNF file in DIMACS format
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the CNF file
    pub fn load_cnf_file<P: AsRef<Path>>(&mut self, path: P) -> Result<(), Error> {
        let path_str = path.as_ref()
            .to_str()
            .ok_or_else(|| Error::InvalidArgument("Invalid path".to_string()))?;

        let c_path = CString::new(path_str)?;

        let result = unsafe {
            sys::yasat_parse_cnf_file(self.inner.as_ptr(), c_path.as_ptr())
        };

        Self::check_error(result)
    }

    /// Load CNF from a string in DIMACS format
    ///
    /// # Arguments
    ///
    /// * `cnf` - CNF formula string
    pub fn load_cnf_string(&mut self, cnf: &str) -> Result<(), Error> {
        let c_cnf = CString::new(cnf)?;

        let result = unsafe {
            sys::yasat_parse_cnf_string(self.inner.as_ptr(), c_cnf.as_ptr())
        };

        Self::check_error(result)
    }

    /// Solve the SAT problem
    ///
    /// Returns `SatResult::Sat` if satisfiable, `SatResult::Unsat` if unsatisfiable.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use yasat::{Solver, SatResult};
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut solver = Solver::new()?;
    /// solver.add_clause(&[1, 2])?;
    ///
    /// match solver.solve()? {
    ///     SatResult::Sat => println!("Satisfiable!"),
    ///     SatResult::Unsat => println!("Unsatisfiable!"),
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn solve(&mut self) -> Result<SatResult, Error> {
        let result = unsafe { sys::yasat_solve(self.inner.as_ptr()) };

        match result {
            sys::yasat_result::YASAT_RESULT_SAT => Ok(SatResult::Sat),
            sys::yasat_result::YASAT_RESULT_UNSAT => Ok(SatResult::Unsat),
            sys::yasat_result::YASAT_RESULT_ERROR => Err(Error::Unknown),
            _ => Err(Error::Unknown),
        }
    }

    /// Get the number of variables in the formula
    pub fn num_variables(&self) -> usize {
        unsafe { sys::yasat_get_num_variables(self.inner.as_ptr()) }
    }

    /// Get the number of clauses in the formula
    pub fn num_clauses(&self) -> usize {
        unsafe { sys::yasat_get_num_clauses(self.inner.as_ptr()) }
    }

    /// Get the assignment for a specific variable
    ///
    /// # Arguments
    ///
    /// * `variable` - Variable number (must be positive)
    ///
    /// # Returns
    ///
    /// * `Some(true)` if assigned true
    /// * `Some(false)` if assigned false
    /// * `None` if unassigned
    pub fn get_assignment(&self, variable: i32) -> Result<Option<bool>, Error> {
        if variable <= 0 {
            return Err(Error::InvalidVariable(variable));
        }

        let result = unsafe {
            sys::yasat_get_assignment(self.inner.as_ptr(), variable)
        };

        match result {
            1 => Ok(Some(true)),
            0 => Ok(Some(false)),
            -1 => Ok(None),
            _ => Err(Error::Unknown),
        }
    }

    /// Get all variable assignments
    ///
    /// Returns a vector where index i contains the assignment for variable i+1.
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// # use yasat::{Solver, SatResult};
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let mut solver = Solver::new()?;
    /// solver.add_clause(&[1, 2])?;
    ///
    /// if solver.solve()? == SatResult::Sat {
    ///     let assignments = solver.get_assignments()?;
    ///     for (i, &value) in assignments.iter().enumerate() {
    ///         println!("x{} = {}", i + 1, value);
    ///     }
    /// }
    /// # Ok(())
    /// # }
    /// ```
    pub fn get_assignments(&self) -> Result<Vec<bool>, Error> {
        let num_vars = self.num_variables();
        if num_vars == 0 {
            return Ok(Vec::new());
        }

        let mut assignments = vec![0i32; num_vars];

        let result = unsafe {
            sys::yasat_get_all_assignments(
                self.inner.as_ptr(),
                assignments.as_mut_ptr(),
                num_vars,
            )
        };

        if result < 0 {
            return Err(Error::Unknown);
        }

        Ok(assignments.into_iter().map(|v| v != 0).collect())
    }

    /// Reset the solver to its initial state
    ///
    /// This clears all clauses and assignments, allowing reuse of the solver.
    pub fn reset(&mut self) {
        unsafe { sys::yasat_reset(self.inner.as_ptr()) }
    }

    fn check_error(error: sys::yasat_error) -> Result<(), Error> {
        match error {
            sys::yasat_error::YASAT_OK => Ok(()),
            sys::yasat_error::YASAT_ERROR_INVALID_INPUT => {
                Err(Error::InvalidArgument("Invalid argument".to_string()))
            }
            sys::yasat_error::YASAT_ERROR_MEMORY => Err(Error::OutOfMemory),
            sys::yasat_error::YASAT_ERROR_PARSE => {
                Err(Error::ParseError("Parse error".to_string()))
            }
            sys::yasat_error::YASAT_ERROR_FILE => Err(Error::IoError("I/O error".to_string())),
            sys::yasat_error::YASAT_ERROR_STATE => Err(Error::NotInitialized),
        }
    }
}

// SAFETY: The C library is thread-safe (each solver instance is independent)
unsafe impl Send for Solver {}
unsafe impl Sync for Solver {}

impl Drop for Solver {
    fn drop(&mut self) {
        unsafe {
            sys::yasat_solver_destroy(self.inner.as_ptr());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_solver() {
        let solver = Solver::new();
        assert!(solver.is_ok());
    }

    #[test]
    fn test_simple_sat() {
        let mut solver = Solver::new().unwrap();

        // (x1 ∨ x2)
        solver.add_clause(&[1, 2]).unwrap();

        let result = solver.solve().unwrap();
        assert_eq!(result, SatResult::Sat);
    }

    #[test]
    fn test_simple_unsat() {
        let mut solver = Solver::new().unwrap();

        // x1 ∧ ¬x1 (contradiction)
        solver.add_clause(&[1]).unwrap();
        solver.add_clause(&[-1]).unwrap();

        let result = solver.solve().unwrap();
        assert_eq!(result, SatResult::Unsat);
    }

    #[test]
    fn test_get_assignments() {
        let mut solver = Solver::new().unwrap();

        // (x1 ∨ x2) ∧ (¬x1 ∨ x3)
        solver.add_clause(&[1, 2]).unwrap();
        solver.add_clause(&[-1, 3]).unwrap();

        let result = solver.solve().unwrap();
        assert_eq!(result, SatResult::Sat);

        let assignments = solver.get_assignments().unwrap();
        assert_eq!(assignments.len(), 3);
    }

    #[test]
    fn test_num_variables_clauses() {
        let mut solver = Solver::new().unwrap();

        solver.add_clause(&[1, 2, 3]).unwrap();
        solver.add_clause(&[-1, 4]).unwrap();

        assert_eq!(solver.num_variables(), 4);
        assert_eq!(solver.num_clauses(), 2);
    }
}
