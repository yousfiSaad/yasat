//! Low-level FFI bindings to the YASAT C API
//!
//! This module contains raw, unsafe bindings to the C API.
//! Users should prefer the safe wrapper in the parent module.

use libc::{c_char, c_int};

/// Opaque handle to a YASAT solver instance
#[repr(C)]
pub struct yasat_solver {
    _private: [u8; 0],
}

/// Result codes returned by yasat_solve()
#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum yasat_result {
    /// Formula is satisfiable
    YASAT_RESULT_SAT = 10,
    /// Formula is unsatisfiable
    YASAT_RESULT_UNSAT = 20,
    /// Result unknown (not yet solved)
    YASAT_RESULT_UNKNOWN = 0,
    /// Error occurred during solving
    YASAT_RESULT_ERROR = -1,
}

/// Error codes for detailed error information
#[repr(C)]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
#[allow(non_camel_case_types)]
pub enum yasat_error {
    /// No error
    YASAT_OK = 0,
    /// Memory allocation failed
    YASAT_ERROR_MEMORY = -1,
    /// Invalid input parameters
    YASAT_ERROR_INVALID_INPUT = -2,
    /// CNF parsing error
    YASAT_ERROR_PARSE = -3,
    /// File I/O error
    YASAT_ERROR_FILE = -4,
    /// Invalid solver state
    YASAT_ERROR_STATE = -5,
}

#[link(name = "yasat")]
extern "C" {
    /// Create a new solver instance
    pub fn yasat_solver_create() -> *mut yasat_solver;

    /// Destroy a solver instance and free resources
    pub fn yasat_solver_destroy(solver: *mut yasat_solver);

    /// Add a clause to the solver
    ///
    /// # Arguments
    /// * `solver` - Solver instance
    /// * `literals` - Array of literals
    /// * `num_literals` - Number of literals
    pub fn yasat_add_clause(
        solver: *mut yasat_solver,
        literals: *const c_int,
        num_literals: c_int,
    ) -> yasat_error;

    /// Parse and load a CNF file
    pub fn yasat_parse_cnf_file(
        solver: *mut yasat_solver,
        filename: *const c_char,
    ) -> yasat_error;

    /// Parse CNF from a string
    pub fn yasat_parse_cnf_string(
        solver: *mut yasat_solver,
        cnf_string: *const c_char,
    ) -> yasat_error;

    /// Solve the SAT problem
    pub fn yasat_solve(solver: *mut yasat_solver) -> yasat_result;

    /// Get the number of variables
    pub fn yasat_get_num_variables(solver: *const yasat_solver) -> c_int;

    /// Get the number of clauses
    pub fn yasat_get_num_clauses(solver: *const yasat_solver) -> c_int;

    /// Get the assignment for a specific variable
    ///
    /// Returns:
    /// * 1 if variable is assigned true
    /// * 0 if variable is assigned false
    /// * -1 if variable is unassigned or error
    pub fn yasat_get_assignment(solver: *const yasat_solver, variable: c_int) -> c_int;

    /// Get all variable assignments
    ///
    /// # Arguments
    /// * `solver` - Solver instance
    /// * `assignments` - Output array (must be pre-allocated)
    /// * `num_vars` - Number of variables
    ///
    /// Returns YASAT_OK on success, error code on failure
    pub fn yasat_get_all_assignments(
        solver: *const yasat_solver,
        assignments: *mut c_int,
        num_vars: c_int,
    ) -> yasat_error;

    /// Reset the solver to initial state
    pub fn yasat_reset(solver: *mut yasat_solver);
}
