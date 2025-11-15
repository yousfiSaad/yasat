/**
 * @file yasat.h
 * @brief C API for YASAT (Yet Another SAT Solver)
 *
 * This header provides a C interface to the YASAT CDCL SAT solver,
 * allowing usage from C, Go, Python, Rust, and other languages via FFI.
 *
 * Basic usage:
 * 1. Create solver: yasat_solver_create()
 * 2. Add clauses: yasat_add_clause() or yasat_parse_cnf_file()
 * 3. Solve: yasat_solve()
 * 4. Get results: yasat_get_assignment()
 * 5. Cleanup: yasat_solver_destroy()
 */

#ifndef YASAT_H
#define YASAT_H

#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Version information */
#define YASAT_VERSION_MAJOR 1
#define YASAT_VERSION_MINOR 0
#define YASAT_VERSION_PATCH 0

/**
 * @brief Opaque handle to a YASAT solver instance
 */
typedef struct yasat_solver yasat_solver;

/**
 * @brief Result codes returned by yasat_solve()
 */
typedef enum {
    YASAT_RESULT_SAT = 10,      /**< Formula is satisfiable */
    YASAT_RESULT_UNSAT = 20,    /**< Formula is unsatisfiable */
    YASAT_RESULT_UNKNOWN = 0,   /**< Result unknown (not yet solved) */
    YASAT_RESULT_ERROR = -1     /**< Error occurred during solving */
} yasat_result;

/**
 * @brief Error codes for detailed error information
 */
typedef enum {
    YASAT_OK = 0,               /**< No error */
    YASAT_ERROR_MEMORY = -1,    /**< Memory allocation failed */
    YASAT_ERROR_INVALID_INPUT = -2,  /**< Invalid input parameters */
    YASAT_ERROR_PARSE = -3,     /**< CNF parsing error */
    YASAT_ERROR_FILE = -4,      /**< File I/O error */
    YASAT_ERROR_STATE = -5      /**< Invalid solver state */
} yasat_error;

/* ========================================================================== */
/* Solver Lifecycle Management                                                */
/* ========================================================================== */

/**
 * @brief Create a new YASAT solver instance
 * @return Pointer to solver instance, or NULL on failure
 */
yasat_solver* yasat_solver_create(void);

/**
 * @brief Destroy a solver instance and free all associated memory
 * @param solver Solver instance to destroy (can be NULL)
 */
void yasat_solver_destroy(yasat_solver* solver);

/* ========================================================================== */
/* Input Methods                                                              */
/* ========================================================================== */

/**
 * @brief Add a single clause to the solver
 * @param solver Solver instance
 * @param literals Array of literals (positive for x_i, negative for ¬x_i)
 * @param num_literals Number of literals in the clause
 * @return YASAT_OK on success, error code on failure
 *
 * Example: To add (x1 ∨ ¬x2 ∨ x3), use literals = {1, -2, 3}, num_literals = 3
 * Variables are numbered starting from 1.
 */
yasat_error yasat_add_clause(yasat_solver* solver,
                             const int* literals,
                             int num_literals);

/**
 * @brief Parse and load a CNF formula from a file
 * @param solver Solver instance
 * @param filename Path to DIMACS CNF file
 * @return YASAT_OK on success, error code on failure
 */
yasat_error yasat_parse_cnf_file(yasat_solver* solver, const char* filename);

/**
 * @brief Parse and load a CNF formula from a string
 * @param solver Solver instance
 * @param cnf_string CNF formula in DIMACS format
 * @return YASAT_OK on success, error code on failure
 */
yasat_error yasat_parse_cnf_string(yasat_solver* solver, const char* cnf_string);

/* ========================================================================== */
/* Solving                                                                    */
/* ========================================================================== */

/**
 * @brief Solve the CNF formula
 * @param solver Solver instance
 * @return YASAT_RESULT_SAT, YASAT_RESULT_UNSAT, or YASAT_RESULT_ERROR
 */
yasat_result yasat_solve(yasat_solver* solver);

/* ========================================================================== */
/* Solution Retrieval                                                         */
/* ========================================================================== */

/**
 * @brief Get the assignment for a specific variable
 * @param solver Solver instance
 * @param var Variable number (1-indexed)
 * @return 1 if var=true, 0 if var=false, -1 if unassigned or error
 *
 * Only valid after yasat_solve() returns YASAT_RESULT_SAT.
 */
int yasat_get_assignment(yasat_solver* solver, int var);

/**
 * @brief Get assignments for all variables
 * @param solver Solver instance
 * @param[out] assignments Array to store assignments (must be pre-allocated)
 * @param num_vars Number of variables (size of assignments array)
 * @return YASAT_OK on success, error code on failure
 *
 * The assignments array will be filled with 1 (true) or 0 (false) for each variable.
 * Array is 0-indexed: assignments[0] corresponds to variable 1.
 * Only valid after yasat_solve() returns YASAT_RESULT_SAT.
 */
yasat_error yasat_get_all_assignments(yasat_solver* solver,
                                      int* assignments,
                                      int num_vars);

/**
 * @brief Get the number of variables in the formula
 * @param solver Solver instance
 * @return Number of variables, or -1 on error
 */
int yasat_get_num_variables(yasat_solver* solver);

/**
 * @brief Get the number of clauses in the formula
 * @param solver Solver instance
 * @return Number of clauses, or -1 on error
 */
int yasat_get_num_clauses(yasat_solver* solver);

/* ========================================================================== */
/* Error Handling                                                             */
/* ========================================================================== */

/**
 * @brief Get the last error message
 * @param solver Solver instance
 * @return Error message string, or NULL if no error
 *
 * The returned string is owned by the solver and should not be freed.
 * It remains valid until the next API call or solver destruction.
 */
const char* yasat_get_error_message(yasat_solver* solver);

/**
 * @brief Get library version string
 * @return Version string in format "major.minor.patch"
 */
const char* yasat_get_version(void);

#ifdef __cplusplus
}
#endif

#endif /* YASAT_H */
