// Package yasat provides Go bindings for the YASAT SAT solver.
//
// YASAT is a CDCL-based SAT solver for boolean satisfiability problems.
// This package wraps the C API to provide an idiomatic Go interface.
//
// Basic usage:
//
//	solver, err := yasat.NewSolver()
//	if err != nil {
//	    log.Fatal(err)
//	}
//	defer solver.Close()
//
//	// Add clauses
//	solver.AddClause(1, -2, 3)  // (x1 ∨ ¬x2 ∨ x3)
//	solver.AddClause(-1, 2)     // (¬x1 ∨ x2)
//
//	// Solve
//	result, err := solver.Solve()
//	if err != nil {
//	    log.Fatal(err)
//	}
//
//	if result == yasat.SAT {
//	    fmt.Println("SAT")
//	    // Get assignments
//	    val, _ := solver.GetAssignment(1)
//	    fmt.Printf("x1 = %v\n", val)
//	} else {
//	    fmt.Println("UNSAT")
//	}
package yasat

/*
#cgo LDFLAGS: -L../../../build -lyasat
#include "../../../src/c_api/yasat.h"
#include <stdlib.h>
*/
import "C"
import (
	"errors"
	"fmt"
	"unsafe"
)

// Result represents the result of solving a CNF formula.
type Result int

const (
	// SAT indicates the formula is satisfiable
	SAT Result = 10
	// UNSAT indicates the formula is unsatisfiable
	UNSAT Result = 20
	// Unknown indicates the result is unknown (solver hasn't been run)
	Unknown Result = 0
	// Error indicates an error occurred during solving
	Error Result = -1
)

// String returns a string representation of the Result.
func (r Result) String() string {
	switch r {
	case SAT:
		return "SAT"
	case UNSAT:
		return "UNSAT"
	case Unknown:
		return "UNKNOWN"
	case Error:
		return "ERROR"
	default:
		return fmt.Sprintf("Result(%d)", r)
	}
}

// Solver represents a YASAT SAT solver instance.
type Solver struct {
	handle *C.yasat_solver
}

// NewSolver creates a new YASAT solver instance.
func NewSolver() (*Solver, error) {
	handle := C.yasat_solver_create()
	if handle == nil {
		return nil, errors.New("failed to create solver")
	}
	return &Solver{handle: handle}, nil
}

// NewSolverFromFile creates a new solver and loads a CNF formula from a file.
func NewSolverFromFile(filename string) (*Solver, error) {
	solver, err := NewSolver()
	if err != nil {
		return nil, err
	}

	if err := solver.ParseCNFFile(filename); err != nil {
		solver.Close()
		return nil, err
	}

	return solver, nil
}

// Close destroys the solver and frees all associated resources.
// After calling Close, the Solver instance should not be used.
func (s *Solver) Close() error {
	if s.handle != nil {
		C.yasat_solver_destroy(s.handle)
		s.handle = nil
	}
	return nil
}

// AddClause adds a single clause to the solver.
// Literals are represented as integers: positive for x_i, negative for ¬x_i.
//
// Example: To add (x1 ∨ ¬x2 ∨ x3), call AddClause(1, -2, 3)
func (s *Solver) AddClause(literals ...int) error {
	if s.handle == nil {
		return errors.New("solver is closed")
	}

	if len(literals) == 0 {
		return nil // Empty clause is valid (though trivially unsatisfiable)
	}

	// Convert Go slice to C array
	cLiterals := make([]C.int, len(literals))
	for i, lit := range literals {
		cLiterals[i] = C.int(lit)
	}

	result := C.yasat_add_clause(s.handle, &cLiterals[0], C.int(len(literals)))
	if result != C.YASAT_OK {
		return s.getError()
	}

	return nil
}

// ParseCNFFile loads a CNF formula from a DIMACS format file.
func (s *Solver) ParseCNFFile(filename string) error {
	if s.handle == nil {
		return errors.New("solver is closed")
	}

	cFilename := C.CString(filename)
	defer C.free(unsafe.Pointer(cFilename))

	result := C.yasat_parse_cnf_file(s.handle, cFilename)
	if result != C.YASAT_OK {
		return s.getError()
	}

	return nil
}

// ParseCNFString loads a CNF formula from a DIMACS format string.
func (s *Solver) ParseCNFString(cnf string) error {
	if s.handle == nil {
		return errors.New("solver is closed")
	}

	cCNF := C.CString(cnf)
	defer C.free(unsafe.Pointer(cCNF))

	result := C.yasat_parse_cnf_string(s.handle, cCNF)
	if result != C.YASAT_OK {
		return s.getError()
	}

	return nil
}

// Solve attempts to find a satisfying assignment for the CNF formula.
// Returns SAT if satisfiable, UNSAT if unsatisfiable, or Error on error.
func (s *Solver) Solve() (Result, error) {
	if s.handle == nil {
		return Error, errors.New("solver is closed")
	}

	result := C.yasat_solve(s.handle)

	switch result {
	case C.YASAT_RESULT_SAT:
		return SAT, nil
	case C.YASAT_RESULT_UNSAT:
		return UNSAT, nil
	case C.YASAT_RESULT_ERROR:
		return Error, s.getError()
	default:
		return Unknown, nil
	}
}

// GetAssignment returns the truth value assigned to a variable in the solution.
// Variables are 1-indexed. Returns true if the variable is true, false if false.
// Only valid after Solve() returns SAT.
func (s *Solver) GetAssignment(variable int) (bool, error) {
	if s.handle == nil {
		return false, errors.New("solver is closed")
	}

	result := C.yasat_get_assignment(s.handle, C.int(variable))
	if result < 0 {
		return false, errors.New("variable not assigned or invalid")
	}

	return result == 1, nil
}

// GetAllAssignments returns the truth values for all variables.
// The returned slice is 0-indexed: result[0] corresponds to variable 1.
// Only valid after Solve() returns SAT.
func (s *Solver) GetAllAssignments() ([]bool, error) {
	if s.handle == nil {
		return nil, errors.New("solver is closed")
	}

	numVars := int(C.yasat_get_num_variables(s.handle))
	if numVars < 0 {
		return nil, errors.New("failed to get number of variables")
	}

	if numVars == 0 {
		return []bool{}, nil
	}

	// Allocate C array
	cAssignments := make([]C.int, numVars)
	result := C.yasat_get_all_assignments(s.handle, &cAssignments[0], C.int(numVars))
	if result != C.YASAT_OK {
		return nil, s.getError()
	}

	// Convert to Go bool slice
	assignments := make([]bool, numVars)
	for i := 0; i < numVars; i++ {
		assignments[i] = cAssignments[i] == 1
	}

	return assignments, nil
}

// NumVariables returns the number of variables in the formula.
func (s *Solver) NumVariables() (int, error) {
	if s.handle == nil {
		return 0, errors.New("solver is closed")
	}

	result := C.yasat_get_num_variables(s.handle)
	if result < 0 {
		return 0, errors.New("failed to get number of variables")
	}

	return int(result), nil
}

// NumClauses returns the number of clauses in the formula.
func (s *Solver) NumClauses() (int, error) {
	if s.handle == nil {
		return 0, errors.New("solver is closed")
	}

	result := C.yasat_get_num_clauses(s.handle)
	if result < 0 {
		return 0, errors.New("failed to get number of clauses")
	}

	return int(result), nil
}

// Version returns the YASAT library version string.
func Version() string {
	return C.GoString(C.yasat_get_version())
}

// getError retrieves the last error message from the C library.
func (s *Solver) getError() error {
	if s.handle == nil {
		return errors.New("solver is closed")
	}

	cMsg := C.yasat_get_error_message(s.handle)
	if cMsg == nil {
		return errors.New("unknown error")
	}

	return errors.New(C.GoString(cMsg))
}
