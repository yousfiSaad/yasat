package yasat

import (
	"os"
	"path/filepath"
	"testing"
)

// TestVersion tests the version string
func TestVersion(t *testing.T) {
	version := Version()
	if version == "" {
		t.Error("Version() returned empty string")
	}
	t.Logf("YASAT version: %s", version)
}

// TestSolverCreate tests solver creation and destruction
func TestSolverCreate(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	if solver == nil {
		t.Fatal("NewSolver() returned nil solver")
	}

	err = solver.Close()
	if err != nil {
		t.Errorf("Close() failed: %v", err)
	}
}

// TestMultipleClose tests that closing a solver multiple times doesn't crash
func TestMultipleClose(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}

	solver.Close()
	solver.Close() // Should not crash
}

// TestSimpleSAT tests a simple satisfiable formula
func TestSimpleSAT(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	// (x1 OR x2) AND (NOT x1 OR x2)
	// Should be SAT with x2=true
	err = solver.AddClause(1, 2)
	if err != nil {
		t.Fatalf("AddClause(1, 2) failed: %v", err)
	}

	err = solver.AddClause(-1, 2)
	if err != nil {
		t.Fatalf("AddClause(-1, 2) failed: %v", err)
	}

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	if result != SAT {
		t.Errorf("Expected SAT, got %v", result)
	}

	// Check that x2 is true
	val, err := solver.GetAssignment(2)
	if err != nil {
		t.Errorf("GetAssignment(2) failed: %v", err)
	}
	if !val {
		t.Error("Expected x2=true, got false")
	}
}

// TestComplexSAT tests a more complex satisfiable formula
func TestComplexSAT(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	// (x1 OR x2) AND (NOT x1 OR x3) AND (NOT x2 OR NOT x3)
	solver.AddClause(1, 2)
	solver.AddClause(-1, 3)
	solver.AddClause(-2, -3)

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	if result != SAT {
		t.Errorf("Expected SAT, got %v", result)
	}

	// Verify we can get assignments individually
	x1, err := solver.GetAssignment(1)
	if err != nil {
		t.Errorf("GetAssignment(1) failed: %v", err)
	}

	x2, err := solver.GetAssignment(2)
	if err != nil {
		t.Errorf("GetAssignment(2) failed: %v", err)
	}

	x3, err := solver.GetAssignment(3)
	if err != nil {
		t.Errorf("GetAssignment(3) failed: %v", err)
	}

	// Verify the solution satisfies the clauses
	if !(x1 || x2) {
		t.Error("Clause (x1 OR x2) not satisfied")
	}
	if !(!x1 || x3) {
		t.Error("Clause (NOT x1 OR x3) not satisfied")
	}
	if !(!x2 || !x3) {
		t.Error("Clause (NOT x2 OR NOT x3) not satisfied")
	}

	// Also test GetAllAssignments
	assignments, err := solver.GetAllAssignments()
	if err != nil {
		t.Fatalf("GetAllAssignments() failed: %v", err)
	}

	numVars, _ := solver.NumVariables()
	if len(assignments) != numVars {
		t.Errorf("Expected %d assignments, got %d", numVars, len(assignments))
	}
}

// TestEmptyClause tests adding an empty clause
func TestEmptyClause(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	// Adding empty clause should not error
	err = solver.AddClause()
	if err != nil {
		t.Errorf("AddClause() with no arguments failed: %v", err)
	}
}

// TestUnitClauses tests unit clauses
func TestUnitClauses(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	// (x1) AND (x2)
	solver.AddClause(1)
	solver.AddClause(2)

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	if result != SAT {
		t.Errorf("Expected SAT, got %v", result)
	}

	x1, _ := solver.GetAssignment(1)
	x2, _ := solver.GetAssignment(2)

	if !x1 {
		t.Error("Expected x1=true")
	}
	if !x2 {
		t.Error("Expected x2=true")
	}
}

// TestParseCNFString tests parsing CNF from a string
func TestParseCNFString(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	cnf := `p cnf 2 2
1 2 0
-1 2 0`

	err = solver.ParseCNFString(cnf)
	if err != nil {
		t.Fatalf("ParseCNFString() failed: %v", err)
	}

	numVars, err := solver.NumVariables()
	if err != nil {
		t.Fatalf("NumVariables() failed: %v", err)
	}
	if numVars != 2 {
		t.Errorf("Expected 2 variables, got %d", numVars)
	}

	numClauses, err := solver.NumClauses()
	if err != nil {
		t.Fatalf("NumClauses() failed: %v", err)
	}
	if numClauses != 2 {
		t.Errorf("Expected 2 clauses, got %d", numClauses)
	}

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	if result != SAT {
		t.Errorf("Expected SAT, got %v", result)
	}
}

// TestParseCNFFile tests parsing CNF from a file
func TestParseCNFFile(t *testing.T) {
	// Create a temporary CNF file
	tmpfile, err := os.CreateTemp("", "test_*.cnf")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	cnfContent := `c Simple test CNF
p cnf 3 3
1 2 0
-1 3 0
-2 -3 0
`
	if _, err := tmpfile.Write([]byte(cnfContent)); err != nil {
		t.Fatal(err)
	}
	if err := tmpfile.Close(); err != nil {
		t.Fatal(err)
	}

	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	err = solver.ParseCNFFile(tmpfile.Name())
	if err != nil {
		t.Fatalf("ParseCNFFile() failed: %v", err)
	}

	numVars, _ := solver.NumVariables()
	if numVars != 3 {
		t.Errorf("Expected 3 variables, got %d", numVars)
	}

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	if result != SAT {
		t.Errorf("Expected SAT, got %v", result)
	}
}

// TestNewSolverFromFile tests the convenience constructor
func TestNewSolverFromFile(t *testing.T) {
	// Create a temporary CNF file
	tmpfile, err := os.CreateTemp("", "test_*.cnf")
	if err != nil {
		t.Fatal(err)
	}
	defer os.Remove(tmpfile.Name())

	cnfContent := `p cnf 2 1
1 2 0
`
	if _, err := tmpfile.Write([]byte(cnfContent)); err != nil {
		t.Fatal(err)
	}
	if err := tmpfile.Close(); err != nil {
		t.Fatal(err)
	}

	solver, err := NewSolverFromFile(tmpfile.Name())
	if err != nil {
		t.Fatalf("NewSolverFromFile() failed: %v", err)
	}
	defer solver.Close()

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	if result != SAT {
		t.Errorf("Expected SAT, got %v", result)
	}
}

// TestInvalidFile tests parsing a non-existent file
func TestInvalidFile(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	err = solver.ParseCNFFile("/nonexistent/file.cnf")
	if err == nil {
		t.Error("Expected error when parsing non-existent file")
	}
}

// TestGetAssignmentBeforeSolve tests getting assignment before solving
func TestGetAssignmentBeforeSolve(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	solver.AddClause(1, 2)

	// Should return error or false when called before solving
	_, err = solver.GetAssignment(1)
	if err == nil {
		t.Log("GetAssignment() before Solve() returned no error (acceptable)")
	}
}

// TestNumVariablesEmpty tests getting number of variables from empty solver
func TestNumVariablesEmpty(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	numVars, err := solver.NumVariables()
	if err != nil {
		t.Fatalf("NumVariables() failed: %v", err)
	}

	if numVars != 0 {
		t.Errorf("Expected 0 variables in empty solver, got %d", numVars)
	}
}

// TestNumClausesAfterAdding tests clause counting
func TestNumClausesAfterAdding(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	solver.AddClause(1, 2)
	solver.AddClause(-1, 3)
	solver.AddClause(2, -3)

	numClauses, err := solver.NumClauses()
	if err != nil {
		t.Fatalf("NumClauses() failed: %v", err)
	}

	if numClauses != 3 {
		t.Errorf("Expected 3 clauses, got %d", numClauses)
	}
}

// TestLargeVariableNumbers tests handling large variable numbers
func TestLargeVariableNumbers(t *testing.T) {
	solver, err := NewSolver()
	if err != nil {
		t.Fatalf("NewSolver() failed: %v", err)
	}
	defer solver.Close()

	// Use variable 100
	solver.AddClause(100, -100)

	numVars, err := solver.NumVariables()
	if err != nil {
		t.Fatalf("NumVariables() failed: %v", err)
	}

	if numVars != 100 {
		t.Errorf("Expected 100 variables, got %d", numVars)
	}
}

// TestResultString tests the Result.String() method
func TestResultString(t *testing.T) {
	tests := []struct {
		result Result
		want   string
	}{
		{SAT, "SAT"},
		{UNSAT, "UNSAT"},
		{Unknown, "UNKNOWN"},
		{Error, "ERROR"},
	}

	for _, tt := range tests {
		got := tt.result.String()
		if got != tt.want {
			t.Errorf("Result(%d).String() = %q, want %q", tt.result, got, tt.want)
		}
	}
}

// TestIntegrationWithRealFile tests with an actual test CNF file if available
func TestIntegrationWithRealFile(t *testing.T) {
	// Try to find a test CNF file
	testFile := "../../../tests/cnf/simple_sat.cnf"

	// Check if file exists
	if _, err := os.Stat(testFile); os.IsNotExist(err) {
		t.Skip("Test CNF file not found, skipping integration test")
	}

	solver, err := NewSolverFromFile(testFile)
	if err != nil {
		t.Fatalf("NewSolverFromFile() failed: %v", err)
	}
	defer solver.Close()

	result, err := solver.Solve()
	if err != nil {
		t.Fatalf("Solve() failed: %v", err)
	}

	// simple_sat.cnf should be satisfiable
	if result != SAT {
		t.Errorf("Expected SAT for simple_sat.cnf, got %v", result)
	}

	t.Logf("Successfully solved %s: %v", filepath.Base(testFile), result)
}

// BenchmarkSolve benchmarks the solve operation
func BenchmarkSolve(b *testing.B) {
	cnf := `p cnf 10 20
1 2 0
-1 3 0
-2 4 0
-3 5 0
-4 6 0
-5 7 0
-6 8 0
-7 9 0
-8 10 0
-9 -10 0
1 -2 0
2 -3 0
3 -4 0
4 -5 0
5 -6 0
6 -7 0
7 -8 0
8 -9 0
9 -10 0
10 -1 0
`

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		solver, _ := NewSolver()
		solver.ParseCNFString(cnf)
		solver.Solve()
		solver.Close()
	}
}
