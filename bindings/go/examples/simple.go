// Simple example demonstrating YASAT Go bindings usage
package main

import (
	"fmt"
	"log"
	"os"

	"github.com/yousfiSaad/yasat/bindings/go/yasat"
)

func main() {
	fmt.Printf("YASAT Go Bindings Example (version %s)\n\n", yasat.Version())

	// Example 1: Manually adding clauses
	fmt.Println("=== Example 1: Manual clause addition ===")
	example1()

	// Example 2: Loading from file
	if len(os.Args) > 1 {
		fmt.Println("\n=== Example 2: Loading from file ===")
		example2(os.Args[1])
	}
}

func example1() {
	// Create solver
	solver, err := yasat.NewSolver()
	if err != nil {
		log.Fatalf("Failed to create solver: %v", err)
	}
	defer solver.Close()

	// Build a simple satisfiable formula:
	// (x1 ∨ x2) ∧ (¬x1 ∨ x3) ∧ (¬x2 ∨ ¬x3)
	fmt.Println("Adding clauses:")
	fmt.Println("  (x1 ∨ x2)")
	solver.AddClause(1, 2)

	fmt.Println("  (¬x1 ∨ x3)")
	solver.AddClause(-1, 3)

	fmt.Println("  (¬x2 ∨ ¬x3)")
	solver.AddClause(-2, -3)

	// Get formula stats
	numVars, _ := solver.NumVariables()
	numClauses, _ := solver.NumClauses()
	fmt.Printf("\nFormula: %d variables, %d clauses\n", numVars, numClauses)

	// Solve
	fmt.Println("\nSolving...")
	result, err := solver.Solve()
	if err != nil {
		log.Fatalf("Solving failed: %v", err)
	}

	fmt.Printf("Result: %s\n", result)

	if result == yasat.SAT {
		fmt.Println("\nSatisfying assignment:")
		assignments, err := solver.GetAllAssignments()
		if err != nil {
			log.Fatalf("Failed to get assignments: %v", err)
		}

		for i, val := range assignments {
			varNum := i + 1
			fmt.Printf("  x%d = %v\n", varNum, val)
		}

		// Verify the solution satisfies the clauses
		fmt.Println("\nVerification:")
		fmt.Printf("  (x1 ∨ x2) = (%v ∨ %v) = %v\n",
			assignments[0], assignments[1], assignments[0] || assignments[1])
		fmt.Printf("  (¬x1 ∨ x3) = (%v ∨ %v) = %v\n",
			!assignments[0], assignments[2], !assignments[0] || assignments[2])
		fmt.Printf("  (¬x2 ∨ ¬x3) = (%v ∨ %v) = %v\n",
			!assignments[1], !assignments[2], !assignments[1] || !assignments[2])
	}
}

func example2(filename string) {
	fmt.Printf("Loading CNF from file: %s\n", filename)

	// Create solver and load file
	solver, err := yasat.NewSolverFromFile(filename)
	if err != nil {
		log.Fatalf("Failed to load file: %v", err)
	}
	defer solver.Close()

	// Get formula stats
	numVars, _ := solver.NumVariables()
	numClauses, _ := solver.NumClauses()
	fmt.Printf("Formula: %d variables, %d clauses\n", numVars, numClauses)

	// Solve
	fmt.Println("\nSolving...")
	result, err := solver.Solve()
	if err != nil {
		log.Fatalf("Solving failed: %v", err)
	}

	fmt.Printf("Result: %s\n", result)

	if result == yasat.SAT {
		fmt.Println("\nSatisfying assignment found:")
		assignments, err := solver.GetAllAssignments()
		if err != nil {
			log.Fatalf("Failed to get assignments: %v", err)
		}

		// Print first 10 variables only (for readability)
		maxToPrint := 10
		if len(assignments) < maxToPrint {
			maxToPrint = len(assignments)
		}

		for i := 0; i < maxToPrint; i++ {
			fmt.Printf("  x%d = %v\n", i+1, assignments[i])
		}

		if len(assignments) > maxToPrint {
			fmt.Printf("  ... and %d more variables\n", len(assignments)-maxToPrint)
		}
	}
}
