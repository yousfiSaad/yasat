#!/usr/bin/env python3
"""
Simple example demonstrating YASAT Python bindings usage.

Run with:
    python simple.py
    python simple.py path/to/problem.cnf
"""

import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent))

from yasat import Solver, Result, solve_file


def example1_manual_clauses():
    """Example 1: Manually adding clauses."""
    print("=== Example 1: Manual clause addition ===")

    # Create solver
    solver = Solver()
    print(f"YASAT version: {Solver.version()}")

    # Build a simple satisfiable formula:
    # (x1 ∨ x2) ∧ (¬x1 ∨ x3) ∧ (¬x2 ∨ ¬x3)
    print("\nAdding clauses:")
    print("  (x1 ∨ x2)")
    solver.add_clause(1, 2)

    print("  (¬x1 ∨ x3)")
    solver.add_clause(-1, 3)

    print("  (¬x2 ∨ ¬x3)")
    solver.add_clause(-2, -3)

    # Get formula stats
    print(f"\nFormula: {solver.num_variables()} variables, {solver.num_clauses()} clauses")

    # Solve
    print("\nSolving...")
    result = solver.solve()
    print(f"Result: {result.name}")

    if result == Result.SAT:
        print("\nSatisfying assignment:")
        assignments = solver.get_all_assignments()

        for i, val in enumerate(assignments, start=1):
            print(f"  x{i} = {val}")

        # Verify the solution
        x1, x2, x3 = assignments
        print("\nVerification:")
        print(f"  (x1 ∨ x2) = ({x1} ∨ {x2}) = {x1 or x2}")
        print(f"  (¬x1 ∨ x3) = ({not x1} ∨ {x3}) = {not x1 or x3}")
        print(f"  (¬x2 ∨ ¬x3) = ({not x2} ∨ {not x3}) = {not x2 or not x3}")


def example2_context_manager():
    """Example 2: Using context manager."""
    print("\n=== Example 2: Context manager usage ===")

    with Solver() as solver:
        # Add a simple unsatisfiable formula
        # (x1) ∧ (x2) ∧ (¬x1 ∨ ¬x2)
        print("Adding UNSAT formula:")
        print("  (x1) ∧ (x2) ∧ (¬x1 ∨ ¬x2)")

        solver.add_clause(1)
        solver.add_clause(2)
        solver.add_clause(-1, -2)

        result = solver.solve()
        print(f"Result: {result.name}")


def example3_from_string():
    """Example 3: Loading from DIMACS string."""
    print("\n=== Example 3: Loading from DIMACS string ===")

    cnf = """c Simple SAT problem
p cnf 3 3
1 2 0
-1 3 0
-2 -3 0"""

    solver = Solver()
    solver.parse_cnf_string(cnf)

    print(f"Loaded {solver.num_variables()} variables, {solver.num_clauses()} clauses")

    result = solver.solve()
    print(f"Result: {result.name}")

    if result == Result.SAT:
        # Get individual assignments
        print("\nVariable assignments:")
        for i in range(1, solver.num_variables() + 1):
            val = solver.get_assignment(i)
            print(f"  x{i} = {val}")


def example4_from_file(filename):
    """Example 4: Loading from file."""
    print(f"\n=== Example 4: Loading from file ===")
    print(f"File: {filename}")

    try:
        result, assignments = solve_file(filename)

        print(f"Result: {result.name}")

        if result == Result.SAT and assignments:
            print("\nSatisfying assignment found:")

            # Print first 10 variables only
            max_to_print = min(10, len(assignments))
            for i in range(max_to_print):
                print(f"  x{i+1} = {assignments[i]}")

            if len(assignments) > max_to_print:
                print(f"  ... and {len(assignments) - max_to_print} more variables")

    except Exception as e:
        print(f"Error: {e}")


def example5_error_handling():
    """Example 5: Error handling."""
    print("\n=== Example 5: Error handling ===")

    from yasat import YasatError

    try:
        solver = Solver()
        # Try to parse invalid file
        solver.parse_cnf_file("nonexistent.cnf")
    except YasatError as e:
        print(f"Caught error: {e}")


def main():
    print("YASAT Python Bindings Examples\n")

    # Run basic examples
    example1_manual_clauses()
    example2_context_manager()
    example3_from_string()
    example5_error_handling()

    # Run file example if filename provided
    if len(sys.argv) > 1:
        example4_from_file(sys.argv[1])
    else:
        # Try to load a test file if available
        test_file = Path(__file__).parent.parent.parent.parent / "tests" / "cnf" / "multi_sat.cnf"
        if test_file.exists():
            example4_from_file(str(test_file))


if __name__ == "__main__":
    main()
