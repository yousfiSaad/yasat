#!/usr/bin/env python3
"""
Simple brute-force SAT solver for verification purposes.
Only works for small formulas (< 20 variables).
"""
import sys
from itertools import product

def parse_cnf(filename):
    """Parse a DIMACS CNF file."""
    clauses = []
    num_vars = 0

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('c'):
                continue
            if line.startswith('p'):
                parts = line.split()
                num_vars = int(parts[2])
                continue

            literals = [int(x) for x in line.split() if x != '0']
            if literals:  # Skip empty clauses in the parse (will handle below)
                clauses.append(literals)
            elif line.strip() == '0':  # Empty clause
                clauses.append([])

    return num_vars, clauses

def evaluate_clause(clause, assignment):
    """Evaluate a clause given an assignment."""
    if not clause:  # Empty clause is always false
        return False

    for lit in clause:
        var = abs(lit)
        val = assignment[var]
        if (lit > 0 and val) or (lit < 0 and not val):
            return True
    return False

def evaluate_formula(clauses, assignment):
    """Check if all clauses are satisfied."""
    return all(evaluate_clause(clause, assignment) for clause in clauses)

def brute_force_sat(num_vars, clauses):
    """Brute force SAT solver."""
    if num_vars == 0:
        return len(clauses) == 0 or all(len(c) > 0 for c in clauses), []

    if num_vars > 20:
        print(f"Too many variables ({num_vars}) for brute force")
        return None, None

    # Try all possible assignments
    for values in product([False, True], repeat=num_vars):
        assignment = {i+1: values[i] for i in range(num_vars)}
        if evaluate_formula(clauses, assignment):
            return True, [1 if v else 0 for v in values]

    return False, None

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: verify_cnf.py <file.cnf>")
        sys.exit(1)

    filename = sys.argv[1]
    num_vars, clauses = parse_cnf(filename)

    sat, solution = brute_force_sat(num_vars, clauses)

    if sat is None:
        print("UNKNOWN (too large)")
    elif sat:
        print(f"SAT {solution}")
    else:
        print("UNSAT")
