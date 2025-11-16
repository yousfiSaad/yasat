# YASAT Test Suite

This directory contains test cases for the YASAT (Yet Another SAT) solver.

## Running Tests

From the repository root, run:
```bash
./tests/run_verification.sh
```

This will test all CNF files against a reference brute-force solver to verify correctness.

## Test Files

### Basic Tests
- `simple_unsat.cnf` - Simple contradiction (x AND ¬x)
- `test_simple_sat.cnf` - Simple satisfiable formula
- `test_single_var_sat.cnf` - Single positive unit clause
- `test_single_var_unsat.cnf` - Single variable contradiction

### Unit Propagation Tests
- `test_unit_propagation.cnf` - Chain of unit propagations
- `test_conflict_after_propagation.cnf` - Unit propagation leading to conflict
- `test_all_false_clause.cnf` - All literals in clause become false (tests the bug fix)

### 3-SAT Tests
- `test_3sat_sat.cnf` - Satisfiable 3-SAT instance
- `test_3sat_unsat.cnf` - Unsatisfiable 3-SAT (all 8 assignments blocked)

### Edge Cases
- `test_empty_clause.cnf` - Empty clause (immediately UNSAT)
- `test_tautology_sat.cnf` - Tautological clauses
- `test_long_clause_sat.cnf` - Long clause with 10 literals
- `test_many_unit_clauses.cnf` - Multiple unit clauses forcing assignments
- `test_chain_implications.cnf` - Chain of implications

## Verification Script

`verify_cnf.py` - Python brute-force SAT solver for small formulas (< 20 variables) used as a reference to verify YASAT's correctness.

## Known Issues

Some graph coloring problems currently fail (pre-existing bug, not related to the recent conflict detection fix).
