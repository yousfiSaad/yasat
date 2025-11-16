# YASAT Test Suite

This directory contains the test suite for the YASAT SAT solver.

## Running Tests

### Automated Verification
From the repository root, run:
```bash
./tests/run_verification.sh
```

This will test all CNF files against a reference brute-force solver to verify correctness.

### Standard Test Suite
```bash
./tests/run_tests.sh
```

This will run all test cases and report which ones pass or fail.

## Test Cases

### Basic SAT Tests (Satisfiable)
- **simple_sat.cnf** - Single variable, single clause
- **test_simple_sat.cnf** - Simple satisfiable formula
- **empty_sat.cnf** - Empty CNF (no clauses, trivially satisfiable)
- **multi_sat.cnf** - Multiple variables with a satisfying assignment
- **test_single_var_sat.cnf** - Single positive unit clause

### Basic UNSAT Tests (Unsatisfiable)
- **simple_unsat.cnf** - Simple contradiction (x AND ¬x)
- **test_single_var_unsat.cnf** - Single variable contradiction
- **small_unsat.cnf** - Small contradictory formula requiring search
- **test_3sat_unsat.cnf** - Unsatisfiable 3-SAT (all 8 assignments blocked)

### Unit Propagation Tests
- **test_unit_propagation.cnf** - Chain of unit propagations
- **test_conflict_after_propagation.cnf** - Unit propagation leading to conflict
- **test_all_false_clause.cnf** - All literals in clause become false (tests the bug fix)

### 3-SAT Tests
- **test_3sat_sat.cnf** - Satisfiable 3-SAT instance
- **test_3sat_unsat.cnf** - Unsatisfiable 3-SAT instance

### Edge Case Tests
Comprehensive edge case coverage to ensure solver robustness:

- **edge_large_clause.cnf** - Large clause with 8 literals
  - Tests handling of clauses with many disjuncts
  - Ensures solver doesn't choke on large clauses

- **edge_unit_clause.cnf** - Unit clause forcing assignment
  - Tests unit propagation at initialization
  - Forces x1=true immediately

- **edge_all_negative.cnf** - All negative literals
  - Tests solver with only negated variables
  - (¬x1 ∨ ¬x2) ∧ (¬x1 ∨ ¬x3)

- **edge_all_positive.cnf** - All positive literals
  - Tests solver with only positive variables
  - (x1 ∨ x2) ∧ (x1 ∨ x3)

- **edge_duplicate_literals.cnf** - Duplicate literals in same clause
  - Tests handling of (x1 ∨ x1 ∨ x2)
  - Ensures duplicates don't break solver

- **edge_tautology.cnf** - Tautological clause
  - Tests handling of (x1 ∨ ¬x1 ∨ x2) - always true clause
  - Solver should handle gracefully

- **test_tautology_sat.cnf** - Tautological clauses

- **edge_unit_chain.cnf** - Chain of unit clauses
  - Multiple unit clauses forcing all assignments
  - Tests full propagation: x1=1, x2=0, x3=1, x4=0

- **edge_horn_sat.cnf** - Horn clauses
  - At most one positive literal per clause
  - Tests special case of polynomial-time SAT

- **edge_pure_literal.cnf** - Pure literal elimination
  - Variable appears only positive (never negated)
  - Tests if solver can exploit pure literals

- **edge_underconstrained.cnf** - Many variables, few constraints
  - 10 variables but only 2 clauses
  - Tests solver on highly underconstrained problems

- **test_empty_clause.cnf** - Empty clause (immediately UNSAT)
- **test_long_clause_sat.cnf** - Long clause with 10 literals
- **test_many_unit_clauses.cnf** - Multiple unit clauses forcing assignments
- **test_chain_implications.cnf** - Chain of implications

### Pigeonhole Principle Tests
- **php_2p_3h.cnf** - 2 pigeons, 3 holes (SAT)
- **php_3p_2h.cnf** - 3 pigeons, 2 holes (UNSAT)
- **test_php_4p_3h.cnf** - 4 pigeons, 3 holes (UNSAT)

### Graph Coloring Tests
- **test_graph_coloring_sat.cnf** - Triangle graph 3-coloring
- **test_graph_simple.cnf** - Simple 2-vertex graph coloring
- **test_3vertices_partial.cnf** - 3 vertices with 2 colors

## Verification Tools

- **verify_cnf.py** - Python brute-force SAT solver for small formulas (< 20 variables) used as a reference to verify YASAT's correctness
- **run_verification.sh** - Automated test runner that compares YASAT output with the reference solver

## Test Statistics

- **Basic SAT**: Multiple tests
- **Basic UNSAT**: Multiple tests
- **Edge cases**: 10+ comprehensive edge case tests
- **Unit propagation**: 3 tests
- **3-SAT**: 2 tests
- **Pigeonhole**: 3 tests
- **Graph coloring**: 3 tests

## Known Issues

### Fixed Issues
- ✅ **Conflict detection bug** - Fixed in recent commit. The solver now correctly detects when all literals in a clause become false during propagation.

### Pre-existing Limitations
- ⚠️ Some graph coloring problems currently fail (pre-existing bug, not related to the recent conflict detection fix)

## Adding New Tests

To add a new test case:

1. Create a CNF file in `tests/cnf/` (or `tests/`) following the DIMACS format
2. Add a test case to the appropriate test script
3. Verify the test passes

## DIMACS CNF Format

Test files use the DIMACS CNF format:
```
c Comments start with 'c'
p cnf <num_variables> <num_clauses>
<clause_1>
<clause_2>
...
```

Each clause is a space-separated list of literals ending with `0`. Positive numbers represent positive literals, negative numbers represent negated literals.

Example:
```
c (a OR b) AND (NOT a OR c)
p cnf 3 2
1 2 0
-1 3 0
```

## Test Categories Explained

### Why These Edge Cases?

1. **Large clauses** - Ensures no buffer overflows or performance issues
2. **Unit clauses** - Tests unit propagation (core SAT technique)
3. **All negative/positive** - Tests polarity handling
4. **Duplicates** - Tests parsing robustness
5. **Tautologies** - Tests trivially satisfied clauses
6. **Unit chains** - Tests propagation cascade
7. **Horn clauses** - Special polynomial case
8. **Pure literals** - Tests optimization opportunities
9. **Underconstrained** - Tests search on easy problems
10. **Pigeonhole** - Classic hard SAT/UNSAT examples
11. **Conflict detection** - Tests the recent bug fix for detecting conflicts when all clause literals are false

These edge cases ensure the solver is robust across a wide variety of CNF formulas.
