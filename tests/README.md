# YASAT Test Suite

This directory contains the test suite for the YASAT SAT solver.

## Running Tests

```bash
./tests/run_tests.sh
```

This will run all test cases and report which ones pass or fail.

## Test Cases

### Basic SAT Tests (Satisfiable)
- **simple_sat.cnf** - Single variable, single clause
- **empty_sat.cnf** - Empty CNF (no clauses, trivially satisfiable)
- **multi_sat.cnf** - Multiple variables with a satisfying assignment

### Basic UNSAT Tests (Unsatisfiable)
- **small_unsat.cnf** - Small contradictory formula requiring search

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

### Pigeonhole Principle Tests
- **php_2p_3h.cnf** - 2 pigeons, 3 holes (SAT)
- **php_3p_2h.cnf** - 3 pigeons, 2 holes (UNSAT)

## Test Statistics

- **Total tests**: 16
- **Basic SAT**: 3 tests
- **Basic UNSAT**: 1 test
- **Edge cases**: 10 tests
- **Pigeonhole**: 2 tests

## Known Limitations

The solver currently has a limitation with direct unit clause contradictions at initialization (e.g., `x` and `¬x` as separate unit clauses). Test cases that require search-based conflict discovery work correctly.

**Example of limitation:**
```
p cnf 1 2
1 0      # x1 must be true
-1 0     # x1 must be false
```
This returns SAT but should return UNSAT. The solver handles contradictions found during search correctly.

## Adding New Tests

To add a new test case:

1. Create a CNF file in `tests/cnf/` following the DIMACS format
2. Add a test case to `run_tests.sh` using the `run_test` function
3. Verify the test passes with `./tests/run_tests.sh`

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

These edge cases ensure the solver is robust across a wide variety of CNF formulas.
