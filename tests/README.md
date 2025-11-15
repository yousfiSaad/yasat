# YASAT Test Suite

This directory contains the test suite for the YASAT SAT solver.

## Running Tests

```bash
./tests/run_tests.sh
```

This will run all test cases and report which ones pass or fail.

## Test Cases

### SAT Tests (Satisfiable)
- **simple_sat.cnf** - Single variable, single clause
- **empty_sat.cnf** - Empty CNF (no clauses, trivially satisfiable)
- **multi_sat.cnf** - Multiple variables with a satisfying assignment
- **php_2p_3h.cnf** - Pigeonhole principle: 2 pigeons, 3 holes (SAT)

### UNSAT Tests (Unsatisfiable)
- **small_unsat.cnf** - Small contradictory formula requiring search
- **php_3p_2h.cnf** - Pigeonhole principle: 3 pigeons, 2 holes (UNSAT)

## Known Limitations

The solver currently has a limitation with direct unit clause contradictions at initialization (e.g., `x` and `¬x` as separate unit clauses). Test cases that require search-based conflict discovery work correctly.

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
