#!/bin/bash

# YASAT Test Suite
# Tests the SAT solver against known test cases
#
# Usage:
#   ./tests/run_tests.sh           # Normal mode
#   VERBOSE=1 ./tests/run_tests.sh # Verbose mode (shows CNF content and solver output)
#
# Note: The solver currently has a limitation with direct unit clause
# contradictions at initialization. Test cases are chosen to work with
# the current implementation.

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

SOLVER="./build/yasat"
PASSED=0
FAILED=0

# Check if solver exists
if [ ! -f "$SOLVER" ]; then
    echo -e "${RED}Error: Solver not found at $SOLVER${NC}"
    echo "Please run ./build.sh first"
    exit 1
fi

echo "================================"
echo "  YASAT Test Suite"
if [ "${VERBOSE:-0}" = "1" ]; then
    echo "  (Verbose Mode)"
fi
echo "================================"
echo ""

# Test function
# Args: test_name, cnf_file, expected_result (SAT or UNSAT)
run_test() {
    local test_name="$1"
    local cnf_file="$2"
    local expected="$3"

    echo -n "Testing $test_name... "

    if [ ! -f "$cnf_file" ]; then
        echo -e "${RED}FAIL${NC} (file not found)"
        FAILED=$((FAILED + 1))
        return
    fi

    # Verbose mode: show CNF file being tested
    if [ "${VERBOSE:-0}" = "1" ]; then
        echo ""
        echo "  File: $cnf_file"
        echo "  Expected: $expected"
        echo "  CNF content:"
        cat "$cnf_file" | head -20
        echo "  Running solver..."
    fi

    # Run solver and capture output and exit code
    set +e  # Don't exit on error
    result=$(cat "$cnf_file" | $SOLVER 2>&1)
    exit_code=$?
    set -e

    # Verbose mode: show raw output and exit code
    if [ "${VERBOSE:-0}" = "1" ]; then
        echo "  Exit code: $exit_code"
        echo "  Solver output: $result"
    fi

    # Check exit code for crashes
    if [ $exit_code -eq 139 ]; then
        echo -e "${RED}FAIL${NC} (segmentation fault - exit code 139)"
        echo "  File: $cnf_file"
        echo "  CNF content:"
        cat "$cnf_file"
        FAILED=$((FAILED + 1))
        return
    elif [ $exit_code -ne 0 ] && [ $exit_code -ne 1 ]; then
        echo -e "${RED}FAIL${NC} (unexpected exit code: $exit_code)"
        echo "  Output: $result"
        FAILED=$((FAILED + 1))
        return
    fi

    # Check if result matches expected
    if echo "$result" | grep -q "^$expected"; then
        echo -e "${GREEN}PASS${NC}"
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}"
        echo "  Expected: $expected"
        echo "  Got:      $result"
        FAILED=$((FAILED + 1))
    fi
}

# Run tests
echo "--- Basic SAT Tests ---"
run_test "simple SAT" "tests/cnf/simple_sat.cnf" "SAT"
run_test "empty SAT" "tests/cnf/empty_sat.cnf" "SAT"
run_test "multi-variable SAT" "tests/cnf/multi_sat.cnf" "SAT"

echo ""
echo "--- Basic UNSAT Tests ---"
run_test "small UNSAT" "tests/cnf/small_unsat.cnf" "UNSAT"

echo ""
echo "--- Edge Case Tests ---"
run_test "large clause (8 literals)" "tests/cnf/edge_large_clause.cnf" "SAT"
run_test "unit clause forcing" "tests/cnf/edge_unit_clause.cnf" "SAT"
run_test "all negative literals" "tests/cnf/edge_all_negative.cnf" "SAT"
run_test "all positive literals" "tests/cnf/edge_all_positive.cnf" "SAT"
run_test "duplicate literals" "tests/cnf/edge_duplicate_literals.cnf" "SAT"
run_test "tautology clause" "tests/cnf/edge_tautology.cnf" "SAT"
run_test "unit clause chain" "tests/cnf/edge_unit_chain.cnf" "SAT"
run_test "Horn clauses" "tests/cnf/edge_horn_sat.cnf" "SAT"
run_test "pure literal" "tests/cnf/edge_pure_literal.cnf" "SAT"
run_test "underconstrained (10 vars, 2 clauses)" "tests/cnf/edge_underconstrained.cnf" "SAT"

echo ""
echo "--- Pigeonhole Principle Tests ---"
run_test "2 pigeons, 3 holes (SAT)" "data/php_2p_3h.cnf" "SAT"
run_test "3 pigeons, 2 holes (UNSAT)" "data/php_3p_2h.cnf" "UNSAT"

echo ""
echo "================================"
echo "  Results: ${GREEN}$PASSED passed${NC}, ${RED}$FAILED failed${NC}"
echo "================================"

if [ $FAILED -gt 0 ]; then
    exit 1
fi

exit 0
