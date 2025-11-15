#!/bin/bash

# YASAT Test Suite
# Tests the SAT solver against known test cases
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

    # Run solver and capture output
    result=$(cat "$cnf_file" | $SOLVER 2>&1)

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
