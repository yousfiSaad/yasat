#!/bin/bash

cd "$(dirname "$0")/.." || exit 1

printf "%-40s %-30s %-30s %s\n" "Test File" "yasat Result" "Reference Result" "Match"
printf "%-40s %-30s %-30s %s\n" "---------------------------------------" "----------------------------" "----------------------------" "-----"

for file in tests/simple_unsat.cnf tests/test_simple_sat.cnf tests/test_unit_propagation.cnf \
            tests/test_conflict_after_propagation.cnf tests/test_all_false_clause.cnf \
            tests/test_3sat_sat.cnf tests/test_3sat_unsat.cnf tests/test_empty_clause.cnf \
            tests/test_tautology_sat.cnf tests/test_single_var_sat.cnf tests/test_single_var_unsat.cnf \
            tests/test_long_clause_sat.cnf tests/test_many_unit_clauses.cnf tests/test_chain_implications.cnf \
            data/php_2p_3h.cnf data/php_3p_2h.cnf; do

    yasat_result=$(cat "$file" | ./build/yasat)
    ref_result=$(python3 tests/verify_cnf.py "$file")

    # Extract just SAT/UNSAT for comparison
    yasat_status=$(echo "$yasat_result" | grep -o "^[A-Z]*")
    ref_status=$(echo "$ref_result" | grep -o "^[A-Z]*")

    if [ "$yasat_status" = "$ref_status" ]; then
        match="✓"
    else
        match="✗ MISMATCH"
    fi

    printf "%-40s %-30s %-30s %s\n" "$file" "$yasat_result" "$ref_result" "$match"
done
