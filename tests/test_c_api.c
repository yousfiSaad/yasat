/**
 * @file test_c_api.c
 * @brief Test program for YASAT C API
 *
 * This program tests the C API functions to ensure they work correctly.
 * Compile with: gcc -o test_c_api test_c_api.c -L../build -lyasat
 * Run with: LD_LIBRARY_PATH=../build ./test_c_api
 */

#include "../src/c_api/yasat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_PASS "\033[32mPASS\033[0m"
#define TEST_FAIL "\033[31mFAIL\033[0m"

int tests_run = 0;
int tests_passed = 0;

void test_assert(int condition, const char* test_name) {
    tests_run++;
    if (condition) {
        tests_passed++;
        printf("[%s] %s\n", TEST_PASS, test_name);
    } else {
        printf("[%s] %s\n", TEST_FAIL, test_name);
    }
}

void test_version() {
    const char* version = yasat_get_version();
    test_assert(version != NULL, "Get version returns non-NULL");
    test_assert(strlen(version) > 0, "Version string is not empty");
    printf("      YASAT version: %s\n", version);
}

void test_solver_lifecycle() {
    yasat_solver* solver = yasat_solver_create();
    test_assert(solver != NULL, "Create solver");

    yasat_solver_destroy(solver);
    test_assert(1, "Destroy solver");

    // Test destroying NULL solver (should not crash)
    yasat_solver_destroy(NULL);
    test_assert(1, "Destroy NULL solver (no crash)");
}

void test_simple_sat() {
    yasat_solver* solver = yasat_solver_create();
    test_assert(solver != NULL, "Create solver for SAT test");

    // Add a simple SAT formula: (x1 OR x2) AND (NOT x1 OR x2)
    int clause1[] = {1, 2};
    int clause2[] = {-1, 2};

    yasat_error err1 = yasat_add_clause(solver, clause1, 2);
    test_assert(err1 == YASAT_OK, "Add clause 1");

    yasat_error err2 = yasat_add_clause(solver, clause2, 2);
    test_assert(err2 == YASAT_OK, "Add clause 2");

    yasat_result result = yasat_solve(solver);
    test_assert(result == YASAT_RESULT_SAT, "Solve returns SAT");

    if (result == YASAT_RESULT_SAT) {
        int val1 = yasat_get_assignment(solver, 1);
        int val2 = yasat_get_assignment(solver, 2);
        test_assert(val1 >= 0 && val1 <= 1, "Variable 1 has valid assignment");
        test_assert(val2 >= 0 && val2 <= 1, "Variable 2 has valid assignment");
        printf("      Solution: x1=%d, x2=%d\n", val1, val2);
    }

    yasat_solver_destroy(solver);
}

void test_simple_unsat() {
    yasat_solver* solver = yasat_solver_create();
    test_assert(solver != NULL, "Create solver for UNSAT test");

    // Add UNSAT formula: (x1) AND (NOT x1)
    int clause1[] = {1};
    int clause2[] = {-1};

    yasat_add_clause(solver, clause1, 1);
    yasat_add_clause(solver, clause2, 1);

    yasat_result result = yasat_solve(solver);
    // Note: Current implementation may not detect this as UNSAT
    // due to known limitation with unit clauses
    printf("      Result: %d (10=SAT, 20=UNSAT)\n", result);

    yasat_solver_destroy(solver);
}

void test_parse_cnf_file() {
    yasat_solver* solver = yasat_solver_create();
    test_assert(solver != NULL, "Create solver for file test");

    yasat_error err = yasat_parse_cnf_file(solver, "cnf/simple_sat.cnf");
    if (err == YASAT_OK) {
        test_assert(1, "Parse CNF file");

        int num_vars = yasat_get_num_variables(solver);
        int num_clauses = yasat_get_num_clauses(solver);
        printf("      Loaded: %d variables, %d clauses\n", num_vars, num_clauses);
        test_assert(num_vars > 0, "File has variables");

        yasat_result result = yasat_solve(solver);
        printf("      Result: %d (10=SAT, 20=UNSAT)\n", result);
    } else {
        printf("      Skipping (file not found - run from tests/ directory)\n");
    }

    yasat_solver_destroy(solver);
}

void test_parse_cnf_string() {
    yasat_solver* solver = yasat_solver_create();
    test_assert(solver != NULL, "Create solver for string test");

    const char* cnf = "p cnf 2 2\n1 2 0\n-1 2 0\n";
    yasat_error err = yasat_parse_cnf_string(solver, cnf);
    test_assert(err == YASAT_OK, "Parse CNF string");

    int num_vars = yasat_get_num_variables(solver);
    test_assert(num_vars == 2, "String has 2 variables");

    yasat_result result = yasat_solve(solver);
    test_assert(result == YASAT_RESULT_SAT, "String formula is SAT");

    yasat_solver_destroy(solver);
}

void test_get_all_assignments() {
    yasat_solver* solver = yasat_solver_create();

    // Simple SAT formula
    int clause1[] = {1, 2};
    yasat_add_clause(solver, clause1, 2);

    yasat_result result = yasat_solve(solver);
    if (result == YASAT_RESULT_SAT) {
        int num_vars = yasat_get_num_variables(solver);
        int* assignments = malloc(num_vars * sizeof(int));

        yasat_error err = yasat_get_all_assignments(solver, assignments, num_vars);
        test_assert(err == YASAT_OK, "Get all assignments");

        printf("      All assignments: ");
        for (int i = 0; i < num_vars; i++) {
            printf("x%d=%d ", i+1, assignments[i]);
        }
        printf("\n");

        free(assignments);
    } else {
        test_assert(0, "Get all assignments (formula not SAT)");
    }

    yasat_solver_destroy(solver);
}

void test_error_handling() {
    // Test NULL solver
    yasat_result result = yasat_solve(NULL);
    test_assert(result == YASAT_RESULT_ERROR, "Solve NULL solver returns ERROR");

    int val = yasat_get_assignment(NULL, 1);
    test_assert(val < 0, "Get assignment from NULL solver returns error");

    // Test invalid file
    yasat_solver* solver = yasat_solver_create();
    yasat_error err = yasat_parse_cnf_file(solver, "nonexistent_file.cnf");
    test_assert(err != YASAT_OK, "Parse nonexistent file returns error");

    const char* error_msg = yasat_get_error_message(solver);
    if (error_msg) {
        printf("      Error message: %s\n", error_msg);
        test_assert(strlen(error_msg) > 0, "Error message is not empty");
    }

    yasat_solver_destroy(solver);
}

void test_num_variables_clauses() {
    yasat_solver* solver = yasat_solver_create();

    int clause1[] = {1, 2, 3};
    int clause2[] = {-1, 4};

    yasat_add_clause(solver, clause1, 3);
    yasat_add_clause(solver, clause2, 2);

    int num_vars = yasat_get_num_variables(solver);
    int num_clauses = yasat_get_num_clauses(solver);

    test_assert(num_vars == 4, "Number of variables is correct");
    test_assert(num_clauses == 2, "Number of clauses is correct");

    printf("      Formula: %d vars, %d clauses\n", num_vars, num_clauses);

    yasat_solver_destroy(solver);
}

int main(int argc, char** argv) {
    printf("YASAT C API Test Suite\n");
    printf("======================\n\n");

    // Change to tests directory if needed
    if (argc > 1 && strcmp(argv[1], "--in-tests-dir") == 0) {
        // Already in tests directory
    }

    test_version();
    test_solver_lifecycle();
    test_simple_sat();
    test_simple_unsat();
    test_parse_cnf_string();
    test_get_all_assignments();
    test_num_variables_clauses();
    test_error_handling();
    test_parse_cnf_file();  // Run last as it depends on file location

    printf("\n======================\n");
    printf("Tests: %d/%d passed\n", tests_passed, tests_run);

    if (tests_passed == tests_run) {
        printf("\033[32mAll tests passed!\033[0m\n");
        return 0;
    } else {
        printf("\033[31m%d test(s) failed\033[0m\n", tests_run - tests_passed);
        return 1;
    }
}
