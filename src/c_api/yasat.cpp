/**
 * @file yasat.cpp
 * @brief C API implementation for YASAT
 *
 * This file implements the C API by wrapping the C++ CDCL_solver class.
 */

#include "yasat.h"
#include "../headers/CDCL_solver.h"
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>

/**
 * @brief Internal solver structure
 *
 * This structure wraps the C++ solver and maintains state needed for the C API.
 */
struct yasat_solver {
    CDCL_solver* cpp_solver;           // Pointer to C++ solver instance
    std::vector<bool> solution;         // Solution vector from last solve
    yasat_result last_result;           // Result of last solve
    std::string error_message;          // Last error message
    std::stringstream cnf_stream;       // Stream for building CNF
    int num_variables;                  // Number of variables declared
    int num_clauses;                    // Number of clauses declared
    int clauses_added;                  // Number of clauses actually added
    bool has_header;                    // Whether DIMACS header was set

    yasat_solver()
        : cpp_solver(nullptr)
        , last_result(YASAT_RESULT_UNKNOWN)
        , num_variables(0)
        , num_clauses(0)
        , clauses_added(0)
        , has_header(false)
    {}

    ~yasat_solver() {
        delete cpp_solver;
    }
};

/* ========================================================================== */
/* Helper Functions                                                           */
/* ========================================================================== */

static void set_error(yasat_solver* solver, const char* message) {
    if (solver) {
        solver->error_message = message;
    }
}

static void clear_error(yasat_solver* solver) {
    if (solver) {
        solver->error_message.clear();
    }
}

/**
 * @brief Finalize the CNF and create the C++ solver
 * @return true on success, false on failure
 */
static bool finalize_cnf(yasat_solver* solver) {
    if (!solver) return false;
    if (solver->cpp_solver) return true; // Already finalized

    try {
        // Create C++ solver from the CNF stream
        solver->cpp_solver = new CDCL_solver(solver->cnf_stream);
        clear_error(solver);
        return true;
    } catch (const std::exception& e) {
        set_error(solver, e.what());
        return false;
    } catch (...) {
        set_error(solver, "Unknown error creating solver");
        return false;
    }
}

/* ========================================================================== */
/* Solver Lifecycle Management                                                */
/* ========================================================================== */

yasat_solver* yasat_solver_create(void) {
    try {
        return new yasat_solver();
    } catch (...) {
        return nullptr;
    }
}

void yasat_solver_destroy(yasat_solver* solver) {
    delete solver;
}

/* ========================================================================== */
/* Input Methods                                                              */
/* ========================================================================== */

yasat_error yasat_add_clause(yasat_solver* solver,
                             const int* literals,
                             int num_literals) {
    if (!solver) return YASAT_ERROR_INVALID_INPUT;
    if (!literals && num_literals > 0) return YASAT_ERROR_INVALID_INPUT;
    if (num_literals < 0) return YASAT_ERROR_INVALID_INPUT;

    // Cannot add clauses after solver is finalized
    if (solver->cpp_solver) {
        set_error(solver, "Cannot add clauses after solving");
        return YASAT_ERROR_STATE;
    }

    try {
        // Track maximum variable number
        int max_var = solver->num_variables;
        for (int i = 0; i < num_literals; ++i) {
            int var = literals[i] > 0 ? literals[i] : -literals[i];
            if (var > max_var) {
                max_var = var;
            }
        }

        // Update variable count
        if (max_var > solver->num_variables) {
            solver->num_variables = max_var;
        }

        // Mark that we're using manual clause addition
        solver->has_header = true;

        // Add the clause to the stream (without header - we'll add it later)
        for (int i = 0; i < num_literals; ++i) {
            solver->cnf_stream << literals[i] << " ";
        }
        solver->cnf_stream << "0\n";
        solver->clauses_added++;

        clear_error(solver);
        return YASAT_OK;
    } catch (const std::exception& e) {
        set_error(solver, e.what());
        return YASAT_ERROR_MEMORY;
    } catch (...) {
        set_error(solver, "Unknown error adding clause");
        return YASAT_ERROR_MEMORY;
    }
}

yasat_error yasat_parse_cnf_file(yasat_solver* solver, const char* filename) {
    if (!solver) return YASAT_ERROR_INVALID_INPUT;
    if (!filename) return YASAT_ERROR_INVALID_INPUT;

    // Cannot parse after solver is created
    if (solver->cpp_solver) {
        set_error(solver, "Cannot parse CNF after solver is created");
        return YASAT_ERROR_STATE;
    }

    try {
        std::ifstream file(filename);
        if (!file.is_open()) {
            set_error(solver, "Failed to open file");
            return YASAT_ERROR_FILE;
        }

        // Read entire file into the CNF stream
        std::string line;
        while (std::getline(file, line)) {
            solver->cnf_stream << line << "\n";

            // Extract variable count from header line
            if (line.length() > 0 && line[0] == 'p') {
                std::istringstream iss(line);
                std::string p, cnf;
                int vars, clauses;
                if (iss >> p >> cnf >> vars >> clauses) {
                    solver->num_variables = vars;
                    solver->num_clauses = clauses;
                }
            }
        }
        file.close();

        // Mark as having header (assuming valid DIMACS file)
        solver->has_header = true;

        clear_error(solver);
        return YASAT_OK;
    } catch (const std::exception& e) {
        set_error(solver, e.what());
        return YASAT_ERROR_FILE;
    } catch (...) {
        set_error(solver, "Unknown error reading file");
        return YASAT_ERROR_FILE;
    }
}

yasat_error yasat_parse_cnf_string(yasat_solver* solver, const char* cnf_string) {
    if (!solver) return YASAT_ERROR_INVALID_INPUT;
    if (!cnf_string) return YASAT_ERROR_INVALID_INPUT;

    // Cannot parse after solver is created
    if (solver->cpp_solver) {
        set_error(solver, "Cannot parse CNF after solver is created");
        return YASAT_ERROR_STATE;
    }

    try {
        // Parse the string line by line to extract header info
        std::istringstream input(cnf_string);
        std::string line;

        while (std::getline(input, line)) {
            solver->cnf_stream << line << "\n";

            // Extract variable count from header line
            if (line.length() > 0 && line[0] == 'p') {
                std::istringstream iss(line);
                std::string p, cnf;
                int vars, clauses;
                if (iss >> p >> cnf >> vars >> clauses) {
                    solver->num_variables = vars;
                    solver->num_clauses = clauses;
                }
            }
        }

        solver->has_header = true; // Assuming valid DIMACS format

        clear_error(solver);
        return YASAT_OK;
    } catch (const std::exception& e) {
        set_error(solver, e.what());
        return YASAT_ERROR_MEMORY;
    } catch (...) {
        set_error(solver, "Unknown error parsing string");
        return YASAT_ERROR_MEMORY;
    }
}

/* ========================================================================== */
/* Solving                                                                    */
/* ========================================================================== */

yasat_result yasat_solve(yasat_solver* solver) {
    if (!solver) return YASAT_RESULT_ERROR;

    try {
        // If we added clauses manually, prepend the DIMACS header
        if (!solver->cpp_solver && solver->has_header && solver->clauses_added > 0) {
            // Get the clauses we've collected
            std::string clauses_content = solver->cnf_stream.str();

            // Clear and rebuild with proper header
            solver->cnf_stream.str("");
            solver->cnf_stream.clear();

            // Write proper DIMACS header
            solver->cnf_stream << "p cnf " << solver->num_variables << " "
                              << solver->clauses_added << "\n";

            // Write all the clauses
            solver->cnf_stream << clauses_content;
        }

        // Finalize and create C++ solver if not already done
        if (!finalize_cnf(solver)) {
            return YASAT_RESULT_ERROR;
        }

        // Solve
        solver->solution.clear();
        bool result = solver->cpp_solver->solve(solver->solution);

        solver->last_result = result ? YASAT_RESULT_SAT : YASAT_RESULT_UNSAT;
        clear_error(solver);
        return solver->last_result;

    } catch (const std::exception& e) {
        set_error(solver, e.what());
        solver->last_result = YASAT_RESULT_ERROR;
        return YASAT_RESULT_ERROR;
    } catch (...) {
        set_error(solver, "Unknown error during solving");
        solver->last_result = YASAT_RESULT_ERROR;
        return YASAT_RESULT_ERROR;
    }
}

/* ========================================================================== */
/* Solution Retrieval                                                         */
/* ========================================================================== */

int yasat_get_assignment(yasat_solver* solver, int var) {
    if (!solver) return -1;
    if (solver->last_result != YASAT_RESULT_SAT) return -1;
    if (var < 1) return -1;

    // Solution vector is 0-indexed, but variables are 1-indexed
    size_t idx = static_cast<size_t>(var - 1);
    if (idx >= solver->solution.size()) return -1;

    return solver->solution[idx] ? 1 : 0;
}

yasat_error yasat_get_all_assignments(yasat_solver* solver,
                                      int* assignments,
                                      int num_vars) {
    if (!solver) return YASAT_ERROR_INVALID_INPUT;
    if (!assignments) return YASAT_ERROR_INVALID_INPUT;
    if (num_vars < 0) return YASAT_ERROR_INVALID_INPUT;
    if (solver->last_result != YASAT_RESULT_SAT) return YASAT_ERROR_STATE;

    try {
        for (int i = 0; i < num_vars && i < static_cast<int>(solver->solution.size()); ++i) {
            assignments[i] = solver->solution[i] ? 1 : 0;
        }
        clear_error(solver);
        return YASAT_OK;
    } catch (...) {
        set_error(solver, "Error retrieving assignments");
        return YASAT_ERROR_MEMORY;
    }
}

int yasat_get_num_variables(yasat_solver* solver) {
    if (!solver) return -1;
    if (!solver->cpp_solver) {
        // If solver not yet created, return the tracked count
        return solver->num_variables;
    }

    try {
        // cnfSize returns 2 * num_variables (for positive and negative literals)
        return static_cast<int>(solver->cpp_solver->cnfSize() / 2);
    } catch (...) {
        return -1;
    }
}

int yasat_get_num_clauses(yasat_solver* solver) {
    if (!solver) return -1;

    if (!solver->cpp_solver) {
        // If solver not yet created, return the tracked count
        // For parsed CNF, use num_clauses; for manually added, use clauses_added
        return solver->num_clauses > 0 ? solver->num_clauses : solver->clauses_added;
    }

    try {
        return static_cast<int>(solver->cpp_solver->getCleanedCnf().size());
    } catch (...) {
        return -1;
    }
}

/* ========================================================================== */
/* Error Handling                                                             */
/* ========================================================================== */

const char* yasat_get_error_message(yasat_solver* solver) {
    if (!solver) return "Invalid solver instance";
    if (solver->error_message.empty()) return nullptr;
    return solver->error_message.c_str();
}

const char* yasat_get_version(void) {
    static char version[32];
    snprintf(version, sizeof(version), "%d.%d.%d",
             YASAT_VERSION_MAJOR,
             YASAT_VERSION_MINOR,
             YASAT_VERSION_PATCH);
    return version;
}
