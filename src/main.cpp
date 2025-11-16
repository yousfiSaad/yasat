
#include "headers/CDCL_solver.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <exception>
#include <string>
#include <cstring>

struct Options {
  bool verbose = false;
  bool show_stats = false;
  bool show_help = false;
  std::string input_file;
};

void print_help(const char* program_name) {
  std::cout << "YASAT - Yet Another SAT Solver\n";
  std::cout << "A CDCL-based SAT solver for CNF formulas\n\n";
  std::cout << "Usage: " << program_name << " [options] [input.cnf]\n\n";
  std::cout << "Options:\n";
  std::cout << "  -h, --help       Show this help message and exit\n";
  std::cout << "  -v, --verbose    Enable verbose output\n";
  std::cout << "  -s, --stats      Show solving statistics\n";
  std::cout << "  input.cnf        Input CNF file in DIMACS format\n";
  std::cout << "                   (if not provided, reads from stdin)\n\n";
  std::cout << "Examples:\n";
  std::cout << "  " << program_name << " problem.cnf          # Solve from file\n";
  std::cout << "  cat problem.cnf | " << program_name << "   # Solve from stdin\n";
  std::cout << "  " << program_name << " -v problem.cnf       # Verbose output\n";
  std::cout << "  " << program_name << " -s problem.cnf       # Show statistics\n\n";
  std::cout << "Exit codes:\n";
  std::cout << "  0  - Success (SAT or UNSAT)\n";
  std::cout << "  1  - Error (invalid input, file not found, etc.)\n";
}

Options parse_args(int argc, char* argv[]) {
  Options opts;

  for (int i = 1; i < argc; i++) {
    std::string arg = argv[i];

    if (arg == "-h" || arg == "--help") {
      opts.show_help = true;
    }
    else if (arg == "-v" || arg == "--verbose") {
      opts.verbose = true;
    }
    else if (arg == "-s" || arg == "--stats") {
      opts.show_stats = true;
    }
    else if (arg[0] == '-') {
      throw std::runtime_error("Unknown option: " + arg + "\nUse -h or --help for usage information.");
    }
    else {
      // Assume it's the input file
      if (!opts.input_file.empty()) {
        throw std::runtime_error("Multiple input files specified. Only one input file is allowed.");
      }
      opts.input_file = arg;
    }
  }

  return opts;
}

int main(int argc, char *argv[]) {
  try {
    Options opts = parse_args(argc, argv);

    // Show help and exit
    if (opts.show_help) {
      print_help(argv[0]);
      return 0;
    }

    // Determine input source
    std::ifstream file_stream;
    std::istream* input = &std::cin;

    if (!opts.input_file.empty()) {
      if (opts.verbose) {
        std::cerr << "Reading from file: " << opts.input_file << std::endl;
      }

      file_stream.open(opts.input_file);
      if (!file_stream.is_open()) {
        throw std::runtime_error("Cannot open file: " + opts.input_file);
      }
      input = &file_stream;
    } else {
      if (opts.verbose) {
        std::cerr << "Reading from stdin..." << std::endl;
      }
    }

    // Parse and solve
    if (opts.verbose) {
      std::cerr << "Parsing CNF formula..." << std::endl;
    }

    CDCL_solver solver{*input};

    if (opts.verbose) {
      std::cerr << "Variables: " << solver.cnfSize() / 2 << std::endl;
      std::cerr << "Clauses: " << solver.cnfSize() << std::endl;
      std::cerr << "Solving..." << std::endl;
    }

    std::vector<bool> solution;
    bool is_sat = solver.solve(solution);

    // Output result
    if (is_sat) {
      std::cout << "SAT";
      if (opts.show_stats || opts.verbose) {
        std::cout << " [";
        if(solution.size() > 0)
          std::cout << (solution[0] ? "1":"0");
        for(size_t i = 1; i < solution.size(); ++i)
          std::cout << ", " << (solution[i] ? "1":"0");
        std::cout << "]";
      }
    } else {
      std::cout << "UNSAT";
    }
    std::cout << std::endl;

    // Show statistics if requested
    if (opts.show_stats) {
#ifdef DEBUG
      std::cerr << "\nStatistics:\n";
      std::cerr << "  Max CNF size: " << solver.maxCnfSize << std::endl;
      std::cerr << "  Max levels: " << solver.maxLevels << std::endl;
      std::cerr << "  Max clause size: " << solver.maxClauseSize << std::endl;
#else
      std::cerr << "Note: Statistics require debug build (make debug)\n";
#endif
    }

    if (opts.verbose) {
      std::cerr << "Done." << std::endl;
    }

    return 0;
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  catch (const std::exception& e) {
    std::cerr << "Unexpected error: " << e.what() << std::endl;
    return 1;
  }
}
