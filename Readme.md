# YASAT - Yet Another SAT Solver

[![CI](https://github.com/yousfiSaad/yasat/actions/workflows/ci.yml/badge.svg)](https://github.com/yousfiSaad/yasat/actions/workflows/ci.yml)
[![Docker](https://github.com/yousfiSaad/yasat/actions/workflows/docker.yml/badge.svg)](https://github.com/yousfiSaad/yasat/actions/workflows/docker.yml)
[![License](https://img.shields.io/badge/license-Proprietary-blue.svg)](license.txt)

A fast, lightweight SAT solver implementing the Conflict-Driven Clause Learning (CDCL) algorithm in modern C++.

## Features

- **CDCL Algorithm**: State-of-the-art conflict-driven clause learning
- **Zero Dependencies**: Uses only the C++ standard library
- **Modern C++**: Built with C++17 standards
- **Robust Input Validation**: Comprehensive error checking for malformed CNF files
- **Flexible I/O**: Read from files or stdin
- **Multiple Build Modes**: Optimized release and debug builds with statistics
- **Comprehensive Testing**: Full test suite with CI/CD integration

## What is SAT?

The Boolean Satisfiability Problem (SAT) asks whether there exists an assignment of true/false values to variables that makes a Boolean formula true. YASAT solves SAT problems expressed in Conjunctive Normal Form (CNF).

### CDCL Algorithm

Conflict-Driven Clause Learning (CDCL) is a modern SAT solving algorithm that:
1. Makes variable assignments (decisions)
2. Propagates implications using unit propagation
3. Analyzes conflicts to learn new clauses
4. Backtracks intelligently using learned information
5. Restarts periodically using Luby sequences

## Installation

### Prerequisites

- C++ compiler with C++17 support (g++ 7+ or clang++ 5+)
- GNU Make
- Git (for cloning the repository)

### Build from Source

```bash
# Clone the repository
git clone https://github.com/yousfiSaad/yasat.git
cd yasat

# Build release version (optimized)
make

# Or build debug version (with statistics)
make debug

# Run tests
make test

# Clean build artifacts
make clean
```

### Using Docker

Pre-built Docker images are available with YASAT and its libraries:

```bash
# Pull latest image
docker pull ghcr.io/yousfisaad/yasat:latest

# Run YASAT
docker run --rm ghcr.io/yousfisaad/yasat:latest --help

# Solve a CNF file
docker run --rm -v $(pwd):/data ghcr.io/yousfisaad/yasat:latest /data/problem.cnf
```

### Pre-built Releases

Download pre-built binaries and libraries from [GitHub Releases](https://github.com/yousfiSaad/yasat/releases):

- Linux (x86_64) - Binary + Shared Library (.so) + Static Library (.a)
- macOS (x86_64) - Binary + Shared Library (.dylib) + Static Library (.a)
- Windows (x86_64) - Binary (.exe) + DLL + Static Library (.a)

Each release includes header files for FFI integration with Go, Python, Rust, etc.

### Building Shared Libraries

Build YASAT as a shared library for use in other languages:

```bash
# Build both shared and static libraries
make lib

# Build only shared library
make shared

# Build only static library
make static

# Install libraries system-wide (requires sudo)
sudo make install
```

Libraries are installed to:
- Linux/macOS: `/usr/local/lib/libyasat.{so,a}`
- Headers: `/usr/local/include/yasat/`

## Usage

### Basic Usage

```bash
# Solve a CNF file
./build/yasat problem.cnf

# Solve from stdin
cat problem.cnf | ./build/yasat

# Show help
./build/yasat --help
```

### Command-Line Options

```
Usage: yasat [options] [input.cnf]

Options:
  -h, --help       Show help message and exit
  -v, --verbose    Enable verbose output (shows solving progress)
  -s, --stats      Show solving statistics (requires debug build)
  input.cnf        Input CNF file in DIMACS format
                   (if not provided, reads from stdin)
```

### Examples

```bash
# Solve with verbose output
./build/yasat -v data/php_2p_3h.cnf

# Show statistics (requires debug build)
make debug
./build/yasat -s data/php_3p_2h.cnf

# Combine options
./build/yasat -v -s problem.cnf
```

### Output

- **SAT**: Formula is satisfiable
- **UNSAT**: Formula is unsatisfiable

In verbose or stats mode, satisfying assignments are shown:
```
SAT [1, 0, 1, 0, 0, 1]
```
Where 1 = true, 0 = false for variables 1, 2, 3, etc.

## CNF File Format

YASAT accepts CNF files in the standard DIMACS format:

```
c This is a comment line
c Comments start with 'c'
p cnf 3 2
1 2 0
-1 3 0
```

Format explanation:
- **Comment lines**: Start with `c`
- **Problem line**: `p cnf <num_variables> <num_clauses>`
- **Clause lines**: Space-separated literals ending with `0`
  - Positive number `n`: Variable n is true
  - Negative number `-n`: Variable n is false (negated)
  - `0`: End of clause marker

Example formula: `(x₁ ∨ x₂) ∧ (¬x₁ ∨ x₃)`

## Testing

YASAT includes a comprehensive test suite:

```bash
# Run all tests
make test

# Run tests manually
./tests/run_tests.sh

# View test documentation
cat tests/README.md
```

### Test Coverage

- ✓ Simple SAT instances
- ✓ Simple UNSAT instances
- ✓ Empty formulas (edge cases)
- ✓ Multi-variable formulas
- ✓ Pigeonhole principle problems
- ✓ Input validation tests

## Development

### Build Modes

```bash
# Release build (optimized, -O3)
make release

# Debug build (symbols, statistics, -g -O0)
make debug

# Run tests
make test

# Clean
make clean

# Show help
make help
```

### Debugging

Build with debug mode to enable:
- Debug symbols for gdb/lldb
- Assertion checks
- Solving statistics (max CNF size, decision levels, clause size)

```bash
make debug
./build/yasat -s problem.cnf
```

### Code Style

- Modern C++17
- No external dependencies
- Comprehensive error handling
- Follows C++ Core Guidelines

## Performance

YASAT is designed for:
- Small to medium SAT instances (< 1000 variables)
- Educational purposes and algorithm understanding
- Baseline SAT solver implementation

For industrial-scale SAT solving, consider:
- MiniSat
- Glucose
- CaDiCaL

## CI/CD

YASAT has a comprehensive CI/CD pipeline powered by GitHub Actions:

### Continuous Integration
- **Multi-platform builds**: Linux, macOS, Windows
- **Matrix testing**: Release and debug builds on all platforms
- **Automated testing**: Full test suite runs on every push and PR
- **Code quality**: Static analysis with cppcheck
- **Build artifacts**: Pre-compiled binaries and libraries

### Automated Releases
- **Version tagging**: Push `vX.Y.Z` tags to trigger releases
- **Multi-platform packages**: Linux, macOS, Windows packages
- **Checksums**: SHA256 verification for all artifacts
- **Changelog**: Auto-generated from commits

### Docker Images
- **GitHub Container Registry**: `ghcr.io/yousfisaad/yasat`
- **Multi-tag support**: `latest`, `dev`, version tags
- **Pre-built libraries**: Shared and static libraries included

**For detailed CI/CD documentation, see [docs/CICD.md](docs/CICD.md)**

## Project Structure

```
yasat/
├── src/
│   ├── main.cpp                    # Entry point with CLI parsing
│   ├── headers/
│   │   ├── CDCL_solver.h           # Solver interface
│   │   └── macros.h                # Utility macros
│   └── implementations/
│       └── CDCL_solver.cpp         # CDCL implementation
├── tests/
│   ├── run_tests.sh                # Test runner script
│   ├── README.md                   # Test documentation
│   └── cnf/                        # Test CNF files
├── data/                           # Example CNF files
├── Makefile                        # Build system
├── build.sh                        # Build script
└── README.md                       # This file
```

## License

Copyright 2024 YOUSFI Saad. All rights reserved.

This software is provided "as is" without warranty of any kind. See [license.txt](license.txt) for full details.

## Contributing

Contributions are welcome! Areas for improvement:
- Additional heuristics (VSIDS, etc.)
- Clause deletion strategies
- Preprocessing techniques
- Performance optimizations
- Additional test cases

## Acknowledgments

- CDCL algorithm based on research by Marques-Silva and Sakallah
- Luby restart strategy from Luby et al.
- Pigeonhole test cases from CNFgen by Massimo Lauria

## References

- [CDCL SAT Solving](https://en.wikipedia.org/wiki/Conflict-driven_clause_learning)
- [DIMACS CNF Format](http://www.satcompetition.org/2009/format-benchmarks2009.html)
- [Modern SAT Solvers](https://www.cs.princeton.edu/~zkincaid/courses/fall18/readings/SATHandbook-CDCL.pdf)

## Contact

For questions, bugs, or suggestions, please open an issue on GitHub.

---

**Happy SAT Solving!** 🎯
