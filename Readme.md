# YASAT - Yet Another SAT Solver

A fast, lightweight SAT solver implementing the Conflict-Driven Clause Learning (CDCL) algorithm in modern C++.

## Features

- **CDCL Algorithm**: State-of-the-art conflict-driven clause learning
- **Shared Library**: Use YASAT from C, Go, Python, and other languages
- **Language Bindings**: Idiomatic bindings for Go and Python included
- **Zero Dependencies**: Core uses only the C++ standard library
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

## Library Usage

YASAT can be used as a shared library from C, Go, Python, and other languages that support FFI.

### Building the Shared Library

```bash
# Build the shared library
make lib

# Build both CLI and library
make all

# Install system-wide (optional)
sudo make install
```

This creates `libyasat.so` (Linux) or `libyasat.dylib` (macOS) in the `build/` directory.

### C API

The C API provides a simple interface for integrating YASAT into C/C++ projects.

**Example:**

```c
#include <yasat.h>
#include <stdio.h>

int main() {
    // Create solver
    yasat_solver* solver = yasat_solver_create();

    // Add clauses: (x1 ∨ x2) ∧ (¬x1 ∨ x3)
    int clause1[] = {1, 2};
    yasat_add_clause(solver, clause1, 2);

    int clause2[] = {-1, 3};
    yasat_add_clause(solver, clause2, 2);

    // Solve
    yasat_result result = yasat_solve(solver);

    if (result == YASAT_RESULT_SAT) {
        printf("SAT\n");
        // Get assignments
        for (int i = 1; i <= 3; i++) {
            int val = yasat_get_assignment(solver, i);
            printf("x%d = %d\n", i, val);
        }
    } else {
        printf("UNSAT\n");
    }

    // Cleanup
    yasat_solver_destroy(solver);
    return 0;
}
```

**Compile:**

```bash
gcc -o myapp myapp.c -L./build -lyasat
LD_LIBRARY_PATH=./build ./myapp
```

**API Documentation:** See [`src/c_api/yasat.h`](src/c_api/yasat.h) for the complete API reference.

### Go Bindings

Idiomatic Go bindings using CGo.

**Example:**

```go
package main

import (
    "fmt"
    "github.com/yousfiSaad/yasat/bindings/go/yasat"
)

func main() {
    solver, _ := yasat.NewSolver()
    defer solver.Close()

    // Add clauses
    solver.AddClause(1, 2)
    solver.AddClause(-1, 3)

    // Solve
    result, _ := solver.Solve()

    if result == yasat.SAT {
        fmt.Println("SAT")
        assignments, _ := solver.GetAllAssignments()
        for i, val := range assignments {
            fmt.Printf("x%d = %v\n", i+1, val)
        }
    } else {
        fmt.Println("UNSAT")
    }
}
```

**Setup:**

```bash
# Build library
make lib

# Set library path
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH

# Run example
cd bindings/go/examples
go run simple.go
```

**Documentation:** See [`bindings/go/README.md`](bindings/go/README.md)

### Python Bindings

Pythonic interface using ctypes.

**Example:**

```python
from yasat import Solver, Result

# Create solver
solver = Solver()

# Add clauses
solver.add_clause(1, 2)
solver.add_clause(-1, 3)

# Solve
result = solver.solve()

if result == Result.SAT:
    print("SAT")
    assignments = solver.get_all_assignments()
    for i, val in enumerate(assignments, start=1):
        print(f"x{i} = {val}")
else:
    print("UNSAT")
```

**Setup:**

```bash
# Build library
make lib

# Install Python package
cd bindings/python
pip install -e .

# Run example
python examples/simple.py
```

**Documentation:** See [`bindings/python/README.md`](bindings/python/README.md)

### Loading from Files

All APIs support loading CNF files:

**C:**
```c
yasat_parse_cnf_file(solver, "problem.cnf");
```

**Go:**
```go
solver, _ := yasat.NewSolverFromFile("problem.cnf")
```

**Python:**
```python
solver.parse_cnf_file("problem.cnf")
# Or use convenience function
result, assignments = solve_file("problem.cnf")
```

### Installation

Install YASAT system-wide to use from any project:

```bash
# Install to /usr/local (requires sudo)
sudo make install

# Or install to user directory
make install INSTALL_PREFIX=~/.local

# Uninstall
sudo make uninstall
```

This installs:
- Library: `/usr/local/lib/libyasat.so`
- Header: `/usr/local/include/yasat.h`
- Binary: `/usr/local/bin/yasat`

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

GitHub Actions automatically:
- Builds both release and debug versions
- Runs full test suite
- Reports build status

## Project Structure

```
yasat/
├── src/
│   ├── main.cpp                    # Entry point with CLI parsing
│   ├── headers/
│   │   ├── CDCL_solver.h           # Solver interface
│   │   └── macros.h                # Utility macros
│   ├── implementations/
│   │   └── CDCL_solver.cpp         # CDCL implementation
│   └── c_api/
│       ├── yasat.h                 # C API header
│       └── yasat.cpp               # C API implementation
├── bindings/
│   ├── go/
│   │   ├── yasat/                  # Go package
│   │   ├── examples/               # Go examples
│   │   └── README.md               # Go documentation
│   └── python/
│       ├── yasat/                  # Python package
│       ├── examples/               # Python examples
│       └── README.md               # Python documentation
├── tests/
│   ├── run_tests.sh                # Test runner script
│   ├── test_c_api.c                # C API tests
│   ├── README.md                   # Test documentation
│   └── cnf/                        # Test CNF files
├── docs/
│   └── FFI_CALLBACKS.md            # FFI callback documentation
├── build/                          # Build artifacts
│   ├── yasat                       # CLI binary
│   └── libyasat.so                 # Shared library
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
