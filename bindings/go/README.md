# YASAT Go Bindings

Go bindings for the YASAT (Yet Another SAT Solver) library.

## Prerequisites

1. Build the YASAT shared library:
   ```bash
   cd ../..  # Go to yasat root directory
   make lib
   ```

2. Ensure the shared library is in the library path:
   ```bash
   # Option 1: Set LD_LIBRARY_PATH (Linux)
   export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH

   # Option 2: Set DYLD_LIBRARY_PATH (macOS)
   export DYLD_LIBRARY_PATH=$PWD/build:$DYLD_LIBRARY_PATH

   # Option 3: Install system-wide
   make install
   ```

## Installation

```bash
go get github.com/yousfiSaad/yasat/bindings/go/yasat
```

Or import directly in your code:

```go
import "github.com/yousfiSaad/yasat/bindings/go/yasat"
```

## Quick Start

```go
package main

import (
    "fmt"
    "log"

    "github.com/yousfiSaad/yasat/bindings/go/yasat"
)

func main() {
    // Create solver
    solver, err := yasat.NewSolver()
    if err != nil {
        log.Fatal(err)
    }
    defer solver.Close()

    // Add clauses: (x1 ∨ x2) ∧ (¬x1 ∨ x3)
    solver.AddClause(1, 2)
    solver.AddClause(-1, 3)

    // Solve
    result, err := solver.Solve()
    if err != nil {
        log.Fatal(err)
    }

    if result == yasat.SAT {
        fmt.Println("SAT")
        // Get assignment for variable 1
        val, _ := solver.GetAssignment(1)
        fmt.Printf("x1 = %v\n", val)
    } else {
        fmt.Println("UNSAT")
    }
}
```

## API Documentation

### Creating a Solver

```go
// Create empty solver
solver, err := yasat.NewSolver()

// Create solver from file
solver, err := yasat.NewSolverFromFile("problem.cnf")

// Always close when done
defer solver.Close()
```

### Adding Clauses

```go
// Variables are 1-indexed integers
// Positive = true, Negative = false

solver.AddClause(1, 2, 3)      // (x1 ∨ x2 ∨ x3)
solver.AddClause(-1, 2)        // (¬x1 ∨ x2)
solver.AddClause(1)            // Unit clause: x1
```

### Loading from File/String

```go
// Load DIMACS CNF file
err := solver.ParseCNFFile("problem.cnf")

// Load DIMACS CNF string
cnf := `p cnf 3 2
1 2 0
-1 3 0`
err := solver.ParseCNFString(cnf)
```

### Solving

```go
result, err := solver.Solve()

switch result {
case yasat.SAT:
    fmt.Println("Satisfiable")
case yasat.UNSAT:
    fmt.Println("Unsatisfiable")
case yasat.Error:
    fmt.Println("Error:", err)
}
```

### Getting Results

```go
// Get single variable assignment (only after SAT)
val, err := solver.GetAssignment(1)  // Get x1
fmt.Printf("x1 = %v\n", val)

// Get all assignments (0-indexed slice)
assignments, err := solver.GetAllAssignments()
for i, val := range assignments {
    fmt.Printf("x%d = %v\n", i+1, val)
}
```

### Solver Information

```go
// Get number of variables
numVars, _ := solver.NumVariables()

// Get number of clauses
numClauses, _ := solver.NumClauses()

// Get library version
fmt.Println("YASAT version:", yasat.Version())
```

## Running the Example

```bash
# Build the library first
cd ../..
make lib

# Set library path
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH

# Run the example
cd bindings/go/examples
go run simple.go

# Or with a CNF file
go run simple.go ../../../tests/cnf/multi_sat.cnf
```

## Running Tests

The Go bindings include a comprehensive test suite:

```bash
# From yasat root directory
make lib
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH

# Run tests
cd bindings/go/yasat
go test -v

# Run specific test
go test -v -run TestSimpleSAT

# Run with coverage
go test -cover

# Run benchmarks
go test -bench=.
```

**Test Coverage:**
- 17 test functions covering all API methods
- Integration tests with real CNF files
- Error handling and edge cases
- Benchmarks for performance testing

## Building Your Application

When building applications that use these bindings:

1. **Development**: Set `LD_LIBRARY_PATH` to point to the build directory
2. **Production**: Install the library system-wide with `make install`
3. **Distribution**: Bundle the shared library with your application

Example build command:
```bash
# Linux
go build -o myapp main.go
LD_LIBRARY_PATH=../../build ./myapp

# With CGO_LDFLAGS to embed rpath
CGO_LDFLAGS="-Wl,-rpath,\$ORIGIN/lib" go build -o myapp main.go
```

## Troubleshooting

### "cannot find -lyasat"

The shared library hasn't been built or isn't in the library path:
```bash
cd ../..
make lib
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH
```

### "yasat.h: No such file or directory"

The CGo directive can't find the header. Make sure you're building from the correct directory or adjust the path in `yasat.go`:
```go
#cgo LDFLAGS: -L/path/to/yasat/build -lyasat
#include "/path/to/yasat/src/c_api/yasat.h"
```

### Runtime: "error while loading shared libraries"

The shared library isn't in the system library path:
```bash
# Temporary (current session)
export LD_LIBRARY_PATH=/path/to/yasat/build:$LD_LIBRARY_PATH

# Permanent (install system-wide)
cd /path/to/yasat
sudo make install
```

## Examples

See the [examples](examples/) directory for complete working examples:

- `simple.go` - Basic usage with manual clause addition and file loading

## License

Same license as YASAT project.
