# YASAT Python Bindings

Python bindings for the YASAT (Yet Another SAT Solver) library.

## Prerequisites

1. Build the YASAT shared library:
   ```bash
   cd ../..  # Go to yasat root directory
   make lib
   ```

2. The Python bindings will automatically search for the library in:
   - `../../build` (development setup)
   - `/usr/local/lib` (system installation)
   - `/usr/lib`
   - `~/.local/lib` (user installation)

   Alternatively, install system-wide:
   ```bash
   make install
   ```

## Installation

### From source (development)

```bash
# Install in development mode
pip install -e .
```

### System installation

```bash
pip install .
```

## Quick Start

```python
from yasat import Solver, Result

# Create solver
solver = Solver()

# Add clauses: (x1 ∨ x2) ∧ (¬x1 ∨ x3)
solver.add_clause(1, 2)
solver.add_clause(-1, 3)

# Solve
result = solver.solve()

if result == Result.SAT:
    print("SAT")
    # Get assignment for variable 1
    val = solver.get_assignment(1)
    print(f"x1 = {val}")
else:
    print("UNSAT")
```

## API Documentation

### Creating a Solver

```python
# Create empty solver
solver = Solver()

# Using context manager (recommended)
with Solver() as solver:
    solver.add_clause(1, 2)
    result = solver.solve()

# Load from file (convenience function)
from yasat import solve_file
result, assignments = solve_file("problem.cnf")
```

### Adding Clauses

Variables are 1-indexed integers. Positive = true, negative = false.

```python
solver.add_clause(1, 2, 3)      # (x1 ∨ x2 ∨ x3)
solver.add_clause(-1, 2)        # (¬x1 ∨ x2)
solver.add_clause(1)            # Unit clause: x1
```

### Loading from File/String

```python
# Load DIMACS CNF file
solver.parse_cnf_file("problem.cnf")

# Load DIMACS CNF string
cnf = """p cnf 3 2
1 2 0
-1 3 0"""
solver.parse_cnf_string(cnf)
```

### Solving

```python
result = solver.solve()

if result == Result.SAT:
    print("Satisfiable")
elif result == Result.UNSAT:
    print("Unsatisfiable")
```

### Getting Results

```python
# Get single variable assignment (only after SAT)
val = solver.get_assignment(1)  # Get x1 (returns bool or None)

# Get all assignments (returns list, 0-indexed)
assignments = solver.get_all_assignments()
for i, val in enumerate(assignments, start=1):
    print(f"x{i} = {val}")
```

### Solver Information

```python
# Get number of variables
num_vars = solver.num_variables()

# Get number of clauses
num_clauses = solver.num_clauses()

# Get library version
print(f"YASAT version: {Solver.version()}")
```

### Error Handling

```python
from yasat import YasatError

try:
    solver = Solver()
    solver.parse_cnf_file("nonexistent.cnf")
except YasatError as e:
    print(f"Error: {e}")
```

## Running the Example

```bash
# Make sure library is built
cd ../..
make lib

# Run the example
cd bindings/python/examples
python simple.py

# Or with a CNF file
python simple.py ../../../tests/cnf/multi_sat.cnf
```

## Running Tests

The Python bindings include a comprehensive test suite:

```bash
# From yasat root directory
make lib
export LD_LIBRARY_PATH=$PWD/build:$LD_LIBRARY_PATH

# Run tests with unittest
cd bindings/python/tests
python test_yasat.py -v

# Or with pytest (if installed)
pytest test_yasat.py -v

# Run specific test class
python test_yasat.py TestSimpleSAT -v

# Run with coverage (requires pytest-cov)
pytest test_yasat.py --cov=yasat --cov-report=html
```

**Test Coverage:**
- 25 test methods across 9 test classes
- All API methods tested
- Integration tests with real CNF files
- Error handling and edge cases
- Context manager behavior
- Convenience functions

## Type Hints

The bindings include full type hints for better IDE support:

```python
from yasat import Solver, Result
from typing import List, Optional

def solve_problem(filename: str) -> Optional[List[bool]]:
    with Solver() as solver:
        solver.parse_cnf_file(filename)
        result: Result = solver.solve()

        if result == Result.SAT:
            assignments: List[bool] = solver.get_all_assignments()
            return assignments
        return None
```

## Examples

The [examples](examples/) directory contains:

- `simple.py` - Comprehensive examples demonstrating all features:
  - Manual clause addition
  - Context manager usage
  - Loading from DIMACS strings
  - Loading from files
  - Error handling

Run with:
```bash
python examples/simple.py
python examples/simple.py path/to/problem.cnf
```

## Troubleshooting

### "Could not find YASAT shared library"

The shared library hasn't been built or isn't in a searchable location:

```bash
# Build the library
cd ../..
make lib

# Or install system-wide
sudo make install
```

### "cannot open shared object file"

On Linux, you may need to update the library cache:

```bash
# After system install
sudo ldconfig

# Or set LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/yasat/build:$LD_LIBRARY_PATH
```

### Import Error

Make sure you've installed the Python package:

```bash
pip install -e .
```

Or add the bindings to your Python path:

```python
import sys
sys.path.insert(0, '/path/to/yasat/bindings/python')
```

## Performance Notes

- The bindings use ctypes, which has some overhead compared to native C
- For performance-critical applications, consider:
  - Batching multiple clause additions
  - Using `parse_cnf_file()` for large formulas instead of adding clauses individually
  - Reusing solver instances when possible

## Requirements

- Python 3.7+
- YASAT shared library (built with `make lib`)

## License

Same license as YASAT project.
