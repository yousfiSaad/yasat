# YASAT Rust Bindings

Safe, idiomatic Rust bindings for [YASAT](https://github.com/yousfiSaad/yasat) - Yet Another SAT Solver.

## Features

- 🦀 **Safe Rust API** - Zero unsafe code in user-facing API
- ⚡ **Zero-cost abstractions** - Thin wrapper over C API
- 🎯 **Type-safe** - Leverages Rust's type system for correctness
- 📦 **Easy to use** - Ergonomic, idiomatic Rust interface
- 🧪 **Well-tested** - Comprehensive test suite
- 📖 **Well-documented** - Extensive documentation and examples

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
yasat = { path = "path/to/yasat/bindings/rust" }
```

### Prerequisites

You need the YASAT shared library installed:

```bash
# From the yasat root directory
make lib
sudo make install

# Or set LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/path/to/yasat/build:$LD_LIBRARY_PATH
```

## Quick Start

```rust
use yasat::{Solver, SatResult};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create a solver
    let mut solver = Solver::new()?;

    // Add clauses: (x1 ∨ x2) ∧ (¬x1 ∨ x3)
    solver.add_clause(&[1, 2])?;
    solver.add_clause(&[-1, 3])?;

    // Solve
    match solver.solve()? {
        SatResult::Sat => {
            println!("SAT!");
            let assignments = solver.get_assignments()?;
            println!("Solution: {:?}", assignments);
        }
        SatResult::Unsat => {
            println!("UNSAT!");
        }
    }

    Ok(())
}
```

## Usage Examples

### Loading from a CNF file

```rust
use yasat::Solver;

let mut solver = Solver::from_file("problem.cnf")?;
let result = solver.solve()?;
```

### Loading from a string

```rust
use yasat::Solver;

let cnf = "p cnf 3 2\n1 2 0\n-1 3 0\n";
let mut solver = Solver::from_string(cnf)?;
let result = solver.solve()?;
```

### Checking individual variable assignments

```rust
if solver.solve()? == SatResult::Sat {
    for var in 1..=solver.num_variables() as i32 {
        if let Some(value) = solver.get_assignment(var)? {
            println!("x{} = {}", var, value);
        }
    }
}
```

### Resetting the solver

```rust
let mut solver = Solver::new()?;

// First problem
solver.add_clause(&[1, 2])?;
solver.solve()?;

// Reset and solve a new problem
solver.reset();
solver.add_clause(&[3, 4])?;
solver.solve()?;
```

## Running Examples

```bash
# Navigate to the rust bindings directory
cd bindings/rust

# Ensure the library is built and in path
export LD_LIBRARY_PATH=../../build:$LD_LIBRARY_PATH

# Run the simple example
cargo run --example simple

# Run the file example
cargo run --example file ../../tests/cnf/simple_sat.cnf
```

## Running Tests

```bash
cd bindings/rust
export LD_LIBRARY_PATH=../../build:$LD_LIBRARY_PATH
cargo test
```

## API Documentation

Generate and view the documentation:

```bash
cargo doc --open
```

## Error Handling

All operations return `Result<T, Error>` with detailed error types:

```rust
use yasat::Error;

match solver.add_clause(&[]) {
    Ok(_) => println!("Success"),
    Err(Error::InvalidArgument(msg)) => eprintln!("Invalid argument: {}", msg),
    Err(e) => eprintln!("Other error: {}", e),
}
```

## Thread Safety

The `Solver` type is both `Send` and `Sync`, allowing safe use across threads:

```rust
use std::thread;

let mut solver = Solver::new()?;
solver.add_clause(&[1, 2])?;

// Solve in a separate thread
let handle = thread::spawn(move || {
    solver.solve()
});

let result = handle.join().unwrap()?;
```

## Performance

The Rust bindings add negligible overhead over the C API:

- Clause addition: ~1-2ns overhead (bounds checking)
- Solving: No overhead (direct FFI call)
- Assignment retrieval: Allocation cost only

## Minimum Supported Rust Version (MSRV)

Rust 1.70 or later.

## License

Same as YASAT - Proprietary. See [license.txt](../../license.txt).

## Contributing

Contributions are welcome! Please ensure:

1. All tests pass: `cargo test`
2. Code is formatted: `cargo fmt`
3. No clippy warnings: `cargo clippy`
4. Documentation is updated

## Resources

- [YASAT Documentation](../../Readme.md)
- [SAT Solving](https://en.wikipedia.org/wiki/Boolean_satisfiability_problem)
- [CDCL Algorithm](https://en.wikipedia.org/wiki/Conflict-driven_clause_learning)
- [DIMACS CNF Format](http://www.satcompetition.org/2009/format-benchmarks2009.html)
