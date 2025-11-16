/// Example demonstrating loading CNF files
///
/// This example shows how to:
/// 1. Load a CNF file in DIMACS format
/// 2. Solve it
/// 3. Display results

use yasat::{Solver, SatResult};
use std::env;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== YASAT Rust Binding - File Example ===\n");

    // Get filename from command line
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: {} <cnf-file>", args[0]);
        eprintln!("\nExample:");
        eprintln!("  {} ../../tests/cnf/simple_sat.cnf", args[0]);
        std::process::exit(1);
    }

    let filename = &args[1];
    println!("Loading CNF file: {}\n", filename);

    // Create solver from file
    let mut solver = Solver::from_file(filename)?;

    println!("Formula statistics:");
    println!("  Variables: {}", solver.num_variables());
    println!("  Clauses: {}\n", solver.num_clauses());

    // Solve
    println!("Solving...");
    let start = std::time::Instant::now();
    let result = solver.solve()?;
    let duration = start.elapsed();

    match result {
        SatResult::Sat => {
            println!("✓ SAT - Formula is satisfiable!");
            println!("Time: {:?}\n", duration);

            let assignments = solver.get_assignments()?;

            if assignments.len() <= 20 {
                // Show full solution for small problems
                println!("Solution:");
                for (i, &value) in assignments.iter().enumerate() {
                    print!("{}", if value { 1 } else { 0 });
                    if (i + 1) % 10 == 0 {
                        println!();
                    } else {
                        print!(" ");
                    }
                }
                println!();
            } else {
                // Show summary for large problems
                println!("Solution found ({} variables)", assignments.len());
                println!("First 10 assignments:");
                for (i, &value) in assignments.iter().take(10).enumerate() {
                    println!("  x{} = {}", i + 1, value);
                }
            }
        }
        SatResult::Unsat => {
            println!("✗ UNSAT - Formula is unsatisfiable!");
            println!("Time: {:?}", duration);
        }
    }

    Ok(())
}
