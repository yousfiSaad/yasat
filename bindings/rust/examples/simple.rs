/// Simple example demonstrating basic YASAT usage
///
/// This example shows how to:
/// 1. Create a solver
/// 2. Add clauses programmatically
/// 3. Solve the formula
/// 4. Retrieve the solution

use yasat::{Solver, SatResult};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== YASAT Rust Binding - Simple Example ===\n");

    // Create a new solver instance
    let mut solver = Solver::new()?;

    println!("Building formula: (x1 ∨ x2) ∧ (¬x1 ∨ x3) ∧ (¬x2 ∨ ¬x3)\n");

    // Add clauses
    // Clause 1: (x1 ∨ x2)
    solver.add_clause(&[1, 2])?;

    // Clause 2: (¬x1 ∨ x3)
    solver.add_clause(&[-1, 3])?;

    // Clause 3: (¬x2 ∨ ¬x3)
    solver.add_clause(&[-2, -3])?;

    println!("Formula statistics:");
    println!("  Variables: {}", solver.num_variables());
    println!("  Clauses: {}\n", solver.num_clauses());

    // Solve the formula
    println!("Solving...");
    let result = solver.solve()?;

    match result {
        SatResult::Sat => {
            println!("✓ SAT - Formula is satisfiable!\n");

            // Get the satisfying assignment
            let assignments = solver.get_assignments()?;

            println!("Solution:");
            for (i, &value) in assignments.iter().enumerate() {
                println!("  x{} = {}", i + 1, value);
            }

            // Verify individual assignments
            println!("\nVerifying solution:");
            for var in 1..=solver.num_variables() as i32 {
                if let Ok(Some(value)) = solver.get_assignment(var) {
                    println!("  Variable {} is {}", var, if value { "true" } else { "false" });
                }
            }
        }
        SatResult::Unsat => {
            println!("✗ UNSAT - Formula is unsatisfiable!");
        }
    }

    println!("\nYASAT version: {}", Solver::version()?);

    Ok(())
}
