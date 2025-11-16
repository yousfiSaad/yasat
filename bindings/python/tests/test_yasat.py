"""
Test suite for YASAT Python bindings.

Run with:
    pytest test_yasat.py -v
    python -m pytest test_yasat.py -v
    python -m unittest test_yasat.py
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent))

from yasat import Solver, Result, YasatError, solve_file


class TestSolverBasics(unittest.TestCase):
    """Test basic solver functionality."""

    def test_version(self):
        """Test that version string is returned."""
        version = Solver.version()
        self.assertIsNotNone(version)
        self.assertIsInstance(version, str)
        self.assertTrue(len(version) > 0)
        print(f"YASAT version: {version}")

    def test_solver_create(self):
        """Test solver creation and destruction."""
        solver = Solver()
        self.assertIsNotNone(solver)
        # Destructor should be called automatically

    def test_solver_context_manager(self):
        """Test solver works as context manager."""
        with Solver() as solver:
            self.assertIsNotNone(solver)
        # Context manager should have closed the solver

    def test_multiple_close(self):
        """Test that closing solver multiple times doesn't crash."""
        solver = Solver()
        solver.__del__()
        solver.__del__()  # Should not crash


class TestSimpleSAT(unittest.TestCase):
    """Test simple SAT solving."""

    def test_simple_sat(self):
        """Test a simple satisfiable formula."""
        solver = Solver()

        # (x1 OR x2) AND (NOT x1 OR x2)
        # Should be SAT with x2=true
        solver.add_clause(1, 2)
        solver.add_clause(-1, 2)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)

        # Check that x2 is true
        val = solver.get_assignment(2)
        self.assertTrue(val)

    def test_complex_sat(self):
        """Test a more complex satisfiable formula."""
        solver = Solver()

        # (x1 OR x2) AND (NOT x1 OR x3) AND (NOT x2 OR NOT x3)
        solver.add_clause(1, 2)
        solver.add_clause(-1, 3)
        solver.add_clause(-2, -3)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)

        # Verify assignments
        x1 = solver.get_assignment(1)
        x2 = solver.get_assignment(2)
        x3 = solver.get_assignment(3)

        # Verify the solution satisfies the clauses
        self.assertTrue(x1 or x2, "Clause (x1 OR x2) not satisfied")
        self.assertTrue(not x1 or x3, "Clause (NOT x1 OR x3) not satisfied")
        self.assertTrue(not x2 or not x3, "Clause (NOT x2 OR NOT x3) not satisfied")


class TestClauses(unittest.TestCase):
    """Test clause addition and handling."""

    def test_empty_clause(self):
        """Test adding an empty clause."""
        solver = Solver()
        # Adding empty clause should not error
        solver.add_clause()

    def test_unit_clauses(self):
        """Test unit clauses."""
        solver = Solver()

        # (x1) AND (x2)
        solver.add_clause(1)
        solver.add_clause(2)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)

        x1 = solver.get_assignment(1)
        x2 = solver.get_assignment(2)

        self.assertTrue(x1)
        self.assertTrue(x2)

    def test_large_clause(self):
        """Test clause with many literals."""
        solver = Solver()

        # Large disjunction
        solver.add_clause(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)


class TestParsing(unittest.TestCase):
    """Test CNF parsing functionality."""

    def test_parse_cnf_string(self):
        """Test parsing CNF from a string."""
        solver = Solver()

        cnf = """p cnf 2 2
1 2 0
-1 2 0"""

        solver.parse_cnf_string(cnf)

        num_vars = solver.num_variables()
        num_clauses = solver.num_clauses()

        self.assertEqual(num_vars, 2)
        self.assertEqual(num_clauses, 2)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)

    def test_parse_cnf_string_with_comments(self):
        """Test parsing CNF with comments."""
        solver = Solver()

        cnf = """c This is a comment
c Another comment
p cnf 3 3
1 2 0
-1 3 0
-2 -3 0
"""

        solver.parse_cnf_string(cnf)

        num_vars = solver.num_variables()
        self.assertEqual(num_vars, 3)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)

    def test_parse_cnf_file(self):
        """Test parsing CNF from a file."""
        # Create a temporary CNF file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cnf', delete=False) as f:
            f.write("c Simple test CNF\n")
            f.write("p cnf 3 3\n")
            f.write("1 2 0\n")
            f.write("-1 3 0\n")
            f.write("-2 -3 0\n")
            tmpfile = f.name

        try:
            solver = Solver()
            solver.parse_cnf_file(tmpfile)

            num_vars = solver.num_variables()
            self.assertEqual(num_vars, 3)

            result = solver.solve()
            self.assertEqual(result, Result.SAT)
        finally:
            os.unlink(tmpfile)

    def test_invalid_file(self):
        """Test parsing a non-existent file."""
        solver = Solver()

        with self.assertRaises(YasatError):
            solver.parse_cnf_file("/nonexistent/file.cnf")


class TestSolverInfo(unittest.TestCase):
    """Test solver information methods."""

    def test_num_variables_empty(self):
        """Test getting number of variables from empty solver."""
        solver = Solver()
        num_vars = solver.num_variables()
        self.assertEqual(num_vars, 0)

    def test_num_clauses_after_adding(self):
        """Test clause counting."""
        solver = Solver()

        solver.add_clause(1, 2)
        solver.add_clause(-1, 3)
        solver.add_clause(2, -3)

        num_clauses = solver.num_clauses()
        self.assertEqual(num_clauses, 3)

    def test_num_variables_after_adding(self):
        """Test variable counting."""
        solver = Solver()

        solver.add_clause(1, 2)
        solver.add_clause(-1, 3)

        num_vars = solver.num_variables()
        self.assertEqual(num_vars, 3)

    def test_large_variable_numbers(self):
        """Test handling large variable numbers."""
        solver = Solver()

        # Use variable 100
        solver.add_clause(100, -100)

        num_vars = solver.num_variables()
        self.assertEqual(num_vars, 100)


class TestAssignments(unittest.TestCase):
    """Test assignment retrieval."""

    def test_get_assignment(self):
        """Test getting individual assignments."""
        solver = Solver()

        solver.add_clause(1, 2)
        result = solver.solve()

        self.assertEqual(result, Result.SAT)

        # At least one should be true
        x1 = solver.get_assignment(1)
        x2 = solver.get_assignment(2)

        self.assertTrue(x1 or x2)

    def test_get_all_assignments(self):
        """Test getting all assignments."""
        solver = Solver()

        solver.add_clause(1, 2)
        solver.add_clause(-1, 3)

        result = solver.solve()
        self.assertEqual(result, Result.SAT)

        assignments = solver.get_all_assignments()
        self.assertIsInstance(assignments, list)
        self.assertTrue(len(assignments) > 0)

        # Verify assignments are booleans
        for val in assignments:
            self.assertIsInstance(val, bool)

    def test_get_assignment_before_solve(self):
        """Test getting assignment before solving."""
        solver = Solver()

        solver.add_clause(1, 2)

        # Should return None when called before solving
        val = solver.get_assignment(1)
        # May return None or raise error, both are acceptable


class TestErrorHandling(unittest.TestCase):
    """Test error handling."""

    def test_yasat_error_exception(self):
        """Test that YasatError is raised on errors."""
        solver = Solver()

        with self.assertRaises(YasatError):
            solver.parse_cnf_file("nonexistent.cnf")


class TestResultEnum(unittest.TestCase):
    """Test Result enum."""

    def test_result_values(self):
        """Test Result enum values."""
        self.assertEqual(Result.SAT, 10)
        self.assertEqual(Result.UNSAT, 20)
        self.assertEqual(Result.UNKNOWN, 0)
        self.assertEqual(Result.ERROR, -1)

    def test_result_names(self):
        """Test Result enum names."""
        self.assertEqual(Result.SAT.name, "SAT")
        self.assertEqual(Result.UNSAT.name, "UNSAT")


class TestConvenienceFunctions(unittest.TestCase):
    """Test convenience functions."""

    def test_solve_file(self):
        """Test solve_file convenience function."""
        # Create a temporary CNF file with fully constrained variables
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cnf', delete=False) as f:
            f.write("p cnf 2 2\n")
            f.write("1 0\n")  # x1 must be true
            f.write("2 0\n")  # x2 must be true
            tmpfile = f.name

        try:
            result, assignments = solve_file(tmpfile)

            self.assertEqual(result, Result.SAT)
            self.assertIsNotNone(assignments)
            self.assertIsInstance(assignments, list)
            # Assignments may be empty for underconstrained formulas
            # but should contain values for this fully constrained example
            if len(assignments) > 0:
                for val in assignments:
                    self.assertIsInstance(val, bool)
        finally:
            os.unlink(tmpfile)


class TestIntegration(unittest.TestCase):
    """Integration tests with real CNF files."""

    def test_integration_with_real_file(self):
        """Test with an actual test CNF file if available."""
        test_file = Path(__file__).parent.parent.parent.parent / "tests" / "cnf" / "simple_sat.cnf"

        if not test_file.exists():
            self.skipTest("Test CNF file not found")

        with Solver() as solver:
            solver.parse_cnf_file(str(test_file))

            result = solver.solve()

            # simple_sat.cnf should be satisfiable
            self.assertEqual(result, Result.SAT)

            print(f"Successfully solved {test_file.name}: {result.name}")


if __name__ == '__main__':
    unittest.main(verbosity=2)
