"""
YASAT Python Bindings

Python interface to the YASAT (Yet Another SAT Solver) library.
Uses ctypes to wrap the C API and provide a Pythonic interface.

Basic usage:
    >>> from yasat import Solver, Result
    >>> solver = Solver()
    >>> solver.add_clause(1, -2, 3)  # (x1 ∨ ¬x2 ∨ x3)
    >>> solver.add_clause(-1, 2)     # (¬x1 ∨ x2)
    >>> result = solver.solve()
    >>> if result == Result.SAT:
    ...     print("SAT")
    ...     print(f"x1 = {solver.get_assignment(1)}")
    ... else:
    ...     print("UNSAT")
"""

from .solver import Solver, Result, YasatError, solve_file

__version__ = "1.0.0"
__all__ = ["Solver", "Result", "YasatError", "solve_file"]
