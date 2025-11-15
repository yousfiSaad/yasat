"""
YASAT Solver Python Bindings

This module provides Python bindings for the YASAT SAT solver using ctypes.
"""

import ctypes
import os
import platform
from enum import IntEnum
from pathlib import Path
from typing import List, Optional, Tuple


class Result(IntEnum):
    """Result of solving a CNF formula."""
    SAT = 10       # Formula is satisfiable
    UNSAT = 20     # Formula is unsatisfiable
    UNKNOWN = 0    # Result unknown (not yet solved)
    ERROR = -1     # Error occurred during solving


class YasatError(Exception):
    """Exception raised for YASAT solver errors."""
    pass


def _find_library():
    """Find the YASAT shared library."""
    # Try common locations
    search_paths = [
        # Relative to this file (development setup)
        Path(__file__).parent.parent.parent.parent / "build",
        # System installation
        Path("/usr/local/lib"),
        Path("/usr/lib"),
        # User installation
        Path.home() / ".local" / "lib",
    ]

    # Determine library name based on platform
    system = platform.system()
    if system == "Linux":
        lib_name = "libyasat.so"
    elif system == "Darwin":
        lib_name = "libyasat.dylib"
    elif system == "Windows":
        lib_name = "yasat.dll"
    else:
        lib_name = "libyasat.so"

    # Search for the library
    for path in search_paths:
        lib_path = path / lib_name
        if lib_path.exists():
            return str(lib_path)

    # Try letting ctypes find it in the system path
    try:
        return ctypes.util.find_library("yasat")
    except:
        pass

    raise YasatError(
        f"Could not find YASAT shared library ({lib_name}). "
        f"Please build it with 'make lib' or install it with 'make install'."
    )


# Load the library
_lib_path = _find_library()
_lib = ctypes.CDLL(_lib_path)

# Define C types
_c_solver_p = ctypes.c_void_p
_c_int = ctypes.c_int
_c_char_p = ctypes.c_char_p

# Define function signatures
_lib.yasat_solver_create.argtypes = []
_lib.yasat_solver_create.restype = _c_solver_p

_lib.yasat_solver_destroy.argtypes = [_c_solver_p]
_lib.yasat_solver_destroy.restype = None

_lib.yasat_add_clause.argtypes = [_c_solver_p, ctypes.POINTER(_c_int), _c_int]
_lib.yasat_add_clause.restype = _c_int

_lib.yasat_parse_cnf_file.argtypes = [_c_solver_p, _c_char_p]
_lib.yasat_parse_cnf_file.restype = _c_int

_lib.yasat_parse_cnf_string.argtypes = [_c_solver_p, _c_char_p]
_lib.yasat_parse_cnf_string.restype = _c_int

_lib.yasat_solve.argtypes = [_c_solver_p]
_lib.yasat_solve.restype = _c_int

_lib.yasat_get_assignment.argtypes = [_c_solver_p, _c_int]
_lib.yasat_get_assignment.restype = _c_int

_lib.yasat_get_all_assignments.argtypes = [_c_solver_p, ctypes.POINTER(_c_int), _c_int]
_lib.yasat_get_all_assignments.restype = _c_int

_lib.yasat_get_num_variables.argtypes = [_c_solver_p]
_lib.yasat_get_num_variables.restype = _c_int

_lib.yasat_get_num_clauses.argtypes = [_c_solver_p]
_lib.yasat_get_num_clauses.restype = _c_int

_lib.yasat_get_error_message.argtypes = [_c_solver_p]
_lib.yasat_get_error_message.restype = _c_char_p

_lib.yasat_get_version.argtypes = []
_lib.yasat_get_version.restype = _c_char_p


class Solver:
    """
    YASAT SAT Solver.

    A CDCL-based SAT solver for boolean satisfiability problems.

    Example:
        >>> solver = Solver()
        >>> solver.add_clause(1, -2, 3)  # (x1 ∨ ¬x2 ∨ x3)
        >>> solver.add_clause(-1, 2)     # (¬x1 ∨ x2)
        >>> result = solver.solve()
        >>> if result == Result.SAT:
        ...     print("SAT")
        ...     assignments = solver.get_all_assignments()
        ...     print(f"Solution: {assignments}")
    """

    def __init__(self):
        """Create a new YASAT solver instance."""
        self._handle = _lib.yasat_solver_create()
        if not self._handle:
            raise YasatError("Failed to create solver")

    def __del__(self):
        """Destroy the solver and free resources."""
        if hasattr(self, '_handle') and self._handle:
            _lib.yasat_solver_destroy(self._handle)
            self._handle = None

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.__del__()
        return False

    def add_clause(self, *literals: int) -> None:
        """
        Add a clause to the solver.

        Args:
            *literals: Variable literals (positive for x_i, negative for ¬x_i)

        Example:
            >>> solver.add_clause(1, -2, 3)  # (x1 ∨ ¬x2 ∨ x3)

        Raises:
            YasatError: If adding the clause fails
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        if not literals:
            return  # Empty clause is valid (though trivially unsatisfiable)

        # Convert to C array
        c_literals = (ctypes.c_int * len(literals))(*literals)
        result = _lib.yasat_add_clause(self._handle, c_literals, len(literals))

        if result != 0:  # YASAT_OK = 0
            self._raise_error("Failed to add clause")

    def parse_cnf_file(self, filename: str) -> None:
        """
        Load a CNF formula from a DIMACS format file.

        Args:
            filename: Path to the CNF file

        Raises:
            YasatError: If loading the file fails
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        c_filename = filename.encode('utf-8')
        result = _lib.yasat_parse_cnf_file(self._handle, c_filename)

        if result != 0:  # YASAT_OK = 0
            self._raise_error(f"Failed to parse CNF file: {filename}")

    def parse_cnf_string(self, cnf: str) -> None:
        """
        Load a CNF formula from a DIMACS format string.

        Args:
            cnf: CNF formula in DIMACS format

        Example:
            >>> cnf = '''p cnf 3 2
            ... 1 2 0
            ... -1 3 0'''
            >>> solver.parse_cnf_string(cnf)

        Raises:
            YasatError: If parsing fails
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        c_cnf = cnf.encode('utf-8')
        result = _lib.yasat_parse_cnf_string(self._handle, c_cnf)

        if result != 0:  # YASAT_OK = 0
            self._raise_error("Failed to parse CNF string")

    def solve(self) -> Result:
        """
        Solve the CNF formula.

        Returns:
            Result.SAT if satisfiable, Result.UNSAT if unsatisfiable

        Raises:
            YasatError: If solving fails
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        result = _lib.yasat_solve(self._handle)

        if result == Result.SAT:
            return Result.SAT
        elif result == Result.UNSAT:
            return Result.UNSAT
        elif result == Result.ERROR:
            self._raise_error("Solving failed")
        else:
            return Result.UNKNOWN

    def get_assignment(self, variable: int) -> Optional[bool]:
        """
        Get the assignment for a specific variable.

        Args:
            variable: Variable number (1-indexed)

        Returns:
            True if variable is true, False if false, None if unassigned/error

        Note:
            Only valid after solve() returns Result.SAT
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        result = _lib.yasat_get_assignment(self._handle, variable)

        if result < 0:
            return None
        return result == 1

    def get_all_assignments(self) -> List[bool]:
        """
        Get assignments for all variables.

        Returns:
            List of boolean values (0-indexed: result[0] is variable 1)

        Raises:
            YasatError: If retrieval fails

        Note:
            Only valid after solve() returns Result.SAT
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        num_vars = self.num_variables()
        if num_vars == 0:
            return []

        # Allocate array for results
        c_assignments = (ctypes.c_int * num_vars)()
        result = _lib.yasat_get_all_assignments(self._handle, c_assignments, num_vars)

        if result != 0:  # YASAT_OK = 0
            self._raise_error("Failed to get assignments")

        # Convert to Python list of bools
        return [bool(c_assignments[i]) for i in range(num_vars)]

    def num_variables(self) -> int:
        """
        Get the number of variables in the formula.

        Returns:
            Number of variables

        Raises:
            YasatError: If retrieval fails
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        result = _lib.yasat_get_num_variables(self._handle)
        if result < 0:
            raise YasatError("Failed to get number of variables")

        return result

    def num_clauses(self) -> int:
        """
        Get the number of clauses in the formula.

        Returns:
            Number of clauses

        Raises:
            YasatError: If retrieval fails
        """
        if not self._handle:
            raise YasatError("Solver is closed")

        result = _lib.yasat_get_num_clauses(self._handle)
        if result < 0:
            raise YasatError("Failed to get number of clauses")

        return result

    @staticmethod
    def version() -> str:
        """
        Get the YASAT library version.

        Returns:
            Version string (e.g., "1.0.0")
        """
        c_version = _lib.yasat_get_version()
        return c_version.decode('utf-8')

    def _raise_error(self, default_msg: str) -> None:
        """
        Raise a YasatError with the last error message from the library.

        Args:
            default_msg: Default message if no error message available
        """
        c_msg = _lib.yasat_get_error_message(self._handle)
        if c_msg:
            msg = c_msg.decode('utf-8')
            raise YasatError(msg)
        else:
            raise YasatError(default_msg)


# Convenience function
def solve_file(filename: str) -> Tuple[Result, Optional[List[bool]]]:
    """
    Convenience function to solve a CNF file.

    Args:
        filename: Path to DIMACS CNF file

    Returns:
        Tuple of (result, assignments)
        - result: Result.SAT or Result.UNSAT
        - assignments: List of boolean values if SAT, None if UNSAT

    Example:
        >>> result, solution = solve_file("problem.cnf")
        >>> if result == Result.SAT:
        ...     print("SAT:", solution)
        ... else:
        ...     print("UNSAT")
    """
    with Solver() as solver:
        solver.parse_cnf_file(filename)
        result = solver.solve()

        if result == Result.SAT:
            return result, solver.get_all_assignments()
        else:
            return result, None
