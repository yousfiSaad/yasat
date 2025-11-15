

#include "headers/CDCL_solver.h"
#include <iostream>
#include <ostream>
#include <vector>
#include <exception>

int main(int argc, char *argv[]) {
  try {
    CDCL_solver solver{std::cin};

    std::vector<bool> solution;
    if (solver.solve(solution))
    {
      std::cout << "SAT [";
      if(solution.size() > 0)
        std::cout << (solution[0] ? "1":"0");
      for(int i = 1; i < solution.size(); ++i)
        std::cout << ", " << (solution[i] ? "1":"0");
      std::cout << "]";
    }
    else
      std::cout << "UNSAT";

    std::cout << std::endl;

    return 0;
  }
  catch (const std::runtime_error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  catch (const std::exception& e) {
    std::cerr << "Unexpected error: " << e.what() << std::endl;
    return 1;
  }

  (void)argc;  // Suppress unused parameter warning
  (void)argv;  // Suppress unused parameter warning
}
