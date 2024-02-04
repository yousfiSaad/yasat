

#include "headers/CDCL_solver.h"
#include <iostream>
#include <ostream>
#include <vector>
int main(int argc, char *argv[]) {
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
