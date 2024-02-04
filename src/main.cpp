

#include "headers/CDCL_solver.h"
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>
// #include <string>
int main(int argc, char *argv[]) {
 
  // std::ifstream fin{"./data/php_10p_9h.cnf"};
  // std::ifstream fin{"./data/php_3p_2h.cnf"};

  // std::string s;
  // while (fin >> s) {
  //   std::cout << s << "\n";
  // }

  CDCL_solver solver{std::cin};
  // CDCL_solver solver{fin};

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
