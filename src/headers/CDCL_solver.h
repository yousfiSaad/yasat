#ifndef CDCL_SOLVER_h
#define CDCL_SOLVER_h

#include <istream>
#include <map>
#include <set>
#include <list>
#include <strstream>
#include <vector>
#include <stack>
#include <cstddef>
#include <iostream>
#include <memory>

// #define NUMBER_OF_VARIABLES 9 * 9 * 9

typedef char state;
typedef std::vector<state> vector_of_states;
typedef int literal;
typedef int idx;
typedef std::vector<literal> vector_of_literals;
typedef std::vector<idx> vector_of_indexes;
typedef std::vector<vector_of_indexes> vector_vector_idx;
typedef std::list<literal> ll;
typedef std::vector<std::stack<idx>> vec_stack_idx ;

struct clause
{
  vector_of_literals literals;
  bool sat;
  clause() : sat{false} {}
  clause(std::istream& in);
};
typedef std::vector<clause> vector_of_clauses;

struct decision
{
  literal lit;
  int level;
  decision(const literal &lit, const int &level);
};
typedef std::stack<decision> stack_of_decisions ;

struct analysisResult
{
  clause newClause;
  int targetLevel;
  int maxLevel;
};

class CDCL_solver
{
  int number_of_variables;
  stack_of_decisions  dStack;
  vector_of_clauses   cnf;
  vector_of_states    values;
  vector_vector_idx   index;
  vector_of_indexes   causes;
  std::vector<char>   levels;
  int               number_of_variables_set;
  vec_stack_idx       literals_set_in_each_level;
  vec_stack_idx       states_satisfied_in_each_level;

  static unsigned luby(unsigned i);
  static inline bool needsRestart(
      unsigned long iterations, int n_prev_restarts);
  bool _setUpLiteralToStartWith(
      const literal &startLit,
      literal &conflictingLiteral,
      bool &valid);
  void clearLevels();
  bool checkCnfSat();
  bool solve_(literal lit = -1,
      unsigned long* pOutN_iterations = nullptr);
  size_t pruneCnf();
  int level();
  void backTrack();
  void addClause(const clause &newClause);
  analysisResult analyze(literal conflictingLiteral);
  bool solved();
  literal findFirstUnsetLiteral(literal startFrom = -1);
  literal chooseUnsetLiteral();
  // literal chooseUnsetLiteral2(char seed = 0);
  literal chooseUnsetLiteral2();
  int undoLevelAndPop();
  literal propagate_(literal lit,
      int level,idx cause); // returns the conflict source (to be used in the analysis)
  bool propagate(literal lit, int level, idx cause);
  bool isLiteralTrue(literal lit) const;
  bool isLiteralFalse(literal lit) const;
  bool isLiteralEmpty(literal) const;
  void setLiteral(literal lit, int level);
  int R;
  void parseCnf(std::istream& in);
  void extractSolution(std::vector<bool>& output);

public://static
  static literal NOT(literal a);
  static literal VAL(idx a);
  static idx IDX(literal a);
  static char SIGN(literal a);

public:
  vector_of_clauses getCleanedCnf() const;
  CDCL_solver(std::istream& in);
  CDCL_solver(CDCL_solver &&) noexcept = default;
  CDCL_solver(const CDCL_solver &) = default;
  CDCL_solver &operator=(CDCL_solver &&) noexcept = default;
  CDCL_solver &operator=(const CDCL_solver &) = default;
  // bool addConstraints(const std::vector<literal>& literals);
  void shrink_to_fit();
  std::ostream &printCnf(std::ostream &out);
  size_t cnfSize();
  size_t satClausesCount();
  bool solve(std::vector<bool>& output);
  bool solve(std::vector<bool>* pOutput);
#ifdef DEBUG 
  std::size_t maxCnfSize;
  int maxLevels;
  std::size_t maxClauseSize;
#endif
};

inline bool CDCL_solver::solved()
{
  return number_of_variables_set == number_of_variables;
}
inline bool CDCL_solver::isLiteralEmpty(literal lit) const
{
  return values[lit] == 0;
}
inline bool CDCL_solver::isLiteralFalse(literal lit) const
{
  return (values[lit] == -1);
}
inline bool CDCL_solver::isLiteralTrue(literal lit) const
{
  return (values[lit] == 1);
}
inline state CDCL_solver::SIGN(literal a)
{
  return a & 1;
}
inline idx CDCL_solver::IDX(literal a)
{
  return a >> 1;
}
inline literal CDCL_solver::VAL(idx a)
{
  return a << 1;
}
inline literal CDCL_solver::NOT(literal a)
{
  return a ^ 1;
}
inline void invFlatIndices(const idx& indx,
    char& i, char&j, char&e)
{
  e = indx%9;
  j = indx/9%9;
  i = indx/81;
}
#endif
