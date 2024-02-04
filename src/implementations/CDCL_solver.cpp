/*
 * Copyright 2024 YOUSFI Saad
 * 
 * This software is protected by copyright laws and international treaties. Unauthorized reproduction or distribution of this software, or any portion of it, may result in severe civil and criminal penalties, and will be prosecuted to the maximum extent possible under the law.
 * 
 * This software is provided "as is" without warranty of any kind, either expressed or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. The entire risk as to the quality and performance of the software is with you. Should the software prove defective, you assume the cost of all necessary servicing, repair, or correction.
 * 
 * YOUSFI Saad reserves the right to make changes to this software at any time without notice.
*/

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <istream>
#include <iterator>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "../headers/CDCL_solver.h"
#include "../headers/macros.h"

// #define LEN_NEEDSRESTART 40
#define MAX_RESTART 800
#define COEFS_NEEDSRESTART 100

decision::decision(const literal &lit, const int &level)
    : lit{lit}, level{level} {}

void CDCL_solver::clearLevels() {
  while (undoLevelAndPop()) {
  }
}
literal CDCL_solver::propagate_(literal nlitetal, int level, idx cause) {
  std::queue<literal> literals_to_propagate;
  literals_to_propagate.push(nlitetal);

  if (isLiteralEmpty(nlitetal))
    causes[nlitetal] = cause;

  while (!literals_to_propagate.empty()) {
    const literal current_literal = literals_to_propagate.front();
    literals_to_propagate.pop();

    if (isLiteralTrue(current_literal))
      continue;
    if (isLiteralFalse(current_literal))
      return current_literal;

    setLiteral(current_literal, level);
    const vector_of_indexes &clausesIdxes_same = index[current_literal];
    for (int i = 0; i < clausesIdxes_same.size(); ++i) {
      const idx clIdx = clausesIdxes_same[i];
      if (!cnf[clIdx].sat) {
        if (level > 0) {
          if (states_satisfied_in_each_level.size() < level + 1)
            states_satisfied_in_each_level.resize(level + 1);
          states_satisfied_in_each_level[level].push(clIdx);
        }
        cnf[clIdx].sat = true;
      }
    }

    const vector_of_indexes &clausesIdxes_compl = index[NOT(current_literal)];
    for (int i = 0; i < clausesIdxes_compl.size(); ++i) {
      const idx clIdx = clausesIdxes_compl[i]; // *rIt;
      const clause &cl = cnf[clIdx];
      if (cl.sat)
        continue; // no conclusions
      int number_of_literals_removed = 0;
      literal remainingLiteral;
      for (const literal &literal_in_clause : cl.literals) {
        if (isLiteralFalse(literal_in_clause))
          number_of_literals_removed++;
        else // et remaining with xor
          remainingLiteral = literal_in_clause;
      }
      if (number_of_literals_removed == cl.literals.size() - 1) {
        if (isLiteralEmpty(remainingLiteral))
          causes[remainingLiteral] = clIdx;
        literals_to_propagate.push(remainingLiteral);
      }
    }
  }
  return -1;
}
// assuming the size is enough
void CDCL_solver::setLiteral(literal lit, int level) {
  values[lit] = 1;
  values[NOT(lit)] = -1;
  levels[IDX(lit)] = level;
  number_of_variables_set++;
  if (level > 0) {
    if (literals_set_in_each_level.size() < level + 1)
      literals_set_in_each_level.resize(level + 1);
    literals_set_in_each_level[level].push(IDX(lit));
  }
}
void CDCL_solver::addClause(const clause &newClause) {
  for (literal nLit : newClause.literals)
    index[nLit].push_back(cnf.size());
  cnf.push_back(newClause);
#ifdef DEBUG
  maxCnfSize = std::max(maxCnfSize, cnf.size());
  maxClauseSize = std::max(maxClauseSize, newClause.literals.size());
#endif
}
bool CDCL_solver::solve(std::vector<bool> &output) { return solve(&output); }
void CDCL_solver::extractSolution(std::vector<bool> &output) {
  output.resize(number_of_variables);
  FOR(i, number_of_variables)
  output[i] = isLiteralTrue(VAL(i));
}
bool CDCL_solver::solve(std::vector<bool> *pOutput) {
  if (!checkCnfSat())
    return false;

  bool solved = solve_();

  if (!solved)
    return false;

  if (pOutput != nullptr)
    extractSolution(*pOutput);

  clearLevels();
  return true;
}

inline bool CDCL_solver::needsRestart(unsigned long iterations,
                                      int n_prev_restarts) {
  return iterations > (luby(n_prev_restarts + 1) * COEFS_NEEDSRESTART);
}

bool CDCL_solver::solve_(literal startLit, unsigned long *pNum_iter) {
  int r = 0;
  // char seed = 0;
  bool valid = true;
  literal conflictingLiteral;
  if (!_setUpLiteralToStartWith(startLit, conflictingLiteral, valid))
    return false;
  unsigned long n_iterations;
  if (pNum_iter == nullptr)
    pNum_iter = &n_iterations;
  (*pNum_iter) = 0;
  while (true) {
    (*pNum_iter)++;
    if (valid) {
      // literal lit = chooseUnsetLiteral2(seed);
      // literal lit = chooseUnsetLiteral2();
      literal lit = findFirstUnsetLiteral();
      if (lit == -1)
        return true;
      int level = dStack.empty() ? 1 : (this->level() + 1);
      dStack.push(decision{NOT(lit), level});
      dStack.push(decision{lit, level});
    } else {
      if (needsRestart(*pNum_iter, r)) {
        // std::cerr<<"Restart " << r << "\n";
        clearLevels();
        if (r > R && *pNum_iter > 100000)
          throw std::runtime_error("takes long time");
        valid = true;
        if (!_setUpLiteralToStartWith(startLit, conflictingLiteral, valid))
          return false;
        r++;
        (*pNum_iter) = 0;
        continue;
      }
      analysisResult analysis = analyze(conflictingLiteral);
      if (analysis.newClause.literals.size() == 1) {
        clearLevels();
        conflictingLiteral =
            propagate_(analysis.newClause.literals.back(), 0, -1);
        valid = conflictingLiteral == -1;
        if (valid) {
          int pruned = pruneCnf();
          if (!_setUpLiteralToStartWith(startLit, conflictingLiteral, valid))
            return false;

          continue;
        }
      } else {
        decision d = dStack.top();
        if (level() == analysis.maxLevel &&
            analysis.maxLevel > analysis.targetLevel) {
          // TODO: review
          while (!dStack.empty() && level() > analysis.targetLevel)
            undoLevelAndPop();
          addClause(analysis.newClause);
          d.level = level() + 1;
          d.lit = NOT(d.lit);
          conflictingLiteral = propagate_(d.lit, level(), cnf.size() - 1);
          valid = conflictingLiteral == -1;
          if (valid)
            dStack.push(d);
          else
            continue;
        } else
          backTrack();
      }
    }

    if (dStack.empty())
      return false;

    conflictingLiteral = propagate_(dStack.top().lit, level(), -1);
    valid = conflictingLiteral == -1;
  }
}
bool CDCL_solver::_setUpLiteralToStartWith(const literal &startLit,
                                           literal &conflictingLiteral,
                                           bool &valid) {
  if (startLit != -1) {
    conflictingLiteral = propagate_(startLit, 1, -1);
    dStack.push(decision{startLit, 1});
    valid = conflictingLiteral == -1;
    if (!valid)
      return false;
  }
  return true;
}
int CDCL_solver::level() { return dStack.empty() ? 0 : dStack.top().level; }
void CDCL_solver::backTrack() {
  bool backToPreviousLevel = false;
  do {
    int lvl = undoLevelAndPop();
    backToPreviousLevel = !dStack.empty() && lvl != level();
  } while (backToPreviousLevel);
}
int CDCL_solver::undoLevelAndPop() {
  if (dStack.empty())
    return 0;
  int level = this->level();
  if (literals_set_in_each_level.size() > level)
    while (!literals_set_in_each_level[level].empty()) {
      idx indx = literals_set_in_each_level[level].top();
      literals_set_in_each_level[level].pop();
      literal lit = VAL(indx);
      values[lit] = values[NOT(lit)] = 0;
      number_of_variables_set--;
    }
  if (states_satisfied_in_each_level.size() > level)
    while (!states_satisfied_in_each_level[level].empty()) {
      idx clIdx = states_satisfied_in_each_level[level].top();
      states_satisfied_in_each_level[level].pop();
      cnf[clIdx].sat = false;
      /* code */
    }
  dStack.pop();
#ifdef DEBUG
  maxLevels = std::max(maxLevels, level);
#endif
  return level;
}
literal CDCL_solver::findFirstUnsetLiteral(literal startFrom) {
  idx nextIdxToTry = 0;
  if (startFrom != -1)
    nextIdxToTry = IDX(startFrom);
  while (nextIdxToTry < number_of_variables &&
         !isLiteralEmpty(VAL(nextIdxToTry)))
    nextIdxToTry++;

  if (nextIdxToTry == number_of_variables)
    return -1;

  return VAL(nextIdxToTry);
}
literal CDCL_solver::chooseUnsetLiteral2() {
  std::array<unsigned int, 81> cells;
  FOR(ij, 81) {
    cells[ij] = 0;
    FOR(e, 9) {
      literal lit = VAL(ij * 9 + e);
      cells[ij] += isLiteralEmpty(lit) ? 1 : 0;
    }
  }
  FOR(ij, 81) {
    cells[ij] <<= 8;
    cells[ij] |= ij;
  }
  std::sort(std::begin(cells), std::end(cells));
  std::array<unsigned int, 81>::iterator it =
      std::lower_bound(std::begin(cells), std::end(cells), 90);
  if (it == std::end(cells))
    return -1;

  char n = (*it) >> 8;
  char n_max = rand() % n;
  char i = 0;
  literal lit = VAL(9 * ((*it) & 0b11111111));
  while (!isLiteralEmpty(lit) || i < n_max) {
    if (isLiteralEmpty(lit))
      i++;
    lit += 2;
  }
  return lit;
}
literal CDCL_solver::chooseUnsetLiteral() {
  int min = __SHRT_MAX__;
  literal retLit = -1;
  // for(int i = 80; i >= 0; i--)
  FOR(i, 81) {
    if (i >= cnf.size())
      return findFirstUnsetLiteral(-1);
    if (cnf[i].sat)
      continue;
    int c = 0;
    literal lit = -1;
    for (literal ilit : cnf[i].literals)
      if (isLiteralEmpty(ilit)) {
        c++;
        lit = ilit;
      }

    if (c == 2 && lit != -1)
      return lit;

    if (c < min && lit != -1) {
      min = c;
      retLit = lit;
    }
  }

  return findFirstUnsetLiteral(-1);
}
analysisResult CDCL_solver::analyze(literal conflictingLiteral) {
  analysisResult res;
  literal nConflictinLiteral = NOT(conflictingLiteral);
  bool changed = true;

  ll lits{conflictingLiteral, nConflictinLiteral};
  std::unordered_set<literal> s{conflictingLiteral, nConflictinLiteral};
  while (changed) {
    changed = false;
    ll::iterator it = lits.begin();
    while (it != lits.end()) {
      ll::iterator bIt = it;
      bIt++;
      if (causes[*it] >= 0) {
        for (literal cLit : cnf[causes[*it]].literals) {
          if (*it == cLit)
            continue;
          literal ncLit = NOT(cLit);
          if (s.find(ncLit) == s.end() && levels[IDX(ncLit)] > 0) {
            s.insert(ncLit);
            // for(const literal prlit: cnf[causes[*it]].literals)
            //   std::cerr << prlit << " ";
            // std::cerr << " : "
            //   << ncLit << " -> " << *it << "\n";
            lits.push_front(ncLit);
            changed = true;
          }
        }
        lits.erase(it);
      }
      it = bIt;
    }
  }
  std::set<int> ss;
  for (literal lit : lits) {
    int lvl = levels[IDX(lit)];
    ss.insert(lvl);
    literal nLit = NOT(lit);
    res.newClause.literals.push_back(nLit);
  }
  std::set<int>::reverse_iterator rIt = ss.rbegin();
  std::set<int>::reverse_iterator rIt2 = ss.rbegin();
  rIt++;
  res.targetLevel = *rIt;
  res.maxLevel = *rIt2;
  return res;
}
size_t CDCL_solver::cnfSize() { return cnf.size(); }
size_t CDCL_solver::satClausesCount() {
  size_t count = 0;
  for (const clause &cl : cnf) {
    count += cl.sat ? 1 : 0;
  }
  return count;
}
size_t CDCL_solver::pruneCnf() {
  // return 0;
  FORS(i, index.size())
  index[i].clear();
  vector_of_clauses newCnf;
  // for (clause &cl : cnf)
  for (int i = 0; i < cnf.size(); i++) {
    clause &cl = cnf[i];
    cl.literals.erase(
        std::remove_if(std::begin(cl.literals), std::end(cl.literals),
                       [&](const literal &lit) { return values[lit] == -1; }),
        std::end(cl.literals));
    if (cl.sat)
      continue;
    newCnf.push_back(cl);
    for (const literal &lit : cl.literals) {
      index[lit].push_back(newCnf.size() - 1);
    }
  }
  FOR(lit, 2 * number_of_variables)
  causes[lit] = -1;
  size_t pruned = cnf.size() - newCnf.size();
  cnf = newCnf;
  return pruned;
}
std::ostream &CDCL_solver::printCnf(std::ostream &out) {
  std::stringstream ss;
  std::set<idx> literals;
  std::map<int, int> literalsmap;
  int clNumber = 0;
  for (const clause &cl : cnf) {
    bool containsEmptyLiteral = false;
    for (const literal &lit : cl.literals) {
      if (isLiteralEmpty(lit)) {
        containsEmptyLiteral = true;
        idx ind = IDX(lit);
        literals.insert(ind);
        literalsmap.try_emplace(ind, literals.size());
        ss << (SIGN(lit) ? "-" : "");
        ss << (literalsmap[IDX(lit)]) << " ";
      }
    }
    if (containsEmptyLiteral) {
      clNumber++;
      // ss << "("<<clNumber<<") ";
      ss << "0\n";
    }
  }
  for (auto &pa : literalsmap) {
    out << "c " << pa.first << " : " << pa.second << "\n";
  }
  out << "p cnf " << literals.size() << " " << clNumber << "\n";
  out << ss.str();
  return out;
}
bool CDCL_solver::checkCnfSat() {
  assert(level() == 0);
  for (const clause &cl : cnf) {
    if (cl.sat)
      continue;
    bool sat = false;
    for (const literal &lit : cl.literals) {
      if (!isLiteralFalse(lit)) {
        sat = true;
        break;
      }
    }
    if (!sat)
      return false;
  }
  return true;
}
void CDCL_solver::shrink_to_fit() {
  pruneCnf();
  cnf.shrink_to_fit();
  FOR(i, index.size())
  index[i].shrink_to_fit();
}
CDCL_solver::CDCL_solver(std::istream &in)
    : number_of_variables_set{0}, R{MAX_RESTART}, number_of_variables{0}
#ifdef DEBUG
      ,
      maxLevels{0}, maxCnfSize{0}, maxClauseSize{0},
#endif
{
  parseCnf(in);
}
unsigned CDCL_solver::luby(unsigned i) {
  for (unsigned k = 1; k < 32; k++)
    if (i == (1 << k) - 1)
      return 1 << (k - 1);
  for (unsigned k = 1;; k++)
    if (1 << (k - 1) <= i && (i < ((1 << k) - 1)))
      return luby(i - (1 << (k - 1)) + 1);
}
vector_of_clauses CDCL_solver::getCleanedCnf() const {
  vector_of_clauses clauses;
  std::copy_if(cnf.begin(), cnf.end(), std::back_inserter(clauses),
               [](const clause &cl) { return !cl.sat; });

  for (clause &cl : clauses)
    cl.literals.erase(
        std::remove_if(cl.literals.begin(), cl.literals.end(),
                       [&](literal lit) { return !isLiteralEmpty(lit); }),
        cl.literals.end());

  clauses.erase(
      std::remove_if(clauses.begin(), clauses.end(),
                     [](const clause &cl) { return cl.literals.empty(); }),
      clauses.end());

  return clauses;
}
void CDCL_solver::parseCnf(std::istream &in) {
  std::string sline;
  int n_clauses = 0;
  while (std::getline(in, sline)) {
    if (sline[0] == 'c')
      continue;
    std::stringstream ss{sline};
    if (sline[0] == 'p') {
      ss >> sline >> sline;
      ss >> number_of_variables >> n_clauses;
      values.resize(2 * number_of_variables);
      index.resize(2 * number_of_variables);
      causes.resize(2 * number_of_variables);
      levels.resize(2 * number_of_variables);
      break;
    }
  }
  FOR(i, n_clauses) {
    clause cl{in};
    addClause(cl);
  }
}
clause::clause(std::istream &in) : sat{false} {
  int ilit;
  while (in >> ilit) {
    if (ilit == 0)
      return;
    literal lit;
    if (ilit < 0) {
      ilit++;
      lit = CDCL_solver::NOT(CDCL_solver::VAL(-ilit));
    } else {
      ilit--;
      lit = CDCL_solver::VAL(ilit);
    }
    literals.push_back(lit);
  }
}
