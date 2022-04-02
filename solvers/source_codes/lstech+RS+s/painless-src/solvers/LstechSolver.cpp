// -----------------------------------------------------------------------------
// Copyright (C) 2017  Ludovic LE FRIOUX
//
// This file is part of PaInleSS.
//
// PaInleSS is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------

// Lstech includes
#include "utils/System.h"
#include "core/Dimacs.h"
#include "simp/SimpSolver.h"

#include "../utils/Logger.h"
#include "../utils/System.h"
#include "../utils/Parameters.h"
#include "../clauses/ClauseManager.h"
#include "../solvers/LstechSolver.h"

using namespace Minisat;

// Macros for minisat literal representation conversion
#define MINI_LIT(lit) lit > 0 ? mkLit(lit - 1, false) : mkLit((-lit) - 1, true)

#define INT_LIT(lit) sign(lit) ? -(var(lit) + 1) : (var(lit) + 1)


static void makeMiniVec(ClauseExchange * cls, vec<Lit> & mcls)
{
   for (size_t i = 0; i < cls->size; i++) {
      mcls.push(MINI_LIT(cls->lits[i]));
   }
}

void cbkLstechExportClause(void * issuer, int lbd, vec<Lit> & cls)
{
	LstechSolver* mp = (LstechSolver*)issuer;

	if (lbd > mp->lbdLimit)
		return;

	ClauseExchange * ncls = ClauseManager::allocClause(cls.size());

	for (int i = 0; i < cls.size(); i++) {
		ncls->lits[i] = INT_LIT(cls[i]);
	}

   ncls->lbd  = lbd;
   ncls->from = mp->id;

   mp->clausesToExport.addClause(ncls);
}

Lit cbkLstechImportUnit(void * issuer)
{
   LstechSolver * mp = (LstechSolver*)issuer;

   Lit l = lit_Undef;

   ClauseExchange * cls = NULL;

   if (mp->unitsToImport.getClause(&cls) == false)
      return l;

   l = MINI_LIT(cls->lits[0]);

   ClauseManager::releaseClause(cls);
   // printf("%d\n", l);
   return l;
}

bool cbkLstechImportClause(void * issuer, int * lbd, vec<Lit> & mcls)
{
   LstechSolver* mp = (LstechSolver*)issuer;

   ClauseExchange * cls = NULL;

   if (mp->clausesToImport.getClause(&cls) == false)
      return false;

   makeMiniVec(cls, mcls);

   *lbd = cls->lbd;
   
   // printf("lbd %d\n", *lbd);
   ClauseManager::releaseClause(cls);

   return true;
}

LstechSolver::LstechSolver(int id) : SolverInterface(id, MAPLE)
{
	lbdLimit = Parameters::getIntParam("lbd-limit", 2);
   solver = new SimpSolver();
	solver->cbkExportClause = cbkLstechExportClause;
	solver->cbkImportClause = cbkLstechImportClause;
	solver->cbkImportUnit   = cbkLstechImportUnit;
	solver->issuer          = this;
}

LstechSolver::~LstechSolver()
{
	delete solver;
}

void 
LstechSolver::addOriginClauses(pakissat *P) {
   solver->varGrowTo(P->vars);
   vec<Lit> lits;
   for (int i = 1; i <= P->clauses; i++) {
      int l = P->clause[i].size();
      lits.clear();
      for (int j = 0; j < l; j++) 
         lits.push(MINI_LIT(P->clause[i][j]));
      solver->addClause_(lits);
   }
}

void
LstechSolver::initshuffle(int id)
{
   if (id) solver->initshuffle = id;
   // if (id) solver->bump_one=rand()%solver->max_var+1;
   // if (id == ID_XOR) {
   //    solver->GE = true;
   // } else {
   //    solver->GE = false;
   // }

   // if (id % 2) {
   //    solver->VSIDS = false;
   // } else {
   //    solver->VSIDS = true;
   // }

   // if (id % 4 >= 2) {
   //    solver->verso = false;
   // } else {
   //    solver->verso = true;
   // }
}


bool
LstechSolver::loadFormula(const char* filename)
{
   assume_var = 0;
   gzFile in = gzopen(filename, "rb");

   parse_DIMACS(in, *solver);

   gzclose(in);

   return true;
}

//Get the number of variables of the formula
int
LstechSolver::getVariablesCount()
{
	return solver->nVars();
}

// Get a variable suitable for search splitting
int
LstechSolver::getDivisionVariable()
{
   return (rand() % getVariablesCount()) + 1;
}

// Set initial phase for a given variable
void
LstechSolver::setPhase(const int var, const bool phase)
{
	solver->setPolarity(var - 1, phase ? true : false);
}

// Bump activity for a given variable
void
LstechSolver::bumpVariableActivity(const int var, const int times)
{
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void
LstechSolver::setSolverInterrupt()
{
   stopSolver = true;

   solver->interrupt();
}

void
LstechSolver::unsetSolverInterrupt()
{
   stopSolver = false;

	solver->clearInterrupt();
}

// Diversify the solver
void
LstechSolver::diversify(int id)
{
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult
LstechSolver::solve(const vector<int> & cube)
{
   unsetSolverInterrupt();

   vector<ClauseExchange *> tmp;

   tmp.clear();
   clausesToAdd.getClauses(tmp);

   for (size_t ind = 0; ind < tmp.size(); ind++) {
      vec<Lit> mcls;
      makeMiniVec(tmp[ind], mcls);

      ClauseManager::releaseClause(tmp[ind]);

      if (solver->addClause(mcls) == false) {
         printf("c unsat when adding cls\n");
         return UNSAT;
      }
   }

   vec<Lit> miniAssumptions;
   for (size_t ind = 0; ind < cube.size(); ind++) {
      miniAssumptions.push(MINI_LIT(cube[ind]));
   }

   lbool res = solver->solveLimited(miniAssumptions);

   if (res == l_True)
      return SAT;

   if (res == l_False)
      return UNSAT;

   return UNKNOWN;
}

void
LstechSolver::addClause(ClauseExchange * clause)
{ 
   printf("wrong call kissat addClause\n");
}

void
LstechSolver::addLearnedClause(ClauseExchange * clause)
{
   if (clause->size == 1) {
      unitsToImport.addClause(clause);
   } else {
      clausesToImport.addClause(clause);
   }
}

void
LstechSolver::addClauses(const vector<ClauseExchange *> & clauses)
{  
   printf("wrong call kissat addClauses\n");
}

void
LstechSolver::addInitialClauses(const vector<ClauseExchange *> & clauses)
{
   printf("wrong call kissat addInitialClauses\n");
}

void
LstechSolver::addLearnedClauses(const vector<ClauseExchange *> & clauses)
{
   for (size_t i = 0; i < clauses.size(); i++) {
      addLearnedClause(clauses[i]);
   }
}

void
LstechSolver::getLearnedClauses(vector<ClauseExchange *> & clauses)
{
   clausesToExport.getClauses(clauses);
}

void
LstechSolver::increaseClauseProduction()
{
   // lbdLimit++;
}

void
LstechSolver::decreaseClauseProduction()
{
   if (lbdLimit > 2) {
      lbdLimit--;
   }
}

SolvingStatistics
LstechSolver::getStatistics()
{
   SolvingStatistics stats;

   stats.conflicts    = solver->conflicts;
   stats.propagations = solver->propagations;
   stats.restarts     = solver->starts;
   stats.decisions    = solver->decisions;
   // stats.memPeak      = memUsedPeak();

   return stats;
}

void
LstechSolver::solverRelease()
{
   delete solver;
}

std::vector<int>
LstechSolver::getModel()
{
   std::vector<int> model;

   for (int i = 0; i < solver->nVars(); i++) {
      if (solver->model[i] != l_Undef) {
         int lit = solver->model[i] == l_True ? i + 1 : -(i + 1);
         model.push_back(lit);
      }
   }

   return model;
}


vector<int>
LstechSolver::getFinalAnalysis()
{
   vector<int> outCls;

   for (int i = 0; i < solver->conflict.size(); i++) {
      outCls.push_back(INT_LIT(solver->conflict[i]));
   }

   return outCls;
}


vector<int>
LstechSolver::getSatAssumptions() {
   vector<int> outCls;
   // vec<Lit> lits;
   // solver->getAssumptions(lits);
   // for (int i = 0; i < lits.size(); i++) {
   //   outCls.push_back(-(INT_LIT(lits[i])));
   // }
   return outCls;
};
