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

#include "../solvers/MapleCOMSPSSolver.h"
#include "../solvers/KissatBonus.h"
#include "../solvers/SolverFactory.h"
#include "../solvers/Reducer.h"
#include "../utils/Parameters.h"
#include "../utils/System.h"

void
SolverFactory::sparseRandomDiversification(
      const std::vector<SolverInterface *> & solvers)
{
   if (solvers.size() == 0)
      return;

   int vars = solvers[0]->getVariablesCount();

   // The first solver of the group (1 LRB/1 VSIDS) keeps polarity = false for all vars
   for (int sid = 1; sid < solvers.size(); sid++) {
      //srand(sid);
      for (int var = 1; var <= vars; var++) {
         if (rand() % solvers.size() == 0) {
            solvers[sid]->setPhase(var, rand() % 2 == 1);
         }
      }
   }
}

void
SolverFactory::nativeDiversification(const std::vector<SolverInterface *> & solvers)
{
   int vars = solvers[0]->getVariablesCount();
   int *id = new int[vars + 1];
   id[0]=-1;
   for (int i = 0; i <= vars; i++) id[i] = i;
   random_shuffle(id + 1, id + vars + 1);
   for (int sid = 0; sid < solvers.size(); sid++) {
      solvers[sid]->setBumpVar(id[sid]);
      solvers[sid]->diversify(sid);
   }
}

void
SolverFactory::initshuffleDiversification(const std::vector<SolverInterface *> & solvers)
{
   for (int sid = 0; sid < solvers.size(); sid++) {
      solvers[sid]->initshuffle(sid);
   }
}

void
SolverFactory::parameterDiversification(const std::vector<SolverInterface *> & solvers)
{
   for (int sid = 0; sid < solvers.size(); sid++) {
      parameter p;
      p.tier1 = 2;
      p.chrono = 1;
      p.stable = 1;
      p.walkinitially = 0;
      p.target = 2;
      p.phase = 1;
      p.heuristic = 0;
      p.margin = 10;

      if (sid == 13 || sid == 14 || sid == 20 || sid == 21)
         p.tier1 = 3;

      if (sid == 3 || sid == 4 || sid == 6 || sid == 8 || sid == 11 || sid == 12 || sid == 13 || sid == 14 || sid == 16 || sid == 18 || sid == 23)
         p.chrono = 0;   

      if (sid == 2 || sid == 23)
         p.stable = 0; 
      else if (sid == 6 || sid == 16)
         p.stable = 2; 

      if (sid == 10 || sid == 22)
         p.walkinitially = 1;

      if (sid == 7 || sid == 8 || sid == 9)
         p.target = 0;
      else if (sid == 0 || sid == 2 || sid == 3 || sid == 4 || sid == 5 || sid == 6 || sid == 10 || sid == 23)
         p.target = 1;
      
      if (sid == 4 || sid == 5 || sid == 8 || sid == 9 || sid == 12 || sid == 13 || sid == 15 || sid == 18 || sid == 19)
         p.phase = 0;
     
      printf("c id=%d\n", sid);
      solvers[sid]->setParameter(p);
   }
}

SolverInterface *
SolverFactory::createMapleCOMSPSSolver()
{
   int id = currentIdSolver.fetch_add(1);

   SolverInterface * solver = new MapleCOMSPSSolver(id);

   //solver->loadFormula(Parameters::getFilename());

   return solver;
}

SolverInterface *
SolverFactory::createKissatBonusSolver()
{
   int id = currentIdSolver.fetch_add(1);
   SolverInterface * solver = new KissatBonusSolver(id);

   //solver->loadFormula(Parameters::getFilename());

   return solver;
}

// SolverInterface *
// SolverFactory::createMapleChronoBTSolver()
// {
//    int id = currentIdSolver.fetch_add(1);

//    SolverInterface * solver = new MapleChronoBTSolver(id);

//    //solver->loadFormula(Parameters::getFilename());

//    return solver;
// }

SolverInterface *
SolverFactory::createReducerSolver(SolverInterface* _solver)
{
   int id = currentIdSolver.fetch_add(1);

   SolverInterface * solver = new Reducer(id, _solver);

   return solver;
}

void
SolverFactory::createMapleCOMSPSSolvers(int maxSolvers,
                                        std::vector<SolverInterface *> & solvers)
{
   solvers.push_back(createMapleCOMSPSSolver());

   double memoryUsed    = getMemoryUsed();
   int maxMemorySolvers = Parameters::getIntParam("max-memory", 240) * 1024 *
                          1024 / memoryUsed;

   if (maxSolvers > maxMemorySolvers) {
      maxSolvers = maxMemorySolvers;
   }

   for (int i = 1; i < maxSolvers; i++) {
      solvers.push_back(cloneSolver(solvers[0]));
   }
}

void
SolverFactory::createKissatBonusSolvers(int maxSolvers,
                                        std::vector<SolverInterface *> & solvers)
{
   for (size_t i = 0; i < maxSolvers; i++) {
      KissatBonusSolver* kissat = (KissatBonusSolver*) createKissatBonusSolver();
      solvers.push_back(kissat);
   }
}


// void
// SolverFactory::createMapleChronoBTSolvers(int nbSolvers,
//                            std::vector<SolverInterface *> & solvers)
// {
//    for (size_t i = 0; i < nbSolvers; i++) {
//       MapleChronoBTSolver* maple = (MapleChronoBTSolver*) createMapleChronoBTSolver();
//       solvers.push_back(maple);
//    }
// }

SolverInterface *
SolverFactory::cloneSolver(SolverInterface * other)
{
   SolverInterface * solver;
   
   int id = currentIdSolver.fetch_add(1);

   switch(other->type) {
      case MAPLE :
	      solver = new MapleCOMSPSSolver((MapleCOMSPSSolver &) *other, id);
      	break;

      default :
         return NULL;
   }

   return solver;
}

void
SolverFactory::printStats(const std::vector<SolverInterface *> & solvers)
{
   printf("c | ID | conflicts  | propagations |  restarts  | decisions  " \
          "| memPeak |\n");

   for (size_t i = 0; i < solvers.size(); i++) {
      SolvingStatistics stats = solvers[i]->getStatistics();

      printf("c | %2zu | %10ld | %12ld | %10ld | %10ld | %7d |\n",
             solvers[i]->id, stats.conflicts, stats.propagations,
             stats.restarts, stats.decisions, (int)stats.memPeak);
   }
}
