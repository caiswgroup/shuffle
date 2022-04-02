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

#include "../solvers/LstechSolver.h"
#include "../solvers/SolverFactory.h"
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
      // solvers[sid]->setBumpVar(id[sid]);
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
}


SolverInterface *
SolverFactory::createLstechSolver()
{
   int id = currentIdSolver.fetch_add(1);
   SolverInterface * solver = new LstechSolver(id);

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


void
SolverFactory::createLstechSolvers(int maxSolvers,
                                        std::vector<SolverInterface *> & solvers)
{
   for (size_t i = 0; i < maxSolvers; i++) {
      LstechSolver* lstech = (LstechSolver*) createLstechSolver();
      solvers.push_back(lstech);
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
