#include <chrono>
#include "pasolve.hpp"
#include <unistd.h>
#include <thread>
#include "simp/SimpSolver.h"
#include "core/Dimacs.h"
using namespace std; 
using namespace Minisat;
pakissat *paS;

void* seqsolve(void *arg) {
    int* pos = (int*) arg;
    int do_pos = *pos, thread_id = do_pos - 1;
    paS->assume_state[do_pos] = 1;
    while (!paS->globalEnding) {
        int do_var = paS->assume_var[do_pos];
        printf("c %d job start: assume var: %d\n", thread_id, do_var);

        vec<Lit> lits;
        paS->T[thread_id].solver = new SimpSolver();
        // while (paS->vars > paS->T[thread_id].solver->nVars()) paS->T[thread_id].solver->newVar();
        // paS->T[thread_id].solver->varGrowTo(paS->vars);
        // if (do_var) {
        //     lits.clear();
        //     lits.push(do_var > 0 ? mkLit(abs(do_var) - 1) : ~mkLit(abs(do_var) - 1));
        //     paS->T[thread_id].solver->addClause_(lits);
        // }
        // for (int i = 1; i <= paS->clauses; i++) {
        //     int l = paS->clause[i].size();
        //     lits.clear();
        //     for (int j = 0; j < l; j++) 
        //         lits.push(paS->clause[i][j] > 0 ? mkLit(abs(paS->clause[i][j]) - 1) : ~mkLit(abs(paS->clause[i][j]) - 1));
        //     paS->T[thread_id].solver->addClause_(lits);
        // }
        
        gzFile in = gzopen(paS->infile.c_str(), "rb");
        parse_DIMACS(in, *paS->T[thread_id].solver);
        assert(abs(do_var) <= paS->T[thread_id].solver->nVars());
        if (do_var) {
            lits.clear();
            lits.push(do_var > 0 ? mkLit(do_var - 1, false) : mkLit((-do_var) - 1, true));
            paS->T[thread_id].solver->addClause(lits);
        }
        printf("c %d job solving\n", thread_id);
        lits.clear();
        paS->T[thread_id].solving = true;
        lbool resM = paS->T[thread_id].solver->solveLimited(lits);
        paS->T[thread_id].solving = false;
        int res = 0;
        if (resM == l_True) res = 10;
        else if (resM == l_False) res = 20;
        printf("c %d job finish: assume var: %d, result: %d\n\n", thread_id, do_var, res);
        if (res) {
            paS->assume_state[do_pos] = res;
            int other_state = paS->assume_state[do_pos ^ 1];
            if (res == 10 && !paS->globalEnding) {
                paS->globalEnding = true;
                paS->finalResult = 10;
                paS->winner = do_var;
            }
            if (res == 20 && !paS->globalEnding && other_state == 20) {
                paS->globalEnding = true;
                paS->finalResult = 20;
                paS->winner = do_var;
            }
        }
        if (!paS->globalEnding) {
            for (int i = paS->assume_stack_head; i < paS->assume_stack_tail; i++)
                if (paS->assume_state[i] == 0) {
                    do_pos = i, paS->assume_state[i] = 1;
                    if (i > paS->assume_stack_head)
                        paS->assume_stack_head = i;
                    break; 
                }
        }
        delete paS->T[thread_id].solver;
    }
    return NULL;
}

void solve(int argc, char **argv) {
    printf("c pakissat: parallel kissat\n\n");
    auto clk_st = std::chrono::high_resolution_clock::now();    
    paS = new pakissat();
    paS->infile = argv[1];
    paS->thread_num = atoi(argv[2]);
    paS->maxtime = atoi(argv[3]);

    paS->read_var_clause_num(paS->infile.c_str());
    paS->alloc_init();
    paS->readfile(paS->infile.c_str());
    printf("c Var num: %d\nc Clause num: %d\n", paS->vars, paS->clauses);
    
    srand(paS->vars);
    paS->sample_var_random();
    printf("c Finish sample assume vars\n");

    printf("c parallel solve start\n");
    int nouse_id[paS->thread_num];
    for (int i = 0; i < paS->thread_num; i++) {
        nouse_id[i] = i + 1;
        pthread_create(&(paS->T[i].thread_ptr), NULL, seqsolve, &nouse_id[i]);
    }

    while (!paS->globalEnding) {
        usleep(100000);        
        auto clk_now = std::chrono::high_resolution_clock::now();
        int solve_time = std::chrono::duration_cast<std::chrono::seconds>(clk_now - clk_st).count();
        if (solve_time >= paS->maxtime) {
            paS->globalEnding = 1;
        }
    }

    for (int i = 0; i < paS->thread_num; i++) {
        if (paS->T[i].solving == true)
            paS->T[i].solver->interrupt();
    }
    for (int i = 0; i < paS->thread_num; i++) {
        pthread_join(paS->T[i].thread_ptr, NULL);
    }

    if (paS->finalResult == 0) printf("s UNKNOWN\n");
    else {
        if (paS->finalResult == 10) printf("s SATISFIABLE\n");
        else printf("s UNSATISFIABLE\n");
        printf("c Winner: %d\n", paS->winner);
        auto clk_now = std::chrono::high_resolution_clock::now();
        int solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(clk_now - clk_st).count();
        printf("c Use times: %d.%d\n", solve_time / 1000, solve_time % 1000);
    }

    delete paS;
    //kissat_print_witness(paS->T[0].solver, paS->vars, false);
}