#include <chrono>
#include "pasolve.hpp"
#include <unistd.h>
#include <thread>
using namespace std; 

extern "C" {
    #include "src/application.h"
    #include "src/internal.h"
    #include "src/witness.h"
}
pakissat *paS;


void* seqsolve(void *arg) {
    int* pos = (int*) arg;
    int do_pos = *pos, thread_id = do_pos - 1;
    paS->assume_state[do_pos] = 1;
    while (!paS->globalEnding) {
        int do_var = paS->assume_var[do_pos];
        printf("c %d job start: assume var: %d\n", thread_id, do_var);

        paS->T[thread_id].solver = kissat_init();
        kissat_mab_parse(paS->T[thread_id].solver);
        kissat_set_option (paS->T[thread_id].solver, "target", 2);
        kissat_reserve(paS->T[thread_id].solver, paS->vars);
        if (do_var) {
            kissat_add(paS->T[thread_id].solver, do_var);
            kissat_add(paS->T[thread_id].solver, 0);
        }
        for (int i = 1; i <= paS->clauses; i++) {
            int l = paS->clause[i].size();
            for (int j = 0; j < l; j++)
                kissat_add(paS->T[thread_id].solver, paS->clause[i][j]);
            kissat_add(paS->T[thread_id].solver, 0);
        }
        // printf("c info: %.2lf %d %d\n", paS->T[thread_id].solver->step_chb, paS->T[thread_id].solver->heuristic, paS->T[thread_id].solver->mab == true ? 1 : 0);
        printf("c %d job solving: unassigned vars %d\n", thread_id, paS->T[thread_id].solver->unassigned);
        paS->T[thread_id].solving = true;
        int res = kissat_solve(paS->T[thread_id].solver);
        paS->T[thread_id].solving = false;
        printf("c %d job finish: assume var: %d, result: %d\n\n", thread_id, do_var, res);
        
        if (res) {
            paS->assume_state[do_pos] = res;
            int other_state = paS->assume_state[do_pos ^ 1];
            if (res == 10 && !paS->globalEnding) {
                paS->globalEnding = true;
                paS->finalResult = 10;
                paS->winner = thread_id;
            }
            if (res == 20 && !paS->globalEnding && other_state == 20) {
                paS->globalEnding = true;
                paS->finalResult = 20;
                paS->winner = thread_id;
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
        // paS->T[thread_id].nconf = paS->T[thread_id].solver->nconflict;
        // paS->T[thread_id].ndec = paS->T[thread_id].solver->nassign;
        kissat_release(paS->T[thread_id].solver);
    }
    return NULL;
}

void solve(int argc, char **argv) {
    printf("c pakissat: parallel kissat_mab\n\n");
    auto clk_st = std::chrono::high_resolution_clock::now();    
    paS = new pakissat();
    paS->infile = argv[1];
    paS->thread_num = atoi(argv[2]);
    paS->maxtime = atoi(argv[3]);

    paS->read_var_clause_num(paS->infile.c_str());
    paS->alloc_init();
    paS->readfile(paS->infile.c_str()); 
    printf("c Var num: %d\nc Clause num: %d\n", paS->vars, paS->clauses);
     
    // auto clk_simp_st = std::chrono::high_resolution_clock::now();
    // paS->simplify_init();
    // paS->simplify_resolution();
    // paS->simplify_binary();
    // auto clk_simp_ed = std::chrono::high_resolution_clock::now();
    // int timest = std::chrono::duration_cast<std::chrono::milliseconds>(clk_simp_ed - clk_simp_st).count();
    // printf("c Finish Simplify, use %d.%d seconds\n", timest / 1000, timest % 1000);

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
            kissat_terminate(paS->T[i].solver);
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

    // for (int i = 0; i < paS->thread_num; i++) {
    //     printf("c thread: %d ,conflict_num: %d ,decide num: %d ,ratio is %.3lf .\n", i, paS->T[i].nconf, paS->T[i].ndec, 1.0 * paS->T[i].nconf / paS->T[i].ndec);
    // }

    delete paS;
    //kissat_print_witness(paS->T[0].solver, paS->vars, false);
}