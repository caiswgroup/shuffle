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
    while (!paS->globalEnding) {
        
        paS->T[thread_id].solver = kissat_init();
        
        
        if (thread_id == 13 || thread_id == 14 || thread_id == 20 || thread_id == 21)
            kissat_set_option (paS->T[thread_id].solver, "tier1", 3);
        else
            kissat_set_option (paS->T[thread_id].solver, "tier1", 2);

        if (thread_id == 3 || thread_id == 4 || thread_id == 6 || thread_id == 8 || thread_id == 11 || thread_id == 12 || thread_id == 13 || thread_id == 14 || thread_id == 16 || thread_id == 18 || thread_id == 23)
            kissat_set_option (paS->T[thread_id].solver, "chrono", 0);
        else
            kissat_set_option (paS->T[thread_id].solver, "chrono", 1);    

        if (thread_id == 2 || thread_id == 23)
            kissat_set_option (paS->T[thread_id].solver, "stable", 0);
        else if (thread_id == 6 || thread_id == 16)
            kissat_set_option (paS->T[thread_id].solver, "stable", 2);
        else
            kissat_set_option (paS->T[thread_id].solver, "stable", 1);

        if (thread_id == 10 || thread_id == 22)
            kissat_set_option (paS->T[thread_id].solver, "walkinitially", 1);
        else
            kissat_set_option (paS->T[thread_id].solver, "walkinitially", 0);

        if (thread_id == 7 || thread_id == 8 || thread_id == 9 || thread_id == 17 || thread_id == 18 || thread_id == 19 || thread_id == 20)
            kissat_set_option (paS->T[thread_id].solver, "target", 0);
        else if (thread_id == 0 || thread_id == 2 || thread_id == 3 || thread_id == 4 || thread_id == 5 || thread_id == 6 || thread_id == 10 || thread_id == 23)
            kissat_set_option (paS->T[thread_id].solver, "target", 1);
        else
            kissat_set_option (paS->T[thread_id].solver, "target", 2);

        
        if (thread_id == 4 || thread_id == 5 || thread_id == 8 || thread_id == 9 || thread_id == 12 || thread_id == 13 || thread_id == 15 || thread_id == 18 || thread_id == 19)
            kissat_set_option (paS->T[thread_id].solver, "phase", 0);
        else
            kissat_set_option (paS->T[thread_id].solver, "phase", 1);

        kissat_reserve(paS->T[thread_id].solver, paS->vars);
        for (int i = 1; i <= paS->clauses; i++) {
            int l = paS->clause[i].size();
            for (int j = 0; j < l; j++)
                kissat_add(paS->T[thread_id].solver, paS->clause[i][j]);
            kissat_add(paS->T[thread_id].solver, 0);
        }
        paS->T[thread_id].solving = true;
        
        int res = kissat_solve(paS->T[thread_id].solver);
        paS->T[thread_id].solving = false;
        printf("c %d job finish: result: %d\n\n", thread_id, res);
        
        if (res) {
            if (res == 10 && !paS->globalEnding) {
                paS->globalEnding = true;
                paS->finalResult = 10;
                paS->winner = thread_id;
            }
            if (res == 20 && !paS->globalEnding) {
                paS->globalEnding = true;
                paS->finalResult = 20;
                paS->winner = thread_id;
            }
        }
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