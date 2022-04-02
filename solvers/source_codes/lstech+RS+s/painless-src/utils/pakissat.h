#ifndef _pakissat_h_INCLUDED
#define _pakissat_h_INCLUDED
#include <atomic>
#include <iostream>
#include "../simplify/pamap.h"
using namespace std;

struct pakissat {
public:
    
    pakissat();

    int vars;
    int clauses;
    int thread_num;
    svec<int> *clause;

    string infile;
    svec<int> assume_var;
    int assume_stack_head, assume_stack_tail;
    atomic<int> *assume_state; // for assume state: 0 is free, 1 is working, 10 is sat, 20 is unsat


    //int *model, *mapto, *mapsign, *clause_state, *resol_var, *thread_time;
    atomic<bool> assumeEnding;

    void alloc_init();
    void sample_var_random();
    void sample_var_frequency();
    void readfile(const char *file);
    void read_var_clause_num(const char *file);
    void print_cnf(const char *file);
};


#endif