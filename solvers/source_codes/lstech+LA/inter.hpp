#include "vect.hpp"
#include <atomic>
#include <iostream>
using namespace std;

namespace Minisat {
	class SimpSolver;
	class Lit;
	template<class T> class vec;
}

struct thread_data {
    Minisat::SimpSolver *solver;
    pthread_t thread_ptr;
    atomic<bool> solving;
};

struct pakissat {
public:
    
    pakissat();

    int vars;
    int clauses;
    int thread_num;
    int simplify_thread_num;
    int orivars;
    int oriclauses;
    vect<int> *clause;
    string infile;
    string simplify_file;
    string simplify_path;

    vect<int> assume_var;
    int assume_stack_head, assume_stack_tail;
    atomic<int> *assume_state; // for assume state: 0 is free, 1 is working, 10 is sat, 20 is unsat

    atomic<int> model_pid; 
    thread_data *T;

    //int *model, *mapto, *mapsign, *clause_state, *resol_var, *thread_time;
    int finalResult;
    int winner;
    int maxtime;
    atomic<bool> globalEnding;

    void alloc_init();
    void sample_var_random();
    void sample_var_frequency();
    void readfile(const char *file);
    void read_var_clause_num(const char *file);
};