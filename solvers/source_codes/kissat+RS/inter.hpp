#include "pamap.hpp"
#include <atomic>
#include <iostream>
using namespace std;
typedef long long ll;

struct kissat;

struct thread_data {
    kissat *solver;
    pthread_t thread_ptr;
    atomic<bool> solving;
    int nconf, ndec;
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
    vec<int> *clause;
    string infile;
    string simplify_file;
    string simplify_path;

    vec<int> assume_var;
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
    void print_cnf(const char *file);
    
    // simplify data structure
    int *f, nlit, *del, *a, *val, *color, *varval, *q, *nxtc, *seen, *resseen;
    int *clean, *neig, wrongF = 0;
    threadsafe_lookup_table<ll, int> C;
    vec<int> *occurp, *occurn;

    ll mapv(int a, int b);
    int find(int x);    
    bool res_is_empty(int var);
    void update_var_clause_label();
    void simplify_init();
    void simplify_resolution();
    void simplify_binary();
};