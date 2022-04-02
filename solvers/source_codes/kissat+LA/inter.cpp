#include "inter.hpp"

pakissat::pakissat():
  vars                  (0),
  clauses               (0),
  orivars               (0),
  oriclauses            (0),
  thread_num            (1),
  simplify_thread_num   (1),
  assume_stack_head     (1),
  assume_stack_tail     (0),
  model_pid             (0),
  finalResult           (0),
  winner                (0),
  maxtime               (5000),
  globalEnding          (false)
{}

void pakissat::alloc_init() {
    clause = new vec<int> [clauses + 1];
    assume_state = new atomic<int> [6010];
    T = new thread_data[thread_num];
    for (int i = 0 ; i < thread_num; i++) {
        T[i].nconf = T[i].ndec = 0;
        T[i].solving = false;
    }
}
