#include "inter.hpp"
#include <algorithm>
int *freq;

void pakissat::sample_var_random() {
    int *id = new int[vars + 1];
    for (int i = 0; i <= vars; i++) id[i] = i;
    random_shuffle(id + 1, id + vars + 1);
    for (int i = 0; i <= min(50 * thread_num, vars); i++) {
        assume_var.push(-id[i]);
        assume_var.push(id[i]);
    }
    delete []id;
    assume_stack_head = thread_num + 1;
    assume_stack_tail = assume_var.size();
    assume_state[0] = 20;
    for (int i = 1; i < assume_stack_tail; i++) assume_state[i] = 0;
}

bool cmp_freq(int x, int y) {
    return freq[x] > freq[y];
}

void pakissat::sample_var_frequency() {
    freq = new int[vars + 1];
    int *id = new int[vars + 1];
    for (int i = 0; i <= vars; i++) freq[i] = 0, id[i] = i;
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size();
        for (int j = 0; j < l; j++)
            freq[abs(clause[i][j])]++;
    }
    sort(id + 1, id + vars + 1, cmp_freq);
    for (int i = 0; i <= min(50 * thread_num, vars); i++) {
        assume_var.push(-id[i]);
        assume_var.push(id[i]);
    }
    delete []id;
    delete []freq;
    assume_stack_head = thread_num + 1;
    assume_stack_tail = assume_var.size();
    assume_state[0] = 20;
    for (int i = 1; i < assume_stack_tail; i++) assume_state[i] = 0;
}