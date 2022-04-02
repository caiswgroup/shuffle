#include "simplify.h"
#include <atomic>
#include <omp.h>
#define TOLIT(x) ((x) > 0 ? (x) : ((-x) + vars))
#define NEG(x) ((x) > vars ? ((x) - vars) : ((x) + vars))

simplify::simplify():
  vars                  (0),
  clauses               (0),
  threads               (1),
  maxlen                (0)
{}

void simplify::alloc_init() {
    clause.growTo(clauses + 1);
}

inline int pnsign(int x) {
    return (x > 0 ? 1 : -1);
}
inline int sign(int x) {
    return x % 2 == 1 ? -1 : 1;
}
inline int tolit(int x) {
    if (x > 0) return x * 2;
    else return (-x) * 2 + 1;
}
inline int toidx(int x) {
    return (x & 1 ? -(x >> 1) : (x >> 1));
}
inline int reverse(int x) {
    if (x % 2 == 1) return x - 1;
    else return x + 1;
}
inline ll simplify::mapv(int a, int b) {
    return 1ll * a * nlit + (ll)b;
}
int simplify::find(int x) {
    if (f[x] == x) return x;
    int fa = f[x];
    f[x] = find(fa);
    val[x] = val[fa] * val[x];
    return f[x];
}

void simplify::simplify_init() {
    f = new int[vars + 10];
    val = new int[vars + 10];
    color = new int[vars + 10];
    varval = new int[vars + 10];
    q = new int[vars + 10];
    clean = new int[vars + 10];
    neig = new int[vars + 10];
    seen = new int[(vars << 1) + 10];
    del.growTo(clauses+1, 0);
    nxtc.growTo(clauses+1, 0);
    abstract = new int[clauses + 10];
    occurp = new svec<int>[vars + 1];
    occurn = new svec<int>[vars + 1];
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size();
        maxlen = max(maxlen, l);
        // del[i] = 0;
        // for (int j = 0; j < l; j++)
        //     if (clause[i][j] > 0) occurp[abs(clause[i][j])].push(i);
        //     else occurn[abs(clause[i][j])].push(i);
    }
    if (threads == 1) {
        resseen = new int[(vars << 1) + 10];
        a = new int[maxlen + 1];
    }
    else {
        // pa_seen = new int*[threads];
        // for (int i = 0; i < threads; i++)
        //     pa_seen[i] = new int[(vars << 1) + 10];
        // pa_a = new int*[threads];
        // for (int i = 0; i < threads; i++)
        //     pa_a[i] = new int[maxlen + 1];
        // val_equal = new int [clauses + 1];
        // val_fail_f = new int [clauses + 1];
        // val_fail_s = new int [clauses + 1];
    }
}


void simplify::release() {
    // for (int i = 0; i <= clauses; i++)
    //     clause[i].clear(true);
    // delete []clause;
    delete []f;
    delete []val;
    delete []color;
    delete []varval;
    // delete []q;
    delete []clean;
    delete []neig;
    delete []seen;
    delete []del;
    delete []nxtc;
    delete []abstract;
    // if (threads == 1) {
    delete []resseen;   
    delete []a;
    // }
    // else {
        // for (int i = 0; i < threads; i++) {
        //     delete []pa_seen[i];
        //     delete []pa_a[i];
        // }
        // delete []pa_seen;
        // delete []pa_a;
    // }
    for (int i = 0; i <= vars; i++)
        occurp[i].clear(true), occurn[i].clear(true);
    delete []occurp;
    delete []occurn;
}


bool simplify::res_is_empty(int x) {
    int op = occurp[x].size(), on = occurn[x].size();
    for (int i = 0; i < op; i++) {
        int o1 = occurp[x][i], l1 = clause[o1].size();
        if (del[o1]) continue;
        for (int j = 0; j < l1; j++)
            if (abs(clause[o1][j]) != x) resseen[abs(clause[o1][j])] = pnsign(clause[o1][j]);
        for (int j = 0; j < on; j++) {
            int o2 = occurn[x][j], l2 = clause[o2].size(), flag = 0;
            if (del[o2]) continue;
            for (int k = 0; k < l2; k++)
                if (abs(clause[o2][k]) != x && resseen[abs(clause[o2][k])] == -pnsign(clause[o2][k])) {
                    flag = 1; break;
                }
            if (!flag) {
                for (int j = 0; j < l1; j++)
                    resseen[abs(clause[o1][j])] = 0;
                return false;
            }
        }
        for (int j = 0; j < l1; j++)
            resseen[abs(clause[o1][j])] = 0;
    }
    return true; 
}

void simplify::simplify_resolution() {
    
    for (int i = 1; i <= vars; i++) {
        occurn[i].clear();
        occurp[i].clear();
        resseen[i] = resseen[i + vars] = clean[i] = 0;
    }
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size();
        del[i] = 0;
        for (int j = 0; j < l; j++)
            if (clause[i][j] > 0) occurp[abs(clause[i][j])].push(i);
            else occurn[abs(clause[i][j])].push(i);
    }
    for (int i = 1; i <= vars; i++) {
        color[i] = 0;
        if (occurn[i].size() == 0 && occurp[i].size() == 0) clean[i] = 1;  
    }

    int l = 1, r = 0;         
    for (int i = 1; i <= vars; i++) {
        int op = occurp[i].size(), on = occurn[i].size();
        if (op * on > op + on || clean[i]) continue;
        if (res_is_empty(i)) {
            q[++r] = i, clean[i] = 1;
        }
    }

    int choose_way = 0, now_turn = 0, seen_flag = 0;
    while (l <= r) {
        ++now_turn;
        for (int j = l; j <= r; j++) {
            int i = q[j];
            int op = occurp[i].size(), on = occurn[i].size();
            for (int j = 0; j < op; j++) del[occurp[i][j]] = 1;
            for (int j = 0; j < on; j++) del[occurn[i][j]] = 1;
        }
        int ll = l; l = r + 1;
        
        svec<int> vars;
        ++seen_flag;
        for (int u = ll; u <= r; u++) {
            int i = q[u];
            int op = occurp[i].size(), on = occurn[i].size();
            for (int j = 0; j < op; j++) {
                int o = occurp[i][j], l = clause[o].size();
                for (int k = 0; k < l; k++) {
                    int v = abs(clause[o][k]);
                    if (!clean[v] && seen[v] != seen_flag)
                        vars.push(v), seen[v] = seen_flag;
                }
            }
            for (int j = 0; j < on; j++) {
                int o = occurn[i][j], l = clause[o].size();
                for (int k = 0; k < l; k++) {
                    int v = abs(clause[o][k]);
                    if (!clean[v] && seen[v] != seen_flag)
                        vars.push(v), seen[v] = seen_flag;
                }
            }
        }
        for (int u = 0; u < vars.size(); u++) {
            int i = vars[u];
            int op = 0, on = 0;
            for (int j = 0; j < occurp[i].size(); j++) op += 1 - del[occurp[i][j]];
            for (int j = 0; j < occurn[i].size(); j++) on += 1 - del[occurn[i][j]];
            if (op * on > op + on) continue;
            if (res_is_empty(i)) {
                q[++r] = i, clean[i] = 1;
            }
        }
    }
    if (r) update_var_clause_label();
}

void simplify::update_var_clause_label() {
    int remain_var = 0;
    for (int i = 1; i <= clauses; i++) {
        if (del[i]) continue;
        int l = clause[i].size();
        for (int j = 0; j < l; j++) {
            if (color[abs(clause[i][j])] == 0) color[abs(clause[i][j])] = ++remain_var;       
        }
    }

    int id = 0;
    for (int i = 1; i <= clauses; i++) {
        if (del[i]) {clause[i].setsize(0); continue;}
        ++id;
        int l = clause[i].size();
        if (i == id) {
            for (int j = 0; j < l; j++)
                clause[id][j] = color[abs(clause[i][j])] * pnsign(clause[i][j]);    
            continue;
        }
        clause[id].setsize(0);
        for (int j = 0; j < l; j++)
            clause[id].push(color[abs(clause[i][j])] * pnsign(clause[i][j]));
    }
    printf("c After simplify: vars: %d -> %d , clauses: %d -> %d ,\n", vars, remain_var, clauses, id);
    for (int i = id + 1; i <= clauses; i++) 
        clause[i].clear(true);
    clause.setsize(id + 1);
    vars = remain_var, clauses = id;
}

void simplify::simplify_binary() {
    auto clk_st = std::chrono::high_resolution_clock::now();    
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size();
        for (int j = 0; j < l; j++)
            clause[i][j] = tolit(clause[i][j]);
    }
    nlit = 2 * vars + 2;
    int simplify = 1, turn = 0;
    for (int i = 1; i <= vars; i++) f[i] = i, val[i] = 1, varval[i] = color[i] = 0;
    for (int i = 1; i <= clauses; i++) del[i] = 0;

    int len = 0;
    for (int i = 1; i <= clauses; i++) {
        if (clause[i].size() != 2) continue;
        nxtc[++len] = i;
        ll id1 = mapv(clause[i][0], clause[i][1]),
           id2 = mapv(clause[i][1], clause[i][0]);
        C.get_bucket(id1).add_or_update_mapping(id1, i);
        C.get_bucket(id2).add_or_update_mapping(id2, i);
    }

    auto clk_s1 = std::chrono::high_resolution_clock::now();
    int timest = std::chrono::duration_cast<std::chrono::milliseconds>(clk_s1 - clk_st).count();
    printf("c add use %d.%d seconds\n", timest / 1000, timest % 1000);

    while (simplify) {
        if (wrongF) break;
        simplify = 0;
        ++turn;        
        for (int k = 1; k <= len; k++) {
            int i = nxtc[k];
            if (clause[i].size() != 2 || del[i]) continue;
            ll id1 = mapv(reverse(clause[i][0]), reverse(clause[i][1])),
               id2 = mapv(clause[i][0], reverse(clause[i][1])),
               id3 = mapv(reverse(clause[i][0]), clause[i][1]);
            int r = C.get_bucket(id1).value_for(id1, 0);
            if (r) {
                simplify = 1;
                del[r] = del[i] = 1;
                int fa = find(clause[i][0] >> 1), fb = find(clause[i][1] >> 1);
                int sig = sign(clause[i][0]) * sign(clause[i][1]) * (-1);
                //sig == 1 : a = b   -1 : a = -b
                if (fa != fb) {
                    if (fa < fb) {
                        f[fa] = fb;
                        val[fa] = sig / (val[clause[i][0] >> 1] * val[clause[i][1] >> 1]);
                        if (varval[fa])
                            varval[fb] = val[fa] * varval[fa];
                    }
                    else if (fa > fb) {
                        f[fb] = fa;
                        val[fb] = sig / (val[clause[i][0] >> 1] * val[clause[i][1] >> 1]);
                        if (varval[fb])
                            varval[fa] = val[fb] * varval[fb];
                    }
                }
                else {
                    assert(sig == val[clause[i][0] >> 1] * val[clause[i][1] >> 1]);
                }
            }
            int d1 = C.get_bucket(id2).value_for(id2, 0);
            if (d1) {
                del[d1] = del[i] = 1;
                simplify = 1;
                varval[clause[i][0] >> 1] = sign(clause[i][0]);
            }
            int d2 = C.get_bucket(id3).value_for(id3, 0);
            if (d2) {
                del[d2] = del[i] = 1;
                simplify = 1;
                varval[clause[i][1] >> 1] = sign(clause[i][1]); 
            }
        }
        auto clk_s2 = std::chrono::high_resolution_clock::now();
        int timest2 = std::chrono::duration_cast<std::chrono::milliseconds>(clk_s2 - clk_st).count();
        printf("c add use %d.%d seconds\n", timest2 / 1000, timest2 % 1000);

        for (int i = 1; i <= vars; i++) {
            int x = find(i);
            assert(f[i] == x);
            if (varval[i] && x != i) {
                if (varval[x]) assert(varval[x] == varval[i] * val[i]);
                else varval[x] = varval[i] * val[i];
            }
        }
        for (int i = 1; i <= vars; i++) 
            if (varval[f[i]]) {
                if (varval[i]) assert(varval[f[i]] == varval[i] * val[i]);
                else varval[i] = varval[f[i]] * val[i];
            }

        len = 0;

        for (int i = 1; i <= clauses; i++) {
            if (del[i]) continue;
            int l = clause[i].size(), oril = l;
            for (int j = 0; j < l; j++) {
                int fa = f[clause[i][j] >> 1];
                a[j] = tolit(sign(clause[i][j]) * val[clause[i][j] >> 1] * fa);
            }
            int t = 0;
            for (int j = 0; j < l; j++) {
                int x = varval[a[j] >> 1];
                if (x) {
                    int k = x * sign(a[j]);
                    if (k == 1) {
                        if (!del[i]) simplify = 1;
                        del[i] = 1, a[t++] = a[j];
                    }
                }
                else a[t++] = a[j];
            }
            if (t == 0) {
                printf("%d: EXSIT WRONG1 %d %d %d\n", i, l, clause[i][0], varval[clause[i][0] >> 1]);
                wrongF = 1;
            }
            if (t != l) simplify = 1, l = t;
            t = 0;
            for (int j = 0; j < l; j++) {
                if (resseen[a[j]] == i) continue;
                resseen[a[j]] = i, a[t++] = a[j];
            }
            if (t != l) simplify = 1, l = t;
            for (int j = 0; j < l; j++)
                if (resseen[reverse(a[j])] == i && !del[i])
                    del[i] = 1, simplify = 1;
            
            for (int j = 0; j < l; j++) resseen[a[j]] = 0;
                
            if (l == 1) {
                if (sign(a[0]) * varval[a[0] >> 1] == -1) {puts("Exsit WRONG2"); wrongF = 1;}
                varval[a[0] >> 1] = sign(a[0]);
                simplify = 1;
            }
            if (!del[i] && l == 2 && oril != 2) {
                nxtc[++len] = i;
                ll id1 = mapv(a[0], a[1]),
                   id2 = mapv(a[1], a[0]);
                C.get_bucket(id1).add_or_update_mapping(id1, i);
                C.get_bucket(id2).add_or_update_mapping(id2, i);
            }
            else if (!del[i] && l == 2 &&  oril == 2) {
                if (a[0] == clause[i][0] && a[1] == clause[i][1]) ;
                else if (a[1] == clause[i][0] && a[0] == clause[i][1]) ;
                else {
                    nxtc[++len] = i;
                    ll id1 = mapv(a[0], a[1]),
                       id2 = mapv(a[1], a[0]);
                    C.get_bucket(id1).add_or_update_mapping(id1, i);
                    C.get_bucket(id2).add_or_update_mapping(id2, i);
                }
            }
            clause[i].clear();
            for (int j = 0; j < l; j++)
                clause[i].push(a[j]);
        }

        for (int i = 1; i <= vars; i++) {
            int x = find(i);
            assert(f[i] == x);
            if (varval[i] && x != i) {
                if (varval[x]) assert(varval[x] == varval[i] * val[i]);
                else varval[x] = varval[i] * val[i];
            }
        }
        for (int i = 1; i <= vars; i++) 
            if (varval[f[i]]) {
                if (varval[i]) assert(varval[f[i]] == varval[i] * val[i]);
                else varval[i] = varval[f[i]] * val[i];
            }
        int outc = len;
    }

    if (wrongF) {puts("Exsit wrong! instance is unsat"); return;} 
    for (int i = 1; i <= clauses; i++) {
        if (del[i]) continue;
        int l = clause[i].size();
        for (int j = 0; j < l; j++) {
            clause[i][j] = sign(clause[i][j]) * (clause[i][j] >> 1);
        }
    }
    update_var_clause_label();
}

void simplify::pa_simplify_binary() {
//     omp_set_num_threads(threads);
//     auto clk_st = std::chrono::high_resolution_clock::now();
// #pragma omp parallel for schedule(guided, 4)      
//     for (int i = 1; i <= clauses; i++) {
//         del[i] = 0;
//         int l = clause[i].size();
//         for (int j = 0; j < l; j++)
//             clause[i][j] = tolit(clause[i][j]);
//     }
//     nlit = 2 * vars + 2;
//     int simplify = 1, turn = 0;
// #pragma omp parallel for schedule(guided, 4)  
//     for (int i = 1; i <= vars; i++) f[i] = i, val[i] = 1, varval[i] = color[i] = 0;
    
    
//     std::atomic<int> len(0);   
//     for (int i = 1; i <= clauses; i++) {
//         if (clause[i].size() != 2) continue;
//         nxtc[++len] = i;
//     }

// #pragma omp parallel for schedule(guided, 4)    
//     for (int j = 1; j <= len; j++) {
//         int i = nxtc[j];
//         ll id1 = mapv(clause[i][0], clause[i][1]),
//            id2 = mapv(clause[i][1], clause[i][0]);
//         C.get_bucket(id1).add_or_update_mapping(id1, i);
//         C.get_bucket(id2).add_or_update_mapping(id2, i);
//     }

//     auto clk_s1 = std::chrono::high_resolution_clock::now();
//     int timest = std::chrono::duration_cast<std::chrono::milliseconds>(clk_s1 - clk_st).count();
//     printf("c add use %d.%d seconds\n", timest / 1000, timest % 1000);

//     while (simplify) {
//         if (wrongF) break;
//         simplify = 0;
//         ++turn;        
        
// #pragma omp parallel
// {
// #pragma omp for schedule(guided, 4)   
//         for (int k = 1; k <= len; k++) {
//             int i = nxtc[k];
//             if (clause[i].size() != 2 || del[i]) continue;
//             ll id1 = mapv(reverse(clause[i][0]), reverse(clause[i][1])),
//                id2 = mapv(clause[i][0], reverse(clause[i][1])),
//                id3 = mapv(reverse(clause[i][0]), clause[i][1]);
//             val_equal[i] = C.get_bucket(id1).value_for(id1, 0);
//             val_fail_f[i] = C.get_bucket(id2).value_for(id2, 0);
//             val_fail_s[i] = C.get_bucket(id3).value_for(id3, 0);
//         }
// #pragma omp barrier     
// #pragma omp master
// {
//         for (int k = 1; k <= len; k++) {
//             int i = nxtc[k];
//             if (clause[i].size() != 2 || del[i]) continue;
//             int r = val_equal[i];
//             if (r) {
//                 simplify = 1;
//                 del[r] = del[i] = 1;
//                 int fa = find(clause[i][0] >> 1), fb = find(clause[i][1] >> 1);
//                 int sig = sign(clause[i][0]) * sign(clause[i][1]) * (-1);
//                 //sig == 1 : a = b   -1 : a = -b
//                 if (fa != fb) {
//                     if (fa < fb) {
//                         f[fa] = fb;
//                         val[fa] = sig / (val[clause[i][0] >> 1] * val[clause[i][1] >> 1]);
//                         if (varval[fa])
//                             varval[fb] = val[fa] * varval[fa];
//                     }
//                     else if (fa > fb) {
//                         f[fb] = fa;
//                         val[fb] = sig / (val[clause[i][0] >> 1] * val[clause[i][1] >> 1]);
//                         if (varval[fb])
//                             varval[fa] = val[fb] * varval[fb];
//                     }
//                 }
//                 else {
//                     assert(sig == val[clause[i][0] >> 1] * val[clause[i][1] >> 1]);
//                 }
//             }
//             if (val_fail_f[i]) {
//                 del[val_fail_f[i]] = del[i] = 1;
//                 simplify = 1;
//                 varval[clause[i][0] >> 1] = sign(clause[i][0]);
//             }
//             if (val_fail_s[i]) {
//                 del[val_fail_s[i]] = del[i] = 1;
//                 simplify = 1;
//                 varval[clause[i][1] >> 1] = sign(clause[i][1]); 
//             }
//         }
//         auto clk_s2 = std::chrono::high_resolution_clock::now();
//         int timest2 = std::chrono::duration_cast<std::chrono::milliseconds>(clk_s2 - clk_st).count();
//         printf("c add use %d.%d seconds\n", timest2 / 1000, timest2 % 1000);

//         for (int i = 1; i <= vars; i++) {
//             int x = find(i);
//             assert(f[i] == x);
//             if (varval[i] && x != i) {
//                 if (varval[x]) assert(varval[x] == varval[i] * val[i]);
//                 else varval[x] = varval[i] * val[i];
//             }
//         }
// }
// #pragma omp barrier 

// #pragma omp for schedule(guided, 4) 
//         for (int i = 1; i <= vars; i++) 
//             if (varval[f[i]]) {
//                 if (varval[i]) assert(varval[f[i]] == varval[i] * val[i]);
//                 else varval[i] = varval[f[i]] * val[i];
//             }

//         len = 0;

// #pragma omp barrier 
// #pragma omp for schedule(guided, 4) 
//         for (int i = 1; i <= clauses; i++) {
//             if (del[i]) continue;
//             int l = clause[i].size(), oril = l, id = omp_get_thread_num();
//             for (int j = 0; j < l; j++) {
//                 int fa = f[clause[i][j] >> 1];
//                 pa_a[id][j] = tolit(sign(clause[i][j]) * val[clause[i][j] >> 1] * fa);
//             }
//             int t = 0;
//             for (int j = 0; j < l; j++) {
//                 int x = varval[pa_a[id][j] >> 1];
//                 if (x) {
//                     int k = x * sign(pa_a[id][j]);
//                     if (k == 1) {
//                         if (!del[i]) simplify = 1;
//                         del[i] = 1, pa_a[id][t++] = pa_a[id][j];
//                     }
//                 }
//                 else pa_a[id][t++] = pa_a[id][j];
//             }
//             if (t == 0) {
//                 printf("%d: EXSIT WRONG1 %d %d %d\n", i, l, clause[i][0], varval[clause[i][0] >> 1]);
//                 wrongF = 1;
//             }
//             if (t != l) simplify = 1, l = t;
//             t = 0;
//             for (int j = 0; j < l; j++) {
//                 if (pa_seen[id][pa_a[id][j]] == i) continue;
//                 pa_seen[id][pa_a[id][j]] = i, pa_a[id][t++] = pa_a[id][j];
//             }
//             if (t != l) simplify = 1, l = t;
//             for (int j = 0; j < l; j++)
//                 if (pa_seen[id][reverse(pa_a[id][j])] == i && !del[i])
//                     del[i] = 1, simplify = 1;
            
//             for (int j = 0; j < l; j++)
//                 pa_seen[id][pa_a[id][j]] = 0;
                
//             if (l == 1) {
//                 if (sign(pa_a[id][0]) * varval[pa_a[id][0] >> 1] == -1) {puts("Exsit WRONG2"); wrongF = 1;}
//                 varval[pa_a[id][0] >> 1] = sign(pa_a[id][0]);
//                 nxtc[len.fetch_add(1) + 1] = i;
//                 simplify = 1;
//             }
//             if (!del[i] && l == 2 && oril != 2) {
//                 nxtc[len.fetch_add(1) + 1] = i;
//             }
//             else if (!del[i] && l == 2 &&  oril == 2) {
//                 if (pa_a[id][0] == clause[i][0] && pa_a[id][1] == clause[i][1]) ;
//                 else if (pa_a[id][1] == clause[i][0] && pa_a[id][0] == clause[i][1]) ;
//                 else nxtc[len.fetch_add(1) + 1] = i;
//             }
//             clause[i].clear();
//             for (int j = 0; j < l; j++)
//                 clause[i].push(pa_a[id][j]);
//         }
// #pragma omp barrier 
// #pragma omp for schedule(guided, 4)
//         for (int j = 1; j <= len; j++) {
//             int i = nxtc[j];
//             assert(clause[i].size() <= 2);
//             if (clause[i].size() == 1) 
//                  varval[clause[i][0] >> 1] = sign(clause[i][0]);
//             else {
//                 ll id1 = mapv(clause[i][0], clause[i][1]),
//                    id2 = mapv(clause[i][1], clause[i][0]);
//                 C.get_bucket(id1).add_or_update_mapping(id1, i);
//                 C.get_bucket(id2).add_or_update_mapping(id2, i);
//             }
//         }
// #pragma omp barrier

// #pragma omp master
// {
//         for (int i = 1; i <= vars; i++) {
//             int x = find(i);
//             assert(f[i] == x);
//             if (varval[i] && x != i) {
//                 if (varval[x]) assert(varval[x] == varval[i] * val[i]);
//                 else varval[x] = varval[i] * val[i];
//             }
//         }
// }
// #pragma omp barrier 
// #pragma omp for schedule(guided, 4) 
//         for (int i = 1; i <= vars; i++) 
//             if (varval[f[i]]) {
//                 if (varval[i]) assert(varval[f[i]] == varval[i] * val[i]);
//                 else varval[i] = varval[f[i]] * val[i];
//             }
// #pragma omp barrier         
// }
//     }
//     if (wrongF) {puts("Exsit wrong! instance is unsat"); return;} 
//     for (int i = 1; i <= clauses; i++) {
//         if (del[i]) continue;
//         int l = clause[i].size();
//         for (int j = 0; j < l; j++) {
//             clause[i][j] = sign(clause[i][j]) * (clause[i][j] >> 1);
//         }
//     }
//     update_var_clause_label();
}

bool simplify::simplify_easy_clause() {
    for (int i = 1; i <= vars; i++) {
        varval[i] = color[i] = 0;
        occurp[i].clear();
        occurn[i].clear();
        resseen[i] = resseen[i + vars] = 0;
    }
    for (int i = 1; i <= clauses; i++) del[i] = 0;
    int head = 1, tail = 0;
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size(), t = 0;
        for (int j = 0; j < l; j++) {
            int lit = TOLIT(clause[i][j]);
            if (resseen[lit] == i) continue;
            if (resseen[NEG(lit)] == i) {del[i] = 1; break;}
            clause[i][t++] = clause[i][j];
            resseen[lit] = i;
        }
        if (del[i]) continue;
        clause[i].setsize(t);
        for (int j = 0; j < t; j++) {
            if (clause[i][j] > 0) occurp[clause[i][j]].push(i);
            else occurn[-clause[i][j]].push(i);
        }
        if (t == 0) return false;
        if (t == 1) {
            int lit = clause[i][0];
            if (varval[abs(lit)]) {
                if (varval[abs(lit)] == pnsign(lit)) continue;
                else return false;
            }
            varval[abs(lit)] = pnsign(lit); 
            q[++tail] = abs(lit); 
            del[i] = 1;
        }
    }
    for (int i = 1; i <= vars + vars; i++) resseen[i] = 0;
    while (head <= tail) {
        int x = q[head++];
        if (varval[x] == 1) {
            for (int i = 0; i < occurp[x].size(); i++)
                del[occurp[x][i]] = 1;
            for (int i = 0; i < occurn[x].size(); i++) {
                int o = occurn[x][i], t = 0;
                if (del[o]) continue;
                for (int j = 0; j < clause[o].size(); j++) {
                    if (varval[abs(clause[o][j])] == pnsign(clause[o][j])) {
                        del[o] = 1; break;
                    }
                    if (varval[abs(clause[o][j])] == -pnsign(clause[o][j])) continue;
                    clause[o][t++] = clause[o][j];
                }
                if (del[o]) continue;
                assert(t < clause[o].size());
                clause[o].setsize(t);
                if (t == 0) return false;
                if (t == 1) {
                    int lit = clause[o][0];
                    if (varval[abs(lit)]) {
                        if (varval[abs(lit)] == pnsign(lit)) continue;
                        else return false;
                    }
                    varval[abs(lit)] = pnsign(lit); 
                    q[++tail] = abs(lit); 
                    del[o] = 1;
                }
            }
        }
        else {
            for (int i = 0; i < occurn[x].size(); i++)
                del[occurn[x][i]] = 1;
            for (int i = 0; i < occurp[x].size(); i++) {
                int o = occurp[x][i], t = 0;
                if (del[o]) continue;
                for (int j = 0; j < clause[o].size(); j++) {
                    if (varval[abs(clause[o][j])] == pnsign(clause[o][j])) {
                        del[o] = 1; break;
                    }
                    if (varval[abs(clause[o][j])] == -pnsign(clause[o][j])) continue;
                    clause[o][t++] = clause[o][j];
                }
                if (del[o]) continue;
                assert(t < clause[o].size());
                clause[o].setsize(t);
                if (t == 0) return false;
                if (t == 1) {
                    int lit = clause[o][0];
                    if (varval[abs(lit)]) {
                        if (varval[abs(lit)] == pnsign(lit)) continue;
                        else return false;
                    }
                    varval[abs(lit)] = pnsign(lit); 
                    q[++tail] = abs(lit); 
                    del[o] = 1;
                }
            }
        }
    }
    update_var_clause_label();
    return true;
}

void simplify::simplify_clause() {
    auto clk_st = std::chrono::high_resolution_clock::now();    
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size();
        for (int j = 0; j < l; j++)
            clause[i][j] = tolit(clause[i][j]);
    }
    nlit = 2 * vars + 2;
    int simplify = 1, turn = 0;
    for (int i = 1; i <= vars; i++) varval[i] = color[i] = 0;
    for (int i = 1; i <= clauses; i++) del[i] = 0;

    while (simplify) {
        if (wrongF) break;
        simplify = 0;
        ++turn;        
        for (int i = 1; i <= clauses; i++) {
            if (del[i]) continue;
            int l = clause[i].size(), oril = l;
            for (int j = 0; j < l; j++) a[j] = clause[i][j];
            int t = 0;
            for (int j = 0; j < l; j++) {
                int x = varval[a[j] >> 1];
                if (x) {
                    int k = x * sign(a[j]);
                    if (k == 1) {
                        if (!del[i]) simplify = 1;
                        del[i] = 1, a[t++] = a[j];
                    }
                }
                else a[t++] = a[j];
            }
            if (t == 0) {
                printf("%d: EXSIT WRONG1 %d %d %d\n", i, l, clause[i][0], varval[clause[i][0] >> 1]);
                wrongF = 1;
            }
            if (t != l) simplify = 1, l = t;
            t = 0;
            for (int j = 0; j < l; j++) {
                if (resseen[a[j]] == i) continue;
                resseen[a[j]] = i, a[t++] = a[j];
            }
            if (t != l) simplify = 1, l = t;
            for (int j = 0; j < l; j++)
                if (resseen[reverse(a[j])] == i && !del[i])
                    del[i] = 1, simplify = 1;
            
            for (int j = 0; j < l; j++) resseen[a[j]] = 0;
                
            if (l == 1) {
                if (sign(a[0]) * varval[a[0] >> 1] == -1) {puts("Exsit WRONG2"); wrongF = 1;}
                varval[a[0] >> 1] = sign(a[0]);
                simplify = 1;
            }
            clause[i].clear();
            for (int j = 0; j < l; j++)
                clause[i].push(a[j]);
        }
    }

    if (wrongF) {puts("Exsit wrong! instance is unsat"); return;} 
    for (int i = 1; i <= clauses; i++) {
        if (del[i]) continue;
        int l = clause[i].size();
        for (int j = 0; j < l; j++) {
            clause[i][j] = sign(clause[i][j]) * (clause[i][j] >> 1);
        }
    }
    if (turn > 1) update_var_clause_label();
}
