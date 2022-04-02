#include "inter.hpp"

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
inline int reverse(int x) {
    if (x % 2 == 1) return x - 1;
    else return x + 1;
}
inline ll pakissat::mapv(int a, int b) {
    return 1ll * a * nlit + (ll)b;
}
int pakissat::find(int x) {
    if (f[x] == x) return x;
    int fa = f[x];
    f[x] = find(fa);
    val[x] = val[fa] * val[x];
    return f[x];
}

void pakissat::simplify_init() {
    f = new int[vars + 10];
    val = new int[vars + 10];
    color = new int[vars + 10];
    varval = new int[vars + 10];
    q = new int[vars + 10];
    clean = new int[vars + 10];
    neig = new int[vars + 10];
    seen = new int[(vars << 1) + 10];
    del = new int[clauses + 10];
    nxtc = new int[clauses + 10];
    resseen = new int[(vars << 1) + 10];
    a = new int[500100];
    occurp = new vec<int>[vars + 1];
    occurn = new vec<int>[vars + 1];

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
}

bool pakissat::res_is_empty(int x) {
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

void pakissat::simplify_resolution() {    
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
        
        vec<int> vars;
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
        
        if (now_turn % 10 == 0) {
            for (int u = 0; u < vars.size(); u++) {
                int i = vars[u];
                int op = occurp[i].size(), on = occurn[i].size();
                assert(!clean[i]);
                int s = 0;
                for (int j = 0; j < op; j++)
                    if (!del[occurp[i][j]]) occurp[i][s++] = occurp[i][j];
                occurp[i].setsize(s);
                s = 0;
                for (int j = 0; j < on; j++)
                    if (!del[occurn[i][j]]) occurn[i][s++] = occurn[i][j];
                occurn[i].setsize(s);
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
    update_var_clause_label();
}

void pakissat::update_var_clause_label() {
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
    vars = remain_var, clauses = id;
}

void pakissat::simplify_binary() {
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