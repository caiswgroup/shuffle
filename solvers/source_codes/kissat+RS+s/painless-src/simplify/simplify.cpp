#include "simplify.h"
#include <m4ri/m4ri.h>
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

int simplify::search_almost_one() {
    nlit = 2 * vars + 2;
    occur = new svec<int>[nlit]; 
    for (int i = 1; i <= clauses; i++) {
        del[i] = 0;
        if (clause[i].size() != 2) continue;
        int x = tolit(clause[i][0]);
        int y = tolit(clause[i][1]);
        ll id1 = mapv(x, y);
        ll id2 = mapv(y, x);
        B.get_bucket(id1).add_or_update_mapping(id1, i);
        B.get_bucket(id2).add_or_update_mapping(id2, i);
        occur[x].push(y);
        occur[y].push(x);
    }
    int *innei = new int[nlit];
    for (int i = 1; i <= vars; i++) {
        innei[i * 2] = innei[i * 2 + 1] = 0;
        seen[i * 2] = seen[i * 2 + 1] = color[i] = 0;
    }
    int t = 0;
    svec<int> ino, nei;
    for (int i = 2; i <= vars * 2 + 1; i++) {
        if (seen[i] || !occur[i].size()) continue;
        seen[i] = 1;
        nei.clear();
        for (int j = 0; j < occur[i].size(); j++)
            if (!seen[occur[i][j]])
                nei.push(occur[i][j]);
        do {
            ino.clear();
            ino.push(i);
            for (int j = 0; j < nei.size(); j++) {
                int v = nei[j], flag = 1;
                for (int k = 0; k < ino.size(); k++) {
                    ll id = mapv(v, ino[k]);
                    int d1 = B.get_bucket(id).value_for(id, 0);
                    if (!d1) {flag = 0; break;}
                    q[k] = d1;
                }
                if (flag) {
                    for (int k = 0; k < ino.size(); k++) {
                        del[q[k]] = 1;
                        ll id1 = mapv(v, ino[k]), id2 = mapv(ino[k], v);
                        B.get_bucket(id1).add_or_update_mapping(id1, 0);
                        B.get_bucket(id2).add_or_update_mapping(id2, 0);
                    }
                    ino.push(v);
            
                }
            }
            if (ino.size() >= 2) {
                card_one.push();
                for (int j = 0; j < ino.size(); j++) {
                    card_one[card_one.size() - 1].push(-toidx(ino[j]));
                }
                // printf("c catch constrain: ");
                // for (int j = 0; j < ino.size(); j++) printf("%d ", toidx(ino[j]));
                // puts("");
            }
        } while (ino.size() != 1);
    }
    return card_one.size();
}

void simplify::upd_occur(int v, int s) {
    int x = abs(v);
    int t = 0;
    if (v > 0) {
        for (int j = 0; j < occurp[x].size(); j++)
            if (occurp[x][j] != s) occurp[x][t++] = occurp[x][j]; 
        assert(t == occurp[x].size() - 1);
        occurp[x].setsize(t);
    }
    else {
        for (int j = 0; j < occurn[x].size(); j++)
            if (occurn[x][j] != s) occurn[x][t++] = occurn[x][j];
        assert(t == occurn[x].size() - 1);
        occurn[x].setsize(t);
    }
}

int simplify::scc_almost_one() {
    for (int i = 1; i <= vars; i++) {
        occurp[i].clear();
        occurn[i].clear();
    }
    cdel.growTo(card_one.size(), 0);
    for (int i = 0; i < card_one.size(); i++) {
        for (int j = 0; j < card_one[i].size(); j++)
            if (card_one[i][j] > 0)
                occurp[card_one[i][j]].push(i);
            else 
                occurn[-card_one[i][j]].push(i);
    }
    int flag = 1;
    do {
        flag = 0;
        for (int i = 1; i <= vars; i++) {
            if (!occurp[i].size() || !occurn[i].size()) continue;
            if (card_one.size() + occurp[i].size() * occurp[i].size() > 1e7 / vars) return 0;
            flag = 1;
            for (int ip = 0; ip < occurp[i].size(); ip++) 
                cdel[occurp[i][ip]] = 1;
            for (int in = 0; in < occurn[i].size(); in++) 
                cdel[occurn[i][in]] = 1;
            for (int ip = 0; ip < occurp[i].size(); ip++) {
                int op = occurp[i][ip];
                for (int in = 0; in < occurn[i].size(); in++) {
                    int on = occurn[i][in];
                    card_one.push();
                    cdel.push(0);
                    int id = card_one.size() - 1;
                    for (int j = 0; j < card_one[op].size(); j++) {
                        int v = card_one[op][j];
                        if (abs(v) != i) {
                            card_one[id].push(v);
                            if (v > 0) occurp[v].push(id);
                            else occurn[-v].push(id);
                        }
                    }
                    for (int j = 0; j < card_one[on].size(); j++) {
                        int v = card_one[on][j];
                        if (abs(v) != i) {
                            card_one[id].push(v);
                            if (v > 0) occurp[v].push(id);
                            else occurn[-v].push(id);
                        }
                    }
                }
            }
            for (int ip = 0; ip < occurp[i].size(); ip++) {
                int op = occurp[i][ip];
                for (int j = 0; j < card_one[op].size(); j++)
                    upd_occur(card_one[op][j], op);
            }
            
            for (int in = 0; in < occurn[i].size(); in++) {
                int on = occurn[i][in];
                for (int j = 0; j < card_one[on].size(); j++)
                    upd_occur(card_one[on][j], on);
            }
        }       
    } while(flag);
    // puts("\n");
    int t = 0;
    for (int i = 0 ; i < card_one.size(); i++) {
        if (cdel[i]) continue;
        ++t;
        // if (card_one[i].size() > 2) {
        //     printf("c catch constrain: ");
        //     for (int j = 0; j < card_one[i].size(); j++) printf("%d ", card_one[i][j]);
        //     puts("");
        // }
    }
    return t;
}


int simplify::card_elimination() {
    //sigma aixi <= b
    for (int i = 0; i < card_one.size(); i++) {
        if (cdel[i]) continue;
        int b = 1;
        mat.push();
        int id = mat.size() - 1;
        mat[id].growTo(vars + 1, 0);
        for (int j = 0; j < card_one[i].size(); j++) {
            if (card_one[i][j] < 0) b--;
            mat[id][abs(card_one[i][j])] += pnsign(card_one[i][j]);
        }
        mat[id][0] = b;
    }
    for (int i = 0; i < card_one.size(); i++)
        card_one[i].clear(true);
    card_one.clear(true);
    cdel.clear(true);
        
    printf("%d\n", mat.size());
    for (int i = 1; i <= clauses; i++) {
        if (del[i]) continue;
        int b = 1;
        for (int j = 0; j < clause[i].size(); j++)
            if (clause[i][j] < 0) b--;
        //convert >= to <=
        mat.push();
        int id = mat.size() - 1;
        mat[id].growTo(vars + 1, 0);
        for (int j = 0; j < clause[i].size(); j++) {
            mat[id][abs(clause[i][j])] += -pnsign(clause[i][j]);
        }
        mat[id][0] = -b;
    }
    printf("%d %d\n", clauses, mat.size());
    int row = mat.size();
    svec<int> mat_del, upp, low;
    mat_del.growTo(row, 0);
    
    for (int v = 1; v <= vars; v++) {
        //elimination v
        upp.clear();
        low.clear();
        // for (int i = 0; i < row; i++) {
        //     if (mat_del[i]) continue;
        //     printf("%d:\t", i);
        //     for (int j = 1; j <= vars; j++)
        //         if (fabs(mat[i][j])>1e-6) printf("%.2lf*%d + ", mat[i][j], j);
        //     printf("<= %.2lf\n", mat[i][0]);
        // }
        // puts("");
        for (int i = 0; i < row; i++) {
            if (mat_del[i]) continue;
            if (fabs(mat[i][v]) < 1e-6) continue;
            mat_del[i] = 1;
            if (mat[i][v] > 0) 
                upp.push(i);
            else 
                low.push(i);
        }
        // printf("%d %d %d\n", v, upp.size(), low.size());
        if (mat.size() + upp.size() * low.size() > 1e7/vars) return 1;
        for (int iu = 0; iu < upp.size(); iu++) {
            int u = upp[iu];
            for (int il = 0; il < low.size(); il++) {
                int l = low[il];
                double b = mat[u][0] / mat[u][v] + mat[l][0] / (-mat[l][v]);
                mat.push();
                mat_del.push(0);
                int id = mat.size() - 1;
                mat[id].growTo(vars + 1, 0);
                for (int j = 0; j <= vars; j++)
                    mat[id][j] = mat[u][j] / mat[u][v] + mat[l][j] / (-mat[l][v]);
                
                //check can get ?
                int ps = 0, ns = 0, t = -1;;
                for (int j = 1; j <= vars; j++) {
                    if (fabs(mat[id][j]) < 1e-6) continue;
                    t = j;
                    if (mat[id][j] > 0) ++ps;
                    if (mat[id][j] < 0) ++ns;
                }
                if (ns + ps == 1) {
                    if (ps) {
                        if (fabs(mat[id][0]) < 1e-6) 
                            printf("c [CE] get unit: %d = 0\n", t);
                        else if (mat[id][0] < 0) {
                            printf("c [CE] get wrong: %d < 0\n", t);
                            return 0;
                        }
                        mat.pop();
                    }
                    if (ns) {
                        double r = mat[id][0] / mat[id][t];
                        if (fabs(r - 1) < 1e-6)
                            printf("c [CE] get unit: %d = 1\n", t);
                        else if (r > 1) {
                            printf("c [CE] get wrong: %d > 1\n", t);
                            return 0;
                        }
                        mat.pop();
                    }
                }
                if (ns + ps == 0) {
                    if (mat[id][0] < -(1e-6)) {
                        printf("c [CE] get empty wrong clause: 0 < %.3lf\n", mat[id][0]);
                        return 0;
                    }
                    mat.pop();    
                }
            }
        }
        row = mat.size();
    }
    return 1;
}

int simplify::simplify_card() {
    int sone = search_almost_one();
    printf("c [CE] almost one cons: %d\n", sone);
    if (!sone) return 1;
    int scc = scc_almost_one();
    if (card_one.size() > 1e7 / vars || !scc) {
        for (int i = 0; i < card_one.size(); i++)
            card_one[i].clear(true);
        card_one.clear(true);
        cdel.clear(true);
        return 1;
    }
    printf("c [CE] scc cons: %d\n", scc);
    int res = card_elimination();
    for (int i = 0; i < mat.size(); i++)
        mat[i].clear(true);
    mat.clear(true);
    return res;
}

bool cmpvar(int x, int y) {
    return abs(x) < abs(y);
}

int simplify::cal_dup_val(int i) {
    for (int j = 0; j < clause[i].size(); j++) a[j] = clause[i][j];
    sort(a, a + clause[i].size(), cmpvar);
    int v = 0;
    for (int j = 0; j < clause[i].size(); j++)
        if (a[j] < 0) v |= (1 << j);
    return v;
}

int simplify::search_xors() {
    svec<int> xorsp;
    svec<bool> dup_table;
    for (int i = 1; i <= clauses; i++) {
        abstract[i] = 0;
        for (int j = 0; j < clause[i].size(); j++)
        abstract[i] |= 1 << (abs(clause[i][j]) & 31);
    }
    for (int i = 1; i <= clauses; i++) {
        if (nxtc[i] || del[i]) continue;
        nxtc[i] = 1;
        int l = clause[i].size();
        if (l <= 2 || l > MAX_XOR) continue;
        int required_num = 1 << (l - 2), skip = 0, mino = clauses + 1, mino_id = 0;
        for (int j = 0; j < l; j++) {
            int idx = abs(clause[i][j]);
            if (occurp[idx].size() < required_num || occurn[idx].size() < required_num) {
                skip = 1; break;
            }
            if (occurp[idx].size() + occurn[idx].size() < mino) {
                mino = occurp[idx].size() + occurn[idx].size();
                mino_id = idx;
            }
        }
        if (skip) continue;
        xorsp.clear();
        xorsp.push(i);
        for (int j = 0; j < occurp[mino_id].size(); j++) {
            int o = occurp[mino_id][j];
            if (!del[o] && !nxtc[o] && clause[o].size() == l && abstract[o] == abstract[i])
                xorsp.push(o);
        }
        for (int j = 0; j < occurn[mino_id].size(); j++) {
            int o = occurn[mino_id][j];
            if (!del[o] && !nxtc[o] && clause[o].size() == l && abstract[o] == abstract[i])
                xorsp.push(o);
        }
        if (xorsp.size() < 2 * required_num) continue;

        int rhs[2] = {0, 0};
        for (int j = 0; j < l; j++)
            seen[abs(clause[i][j])] = i;
        dup_table.clear();
        dup_table.growTo(1 << MAX_XOR, false);
        
        for (int j = 0; j < xorsp.size(); j++) {
            int o = xorsp[j], dup_v;
            bool xor_sign = true;
            for (int k = 0; k < clause[o].size(); k++) {
                if (!seen[abs(clause[o][k])]) goto Next;
                if (clause[o][k] < 0) xor_sign = !xor_sign;
            }
            dup_v = cal_dup_val(o);
            if (dup_table[dup_v]) continue;
            dup_table[dup_v] = true;
            rhs[xor_sign]++;
            if (j) nxtc[o] = 1;
        Next:;
        }
        assert(rhs[0] <= 2 * required_num);
        assert(rhs[1] <= 2 * required_num);
        if (rhs[0] == 2 * required_num)
            xors.push(xorgate(i, 0, l));
        if (rhs[1] == 2 * required_num)
            xors.push(xorgate(i, 1, l));
    }
    return xors.size();
}

bool cmpxorgate(xorgate a, xorgate b) {
    return a.sz > b.sz;
}

int simplify::ecc_var() {
    scc_id.clear();
    scc_id.growTo(vars + 1, -1);
    scc.clear();
    //sort(xors.data, xors.data + xors.sz, cmpxorgate);
    svec<int> xids;

    for (int i = 0; i < xors.size(); i++) {
        int x = xors[i].c;
        xids.clear();
        for (int j = 0; j < clause[x].size(); j++)
            if (scc_id[abs(clause[x][j])] != -1)
                xids.push(scc_id[abs(clause[x][j])]);

        if (xids.size() == 0) {
            scc.push();
            for (int j = 0; j < clause[x].size(); j++) {
                scc_id[abs(clause[x][j])] = scc.size() - 1;
                scc[scc.size() - 1].push(abs(clause[x][j]));
            }
        }
        else if (xids.size() == 1) {
            int id = xids[0];
            for (int j = 0; j < clause[x].size(); j++) {
                if (scc_id[abs(clause[x][j])] == -1) {
                    scc_id[abs(clause[x][j])] = id;
                    scc[id].push(abs(clause[x][j]));
                }
            }
        }
        else {
            int id_max = -1; int sz_max = 0;
            for (int j = 0; j < xids.size(); j++)
                if (scc[xids[j]].size() > sz_max)
                    sz_max = scc[xids[j]].size(), id_max = xids[j];
            for (int j = 0; j < xids.size(); j++) {
                if (xids[j] != id_max) {
                    svec<int>& v = scc[xids[j]];
                    for (int k = 0; k < v.size(); k++) {
                        scc_id[v[k]] = id_max;
                        scc[id_max].push(v[k]);
                    }
                    v.clear();
                }
            }
            for (int j = 0; j < clause[x].size(); j++) {
                if (scc_id[abs(clause[x][j])] == -1) {
                    scc_id[abs(clause[x][j])] = id_max;
                    scc[id_max].push(abs(clause[x][j]));
                }
            }
        }
    }
    return scc.size();
}

int simplify::ecc_xor() {
    for (int i = 0; i < scc.size(); i++) seen[i] = -1;
    for (int i = 0; i < xors.size(); i++) {
        int id = scc_id[abs(clause[xors[i].c][0])];
        if (seen[id] == -1) xor_scc.push(), seen[id] = xor_scc.size() - 1;
        int id2 = seen[id];
        xor_scc[id2].push(i);
    }
    return xor_scc.size();
}

int simplify::gauss_elimination() {
    gauss_eli_unit = gauss_eli_binary = 0;
    svec<int> v2mzd(vars + 1, -1);
    svec<int> mzd2v;
    for (int i = 0; i < xor_scc.size(); i++) {
        if (xor_scc[i].size() == 1) continue;
        int id = scc_id[abs(clause[xors[xor_scc[i][0]].c][0])];
        assert(scc[id].size() > 3);
        if (scc[id].size() > 1e7 / xor_scc[i].size()) continue;
        mzd2v.clear();
        sort(scc[id].data, scc[id].data + scc[id].size(), cmpvar);
        for (int j = 0; j < scc[id].size(); j++) {
            assert(scc[id][j] > 0);
            assert(scc[id][j] <= vars);
            v2mzd[scc[id][j]] = j;
            mzd2v.push(scc[id][j]);
        }
        int cols = scc[id].size() + 1;
        mzd_t* mat = mzd_init(xor_scc[i].size(), cols);
        for (int row = 0; row < xor_scc[i].size(); row++) {
            int k = xors[xor_scc[i][row]].c;
            for (int j = 0; j < clause[k].size(); j++) 
                mzd_write_bit(mat, row, v2mzd[abs(clause[k][j])], 1);
            if (xors[xor_scc[i][row]].rhs) 
                mzd_write_bit(mat, row, cols - 1, 1); 
        }
        mzd_echelonize(mat, true);
        mzd_free(mat);
        for (int row = 0, rhs; row < xor_scc[i].size(); row++) {
            svec<int> ones;
            for (int col = 0; col < cols - 1; col++) 
                if (mzd_read_bit(mat, row, col)) {
                    if (ones.size() == 2) goto NextRow;
                    ones.push(mzd2v[col]);
                }
            
            rhs = mzd_read_bit(mat, row, cols - 1);
            if (ones.size() == 1) {
                ++gauss_eli_unit;
                if (varval[ones[0]]) assert(varval[ones[0]] == (rhs ? 1 : -1));
                else {
                    varval[ones[0]] = rhs ? 1 : -1;
                    q[++tail] = ones[0];
                }
                printf("gauss: unit: %d = %d\n", ones[0], !rhs);
            }
            else if (ones.size() == 2) {
                ++gauss_eli_binary;
                assert(clauses == clause.size() - 1);
                int p = ones[0], q = rhs ? ones[1] : -ones[1];
                clause.push();
                ++clauses;
                clause[clauses].push(p);
                clause[clauses].push(q);
                clause.push();
                ++clauses;
                clause[clauses].push(-p);
                clause[clauses].push(-q);
                //printf("gauss: binary: %d %d = %d\n", ones[0], ones[1], !rhs);
            }
            else if (rhs) {printf("gauss: get empty clause\n"); return false;}
        NextRow:;
        }
    }
    return true;
}

bool simplify::simplify_gauss() {
    for (int i = 1; i <= vars; i++) {
        varval[i] = color[i] = seen[i] = 0;
        occurp[i].clear();
        occurn[i].clear();
        resseen[i] = resseen[i + vars] = 0;
    }
    for (int i = 1; i <= clauses; i++) del[i] = nxtc[i] = 0;

    // for (int i = 1; i <= clauses; i++) {
    //     int l = clause[i].size();
    //     for (int j = 0; j < l; j++) {
    //         if (clause[i][j] > 0) occurp[clause[i][j]].push(i);
    //         else occurn[-clause[i][j]].push(i);
    //     }
    // }

    head = 1, tail = 0;
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
    
    int nxors = search_xors();
    printf("c [GE] XORs: %d (time: 0.00)\n", nxors);
    if (!nxors) return true;
    int nvarscc = ecc_var();
    printf("c [GE] VAR SCC: %d\n", nvarscc);
    int nxorscc = ecc_xor();
    printf("c [GE] XOR SCCs: %d (time: 0.00)\n", nxorscc);
    int res = gauss_elimination();
    printf("c [GE] unary xor: %d, bin xor: %d, bin added\n", gauss_eli_unit, gauss_eli_binary);
    if (!res) {printf("c [GE] UNSAT\n"); }
    xors.clear(true);
    scc_id.clear(true);
    for (int i = 0; i < scc.size(); i++)
        scc[i].clear(true);
    scc.clear(true);
    for (int i = 0; i < xor_scc.size(); i++)
        xor_scc[i].clear(true);
    xor_scc.clear(true);
    if (!res) return false;
    if (gauss_eli_binary) {
        del.growTo(clauses + 1, 0);
        nxtc.growTo(clauses + 1, 0);
    }
    if (!gauss_eli_unit) {
        update_var_clause_label();
        return true;
    }

    for (int i = clauses - 2 * gauss_eli_binary + 1; i <= clauses; i++) {
        int l = clause[i].size();
        for (int j = 0; j < l; j++) {
            if (clause[i][j] > 0) occurp[clause[i][j]].push(i);
            else occurn[-clause[i][j]].push(i);
        }
    }
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
