#include "inter.hpp"
#include <cstring>
#include <fstream>
using namespace std;

void pakissat::readfile(const char *file) {
    string infile(file);
	ifstream fin(infile);
    fin.seekg(0, fin.end);
	size_t file_len = fin.tellg();
	fin.seekg(0, fin.beg);
	char *data = new char[file_len + 1];
	fin.read(data, file_len);
	fin.close();
	data[file_len] = '\0';
    char *p = data;
    do {
        while (*p != 'p') p++;
        if (*(p + 2) == 'c' && *(p + 3) == 'n' && *(p + 4) == 'f')
            {p += 5; break;}
        else p++;
    } while (true);

    int n = 0, m = 0;
    while (*p < '0' || *p > '9') ++p;
	while (*p >= '0' && *p <= '9')
		n = n * 10 + *p - '0', ++p;

    while (*p < '0' || *p > '9') ++p;
	while (*p >= '0' && *p <= '9')
		m = m * 10 + *p - '0', ++p;
    
    int v = 0, sign = 0;
    for (int i = 1; i <= m; i++) {
        do {
            v = 0, sign = 1;
            while ((*p < '0' || *p > '9') && (*p != '-')) ++p;
            if (*p == '-') sign = -1, p++;
	        while (*p >= '0' && *p <= '9')
		        v = v * 10 + *p - '0', ++p;
            if (v) {
                clause[i].push(v * sign);
                // if (sign == 1) occurp[v].push_back(i);
                // else occurn[v].push_back(i);
            }
        } while (v);       
    }
    delete []data;
}

void pakissat::read_var_clause_num(const char *file) {
    freopen(file, "r", stdin);
    char s[500];
    int lastp = 0;
    while (scanf("%s", s) != -1) {
        if (strlen(s) == 3 && s[0] == 'c' && s[1] == 'n' && s[2] == 'f' && lastp == 1) {
            break;
        }
        if (strlen(s) == 1 && s[0] == 'p') lastp = 1;
        else lastp = 0;
    }
    scanf("%d %d", &vars, &clauses);
    fclose(stdin);
}

void pakissat::print_cnf(const char *file) {
    freopen(file, "w", stdout);
    printf("p cnf %d %d\n", vars, clauses);
    for (int i = 1; i <= clauses; i++) {
        int l = clause[i].size();
        for (int j = 0; j < l; j++)
            printf("%d ", clause[i][j]);
        puts("0");
    }
    fclose(stdout);
}