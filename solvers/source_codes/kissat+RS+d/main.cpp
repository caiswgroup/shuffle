#include "pasolve.hpp"
#include <chrono>
#include <thread> 
using namespace std;

int main(int argc, char **argv) { 
    //1:file, 2:thread_num, 3:simplify_file 4:maxtime
    auto clk_st = std::chrono::high_resolution_clock::now(); 
    solve(argc, argv);
    auto clk_now = std::chrono::high_resolution_clock::now();
    int solve_time = std::chrono::duration_cast<std::chrono::milliseconds>(clk_now - clk_st).count();
    printf("c Use times: %d.%d\n", solve_time / 1000, solve_time % 1000);
    return 0;
}