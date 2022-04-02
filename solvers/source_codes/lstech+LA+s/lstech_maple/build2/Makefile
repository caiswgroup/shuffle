EXEC 	= g++
CFLAGS 	= -std=c++11 -O3 -Wall -g
ZLIB 	= -lz
heads = ../lstech/mtl/XAlloc.h ../lstech/mtl/Vec.h ../lstech/mtl/Sort.h ../lstech/mtl/Queue.h \
../lstech/mtl/Map.h ../lstech/mtl/IntTypes.h ../lstech/mtl/Heap.h ../lstech/mtl/Alloc.h \
../lstech/mtl/Alg.h ../lstech/core/SolverTypes.h ../lstech/core/Dimacs.h \
../lstech/utils/ParseUtils.h
deps = System.o SimpSolver.o ccnr.o Solver.o Options.o
maino = Main.o
includes = -I../lstech/
all: Relaxed libRelaxed.a


libRelaxed.a : $(deps)
	ar r libRelaxed.a $(deps)

Relaxed : $(deps) $(heads) $(maino)
	$(EXEC) -o Relaxed $(deps) $(maino) $(CFLAGS) $(ZLIB)

Main.o : ../lstech/simp/Main.cc SimpSolver.o
	$(EXEC) $(includes) -c $< -o Main.o $(CFLAGS) $(ZLIB)

System.o : ../lstech/utils/System.cc ../lstech/utils/System.h
	$(EXEC) $(includes) -c $< $(CFLAGS)  

SimpSolver.o : ../lstech/simp/SimpSolver.cc ../lstech/simp/SimpSolver.h 
	$(EXEC) $(includes) -c $< $(CFLAGS)

Solver.o : ../lstech/core/Solver.cc ../lstech/core/Solver.h
	$(EXEC) $(includes) -c $< $(CFLAGS)

ccnr.o : ../lstech/utils/ccnr.cc ../lstech/utils/ccnr.h 
	$(EXEC) $(includes) -c $< $(CFLAGS)

Options.o : ../lstech/utils/Options.cc ../lstech/utils/Options.h
	$(EXEC) $(includes) -c $< $(CFLAGS)


clean:
	-rm *.o *.a Relaxed

.FORCE:
.PHONY: clean
