SRCS = $(shell find . -name "*.cpp") 

OBJS = $(addsuffix .o, $(basename $(SRCS)))

EXEC = pakissat

LIBS = -lkissat -L../kissat_sat/build/kissat/build/ \
       -lpthread -lz -lm -static

CXXFLAGS = -I../kissat_sat/build/kissat -std=c++17 -O3

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC)