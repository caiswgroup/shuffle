SRCS = $(shell find . -name "*.cpp") 

OBJS = $(addsuffix .o, $(basename $(SRCS)))

EXEC = pakissat

LIBS = -llstech -L../lstech_maple/build/release/lib/ \
       -lpthread -lz -lm -static

CXXFLAGS = -I../lstech_maple/lstech -std=c++17 -O3

$(EXEC): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(EXEC)