CXXFLAGS	= -march=native -O3
CXXFLAGS	+= -std=c++11
CXXFLAGS	+= -isystem ./include
CXXFLAGS	+= -MMD
CC	= g++

TARGET	=  EvalLossLocal

all: $(TARGET)

EvalLossLocal: EvalLossLocal.o

clean:
	rm -f *~ $(TARGET) *.d *.o

-include *.d
