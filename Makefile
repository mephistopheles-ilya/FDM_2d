ifeq ($(origin CXX),default)
  CXX = g++
endif

CXXFLAGS ?=  -O3 -march=native -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align\
			-Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security\
			-Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long\
			-Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format\
			#-fsanitize=leak,undefined,address

CSRC = main.cpp fill_matrix.cpp matrix_pattern.cpp 
COBJ = main.o fill_matrix.o matrix_pattern.o 

#LDFLAGS ?= -fsanitize=leak,undefined,address

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

.PHONY: all
all: a.out

a.out: $(COBJ)
	$(CXX) $^ -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf *.o
