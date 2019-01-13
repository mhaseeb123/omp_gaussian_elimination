# Makefile for MPIC programs

BUILD = RELEASE

LINKCC = $(CXX)

# CXXFLAGS denotes flags for the C++ compiler
CXX = g++

ifeq (${BUILD}, DEBUG)
 BUILD_FLAGS = -fopenmp -O0 -g
else
 BUILD_FLAGS = -fopenmp -O3
endif

CXXFLAGS = ${BUILD_FLAGS} -Wall

SRCS = $(wildcard src/*.cpp)
EXEC := omp.exe

all: $(EXEC)

omp.exe:
	$(LINKCC) $(BUILD_FLAGS) -o $@ $(SRCS)

clean:
	rm -rf $(EXEC)

.PHONY: all clean
