SRC=$(wildcard src/*.cpp)
INC=$(wildcard src/*.h)
OBJS=$(patsubst src/%.cpp, bld/%.o, $(SRC))
CXXFLAGS=-std=c++11 -g -Wall -Werror -fopenmp

CUDA_SRC=$(wildcard src/*.cu)
CUDA_INC=$(wildcard src/*.cuh)
CUDA_OBJS=$(patsubst src/%.cu, bld/cuda/%.o, $(CUDA_SRC))

isofinder: $(CUDA_OBJS) $(CUDA_INC) $(OBJS) $(INC)
	nvcc -o $@ $(OBJS) $(CUDA_OBJS) -lseqio -lcuda -lgomp

# Compile a source file.
bld/%.o: src/%.cpp Makefile
	@mkdir -p bld
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile a cuda file.
bld/cuda/%.o: src/%.cu Makefile
	@mkdir -p bld/cuda
	nvcc -arch=sm_13 -c $< -o $@

clean:
	rm -rf bld
	rm -f isofinder
