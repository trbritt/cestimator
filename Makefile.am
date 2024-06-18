CXX=g++
RM=rm -f
CPPFLAGS=-O3 -fPIC -fopenmp -Wall -std=c++17 -march=native -fno-math-errno
LDFLAGS=-lm -lz @ldflags@
DFLAGS=@dflags@

MODULES   := utils estimators evaluators optimizers SDE
SRC_DIR   := $(addprefix ./src/,$(MODULES))
BUILD_DIR := $(addprefix ./bin/,$(MODULES))
EXEDIR    := ./examples

SRC	      := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp)) 
EXES      := $(wildcard $(EXEDIR)/*.cpp)

OBJ	      := $(patsubst ./src/%.cpp,./bin/%.o,$(SRC)) 
OBJ_EXES  := $(patsubst ./examples/%.cpp,./bin/examples/%.x,$(EXES)) 
INCL      := $(addprefix -I,$(SRC_DIR)) -I/usr/include/eigen3 @iflags@

vpath %.cpp $(SRC_DIR)

define make-goal
$1/%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(DFLAGS) $(INCL) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs lib/libcestimator.so 

lib/libcestimator.so: $(OBJ)
	$(CXX) $^ -shared -o $@ $(LDFLAGS) $(LDLIBS)

checkdirs: $(BUILD_DIR)

examples: $(OBJ_EXES) lib/libcestimator.so

format:
	@uncrustify -c ucfy.cfg --no-backup ./src/**/*.hpp ./src/**/*.cpp ./examples/*.cpp

tests: tests/run_tests.x

tests/run_tests.x:
	@cmake tests -B./tests
	@make -C tests
	tests/run_tests.x
	
bin/examples/%.x: examples/%.cpp
	$(CXX) $(CPPFLAGS) $(DFLAGS) $(INCL) -I./src $^ -o $@ $(LDFLAGS) -L./lib -lcestimator $(LDLIBS)

$(BUILD_DIR):
	@mkdir -p $@
	@mkdir -p bin/examples
	@mkdir -p lib

clean:
	@rm -rf $(BUILD_DIR) bin/examples/*

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))

