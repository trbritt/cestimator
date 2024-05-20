CXX=g++
RM=rm -f
CPPFLAGS=-O3 -Wall -std=c++17 -march=native -fno-math-errno
LDFLAGS=-lm -lz @ldflags@
DFLAGS=@dflags@

MODULES   := utils estimators
SRC_DIR   := $(addprefix ./src/,$(MODULES))
BUILD_DIR := $(addprefix ./bin/,$(MODULES))

SRC_MAIN  := ./src/cestimator.cpp
SRC	      := $(SRC_MAIN) $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp)) 
OBJ_MAIN  := ./bin/cestimator.o
OBJ	      := $(OBJ_MAIN) $(patsubst ./src/%.cpp,./bin/%.o,$(SRC)) 
INCL      := $(addprefix -I,$(SRC_DIR)) -I/usr/include/eigen3

vpath %.cpp $(SRC_DIR)

define make-goal
$1/%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(DFLAGS) $(INCL) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs bin/cestimator.x

bin/cestimator.x: $(OBJ) 
	$(CXX) $^ -o $@ $(LDFLAGS) $(LDLIBS)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) ./bin/cestimator*

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))

./bin/cestimator.o: $(SRC_MAIN)
	$(CXX) $(CPPFLAGS) $(DFLAGS) $(INCL) -c $< -o $@
