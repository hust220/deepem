SRC_DIR   := src/pp src/cppsci
TEST_DIR  := src/test
APPS_DIR  := apps
BIN_DIR   := bin

ifeq (,${ROME_CC})
	ROME_CC := mpiicpc
endif

ifneq (,$(findstring mpi,${ROME_CC}))
	MACROS  := -DUSEMPI
else
	MACROS  :=
endif

#the complier and build option
CC := $(ROME_CC)
LD := $(ROME_CC)

# -O3 optimization may get some error
FLAGS := -Isrc/cppsci -mkl -fopenmp -parallel-source-info=2 -xHost $(TBB) -std=c++11 -static-intel $(FLOAT_FLAGS)

ifeq (true,${RAND})
	FLAGS        := $(FLAGS) -DRAND
endif

ifeq (true,${MPI})
	FLAGS        := $(FLAGS) -DUSEMPI
endif

ifeq (1,${INFO})
	FLAGS        := $(FLAGS) -DJN_SHOW_INFO
endif
ifeq (2,${INFO})
	FLAGS        := $(FLAGS) -DJN_SHOW_INFO -DMPI_INFO_ALL
endif

# build on KNL
ifeq (true,${FLOAT})
	FLOAT_FLAGS := -DFLOAT_PRECISION
else
	FLOAT_FLAGS := 
endif

#debug for vtune and some vectorization report
ifeq (true,${DEBUG})
	BUILD_PREFIX := build/debug
	#DEBUG     := -g -debug inline-debug-info -qopt-report-phase=vec -qopt-report=5
	FLAGS        := $(FLAGS) -g -gdwarf-2
else
	BUILD_PREFIX := build/release
	FLAGS        := $(FLAGS) -DNDEBUG -O3
endif

#The following line disables offload.  By default we not offload, no flags needed
ifneq (true,${ROME_OFFLOAD})
	OFFLOAD      := -qno-offload
endif
#This line allows offload debugging
#OFFLOAD   := -DSERIAL_OFFLOAD
BUILD_DIR := $(addprefix $(BUILD_PREFIX)/, $(SRC_DIR) $(TEST_DIR) $(APPS_DIR))

SRC_CPP   := $(foreach sdir,$(SRC_DIR),  $(wildcard $(sdir)/*.cpp))
TEST_CPP  := $(foreach sdir,$(TEST_DIR), $(wildcard $(sdir)/*.cpp))
APPS_CPP  := $(foreach sdir,$(APPS_DIR), $(wildcard $(sdir)/*.cpp))

SRC_OBJ   := $(patsubst %.cpp, $(BUILD_PREFIX)/%.o, $(SRC_CPP))
TEST_OBJ  := $(patsubst %.cpp, $(BUILD_PREFIX)/%.o, $(TEST_CPP))
APPS_OBJ  := $(patsubst %.cpp, $(BUILD_PREFIX)/%.o, $(APPS_CPP))

#ALL_APPS  := $(notdir $(patsubst %.cpp, %, $(APPS_CPP))) test
ALL_APPS  := rome_picker rome_picker_test

vpath %.cpp $(SRC_DIR) $(TEST_DIR) $(APPS_DIR)

define make-goal

$1/%.o: %.cpp
	$(CC) $(FLAGS) $(MACROS) $(OFFLOAD) -c $$< -o $$@

$1/%.o: %.cc
	$(CC) $(FLAGS) $(MACROS) $(OFFLOAD) -c $$< -o $$@

endef

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))

.PHONY: all checkdirs clean $(ALL_APPS)

all: $(ALL_APPS)

#rome_pp: checkdirs $(SRC_OBJ) $(BUILD_PREFIX)/$(APPS_DIR)/rome_pp.o
#	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) $(BUILD_PREFIX)/$(APPS_DIR)/$@.o $(SRC_OBJ) -o bin/$@

rome_picker: checkdirs $(APPS_OBJ) $(SRC_OBJ)
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) $(APPS_OBJ) $(SRC_OBJ) -o bin/$@

rome_picker_test: checkdirs $(TEST_OBJ) $(SRC_OBJ)
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) $(TEST_OBJ) $(SRC_OBJ) -o bin/$@

checkdirs: $(BUILD_DIR) $(BIN_DIR)

$(BUILD_DIR):
	@mkdir -p $@

$(BIN_DIR):
	@mkdir bin

clean:
	@rm -rf bin build


