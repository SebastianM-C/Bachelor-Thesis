# Folders
BIN = ../../../Bin
OBJ = $(BIN)/Obj

FC = gfortran
FCFLAGS = -march=native -O2 -g -fopt-info -Wall -Wextra -Wconversion \
        	-fbacktrace #-ffpe-trap=zero,overflow,underflow,denormal
#FLFLAGS =

EXECUTABLE_NAME = hamilt

# Files
SOURCE_FILES := $(wildcard *.for)

INPUT_FILES = hamilt.inp
EXECUTABLE_FILES = $(EXECUTABLE_NAME:%=$(BIN)/%)
OBJECT_FILES     = $(SOURCE_FILES:%.for=$(OBJ)/%.o)
# ^^^ A more succinct expression for $(OBJECT_FILES), using
#     http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

COMPILE = $(FC) $(DEPFLAGS) $(FCFLAGS) -c

build: $(EXECUTABLE_FILES) input

clean:
	rm -rf $(BIN)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean

$(EXECUTABLE_FILES) : $(OBJECT_FILES)
	@$(FC) $(FLFLAGS) -o $@ $^
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
	@echo "Build successful!"

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJ)/%.o : %.for
	@echo Compiling $<
	@mkdir -p $(@D)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029
	@$(COMPILE) -o $@ $<

input : $(BIN)/$(INPUT_FILES)
$(BIN)/$(INPUT_FILES) : $(INPUT_FILES)
	@cp $(INPUT_FILES) $(BIN)
	@echo "Copied input"
