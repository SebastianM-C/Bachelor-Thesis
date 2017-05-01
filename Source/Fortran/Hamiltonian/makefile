# Folders
BIN = ../../../Bin
OBJ = $(BIN)/Obj

FC = gfortran
FCFLAGS = -march=native -O2 -g \
        	-fopt-info -Wall -Wextra -Wconversion -fbacktrace \
        	#-ffpe-trap=zero,overflow,underflow,denormal
FLFLAGS =

MAKEDEPEND=fort_depend.py
DEP_FILE = makefile.dep

EXECUTABLE_NAME = hamilt

# Files
SOURCE_FILES := $(wildcard *.f90)

INPUT_FILES = hamilt.inp
EXECUTABLE_FILES = $(EXECUTABLE_NAME:%=$(BIN)/%)
OBJECT_FILES     = $(SOURCE_FILES:%.f90=$(OBJ)/%.o)
# ^^^ A more succinct expression for $(OBJECT_FILES), using
#     http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

COMPILE = $(FC) $(FCFLAGS)

build: $(EXECUTABLE_FILES) input

clean:
	rm -rf $(BIN)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: depend build clean

$(EXECUTABLE_FILES) : $(OBJECT_FILES)
	@$(COMPILE) -o $@ $^ $(FLFLAGS)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
	@echo "Build successful!"

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJ)/%.o : %.f90
	@echo Compiling $<
	@mkdir -p $(@D)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029
	@$(COMPILE) -c -o $@ $<

all: $(EXECUTABLE_NAME) $(DEP_FILE)

depend: $(DEP_FILE)

$(DEP_FILE): $(OBJECTS)
	@echo "Making dependencies!"
	python2 $(MAKEDEPEND) -w -b $(OBJ) -f $(SOURCE_FILES)

input : $(BIN)/$(INPUT_FILES)
$(BIN)/$(INPUT_FILES) : $(INPUT_FILES)
	@cp $(INPUT_FILES) $(BIN)
	@echo "Copied input"

# Dependencies
include $(DEP_FILE)
