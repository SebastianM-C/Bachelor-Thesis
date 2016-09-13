FC = gfortran
FCFLAGS = -march=native -O2 -g -Wall -Wextra -Wconversion \
        	-fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
#FLFLAGS =

EXECUTABLE_NAME=hamilt

# Folders
BIN = ../../../Bin
OBJ = $(BIN)/Obj

# Files
SOURCE_FILES=\
    hamiltv1.for \
    elem.for \
    sdiago.for

EXECUTABLE_FILES = $(EXECUTABLE_NAME:%=$(BIN)/%)
OBJECT_FILES     = $(SOURCE_FILES:%.for=$(OBJ)/%.o)
# ^^^ A more succinct expression for $(OBJECT_FILES), using
#     http://www.gnu.org/software/make/manual/make.html#Substitution-Refs

build: $(EXECUTABLE_FILES)

clean:
	rm -r -f $(BIN)

# http://www.gnu.org/software/make/manual/make.html#Phony-Targets
.PHONY: build clean

$(EXECUTABLE_FILES): $(OBJECT_FILES)
	@$(FC) $(FLFLAGS) -o $@ $^
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#Automatic-Variables
	@cp hamilt.inp $(BIN)
	@echo "Build successful!"

# http://www.gnu.org/software/make/manual/make.html#Static-Pattern
$(OBJECT_FILES): $(OBJ)/%.o: %.for
	@echo Compiling $<
	@mkdir -p $(@D)
	@# ^^^ http://www.gnu.org/software/make/manual/make.html#index-_0024_0028_0040D_0029
	@$(FC) -c $(FCFLAGS) $< -o $@
