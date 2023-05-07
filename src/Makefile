# Paths definition
IDIR =../include
BDIR =../bin
ODIR =./obj
MDIR =./mod
SDIR =.
DIR  = ~/mis_bin
#libpath=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so \
	/usr/lib/x86_64-linux-gnu/blas/libblas.so
libpath = -lblas -llapack

# Compiler/Linker settings
FC = gfortran -I$(IDIR) -J$(MDIR)  # -Wno-argument-mismatch 
FLFLAGS = -g
FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -fcheck=bounds -fmax-errors=5
PROGRAM = computherm
SRCF90  = g16_computherm.f90
SRCMOD  = parameters.f90   \
          lineshape.f90    \
          utils.f90        \
          sorting.f90      \
          mathtools.f90    \
          wrap_lapack.f90  \
          inertia.f90      \
          geometry.f90     \
          datatypes.f90    \
          g16files.f90     \
          input.f90        \
          statmech.f90     \
          superpose.f90    \
          screening.f90

# Find all source files, create a list of corresponding object files
OBJ     = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90))
OBJ    := $(OBJ) $(patsubst %.F90,$(ODIR)/%.o,$(SRC))
MODOBJ  = $(patsubst %.f90,$(ODIR)/%.o,$(SRCMOD)) 
MODOBJ := $(patsubst %.F90,$(ODIR)/%.o,$(MODOBJ))
MODMOD  = $(patsubst $(ODIR)/%.o,$(MDIR)/%.mod,$(MODOBJ))

dir_guard=@mkdir -p $(@D)

default: $(MODOBJ) $(OBJ)
	@mkdir -p $(BDIR)
	$(FC) $(MODOBJ) $(OBJ) -o $(BDIR)/$(PROGRAM) $(libpath)

# Compiler steps for all objects
$(ODIR)/%.o: %.F90
	@mkdir -p $(MDIR)
	$(dir_guard)
	$(FC) $(FCFLAGS) -o $@ $^ 

$(ODIR)/%.o: %.f90
	@mkdir -p $(MDIR)
	$(dir_guard)
	$(FC) $(FCFLAGS) -o $@ $^ 

# Linker
$(BDIR)/$(PROGRAM) : $(OBJ)
	$(dir_guard)
	$(FC) $(FLFLAGS) -o $@ $< $(libpath) 

# Dependencies : Main program depends on all modules
$(ODIR)/$(PRGOBJ) : $(MODOBJ)

debug: 

	@echo "SRCF90  = $(SRCF90)"
	@echo "OBJ     = $(OBJ)"
	@echo "SRCMOD  = $(MOD)"
	@echo "MODOBJ  = $(MODOBJ)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRGOBJ  = $(PRGOBJ)"
	@echo "MODMOD  = $(MODMOD)"
#	DEBUG=1 make clean $(BDIR)/$(PROGRAM)

install:
	cp $(BDIR)/$(PROGRAM) $(DIR)

clean:
	rm -f $(ODIR)/*.o $(MDIR)/*mod $(BDIR)/$(PROGRAM) $(BDIR)/debug

