FC = gfortran
SRCDIR = src
OBJDIR = obj

TARGET = a.out

FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -fcheck=bounds -fallow-argument-mismatch -funroll-loops -std=f2008 -Wall -Wextra -Wpedantic -fopenmp -I$(OBJDIR) -J$(OBJDIR)

#----FFTW library------------------- (adjust prefix if needed)--------------
# For personal Linux machine
FFTW_PREFIX := $(shell pkg-config --variable=prefix fftw3)

# For personal MAC machine
#FFTW_PREFIX := $(shell brew --prefix fftw)

#--- End prefix modification -----------------------------------------------

FFLAGS += -I$(FFTW_PREFIX)/include 
LIBS = -L$(FFTW_PREFIX)/lib  -lfftw3_omp -lfftw3 -lm

#----- END FFTW library-------------------

# Source files
SRC =	$(SRCDIR)/main.f90 \
		$(SRCDIR)/fileIO.f90 \

MOD_FILES = $(SRCDIR)/geom.f90 \
			$(SRCDIR)/delta.f90 \
			$(SRCDIR)/my_fftw3.f90 \

# Object files
OBJ = $(addprefix $(OBJDIR)/, $(notdir $(SRC:.f90=.o)))
MOBJ = $(addprefix $(OBJDIR)/, $(notdir $(MOD_FILES:.f90=.o)))

# Default target
all: $(OBJDIR) $(TARGET)

$(OBJDIR):
	mkdir -p $(OBJDIR)

# Pattern rule to compile .f90 files into .o files in SRCDIR
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Linking the final executable
$(TARGET): $(MOBJ) $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(MOBJ) $(OBJ) $(LIBS)

# Clean up generated files
clean:
	@rm -f $(OBJDIR)/*.o $(TARGET)
	@rm -f $(OBJDIR)/*.mod

veryclean: clean
	@rm -rf $(OBJDIR)
	@rm  output*

.PHONY : clean veryclean