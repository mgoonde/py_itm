# -*- Makefile -*-


F90 := gfortran
LIB   := libitm

IRALIB := ../IRA/src/libira.so


ifeq ($(F90),)
  $(error " FORTRAN COMPILER NOT DEFINED!!")
endif

# -- gfortran flags
DL := -ldl -fPIC
FFLAGS := -Wall -O0 -cpp $(DL) -fdefault-real-8 -fopenmp
FDEBUG := -g -fcheck=all -fbacktrace  #-DDEBUG


FLAGS  := $(FFLAGS) $(FDEBUG)

# ---------------------------------------------------------------------------

SRC  := .
OBJ  := $(SRC)/Obj
LIBD := $(SRC)/../lib

#..............................................


MOD := -I$(OBJ) -J$(OBJ)



module := $(sort $(wildcard $(SRC)/m_*.f90))
modulo := $(module:$(SRC)%.f90=$(OBJ)%.o)

apifile := $(SRC)/itm.f90
apifile.o := $(apifile:$(SRC)/%.f90=$(OBJ)/%.o)

# ---------------------------------------------------------------------------

default: lib
pp:
	@echo $(routine)

#gencat: file test $(OBJ) $(BIN) $(EXE)
lib: $(OBJ) $(LIBD) $(LIB)



# ----------------------------------------------------------------------------
# create folders

$(OBJ):
	@if [ ! -d $(OBJ) ]; then mkdir $(OBJ) ; fi
$(LIBD):
	@if [ ! -d $(LIBD) ]; then mkdir $(LIBD) ; fi
$(BIN):
	@if [ ! -d $(BIN) ]; then mkdir $(BIN) ; fi


# ----------------------------------------------------------------------------
# rules for each

$(LIB): $(modulo) $(apifile.o) $(nautyfile.o)
	${F90} -o $(LIBD)/$@.so -shared -fPIC $^ -ldl $(IRALIB)



# ...............................................................................
# object files

$(OBJ)/%.o: $(SRC)/%.f90
	$(F90) $(MOD)  -c $(FLAGS)  $< -o $@


# --- SYSTEM INSTALLATION::
test:
	@echo "*module    ::" $(module)
	@echo "*modulo    ::" $(modulo)
	@echo "*header    ::" $(header)
	@echo "*header.o  ::" $(header.o)
	@echo "*style     ::" $(style)
	@echo "*style.o   ::" $(style.o)
	@echo "*routine   ::" $(routine)
	@echo "*routine.o ::" $(routine.o)


# --- CLEANUP
clean:
	rm -rf *.mod $(OBJ) $(EXE) $(BIN) $(LIB) $(LIBD)


# --- DEPENDENCIES
include ./make.depend





