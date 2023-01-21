
FC = gfortran #COMPILER
FFLAGS = -g  
SRC = general_mod.f90 LJ_mod.f90 ST_mod.f90 main.f90

OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

program.exe: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.mod *.o *.xyz *.log program.exe

