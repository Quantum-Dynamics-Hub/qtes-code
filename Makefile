# makefile
objects = qtes.o cqf.o dumpef.o read_input_mpi.o set_traj.o ./dftb/libdftb.a
FC = mpif90
switch = -mkl=sequential
EXEC = qtes.x
LIB = -openmp

$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects)
%.o: %.f
	$(FC) -c $(LIB) $(switch) $< 

allclean:	  
	rm -rf *.o qtes.x *.log core.* fort.* data/*
