FC     = mpiifort  -fc=ifort  -r8 -i4  -O3 -xAVX  -mcmodel=medium -qmkl=sequential

OBJ     =  main.o solar.o biofouling.o temperature.o salinity.o chlorophyll.o gauss_rng2.o
main: $(OBJ)
	$(FC)  $(OBJ) $(FFLAGS) -o main

%.o      : %.f90
	$(FC)  -c $*.f90 $(FFLAGS)


clean    : 
	rm -rf main *.o *.mod 
