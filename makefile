#This is a makefile
FC = gfortran
FFLAGS = -O2 -fopenmp

ExpoAction: Expo.o Nullarbor.o
	$(FC) $(FFLAGS) Nullarbor.o Expo.o -o ../ExpoAction

Expo.o: Nullarbor.o ExpoAction_1.0.1.f90
	$(FC) $(FFLAGS) -c ExpoAction_1.0.1.f90 -o Expo.o

Nullarbor.o: Nullarbor_1.0.1.f90
	$(FC) $(FFLAGS) -c Nullarbor_1.0.1.f90 -o Nullarbor.o

.PHONY: clean 

clean:
	rm *.o *.mod
