#This is a makefile
FC = gfortran
FFLAGS = -O2 -fopenmp

OBJECTS = Nullarbor.o ExpoAction.o
MODS = nullarbor.mod expo.mod

PARENTDIRS = Input Output 
CHILDDIRS = Input/Matrices Input/Vectors Output/LogFiles

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ 

ExpoAction: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o ExpoAction

ExpoAction.o: Nullarbor.o ExpoAction.f90

Nullarbor.o: Nullarbor.f90

$(shell mkdir -p $(PARENTDIRS) $(CHILDDIRS))

.PHONY: test clean veryclean

test:
	@echo "Eventually this will test things."
	

clean:
	rm ExpoAction $(OBJECTS) $(MODS)
	
veryclean:
	rm ExpoAction $(OBJECTS) $(MODS)
	rm -rf $(PARENTDIRS)
