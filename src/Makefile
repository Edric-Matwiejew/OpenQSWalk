#ExpoAction Makefile

#Set install directory.
INSTALL_DIR = ..

#Compiler flags.
FC = ifort
FFLAGS = -g -O3 -qopt-report=5 -vec -simd -qopenmp -shared-intel -ldl -mavx2 -xCORE-AVX2
FFLAGS2 = -qoverride-limits -I$/opt/intel/parallel_studio_xe_2019.0.045/include/ia64 
FFLAGS3 = -L$/opt/intel/parallel_studio_xe_2019.0.045/lib64
FFLAGS4 = -Iadvisor

#Created objects and modules.
OBJECTS = Nullarbor.o ExpoAction.o hb_read.o
MODS = nullarbor.mod expo.mod hb_read.mod

#Creating file structure and copying data files to directories.
PROGRAM_DIR := $(INSTALL_DIR)/ExpoAction

INPUT_DIR := $(PROGRAM_DIR)/Input
OUTPUT_DIR := $(PROGRAM_DIR)/Output
OPTIONS_DIR := $(PROGRAM_DIR)/Options
PARENT_DIRECTORIES := $(INPUT_DIR), $(OUTPUT_DIR), $(OPTIONS_DIR)

MATRICES_DIR := $(INPUT_DIR)/Matrices 
VECTOR_DIR := $(INPUT_DIR)/Vectors
LOGFILES_DIR := $(OUTPUT_DIR)/LogFiles
DATA_DIR := $(OPTIONS_DIR)/data
CHILD_DIRECTORIES := $(MATRICES_DIR) $(VECTOR_DIR) $(LOGFILES_DIR) \
											$(DATA_DIR) 
											
INSTALL_DATA = ../data

$(shell mkdir -p $(CHILD_DIRECTORIES))

$(shell cp $(INSTALL_DATA)/ExpoOptions.txt $(OPTIONS_DIR))
$(shell cp $(INSTALL_DATA)/InputList $(INPUT_DIR))
$(shell cp $(INSTALL_DATA)/Test10.dat $(VECTOR_DIR))
$(shell cp $(INSTALL_DATA)/Test10by10.dat $(MATRICES_DIR))
$(shell cp $(INSTALL_DATA)/Test800.dat $(VECTOR_DIR))
$(shell cp $(INSTALL_DATA)/Test800by800.dat $(MATRICES_DIR))
$(shell cp $(INSTALL_DATA)/theta.dat $(DATA_DIR))

#Compile ExpoAction

%.o: %.f90
	$(FC) $(FFLAGS) $(FFLAGS2) $(FFLAGS3) $(FFLAGS4) -c $< -o $@ 




ExpoAction: $(OBJECTS)
	$(FC) $(FFLAGS) $(FFLAGS2) $(FFLAGS3) $(FFLAGS4) $(OBJECTS) -o $(PROGRAM_DIR)/ExpoAction

ExpoAction.o: hb_read.o Nullarbor.o ExpoAction.f90

Nullarbor.o: hb_read.f90 Nullarbor.f90
	$(FC) $(FFLAGS) $(FFLAGS2) $(FFLAGS3) $(FFLAGS4) -c hb_read.f90 Nullarbor.f90


.PHONY: test clean veryclean

#Test, not yet functioning 
test:
	@echo "Eventually this will test things."

#Clean .o and .mod from src	
clean:
	rm ExpoAction $(OBJECTS) $(MODS)

#Uninstall ExpoAction, removing ALL files.	
veryclean:
	rm $(OBJECTS) $(MODS)
	rm -rf $(PROGRAM_DIR)
