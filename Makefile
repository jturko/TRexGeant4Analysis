
# ROOT
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs) -lMinuit
ROOTINC      := -I$(shell root-config --incdir)

LIB_DIR = $(HOME)/lib

#CLHEP_LIBS = -I /usr/local/opt/geant/CLHEP/include -l/usr/local/opt/geant/CLHEP/lib/libCLHEP.a

ALLIBS = $(ROOTLIBS) $(ROOTGLIBS) $(CLHEP_LIBS) -lSpectrum

CPP         = g++
CFLAGS		= -g -fPIC $(ROOTCFLAGS) -Wall  -Wno-write-strings -DSIMULATION_PATH=\"$(PWD)\" -std=c++11

LFLAGS		= -g -fPIC
LIBS 		= $(ALLIBS)

COMM_DIR = $(HOME)/CommandLineInterface 
SIM_DIR = ../include

INCLUDE_DIRS = $(SIM_DIR) $(COMM_DIR) $(PWD) $(PWD)/include $(G4INCLUDE)

FILES = Reconstruction Nucleus Compound HitSim Particle Kinematics ParticleMC ReconstructSimDictionary TRexSettings Germanium

O_FILES = $(addsuffix .o, $(addprefix .build/, $(FILES) ) )

DEPENDENCIES = include/Particle.hh \
	../include/ParticleMC.hh \
	../include/Germanium.hh \
	../include/TRexSettings.hh \
	include/RootLinkDef.h

TOTAL_LIBS = $(ALLIBS) -L$(LIB_DIR) -L$(G4TMP)/$(G4SYSTEM)/TRexGeant4 -lCommandLineInterface

INCLUDES =  $(addprefix -I,$(INCLUDE_DIRS))

.PHONY: all clean

all: $(O_FILES) ReconstructSim libReconstructSim.so 

#===================================================================
#generation of root dictionary
#===================================================================

.build/ReconstructSimDictionary.cc: $(DEPENDENCIES) 
	rootcint -f $@ -c $(DEPENDENCIES) 

.build/ReconstructSimDictionary.o: .build/ReconstructSimDictionary.cc
	$(CPP) -fPIC $(CFLAGS) $(INCLUDES) -c $< -o $@

#==================================================================
# generation of object files
#==================================================================

.build/%.o : src/%.cc include/%.hh
	@mkdir -p $(dir $@)
	$(CPP) -fPIC $(CFLAGS) $(INCLUDES)  -c $< -o $@

.build/TRexSettings.o: ../src/TRexSettings.cc ../include/TRexSettings.hh
	@mkdir -p $(dir $@)
	$(CPP) -fPIC $(CFLAGS) $(INCLUDES)  -c $< -o $@

.build/ParticleMC.o: ../src/ParticleMC.cc ../include/ParticleMC.hh
	@mkdir -p $(dir $@)
	$(CPP) -fPIC $(CFLAGS) $(INCLUDES)  -c $< -o $@

.build/Germanium.o: ../src/Germanium.cc ../include/Germanium.hh
	@mkdir -p $(dir $@)
	$(CPP) -fPIC $(CFLAGS) $(INCLUDES) -c $< -o $@

#==================================================================
# generation of executables
#==================================================================
ReconstructSim : ReconstructSim.cc $(O_FILES)
	@echo building reco class object
	$(CPP) -fPIC $(CFLAGS) $(INCLUDES)  -c $< -o .build/ReconstructSim.o
	@echo building executable
	$(CPP) $(CFLAGS) $(INCLUDES)  $(TOTAL_LIBS) $(O_FILES) .build/ReconstructSim.o -o $@
	cp ReconstructSim $(HOME)/bin
	@echo Done

#=================================================================
# make .so file
#=================================================================
libReconstructSim.so : $(O_FILES)
	@echo Creating shared library $@
	@$(CXX) -shared -o $(LIB_DIR)/$@  $(O_FILES)   $(CFLAGS) $(INCLUDES)  $(LIBS)


#==================================================================
# clean
#==================================================================

clean : 
	@rm -f .build/* *.o DictionarySim.* libReconstructSim.so ReconstructSim 
