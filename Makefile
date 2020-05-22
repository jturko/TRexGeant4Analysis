
# ROOT
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs) -lMinuit
ROOTINC      := $(shell root-config --incdir)

LIB_DIR = $(HOME)/lib

#CLHEP_LIBS = -I /usr/local/opt/geant/CLHEP/include -l/usr/local/opt/geant/CLHEP/lib/libCLHEP.a

ALLIBS = $(CLHEP_LIBS) $(ROOTGLIBS) -lSpectrum 

CPP         = g++
CFLAGS		= -g -fPIC -Wall -Wno-write-strings -Wl,--no-as-needed -DSIMULATION_PATH=\"$(PWD)\" -std=c++11 $(ROOTCFLAGS)

LFLAGS		= -g -fPIC
LIBS 		= $(ALLIBS)

COMM_DIR = $(HOME)/programs/CommandLineInterface 
SIM_DIR = ../include

INCLUDE_DIRS = $(COMM_DIR) $(PWD) $(PWD)/include $(SIM_DIR) $(G4INCLUDE)

FILES = Reconstruction Nucleus Compound HitSim Particle Kinematics ParticleMC ReconstructSimDictionary TRexSettings Germanium

O_FILES = $(addsuffix .o, $(addprefix .build/, $(FILES) ) )

DEPENDENCIES = include/Particle.hh \
	../include/ParticleMC.hh \
	../include/Germanium.hh \
	../include/TRexSettings.hh \
	include/RootLinkDef.h

TOTAL_LIBS = -L$(G4WORKDIR)/bin/$(G4SYSTEM)/TRexGeant4 -L$(LIB_DIR) -lCommandLineInterface $(ALLIBS)

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
