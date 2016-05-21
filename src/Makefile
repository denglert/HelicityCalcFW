#####################################################
## -- Objects/libraries needed for main programs-- ##
#####################################################

OBJ_rh_tautau = $(CMN) ../lib/test_rh_tautau.cpp.o ../lib/HelicityTools.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o 
OBJ_rh_6f_k0_1001     = $(CMN) ../lib/HelicityTools.cpp.o ../lib/FortranUtils.f.o ../lib/test_rh_6f_k0_1001.cpp.o     ../lib/rh_6f_k0_1001.f.o
OBJ_rh_6f_k0_physical = $(CMN) ../lib/HelicityTools.cpp.o ../lib/FortranUtils.f.o ../lib/test_rh_6f_k0_physical.cpp.o ../lib/rh_6f_k0_physical.f.o

OBJ_test = $(CMN) ../lib/HelicityTools.cpp.o ../lib/main_test.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o ../lib/rh_6f_k0_1001.f.o ../lib/FortranUtils.f.o

OBJ_CUBA = $(CMN) ../lib/test_cuba.cpp.o ../lib/vegas.f.o

OBJ_test_HelicityTools = $(CMN) ../lib/test_HelicityTools.cpp.o ../lib/HelicityTools.cpp.o ../lib/rh_tautau.f.o

OBJ_test_FortranUtils = ../lib/TestFortranUtils.f.o ../lib/FortranUtils.f.o

OBJ_test_PhaseSpaceTools = $(CMN) ../lib/test_PhaseSpaceTools.cpp.o ../lib/PhaseSpaceTools.cpp.o ../lib/UtilFunctions.cpp.o
OBJ_test_PhaseSpaceTools2 = $(CMN) ../lib/test_PhaseSpaceTools2.cpp.o ../lib/PhaseSpaceTools.cpp.o ../lib/UtilFunctions.cpp.o
OBJ_test_PhaseSpaceTools3 = $(CMN) ../lib/test_PhaseSpaceTools3.cpp.o ../lib/PhaseSpaceTools.cpp.o ../lib/UtilFunctions.cpp.o

OBJ_test_PhaseSpaceIntegration = $(CMN) ../lib/test_PhaseSpaceIntegration.cpp.o ../lib/PhaseSpaceTools.cpp.o
OBJ_test_PhaseSpaceIntegration2 = $(CMN) ../lib/test_PhaseSpaceIntegration2.cpp.o ../lib/PhaseSpaceTools.cpp.o
OBJ_test_PhaseSpaceIntegration3 = $(CMN) ../lib/test_PhaseSpaceIntegration3.cpp.o ../lib/PhaseSpaceTools.cpp.o
OBJ_test_PhaseSpaceIntegration4 = $(CMN) ../lib/test_PhaseSpaceIntegration4.cpp.o ../lib/PhaseSpaceTools.cpp.o

OBJ_Scalar_2Fermion_6Fermion_SingleBranchIntegration = $(CMN)  ../lib/PhaseSpaceTools.cpp.o ../lib/rh_6f_k0_1001.f.o ../lib/HelicityTools.cpp.o ../lib/FortranUtils.f.o ../lib/Scalar_2Fermion_6Fermion_SingleBranchIntegration.cpp.o

OBJ_Scalar_2Fermion_6Fermion_FullIntegration = $(CMN)  ../lib/PhaseSpaceTools.cpp.o ../lib/rh_6f_k0_1001.f.o ../lib/HelicityTools.cpp.o ../lib/FortranUtils.f.o ../lib/Scalar_2Fermion_6Fermion_FullIntegration.cpp.o

OBJ_Analyzer-h-tautau-6f = $(CMN) ../lib/PhaseSpaceTools.cpp.o ../lib/rh_6f_k0_1001.f.o ../lib/HelicityTools.cpp.o ../lib/FortranUtils.f.o ../lib/Analyzer-h-tautau-6f.cpp.o

OBJ_LifeTime_Muon_Polarized = $(CMN) ../lib/LifeTime_Muon_Polarized.cpp.o ../lib/PhaseSpaceTools.cpp.o

OBJ_SpinPolarization = $(CMN) ../lib/SpinPolarization.cpp.o ../lib/PhaseSpaceTools.cpp.o

CMN = ../lib/UtilFunctions.cpp.o 


#########################
## -- Main programs -- ##
#########################
## Linking
LINKER = ifort

binaries : ../bin/test_rh_tautau ../bin/test_rh_6f_k0_1001 ../bin/test_rh_6f_k0_physical ../bin/test_cuba ../bin/main_test ../bin/test_HelicityTools ../bin/test_FortranUtils ../bin/test_PhaseSpaceTools ../bin/test_PhaseSpaceIntegration ../bin/test_PhaseSpaceIntegration2 ../bin/test_PhaseSpaceIntegration3 ../bin/test_PhaseSpaceIntegration4 ../bin/Scalar_2Fermion_6Fermion_SingleBranchIntegration ../bin/LifeTime_Muon_Polarized ../bin/SpinPolarization ../bin/test_PhaseSpaceTools2 ../bin/test_PhaseSpaceTools3  ../bin/Scalar_2Fermion_6Fermion_FullIntegration ../bin/Analyzer-h-tautau-6f

../bin/test_rh_tautau : $(OBJ_rh_tautau)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_rh_6f_k0_physical : $(OBJ_rh_6f_k0_physical)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_rh_6f_k0_1001 : $(OBJ_rh_6f_k0_1001)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_HelicityTools: $(OBJ_test_HelicityTools)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceTools: $(OBJ_test_PhaseSpaceTools)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceTools2: $(OBJ_test_PhaseSpaceTools2)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceTools3: $(OBJ_test_PhaseSpaceTools3)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/LifeTime_Muon_Polarized: $(OBJ_LifeTime_Muon_Polarized)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/SpinPolarization: $(OBJ_SpinPolarization)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceIntegration: $(OBJ_test_PhaseSpaceIntegration)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceIntegration2: $(OBJ_test_PhaseSpaceIntegration2)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceIntegration3: $(OBJ_test_PhaseSpaceIntegration3)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/test_PhaseSpaceIntegration4: $(OBJ_test_PhaseSpaceIntegration4)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/Scalar_2Fermion_6Fermion_SingleBranchIntegration: $(OBJ_Scalar_2Fermion_6Fermion_SingleBranchIntegration)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/Scalar_2Fermion_6Fermion_FullIntegration: $(OBJ_Scalar_2Fermion_6Fermion_FullIntegration)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/Analyzer-h-tautau-6f: $(OBJ_Analyzer-h-tautau-6f)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

../bin/test_FortranUtils: $(OBJ_test_FortranUtils)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/main_test : $(OBJ_test)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_cuba : $(OBJ_CUBA)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) $(CUBALFLAGS)
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) $(CUBALFLAGS) -o $@ $^
endif

#########################################################

####################
## -- Fortran -- ###
####################

#FC    = gfortran
FC     = ifort

FCFLAGS = -c

#FLFLAGS = -cxxlib -nofor_main -\#
FLFLAGS = -cxxlib -nofor_main $(ROOTLFLAGS)

../lib/%.f.o : %.f
	$(FC) $(FCFLAGS) $< -o $@ 

#########################################################

#################
### -- C++ -- ###
#################

CC     = g++
#CC     = icc

DEPDIR := .d
$(shell mkdir -p $(DEPDIR) >/dev/null)
CCCUSTOMFLAGS = -I../inc/
CCGENERALFLAGS = -I../inc/  $(DEPFLAGS)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
CLFLAGS = $(ROOTLFLAGS)
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

# ROOT flags can also be set to environment variables
ROOTCFLAGS = `root-config --cflags`
ROOTLFLAGS = -L/home/de3u14/lib/build/hep/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -limf -lm -ldl

# CUBAFLAGS flags can also be set to environment variables
CUBACFLAGS = -I/home/de3u14/lib/build/hep/Cuba/Cuba-4.2
CUBALFLAGS = -L/home/de3u14/lib/build/hep/Cuba/Cuba-4.2 -lcuba

#####################
### Compilation rule for specific .cpp sources

../lib/test_cuba.cpp.o : test_cuba.cpp
	$(CC) $(ROOTCFLAGS) $(CCCUSTOMFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/test_PhaseSpaceIntegration.cpp.o : test_PhaseSpaceIntegration.cpp
	$(CC) $(ROOTCFLAGS) $(CCCUSTOMFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/test_PhaseSpaceIntegration2.cpp.o : test_PhaseSpaceIntegration2.cpp
	$(CC) $(ROOTCFLAGS) $(CCCUSTOMFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/test_PhaseSpaceIntegration3.cpp.o : test_PhaseSpaceIntegration3.cpp
	$(CC) $(ROOTCFLAGS) $(CCCUSTOMFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/test_PhaseSpaceIntegration4.cpp.o : test_PhaseSpaceIntegration4.cpp
	$(CC) $(ROOTCFLAGS) $(CCCUSTOMFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/Scalar_2Fermion_6Fermion_SingleBranchIntegration.cpp.o : Scalar_2Fermion_6Fermion_SingleBranchIntegration.cpp
	$(CC) $(ROOTCFLAGS) $(CCCUSTOMFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/Scalar_2Fermion_6Fermion_FullIntegration.cpp.o : Scalar_2Fermion_6Fermion_FullIntegration.cpp
	$(CC) $(ROOTCFLAGS) $(CCFLAGS) $(CUBACFLAGS) -c $< -o $@

../lib/LifeTime_Muon_Polarized.cpp.o : LifeTime_Muon_Polarized.cpp
	$(CC) $(ROOTCFLAGS) $(CCFLAGS) $(CUBACFLAGS) -c $< -o $@

### General compilation rule for .cpp sources
#%.cpp.o : %.cpp
../lib/%.cpp.o : %.cpp
	$(CC) $(ROOTCFLAGS) $(CCGENERALFLAGS) -c $< -o $@
	$(POSTCOMPILE)

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

SRCS = Analyzer-h-tautau-6f.cpp

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))$(CC) $(ROOTCFLAGS) $(CCFLAGS) -c $< -o $@
