OBJ_rh_tautau = $(CMN) ../lib/test_rh_tautau.cpp.o ../lib/HelicityTools.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o ../lib/rh_6f.f.o

OBJ_test = $(CMN) ../lib/HelicityTools.cpp.o ../lib/main_test.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o ../lib/rh_6f.f.o

OBJ_rh_6f = $(CMN) ../lib/HelicityTools.cpp.o ../lib/test_rh_6f.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o ../lib/rh_6f.f.o

OBJ_CUBA = $(CMN) ../lib/main_cuba_test.cpp.o ../lib/vegas.f.o

CMN = ../lib/UtilFunctions.cpp.o 

#########################
## -- Main programs -- ##
#########################
LINKER = ifort

binaries : ../bin/test_rh_tautau ../bin/test_rh_6f ../bin/main_cuba_test ../bin/main_test

../bin/test_rh_tautau : $(OBJ_rh_tautau)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

../bin/test_rh_6f : $(OBJ_rh_6f)
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

../bin/main_cuba_test : $(OBJ_CUBA)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $^ $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $^
endif

#########################################################

####################
## -- Fortran -- ###
####################

#FC    = gfortran
FC     = ifort

FCFLAGS = -c

#FLFLAGS = -cxxlib -nofor_main -\#
FLFLAGS = -cxxlib -nofor_main $(ROOTLFLAGS) $(CUBALFLAGS)

../lib/%.f.o : %.f
	$(FC) $(FCFLAGS) $< -o $@ 

#################
### -- C++ -- ###
#################

CC     = g++
#CC     = icc

CCFLAGS = $(CUBACFLAGS) -I../inc/ -MMD -MF .depend_cpp
CLFLAGS = $(CUBALFLAGS) $(ROOTLFLAGS)

# ROOT flags can also be set to environment variables
ROOTCFLAGS = `root-config --cflags`
ROOTLFLAGS = -L/home/de3u14/lib/build/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -limf -lm -ldl

# CUBAFLAGS flags can also be set to environment variables
CUBACFLAGS = -I/home/de3u14/lib/build/Cuba/Cuba-4.2
CUBALFLAGS = -L/home/de3u14/lib/build/Cuba/Cuba-4.2 -lcuba

../lib/%.cpp.o : %.cpp
	$(CC) $(ROOTCFLAGS) $(CCFLAGS) -c $< -o $@

include .depend_cpp
