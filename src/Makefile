OBJ = ../lib/main.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o ../lib/UtilFunctions.cpp.o ../lib/rh_6f.f.o

OBJ02 = ../lib/main02.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o ../lib/vegas.f.o ../lib/UtilFunctions.cpp.o ../lib/rh_6f.f.o

##########################################################
LINKER = ifort

../bin/main : $(OBJ)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $(OBJ) $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $(OBJ) 
endif

../bin/main02 : $(OBJ02)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $(OBJ02) $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $(OBJ02) 
endif


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

#################
### -- C++ -- ###
#################

CC     = g++
#CC     = icc

CCFLAGS = -I../inc/ -MMD -MF .depend_cpp
CLFLAGS = $(ROOTLFLAGS)

# ROOT flags can also be set to environment variables
ROOTCFLAGS = `root-config --cflags`
ROOTLFLAGS = -L/home/de3u14/lib/build/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -limf -lm -ldl

# CUBAFLAGS flags can also be set to environment variables
CUBACFLAGS = -I/home/de3u14/lib/build/Cuba/Cuba-4.2
CUBALFLAGS = -L/home/de3u14/lib/build/Cuba/Cuba-4.2 -lcuba -m

../lib/%.cpp.o : %.cpp
	$(CC) $(ROOTCFLAGS) $(CCFLAGS) -c $< -o $@
include .depend_cpp
