LINKER = ifort

#FC    = gfortran
FC     = ifort
FCFLAGS = -c
#FLFLAGS = -cxxlib -nofor_main -\#
FLFLAGS = -cxxlib -nofor_main $(ROOTLFLAGS) 
CC     = g++
#CC     = icc
CCFLAGS =
CLFLAGS = $(ROOTLFLAGS)

ROOTCFLAGS = `root-config --cflags`
ROOTLFLAGS = -L/home/de3u14/lib/build/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -limf -lm -ldl

OBJ = ../lib/main.cpp.o ../lib/rh_tautau.f.o ../lib/rambo.f.o ../lib/ranmar.f.o

##########################################################

../bin/main : $(OBJ)
ifeq ($(LINKER),ifort)
	$(FC) -o $@ $(OBJ) $(FLFLAGS) 
endif
ifeq ($(LINKER),g++)
	$(CC) $(CLFLAGS) -o $@ $(OBJ) 
endif


../lib/%.cpp.o : %.cpp
	$(CC) $(ROOTCFLAGS) $(CCFLAGS) -c $< -o $@ 


../lib/%.f.o : %.f
	$(FC) $(FCFLAGS) $< -o $@ 
