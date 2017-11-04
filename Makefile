SHELL = /bin/sh
# CXX= icpc -c 
# We recommend using "make -D DEBUG" and "make -D LOCALFENK" instead of following options
# Compile sfbox with additional debug code in "#if DEBUG" blocks
# DEBUG = 1               
# Use ~/lib and ~/include instead of normal fenklib locations
# LOCALFENK = 1  

#CXX= icpc -c # if this line is enabled the intel c++ compiler is used instead of Gcc.

Fabien = 0
LAPACK = 1 

NVCC := nvcc 

DEBUG_STR =              
ifdef DEBUG
    DEBUG_STR = -DDEBUG
endif

ifndef LOCALFENK
	INC_DIR = /usr/local/include
	LIB_DIR = /usr/local/lib
else
    INC_DIR = $(HOME)/include
    LIB_DIR = $(HOME)/lib
endif

ifdef Johan
	INC_DIR = $(HOME)/sfbox/trunk/src
	LIB_DIR = $(HOME)/sfbox/trunk/src
endif

ifdef Fabien
	#INC_DIR = /usr/local/clapack-3.2.1/include
	#LIB_DIR = /usr/local/clapack-3.2.1/lib
	MATRIXLIB = -llapack -lblas -lf2c -lm -lpthread
endif

ifndef Fabien
	MATRIXLIB = -llapack -lblas -lf2c -lm -lpthread
endif

LIBS=-lfenk  

CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -Wall
ifndef DEBUG
CPPFLAGS+= -O3
endif

ifdef OMP
    LIBS = -lfenk_omp -lgomp -lpthread  
    CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -g -fopenmp -Wall
    ifndef DEBUG
	CPPFLAGS+= -O2
    endif
endif

ifdef MAC
    LIBS = -lfenk
    CPPFLAGS = -Dunix -Daux -I$(INC_DIR) -g $(DEBUG_STR)
    ifndef DEBUG
        CPPFLAGS+= -O3
    endif
endif

ifdef i7
    LIBS = -lfenk_omp -lgomp -lpthread
    CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -g -fopenmp
    ifndef DEBUG
        CPPFLAGS+= -O3
    endif
endif

ifdef BLAS
    LIBS = -lfenk_omp -lgomp -lpthread -lblas 
    CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -DHAS_BLAS -g -fopenmp
    ifndef DEBUG
        CPPFLAGS+= -O3
    endif
endif

ifdef NoOpenMP
    LIBS = -lfenk_omp -lgomp -lpthread -lblas
    CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -DHAS_BLAS  
    ifndef DEBUG
        CPPFLAGS+= -O3
    endif
endif

ifdef NoOpenMPNoBLAS
    LIBS = -lfenk_omp -lgomp -lpthread 
    CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR)  
    ifndef DEBUG
        CPPFLAGS+= -O3
    endif
endif

ifdef LAPACK
#	CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -O3 -Wall -DHAS_LAPACK -DNO_BLAS_WRAP 
	CPPFLAGS = -I$(INC_DIR) $(DEBUG_STR) -g -O3 -Wall -DHAS_LAPACK -DNO_BLAS_WRAP
endif
CFLAGS = $(CPPFLAGS)

OBJECTS = sfnewton.o  sf_vectoriterations.o newttool.o misc.o Vary.o SegmentParam.o \
tools.o SF_MolList2ndO.o SF_Homopol2ndO.o SF_Monomer2ndO.o Lat3D.o Lat3DFCC1stO.o\
Lat3DHEX1stO.o Lat3D1stO.o  LatticeRange.o LatticeRange3D.o Lat3D2ndO.o Lat2D2ndO.o Lat1D2ndO.o\
SF_Copol2ndO.o SF_Copol2ndO_stiff_range.o Lat1DSphere2ndO.o Lat1DCyl2ndO.o Lat2DCyl2ndO.o\
SF_SystemFraaije.o SF_System3rdGen.o SF_System.o SF_SystemMC.o SF_SystemLD.o\
SF_StateRatioList.o SF_StateRatio.o SF_State.o SF_SolveFraaije.o SF_SolveCopel.o\
SF_SolvePikar.o SF_SolveCG.o SF_SolvePhiU.o SF_Solve_Johan.o\
SF_Solve3rdGen.o SF_Solve.o SF_SegmentList.o SF_SegmentBlock.o \
SF_Segment.o SF_ReactionList.o SF_Reaction.o SF_Monomer1stO.o \
SF_MoleculeList.o SF_Molecule.o SF_MolStructure.o SF_MolState.o \
SF_MolSegment.o SF_MolList1stO.o SF_Homopol1stO.o SF_HomoRing1stO.o SF_Branched1stO.o SF_Dend1stO.o \
SF_AsymDend1stO.o SF_Comb1stO.o SF_Comb2ndO.o\
SF_Copol1stO.o SF_Box.o OutputLine.o Output.o NoArtefact.o Message.o \
LatticeRange2D.o LatticeRange1D.o LatticeRangeFile.o Lat2DFlat1stOSten.o \
Lat2DFlat1stOS.o Lat2DFlat1stO.o Lat2DFlat.o Lat2DCylinderSten.o Lat2DCylinder.o \
Lat2DCyl1stOSten.o Lat2DCyl1stOStenS.o Lat2DCyl1stOS.o Lat2DCyl1stO.o \
Lat1DSphere1stOS.o Lat1DSphere1stO.o Lat1DSphere.o Lat1DFlat1stOS.o \
Lat1DFlat1stO.o Lat1DFlat.o Lat1DCylinder.o Lat1DCyl1stOS.o \
Lat1DCyl1stO.o Input.o SF_LinkNode.o SF_ParserBranched.o   

ifdef CUDA
    LIBS += -L/usr/include -lcuda -lcudart -lcublas
    CPPFLAGS += -DCUDA
    NVCCFLAGS += -DCUDA -arch sm_13  
    OBJECTS +=Cuda_tools.o
endif

%.o: %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

ifdef LAPACK
     LIBS +=  $(MATRIXLIB)
endif

sfbox: $(OBJECTS)
#	icpc  -o sfbox $(OBJECTS) -L$(LIB_DIR) $(LIBS) 
	g++ -o sfbox $(OBJECTS) -L$(LIB_DIR) $(LIBS)  
plugin: $(OBJECTS)
	gcc -L$(LIB_DIR) -o "c:/Program Files/Dullware Houston/Plugins/sfbox" $(OBJECTS) $(LIBS)

solaris: $(OBJECTS)
	gcc -L/usr/local/lib -o sfbox $(OBJECTS) -lfenk -lnsl -lsocket
	
Slope: Slope.o
	g++ -o Slope Slope.o -L$(LIB_DIR) $(LIBS)

Box : Box.o
	g++ -o Box Box.o -L$(LIB_DIR) $(LIBS)

Test : Test.o
	g++ -o Test Test.o -L$(LIB_DIR) $(LIBS)

	
InGuess: InGuess.o
	g++ -o InGuess InGuess.o -L$(LIB_DIR) $(LIBS)
	
Duplicate: Duplicate.o
	g++ -o Duplicate Duplicate.o -L$(LIB_DIR) $(LIBS)
	
Fine: Fine1D.o
	g++ -o Fine1D Fine1D.o -L$(LIB_DIR) $(LIBS)
	
Tubes: Tubes.o
	g++ -o Tubes Tubes.o -L$(LIB_DIR) $(LIBS)
	
clean:
	rm -f $(OBJECTS)

Clean:
	rm -f $(OBJECTS)
	rm *.o
