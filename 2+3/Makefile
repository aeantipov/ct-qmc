# ct-qmc Makefile example
# Patch the ct-qmc.cpp file and compile it with a 'make' (MPI-version) or 'make ct-qmc' (1-proc version)
# this directory can be placed elsewhere, check that QMCDIR points to the directory where the sources of ct-qmc exist, e.g. where the main_mpi.cpp / main.cpp .

CPP= g++ 
MPICPP = mpic++ 
QMCDIR=..

PART= 2
#This flag determines how many particles are to include in QMC calc. Affects the size of all matrices and henceworth its increasing leads to substantial growing of computation time

ZONE= 2
#This flag determines how many zones are to include in QMC calc. A zone differs from an orbital by a definition that all non-diagonal correlators are equal to 0

MPI = -Duse_mpi
# This flag determines whether to use mpi or not

PARTFLAG = -Dorb${PART}
ZONEFLAG= -Dzone${ZONE}
CFLAGS = $(PARTFLAG) $(ZONEFLAG) -I${QMCDIR}
CXXFLAGS = ${CFLAGS} -O3 -Wno-non-template-friend -Wno-write-strings

BINNAME=ct-qmc

default: ct-qmc-mpi

minictqmcworld.o : ${QMCDIR}/minictqmcworld.cpp ${QMCDIR}/minictqmcworld.h
	$(MPICPP) $(CXXFLAGS) $(INCLUDE) -c ${QMCDIR}/minictqmcworld.cpp 

ct-qmc: ct-qmc.cpp ${QMCDIR}/main.cpp ${QMCDIR}/ini.cpp
	$(CPP) $(CXXFLAGS) ct-qmc.cpp -o ${BINNAME}${PART}${ZONE}

ct-qmc-mpi: ct-qmc.cpp ${QMCDIR}/main_mpi.cpp ${QMCDIR}/ini.cpp ${QMCDIR}/input.cpp minictqmcworld.o 
	$(MPICPP) $(CXXFLAGS) $(MPI) ct-qmc.cpp minictqmcworld.o -o ${BINNAME}-mpi${PART}${ZONE}

change : ${QMCDIR}/change.cpp
	$(CPP) ${QMCDIR}/change.cpp -o change

clean:
	${RM} ct-qmc *.o

cleandata:
	cp _input.dat _input.dat.backup
	rm -f *.dat
	rm -f *.it*
	mv _input.dat.backup _input.dat

