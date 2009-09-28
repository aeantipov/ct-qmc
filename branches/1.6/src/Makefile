# Project: ctqmc

#-Wno-non-template-friend
CPP= g++ -O3 -Wno-non-template-friend -Wno-write-strings
#CPP = icpc -O3 -Wall -Wno-non-virtual-dtor -ansi -pedantic
#NORB= -Dorb1
#NORB= -Dorb2

CTQMC_PATH=./ct-qmc/splitting10

PART= -Dorb1
#Determines how many particles do we have

ZONE= -Dzone4
#Determines how many zones do we have

MPI = -Dmpi -I/usr/include/eigen2
#Determines whether we use mpi or not

CFLAGS = $(PART) $(ZONE) $(MPI)

default: ct-qmc

ct-qmc: ct-qmc.cpp ../main.cpp ../ini.cpp
	$(CPP) $(CFLAGS) ct-qmc.cpp -o ct-qmc

change : ../change.cpp
	$(CPP) ../change.cpp -o ../change

clean:
	${RM} ct-qmc *.o

cleandata:
	cp _input.dat _input.dat.backup
	rm -f *.dat
	rm -f *.it*
	mv _input.dat.backup _input.dat





