# Project: ctqmc

#CPP = g++ -O3 -Wall -ansi -pedantic
CC = g++ -O3 -Wall -Wno-non-virtual-dtor -ansi -pedantic 
#CC = g++-4.1.2 -O3 
#-Wno-non-template-friend
CPP= g++-4.1.2 -O3 -Wno-non-template-friend
#CPP = icpc -O3 -Wall -Wno-non-virtual-dtor -ansi -pedantic
#NORB= -Dorb1
#NORB= -Dorb2
CFLAGS = $(INCS)

CTQMC_PATH=./ct-qmc/splitting10

CXXFLAGS = $(CXXINCS)

PART= -Dorb1
#Determines how many particles do we have

ZONE= -Dzone4
#Determines how many zones do we have

all: ct-qmc

ct-qmc: ct-qmc.cpp
	$(CPP) ct-qmc.cpp $(PART) $(ZONE)  -o ct-qmc

clean:
	${RM} $(OBJ) $(BIN) $(DBIN) $(BBIN) ct-qmc ../change *.o

cleandata:
	cp _input.dat _input.dat.backup
	rm -f *.dat
	rm -f *.it*
	mv _input.dat.backup _input.dat





