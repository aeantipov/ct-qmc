# Project: ctqmc

#-Wno-non-template-friend
CPP= g++ -O3 -Wno-non-template-friend -Wno-write-strings

PART= -Dorb1
#Determines how many particles do we have

ZONE= -Dzone1
#Determines how many zones do we have

CFLAGS = $(PART) $(ZONE)

default: ct-qmc

ct-qmc: ct-qmc.cpp ../main.cpp ../ini.cpp ../default.h
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





