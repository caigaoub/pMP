CC=g++
CFLAGS=-W -Wall -ansi -pedantic -std=c++11

# Cplex configuration
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR=/opt/ibm/ILOG/CPLEX_Studio201/cplex/
CONCERTDIR=/opt/ibm/ILOG/CPLEX_Studio201/concert/
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl
CONCERTINCDIR = $(CONCERTDIR)include/
CPLEXINCDIR   = $(CPLEXDIR)include/
CCFLAGS = $(CFLAGS) -DIL_STD -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) -w
SRC = ./src/

main:  $(SRC)main.cpp  $(SRC)DataHandler.cpp $(SRC)pMedian.cpp  $(SRC)Supporter.cpp
	$(CC) $(CCFLAGS) -o $@ $^ $(CCLNFLAGS)


clean:
	rm -rf *.o

.PHONY: all main clean