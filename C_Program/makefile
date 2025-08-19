CC=gcc
#OT=-O3 -tpp7 -axW -xW -unroll0
OT=-ffast-math -O3
OP=$(OT)
LOP=-lm -o

# For debugging with x86
#CC=gcc
#OT=-g -Wall -rdynamic
#OP= $(OT)
#LOP=-lm -o

CHDR=./
CSRC=./
MSRC=./

MAIN=max_modularity
VER=


EXE=$(MAIN)$(VER)$(SIZE)
TAG=$(EXE).out

OBJ=$(MAIN).o mt19937-64.o modularity_functions.o

$(TAG)	: $(OBJ)
	$(CC) $(OBJ) $(LOP) $(TAG)
	rm *.o

$(MAIN).o :
	$(CC) $(OP) -I$(CHDR) -c $(MSRC)$(MAIN).c

mt19937-64.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)mt19937-64.c

network_functions.o	:
	$(CC) $(OP) -I$(CHDR) -c $(CSRC)modularity_functions.c

clear	:
	rm data/*; rm index.dat; cp work/words.lst ./

clean	:
	rm *.o
