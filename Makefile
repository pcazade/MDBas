#################################################################
########################   MakeVars   ###########################
#################################################################

#use : make DEBUG=ON for a debug build
DEBUG=OFF
COMP=gcc
OMP=OFF

ifeq ($(COMP),icc)

ifeq ($(OMP),ON)
CC=icc -openmp
else
CC=icc
endif

else

ifeq ($(OMP),ON)
CC=gcc -fopenmp
else
CC=gcc
endif

endif

ifeq ($(DEBUG),OFF) 
CC_OPT=-I"./dSFMT" -I"./include" -std=c99 -Wall -O2 -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=19937
else
CC_OPT=-I"./dSFMT" -I"./include" -std=c99 -Wall -Wextra -O0 -g -msse2 -DHAVE_SSE2 -DDSFMT_MEXP=19937
endif

CC_SFMT_OPT=-I"./dSFMT" -std=c99 -O2 -msse2 -fno-strict-aliasing -DHAVE_SSE2 -DDSFMT_MEXP=19937

LD_OPT=-lm 

MKDIR=mkdir -p ./obj/dSFMT
 
CIBLE=MDBas
 
SRC=$(wildcard ./src/*.c)
dSRC=$(wildcard ./dSFMT/*.c)
 
OBJ=$(patsubst ./src/%.c,./obj/%.o,$(SRC))
dOBJ=$(patsubst ./dSFMT/%.c,./obj/dSFMT/%.o,$(dSRC)) 
 
#################################################################
########################   Makefile   ###########################
#################################################################
 
all:$(CIBLE)
	@echo "Compilation Success"
 
./obj/%.o:./src/%.c 
	$(CC) $(CC_OPT) -c $< -o $@ 

./obj/dSFMT/%.o:./dSFMT/%.c
	@$(MKDIR)
	$(CC) $(CC_SFMT_OPT) -c $< -o $@
 
$(CIBLE):$(dOBJ) $(fOBJ) $(OBJ)
	$(CC) $(dOBJ) $(fOBJ) $(OBJ) -o $@ $(LD_OPT)

clean:
	rm -f $(CIBLE) ./obj/*.o

clean_all:
	rm -f $(CIBLE) ./obj/*.o ./obj/dSFMT/*.o

