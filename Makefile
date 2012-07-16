CC=gcc
EXEC=MDBAS

DPAR=-g -O0 -Wall
RPAR=-O2
LIB=-lm
DEF=-DDSFMT_MEXP=19937

all:
	make release
	
debug:
	$(CC) $(DPAR) $(DEF) *.c dSFMT/*.c -o $(EXEC) $(LIB)	

release:
	$(CC) $(RPAR) $(DEF) *.c dSFMT/*.c -o $(EXEC) $(LIB)

clean:
	rm -f *.o $(EXEC)

