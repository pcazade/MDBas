#################################################################
########################   MakeVars   ###########################
#################################################################

CC=gcc

CIBLE=user_functions.so

INC_DIR=../include

CC_OPT=-I$(INC_DIR) -std=gnu99 -Wall -O2 -msse2 -fPIC
SO_OPT=-shared -Wl,-soname,$(CIBLE)
 
SRC=$(wildcard ./*.c)
 
OBJ=$(patsubst ./%.c,./%.o,$(SRC))
 
#################################################################
########################   Makefile   ###########################
#################################################################
 
all:$(CIBLE)
	@echo "Compilation Success"
 
./%.o:./%.c 
	$(CC) $(CC_OPT) -c $< -o $@ 

$(CIBLE):$(OBJ)
	$(CC) $(OBJ) $(SO_OPT) -o $@

clean:
	rm -f $(CIBLE) ./*.o

