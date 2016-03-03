# Compiler 
CC = mpic++
# Include
IFLAGS = /opt/intel/impi/include64/
#libraries
#LFLAGS = -shared-intel
LFLAGS = 
all: alg
alg: alg.o 
	$(CC) $(LFLAGS) alg.o  -o alg -lm 
alg.o: alg.c
	$(CC) $(LFLAGS) -c  alg.c 
clean:
	find . -exec touch {} \;
	rm -rf *o alg
