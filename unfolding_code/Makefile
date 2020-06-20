#PREDEFINED VARIABLES
CC = icc
FF = gfortran

CFLAGS = -g -c -O3 -w
LFLAGS = -O3 -ldl -lm -lgsl -lgslcblas

PROGRAM = new

MODULOSC = $(PROGRAM).c lib.c 
MODULOSO = $(PROGRAM).o lib.o 

compile:	clean
	@$(CC) $(CFLAGS) $(MODULOSC)
	@$(CC) -o $(PROGRAM) $(MODULOSO) $(LFLAGS)

run:    compile
	@./$(PROGRAM)

clean:
	@rm -rf *.o $(PROGRAM) *.c~ *.dat~ *~ *.log *.mod
