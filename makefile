NAME = my_test.exe

BASE = $(HOME)/my/path/to/fortran_VSH
EXEC = $(BASE)/$(NAME)
OBJD = $(BASE)/obj
SRCD = $(BASE)/src

F77 = gfortran -O -I $(SRCD)

FFLAGS = 

VPATH = $(OBJD):$(SRCD)

OBJS = main.o globals.o kinds.o vsh.o

$(NAME): $(OBJS)
	@echo Building executable
	cd $(OBJD); $(F77) -o $(EXEC) $(OBJS)

.f.o:
	@echo Building $*.o 
	cd $(OBJD); $(F77) -c $<

clean:
	rm -f terminal* makefile~ *.dat~ obj/*.o $(EXEC)
