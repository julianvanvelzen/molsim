CC          = mpicc
CLINKER     = mpicc

CFLAGS      = -Wall -O4 -g -w
LIBS        = -lm 
DEPEND		= makedepend

SRC     = pdm.c ran_uniform.c initialize.c helpers.c mdloop.c
OBJS    = pdm.o ran_uniform.o initialize.o helpers.o mdloop.o
EXECS   = pdm

default: pdm

all: $(EXECS)

pdm:$(OBJS)
	$(CLINKER) $(OPTFLAGS) -o pdm $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o *~ $(EXECS)

.c.o:
	$(CC) $(CFLAGS) -c $*.c

pdm.o: system.h 
initialize.o: system.h
ran_uniform.o: system.h
helpers.o: system.h
mdloop.o: system.h