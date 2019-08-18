EXE    = qtree worms ballstick cpk
CC     = gcc
OFILES = qtree.o graphics.o commands.o
COPT   = -I$(HOME)/include -ansi -Wall -O3
LOPT   = -L$(HOME)/lib

LIBS   = -lbiop -lgen -lm -lxml2

all : $(EXE)

qtree :  $(OFILES) $(LFILES)
	$(CC) $(COPT) $(LOPT) -o $@ $(OFILES) $(LIBS)

worms :  worms.o $(UFILES)
	$(CC) $(COPT) $(LOPT) -o $@ worms.o $(LIBS)

ballstick : ballstick.o $(UFILES)
	$(CC) $(COPT) $(LOPT) -o $@ ballstick.o $(LIBS)

cpk : cpk.o $(UFILES)
	$(CC) $(COPT) $(LOPT) -o $@ cpk.o $(LIBS)

.c.o  :
	$(CC) $(COPT) -o $@ -c $<

clean :
	\rm -f *.o

distclean : clean
	\rm $(EXE)


