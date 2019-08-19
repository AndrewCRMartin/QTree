EXE    = qtree worms ballstick cpk
CC     = gcc
OFILES = qtree.o graphics.o commands.o
COPT   = -I$(HOME)/include -ansi -Wall -O3
LOPT   = -L$(HOME)/lib
LIBS   = -lbiop -lgen -lm -lxml2

# If using PNG
GLIBS   = -lpng
GOFILES = writepng.o
GSUPP   = -DSUPPORT_PNG

all : $(EXE)

qtree :  $(OFILES) $(LFILES) $(GOFILES)
	$(CC) $(COPT) $(LOPT) -o $@ $(OFILES) $(GOFILES) $(GLIBS) $(LIBS)

worms :  worms.o $(UFILES)
	$(CC) $(COPT) $(LOPT) -o $@ worms.o $(LIBS)

ballstick : ballstick.o $(UFILES)
	$(CC) $(COPT) $(LOPT) -o $@ ballstick.o $(LIBS)

cpk : cpk.o $(UFILES)
	$(CC) $(COPT) $(LOPT) -o $@ cpk.o $(LIBS)

.c.o  :
	$(CC) $(COPT) $(GSUPP) -o $@ -c $<

clean :
	\rm -f *.o

distclean : clean
	\rm $(EXE)


