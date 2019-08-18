CC     = cc -ansi -Wall -pedantic
OFILES = qtree.o graphics.o commands.o 
#COPT   = -I/home/bsm/martin/include -L/home/bsm/martin/lib 
COPT   = -I/home/amartin/include -L/home/amartin/lib 

all : qtree worms ballstick cpk


qtree :  $(OFILES) 
	$(CC) $(COPT) -o qtree $(OFILES) -lbiop -lgen -lm

worms :  worms.o 
	$(CC) $(COPT) -o worms worms.o -lbiop -lgen -lm 

ballstick : ballstick.o 
	$(CC) $(COPT) -o ballstick ballstick.o -lbiop -lgen -lm

cpk : cpk.o 
	$(CC) $(COPT) -o cpk cpk.o -lbiop -lgen -lm

.c.o  :
	$(CC) $(COPT) -o $@ -c $<

clean :
	\rm *.o
