OFILES = qtree.o graphics.o commands.o 
COPT   = STRICT ANSI MATH=8 IGNORE=51 IGNORE=149 OPT
CC     = sc
LOPT   = LIB lib:bookFPU.lib LIB:scm881.lib LIB:sc.lib LIB:amiga.lib LINK

qtree :  $(OFILES) 
	$(CC) $(COPT) $(OFILES) $(LOPT)

worms :  worms.o 
	$(CC) $(COPT) worms.o $(LOPT)

BallStick : BallStick.o 
	$(CC) $(COPT) BallStick.o $(LOPT)

.c.o  :
	$(CC) $(COPT) $<
   
   
