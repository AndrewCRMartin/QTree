CC     = cc
OFILES = qtree.o graphics.o commands.o 
LFILES = bioplib/RotPDB.o bioplib/ReadPDB.o bioplib/help.o \
bioplib/matrix.o bioplib/parse.o bioplib/general.o \
bioplib/fsscanf.o bioplib/ModPDB.o bioplib/ApMatPDB.o \
bioplib/OrigPDB.o bioplib/GetCGPDB.o bioplib/WindIO.o \
bioplib/WritePDB.o bioplib/angle.o bioplib/array.o \
bioplib/ParseRes.o bioplib/PDBList.o bioplib/CalcPDB.o \
bioplib/ZonePDB.o


UFILES = bioplib/ReadPDB.o bioplib/angle.o \
bioplib/array.o bioplib/fsscanf.o bioplib/PDBList.o \
bioplib/WritePDB.o bioplib/general.o bioplib/WindIO.o

all : qtree worms ballstick cpk

clean :
	\rm -f $(OFILES) $(LFILES) $(UFILES) ballstick.o cpk.o worms.o

qtree :  $(OFILES) $(LFILES)
	$(CC) -o qtree $(OFILES) $(LFILES) -lm

worms :  worms.o $(UFILES)
	$(CC) -o worms worms.o $(UFILES) -lm

ballstick : ballstick.o $(UFILES)
	$(CC) -o ballstick ballstick.o $(UFILES) -lm

cpk : cpk.o $(UFILES)
	$(CC) -o cpk cpk.o $(UFILES) -lm

.c.o  :
	$(CC) -o $@ -c $<
   
   