CC     = cc
OFILES = qtree.o graphics.o commands.o 
LFILES = bioplib/RotPDB.o bioplib/ReadPDB.o bioplib/help.o \
bioplib/parse.o  \
bioplib/fsscanf.o bioplib/ApMatPDB.o \
 bioplib/GetCGPDB.o bioplib/WindIO.o \
bioplib/WritePDB.o bioplib/angle.o  \
bioplib/ParseRes.o \
bioplib/array2.o \
bioplib/padterm.o \
bioplib/OriginPDB.o \
bioplib/InPDBZone.o \
bioplib/KillLeadSpaces.o \
bioplib/StringToUpper.o \
bioplib/OpenFile.o \
bioplib/MatMult3_33.o \
bioplib/FindNextResidue.o \
bioplib/chindex.o \
bioplib/TranslatePDB.o \
bioplib/CreateRotMat.o 







UFILES = bioplib/ReadPDB.o bioplib/angle.o \
bioplib/fsscanf.o \
bioplib/WritePDB.o \
bioplib/OpenStdFiles.o \
bioplib/CopyPDB.o \
bioplib/chindex.o \
bioplib/FindNextResidue.o \
bioplib/padterm.o






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
   
   