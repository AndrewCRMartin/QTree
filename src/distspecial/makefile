CC     = cc
EXE = qtree worms ballstick cpk
OFILES = qtree.o graphics.o commands.o 
LFILES = bioplib/RotPDB.o bioplib/ReadPDB.o bioplib/help.o \
bioplib/parse.o bioplib/throne.o bioplib/strcatalloc.o \
bioplib/fsscanf.o bioplib/ApMatPDB.o \
bioplib/WindIO.o bioplib/GetCGPDB.o \
bioplib/WritePDB.o bioplib/angle.o \
bioplib/ParseRes.o bioplib/StoreString.o \
bioplib/MatMult3_33.o bioplib/hash.o bioplib/prime.o \
bioplib/stringutil.o bioplib/PDBHeaderInfo.o \
bioplib/GetWord.o bioplib/array2.o bioplib/stringcat.o \
bioplib/GetPDBChainLabels.o bioplib/BuildConect.o \
bioplib/IndexPDB.o bioplib/FindNextResidue.o \
bioplib/FindResidue.o bioplib/padterm.o bioplib/KillLeadSpaces.o \
bioplib/StringToUpper.o bioplib/OpenFile.o bioplib/chindex.o \
bioplib/FreeStringList.o bioplib/TranslatePDB.o \
bioplib/OriginPDB.o \
bioplib/CreateRotMat.o  bioplib/InPDBZone.o

UFILES = bioplib/ReadPDB.o bioplib/angle.o \
bioplib/array2.o bioplib/fsscanf.o \
bioplib/WritePDB.o bioplib/WindIO.o bioplib/StoreString.o \
bioplib/hash.o bioplib/prime.o bioplib/stringutil.o \
bioplib/BuildConect.o bioplib/IndexPDB.o  bioplib/FindNextResidue.o \
bioplib/FindResidue.o bioplib/FreeStringList.o  bioplib/PDBHeaderInfo.o \
bioplib/throne.o  bioplib/strcatalloc.o bioplib/GetWord.o \
bioplib/stringcat.o bioplib/GetPDBChainLabels.o bioplib/padterm.o \
bioplib/chindex.o bioplib/CopyPDB.o bioplib/OpenStdFiles.o


all : $(EXE)


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

clean :
	\rm -f *.o biopib/*.o

distclean : clean
	\rm -f $(EXE)
