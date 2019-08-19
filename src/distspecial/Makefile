EXE    = qtree worms ballstick cpk
CC     = cc 
COPT   = -ansi -Wall -O3 -Wno-unused-function
OFILES = qtree.o graphics.o commands.o 
LIBS   = -lm

# If using PNG - You need the libpng development library to be installed
# comment these out if you don't want this support
GLIBS   = -lpng
GOFILES = writepng.o
GSUPP   = -DSUPPORT_PNG

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
	$(CC) $(COPT) -o $@ $(OFILES) $(LFILES) $(LIBS) $(GLIBS) $(LIBS)

worms :  worms.o $(UFILES)
	$(CC) $(COPT) -o $@ worms.o $(UFILES) $(LIBS)

ballstick : ballstick.o $(UFILES)
	$(CC) $(COPT) -o $@ ballstick.o $(UFILES) $(LIBS)

cpk : cpk.o $(UFILES)
	$(CC) $(COPT) -o $@ cpk.o $(UFILES) $(LIBS)

.c.o  :
	$(CC) $(COPT) $(GSUPP) -o $@ -c $<

clean :
	\rm -f *.o bioplib/*.o

distclean : clean
	\rm -f $(EXE)
