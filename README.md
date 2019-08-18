



                              QTree V1.7
                              ==========


                           SciTech Software

                             23 Stag Leys,
                               Ashtead,
                                Surrey.
                               KT21 2TD.

                        Tel.: +44 (0)372 275775


   =====================================================================

      EMail: amartin@scitec.adsp.sub.org
             cbmehq!cbmuk!scitec!amartin@cbmvax.commodore.com
             {uunet|rutgers}...!cbmvax!cbmehq!cbmuk!scitec!amartin
      -or-
             andrew@uk.ac.ox.biop

   =====================================================================


   QTree is Copyright (c) 1993-4, Dr. Andrew C.R. Martin


   The executable program and source code are freely distributable,
   providing all files are distributed intact. Further versions of 
   the file graphics.c may also be distributed with the software to
   provide supprt for alternative hardware platforms and display
   devices.

   With the exception of additional versions of graphics.c, you may not
   distribute any changes you make to QTree. Please send any such 
   changes to SciTech Software so that they can be incorporated into
   future releases of the software.

   Please note that the source code may NOT be used for any other 
   purposes without the prior permission of SciTech Software. This is
   because it is intended to produce a library of such routines which
   will be made available (probably in the form of a book).
                              
   =====================================================================




Introduction
============

   QTree is a program for generating space-filling pictures of molecules
from PDB files. The code is written to be as portable as possible. Only 
the routines in graphics.c need to be changed to support the display 
device being used. 

   QTree and its support programs, Worms and BallStick have been
compiled and tested on Commodore Amiga, Evans & Sutherland ESV running
Unix (native compiler), and on a PC running Linux (Gnu C compiler).
It should run on anything else with no problems!

   As of V1.5 of QTree, graphics.c has been changed simply to create
a 24bit RGB image in MTV raytracer format. No direct graphics device
support is provided. You are encouraged to maintain this form and
thus to write a separate utility to display MTV files on your graphics
device. Source code for reading an MTV file is provided. If you
really must, you may modify graphics.c to drive your graphics
directly, but the other route is preferable for ease of portability
of the code.

   Worms is a preprocessor program for QTree which allows the generation 
of worm style images showing a B-spline smoothed image of the C-alpha 
trace. It reads a PDB file and writes a PDB file containing the 
interpolated smoothed atom positions necessary to generate the worms 
image.

   BallStick is a preprocessor program for QTree which allows the 
generation of ball and stick style images. It reads a PDB file and writes 
a PDB file containing interpolated sphere positions along each bond.


   In the explanations which follow, items enclosed in square brackets 
are optional; items in angle brackets should be replaced by suitable 
values.



Running QTree
=============

   Required files
   --------------
   
   qtree                The executable program
   qtree.hlp            The help file. This should be in the current
                        directory, or a directory with the logical name
                        HELP: (not Unix)
   
   
   Running the program
   -------------------

   Simply type:
   
      qtree <file.pdb> <file.mtv>         (CPK & Worms images)

   A control file (read the help information with qtree -h) may be
specified. For example:
   
      qtree -c <ctrl.dat> <file.pdb> <file.mtv>

   If generating ball and stick images after using the BallStick 
preprocessor, you should use the -b flag. For example:
   
      qtree -b -c <ctrl.dat> <file.stk> <file.mtv>

   The size of the square used to perform the calculations may be 
specified using the -r option. (Defaults to 512). For example:

      qtree -r 256 <file.pdb> <file.mtv>

   The size of the screen used (i.e. the dimensions of the output file)
are specified with the -s option. (Defaults to 800x600.) For example:

      qtree -r 256 -s 256 256 <file.pdb> <file.mtv>





   24-bit output is created in a format compatible with the MTV ray 
tracer.
      



   To get further help, type
   
      qtree -h
      
      

      
Running Worms
=============

   Required files
   --------------

   worms                The executable program

   Running the program
   -------------------

   Simply type:
   
      worms [-n <n>] [-s <n>] [-d] <in.pdb> <out.pdb>
   
   The -n switch is followed by an integer which specifies the
   number of sphere to be placed between smoothed atom positions
   (default: 30).

   The -d switch causes division smoothing to be used instead of
   B-spline smoothing.
   
   The -s switch is followed by an integer division smoothing factor 
   (default: 4). This option is ignored if -d has not been specified.
   
   
   
Running BallStick
=================

   Required files
   --------------
   
   BallStick            The executable program
   
   Running the program
   -------------------
      
   Simply type:
   
      BallStick [-n <n>] [-b <b>] [-s <s>] [-d] <in.pdb> <out.pdb>
      
   The -n switch is followed by an integer specifying the number of small
   spheres to be placed along each bond (default: 30). With large 
   structures it may be necessary to reduce this number to reduce memory 
   usage by qtree.
   
   The -b switch is followed by a floating point number which specifies
   the ball radius (default 0.4).
   
   The -s switch is followed by a floating point number which specifies
   the stick radius (default 0.2).

   The -d switch stops sticks from being created between disulphides.
   
      
      
      
Compiling QTree, Worms and BallStick
====================================

   Required files
   --------------
   
      qtree.c        The QTree program
      qtree.h        Structure definitions, globals and flags
      qtree.p        Prototypes for qtree.c
      graphics.c     Graphics support routines
      graphics.p     Prototypes for graphics.c
      commands.c     Command parser setup and handling
      commands.p     Prototypes for commands.c

      worms.c        The Worms program
      
      ballstick.c    The Ball and Stick program

      
   Required library files and includes
   -----------------------------------

      These should be placed in a sub-directory of the current directory
      called `bioplib'.
   
      macros.h          Generally useful C programming macros
      MathType.h        Maths type definitions
      SysDefs.h         General system type definitions
      pdb.h             PDB structure definitions and library prototypes
      help.h            Help interface prototype include
      parse.h           Command line parser defines and prototype include
      general.h         General C library routine prototype include
      matrix.h          Matrix handling prototype include
      WindIO.h          Window I/O routine prototype include
      fsscanf.h         Hard formatted string scan include
      CursWind.h        Curses/windowing prototype include
      angle.h           Angle calculation include
      array.h           2D array creation include

      ReadPDB.p         Prototype files
      WritePDB.p            "       "
      getcofgPDB.p          "       "
      originpdb.p           "       "
      rotatepdb.o           "       "
      translatepdb.p        "       "
      FitPDB.p              "       "
      FindZonePDB.p         "       "
      HAddPDB.p             "       "
      FixPDB.p              "       "
      SelectAtomsPDB.p      "       "
      RotatePDB.p           "       "
      ReadSecPDB.p          "       "
      RenumAtomsPDB.p       "       "
      PDBUtil.p             "       "
      ApplyMatrixPDB.p      "       "
      parse.p               "       "
      general.p             "       "
      matrix.p              "       "
      WindIO.p              "       "
      fsscanf.p             "       "
      CursWind.p            "       "
      angle.p               "       "
      array.p               "       "

      RotatePDB.c       Rotate a PDB linked list
      ReadPDB.c         Read a PDB file
      WritePDB.c        Write a PDB file
      Help.c            Help file handler
      Matrix.c          Matrix manipulation routines
      Parse.c           Comman line parser
      General.c         General routines
      fsscanf.c         Hard formatted string scan
      WindIO.c          Window I/O routines
      TranslatePDB.c    Move a PDB file
      ApplyMatrixPDB.c  Apply a rotation matrix to a PDB linked list
      OriginPDB.c       Move a PDB linked list
      GetCofGPDB.c      Get centre of geometry from a PDB linked list
      PDBUtil.c         Misc PDB handling utilities
      array.c           2D array generation code


   Compilation procedure
   ---------------------
   1. Create a directory called `QTree',
   2. Install the required files listed above in this directory,
   3. Modify graphics.c as required to support the output device you
      are using,
   4. Create a subdirectory of QTree called `bioplib',
   5. Install the library files and includes listed above in this 
      directory,
   6. Enter the bioplib directory and compile each of the .c files,
   7. Enter the QTree directory and compile each of the .c files,
   8. Link the object files in the QTree directory with those in the
      bioplib subdirectory.


   There are 3 compile time options (controlled by symbols defined in
   qtree.h) which may be used to switch off features of qtree thus 
   increasing its running speed. These are documented in qtree.c

   N.B. Compiling on MS-DOS machines will require modification of the
   path names for include files to reflect the backslash path separator
   (i.e. change bioplib/... to bioplib\...). The VAX/VMS C compiler is 
   intelligent enough to cope with slash separted paths even though VMS 
   uses a dot (.) to separate paths.  A number of the filenames are too 
   long for use with MS-DOS and these will also have to be changed. 
   (Suggestion: use a decent operating system!!)


Possible Future Enhancements
============================

   1. Progress indicator
   2. Non-square image support
   
