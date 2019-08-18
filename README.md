



                              QTree V2.3
                              ==========

                        Dr. Andrew C.R. Martin
                           SciTech Software


                            Now working at:

              Biomolecular Structure and Modelling Unit,
          Department of Biochemistry and Molecular Biology,
                      University College London,
                            Gower Street,
                           LONDON WC1E 6BT

              Some modifications made while working at:

               School of Animal and Microbial Sciences,
                        University of Reading,
                            Whiteknights,
                             P.O.Box 228,
                           Reading RG6 6AJ.


   =====================================================================

                     EMail: andrew@bioinf.org.uk

   =====================================================================


   QTree is Copyright (c) 1993-2007, Dr. Andrew C.R. Martin


   This program is not in the public domain.

   It may not be copied or made available to third parties, but may be
   freely used by non-profit-making organisations who have obtained it
   directly from the author or by FTP.
   
   You are requested to send EMail to the author to say that you are
   using this code so that you may be informed of future updates.
   
   The code may not be made available on other FTP sites without express
   permission from the author.
   
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If
   someone else breaks this code, the author doesn't want to be blamed
   for code that does not work! You may not distribute any
   modifications, but are encouraged to send them to the author so
   that they may be incorporated into future versions of the code.
   
   The code may not be sold commercially or used for commercial purposes
   without prior permission from the author.
                                                
   =====================================================================




Introduction
============

   QTree is a program for generating CPK (Corey, Pauling, Koltun)
space-filling pictures of molecules from PDB files.  It can also
generate worms and ball-and-stick images using preprocessing filters.
The code is written to be as portable as possible. Only the routines
in graphics.c need to be changed to support the display device being
used.

   QTree and its support programs, CPK, Worms and BallStick have been
compiled and tested on Commodore Amiga, Silicon Graphics and Evans &
Sutherland ESV running Unix (native compilers), and on a PC running
Linux (Gnu C compiler).  It should run on anything else with no
problems!

   As of V1.5 of QTree, graphics.c has been changed simply to create
a 24bit RGB image in MTV raytracer format. No direct graphics device
support is provided. You are encouraged to maintain this form and
thus to write a separate utility to display MTV files on your graphics
device. Source code for reading an MTV file is provided. If you
really must, you may modify graphics.c to drive your graphics
directly, but the other route is preferable for ease of portability
of the code.

   Under UNIX/X-Windows, MTV graphics files may be displayed using the
ImageMagick display program. The convert program from the ImageMagick
package may be used in a pipeline to generate GIF output (or whatever
other form you prefer).

   Worms is a preprocessor program for QTree which allows the generation 
of worm style images showing a B-spline smoothed image of the C-alpha 
trace. It reads a PDB file and writes a PDB file containing the 
interpolated smoothed atom positions necessary to generate the worms 
image.

   BallStick is a preprocessor program for QTree which allows the 
generation of ball and stick style images. It reads a PDB file and writes 
a PDB file containing interpolated sphere positions along each bond.

   CPK is a preprocessor program for QTree which places the default
atom radii in the occupancy column. This allows BallStick and CPK
representations to be combined in one image.

   As of V1.11, the files generated by BallStick are different: the radii
are stored in the occupancy column rather than the temperature factor
column. The allows the colouring on temperature factors (introduced in
V1.10) to be performed with ball & stick representations. This means
that BallStick output files generated with versions prior to V1.11
will not be compatible with the latest version of QTree. A simple
Perl or awk script may be used to swap the B-value and occupancy columns
of an old BallStick output file if required.

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
                        HELP: (VMS/AmigaDOS), or in a directory whose
                        name is placed in the environment variable
                        HELPDIR (Unix)
   
   
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

   There is also a -q switch which causes the program to run quietly
without gerenating any copyright of informational messages.



   24-bit output is created in a format compatible with the MTV ray 
tracer. This may be displayed using the ImageMagick package under Unix.
      


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
   
      worms [-n <n>] [-s <n>] [-d] [-q] <in.pdb> <out.pdb>
   
   The -n switch is followed by an integer which specifies the
   number of sphere to be placed between smoothed atom positions
   (default: 30).

   The -d switch causes division smoothing to be used instead of
   B-spline smoothing.
   
   The -s switch is followed by an integer division smoothing factor 
   (default: 4). This option is ignored if -d has not been specified.
   
   The -q switch causes the program to run quietly without gerenating
   any copyright of informational messages.
   
   
Running BallStick
=================

   Required files
   --------------
   
   ballstick            The executable program
   
   Running the program
   -------------------
      
   Simply type:
   
      ballstick [-n <n>] [-b <b>] [-s <s>] [-d] [-q] <in.pdb> <out.pdb>
      
   The -n switch is followed by an integer specifying the number of small
   spheres to be placed along each bond (default: 30). With large 
   structures it may be necessary to reduce this number to reduce memory 
   usage by qtree. With small structures, you may wish to increase this
   value to get smoother sticks.
   
   The -b switch is followed by a floating point number which specifies
   the ball radius (default 0.4).
   
   The -s switch is followed by a floating point number which specifies
   the stick radius (default 0.2).

   The -d switch stops sticks from being created between disulphides.

   The -q switch causes the program to run quietly without gerenating
   any copyright of informational messages.



Running CPK
===========
   
   Required files
   --------------
   
   cpk                  The executable program
   
   Running the program
   -------------------
      
   Simply type:
   
      cpk [-q] <in.pdb> <out.pdb>

   The -q switch causes the program to run quietly without gerenating
   any copyright of informational messages.
      




Pipes Under Unix
================

   As of V2.0, all programs in the QTree package support I/O through
pipes. Therefore, to generate a worms image and display it using the
ImageMagick display program, you may use a command line like:

      worms <in.pdb> | qtree | display mtv:-

   To generate a combined image where file1.pdb is to be rendered as
CPK while file2.pdb is to be rendered as ball and stick, one would use
the following:

      cpk       file1.pdb    image.bst
      ballstick file2.pdb  >>image.bst
      qtree -b image.bst | display mtv:-




Documentation
=============

   There is no separate documentation other than what is supplied here
since there is an extensive help facility built into the QTree program.
Type qtree -h to access this help facility. (N.B. You must either have
qtree.hlp in your current directory, or in the directory pointed to
by the environment variable HELPDIR to access help.)

   A simple command file for QTree, restype.dat, is supplied which 
colours on residue type with a shaded background and specular reflections.

      
      
Compiling QTree, Worms, CPK and BallStick
=========================================

   When you unpack the tar file, you will obtain a QTreeVx.xx directory
with a subdirectory called bioplib.

   Enter the QTreeVx.xx directory and modify the makefile to change the
name of the C compiler for your system if required. You may also modify
qtree.h to change some default values and switch off compilation of some
optional sections. The following values may be changed:

SIZE            This is the actual image size and must be a power of 2
                (Overridden with -r flag)
XSIZE, YSIZE    This is the size of the background on which the image will
                be placed. Both values must be larger than SIZE
                (Overridden with -s flag)
DEPTHCUE        If this is not defined, code to handle depth cueing will
                not be compiled
SPEC            If this is not defined, code to handle specular reflections
                will not be compiled
SHOW_INFO       If this defined, the run time and pixel saturation will be
                shown at the end of the run. ESVs don't seem to have the
                required system calls...
OVERLAP_SLAB    If defined, Z-slab regions will include any atoms which
                partially overlap the slab region; if not the centre of
                the atom must be within the region.


Below are listed the files contained in this distribution

   Source files
   ------------
      qtree.c        The QTree program
      qtree.h        Structure definitions, globals and flags
      qtree.p        Prototypes for qtree.c
      graphics.c     Graphics support routines
      graphics.p     Prototypes for graphics.c
      commands.c     Command parser setup and handling
      commands.p     Prototypes for commands.c

      worms.c        The Worms program
      
      BallStick.c    The Ball and Stick program
      cpk.c          The CPK program

   Additional source files
   -----------------------
      ColourGraphics.c  An example of driving colour graphics directly
                        (replaces graphics.c)
      GreyGraphics.c    An example of driving grey-scale graphics directly
                        (replaces graphics.c)
      mtvham.c          Code to display an MTV file in HAM mode on the 
                        Commodore Amiga
      ReadMTV.c         Code to read an MTV graphics file

      
   Required library files and includes
   -----------------------------------
      These should be placed in a sub-directory of the QTree directory
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

      RotPDB.c          Rotate a PDB linked list
      ReadPDB.c         Read a PDB file
      WritePDB.c        Write a PDB file
      help.c            Help file handler
      matrix.c          Matrix manipulation routines
      parse.c           Comman line parser
      general.c         General routines
      fsscanf.c         Hard formatted string scan
      WindIO.c          Window I/O routines
      ModPDB.c          Modify coordinates in a PDB linked list
      ApMatPDB.c        Apply a rotation matrix to a PDB linked list
      OrigPDB.c         Move a PDB linked list to origin
      GetCGPDB.c        Get centre of geometry from a PDB linked list
      ParseRes.c        Parse a residue specification
      CalcPDB.c         Perform calculations on PDB linked list
      PDBList.c         Manipulate PDB linked list
      array.c           2D array generation code
      angle.c           Angle and torsion calculation code


   N.B. Compiling on MS-DOS machines will require modification of the
   path names for include files to reflect the backslash path separator
   (i.e. change bioplib/... to bioplib\...). The VAX/VMS C compiler is 
   intelligent enough to cope with slash separted paths even though VMS 
   uses a dot (.) to separate paths.  A number of the filenames are too 
   long for use with MS-DOS and these will also have to be changed. 
   (Suggestion: use a decent operating system like Linux!!)


Possible Future Enhancements
============================

   1. Progress indicator
   2. Non-square image support
   3. Colour by CHAIN
   

Revision History
================

   V1.0  19.07.93 Original
   V1.1  28.07.93 Added support for ball and stick
                  Added B-spline smoothing
   V1.2  29.07.93 Corrected bug in failure to open file
                  Added MTV file output support and background colouring
                  Added support for disulphide bonds
   V1.3  11.08.93 Corrected usage message
                  Fixed memory leak
   V1.4  07.10.93 Fixed bug in FindChainPDB()
                  Added option to create output file of specifed
                  dimensions
                  Removed direct graphic display
   V1.5  14.09.93 Added sphere scaling option
   V1.6  04.01.94 Fixed bug in argument parsing
                  Some externals made static static and added casts for 
                  GCC
                  Fixed bug in BSplineSmoothPDB()
   V1.7  24.03.94 Minimal tidying up
                  Added SLAB option and warning messages.
                  Removed ParseResSpec() as this is now in the library
                  Applies sphere scaling when radius comes from B-value
   V1.8  09.05.94 Fixed rounding error bug in B-spline smoothing
   V1.9  13.05.94 Various fixes to B-spline smoothing which was getting
                  coords out of sync with labels...
   V1.10 24.06.94 Default colouring may now be done on temperature
                  factor using command TEMPERATURE
   V1.11 04.10.94 With -b, reads radii from occ rather than bval
   V1.12 21.12.94 Improved Usage message
   V2.0  28.03.95 Modified to allow I/O through pipes
                  Fixed DoZone() which wasn't working properly for
                  Ball and Stick images which split the PDB data into
                  two sections.
                  Ctrl-C handled properly on non-AmigaDOS systems.
   V2.1  23.10.95 Warnings and errors all go to stderr
   V2.2  14.10.03 Added BOUNDS and RADIUS commands
   V2.2a 18.10.07 Clean compile with -ansi -Wall -pedantic
   V2.3  18.10.07 Added HIGHIGHT / BORDERWIDTH
