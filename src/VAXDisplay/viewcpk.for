C  N.B. This was the old software to read from the last attempt using
C  `CPK' rather than `QTree'.

      PROGRAM VIEWCPK
C***********************************************************************
C
C  Program:    VIEWCPK
C  File:       VIEWCPK.FOR
C  
C  Version:    V2.0
C  Date:       08.09.88
C  Function:   Display a 256 colour datafile on a colour VaxStation
C  
C  Copyright:  SciTech Software 1992
C  Author:     Dr. Andrew C. R. Martin
C  Address:    SciTech Software
C              23, Stag Leys,
C              Ashtead,
C              Surrey,
C              KT21 2TD.
C  Phone:      +44 (0372) 275775
C  EMail:      UUCP: cbmuk!cbmuka!scitec!amartin
C              JANET: andrew@uk.ac.ox.biop
C              
C***********************************************************************
C
C  This program is not in the public domain, but it may be freely copied
C  and distributed for no charge providing this header is included.
C  The code may be modified as required, but any modifications must be
C  documented so that the person responsible can be identified. If
C  someone else breaks this code, I don't want to be blamed for code 
C  that does not work! The code may not be sold commercially without 
C  prior permission from the author, although it may be given away free
C  with commercial products, providing it is made clear that this 
C  program is free and that the source code is provided with the program
C
C***********************************************************************
C
C  Description:
C  ============
C
C  Reads a VaxCPK image file and displays it on a VAX colour 
C  workstation.
C  This version clips the bottom 85 lines to center the 1024 line image
C  on the 854 lines of the workstation.
C
C
C***********************************************************************
C
C  Usage:
C  ======
C
C***********************************************************************
C
C  Revision History:
C  =================
C
C  V2.0   18.02.92
C  Changed to read the datafile in C for use with my VaxCPK
C  ray-tracer. The VLT file is assumed to be AMDATA:SHADE.VLT; the
C  name of the image file is prompted for.
C
C***********************************************************************

      parameter nx=1024,ny=1024

      byte   fbuf
      COMMON /IMGCMN/ fbuf(nx,ny)

      integer vcm_size
      real red(256),green(256),blue(256)
      byte image(nx)
      character*80 junk
c
      include 'sys$library:uisentry'
      include 'sys$library:uisusrdef'
c
c
      data vcm_size/256/
      data icolunit,imageunit/44,1/
      data xmax,ymax/35,30/
      data nreject/85/
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Set up an attributes structure to enable a borderless window
c
      structure/place/
      integer*4 code_5
      integer*4 pos_x
      integer*4 code_6
      integer*4 pos_y
      integer*4 code_7
      integer*4 attr
      integer*4 end_of_list
      end structure
c
      record /place/w_attributes
c
      w_attributes.code_5 = wdpl$c_abs_pos_x
      w_attributes.pos_x  = 0.0
      w_attributes.code_6 = wdpl$c_abs_pos_y
      w_attributes.pos_y  = 0.0
      w_attributes.code_7 = wdpl$c_attributes
      w_attributes.attr   = wdpl$m_noborder
      w_attributes.end_of_list = wdpl$c_end_of_list
c
c     Attribute record complete. 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Initialise the display
c
      id_vcm=uis$create_color_map(vcm_size)
c
      id_vd =uis$create_display(0.0,0.0,xmax,ymax,xmax,ymax,id_vcm)
c
c     Read the VLT and set up the colours
c
      open(unit=icolunit,file='amdata:shade.vlt',status='old',
     +     form='formatted',readonly,shared)
      call readvlt(icolunit,id_vd)
      print *,'NOTE: The workstation clips the top and bottom 85 lines'
      print *,'      of the image. Use the SMALL option in CPK if this'
      print *,'      is a problem! (Use the version of CPK in ',
     +        ' $8$dra0:[andrew.progs])'
      print *
      print *,'Colour file read, reading image data....'
c
c     Read the image file
c
C      open(unit=imageunit,status='old',readonly,shared,
C     +form='unformatted')
c
      call readimage(nreject)
c
c     Display the image and open the window
c
   99 call uis$image(id_vd,1,0.0,0.0,xmax,ymax,nx,ny-nreject,8,fbuf)
c
      id_wd =uis$create_window(id_vd,'sys$workstation',
     +'*** View CPK ***',0.0,0.0,xmax,ymax,xmax,ymax,w_attributes)
c
      read(5,1099) junk
 1099 format(a80)
c
      call uis$delete_display(id_vd)
c
      stop
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE READVLT(COLUNIT,ID_VD)
C
C     Reads a VLT from disc, then loads to a Workstation colour map &
C     initialises the workstation for multi-colours.
C
      INCLUDE 'sys$library:uisentry'
      INCLUDE 'SYS$LIBRARY:UISUSRDEF'
C
      INTEGER COLUNIT,ID_VD
      INTEGER*4 LEV(256),IGB(256),IR(256)
      REAL RED(256),GREEN(256),BLUE(256)
      CHARACTER*80 TITLE
C
C     Open the VLT and read in the data
C
      OPEN(UNIT=COLUNIT,STATUS='OLD',READONLY,SHARED)
C
      READ(COLUNIT,10) TITLE
   10 FORMAT(A80)
      READ(COLUNIT,20) NENT
   20 FORMAT(8I10)
      READ(COLUNIT,20,ERR=90,END=91) (IGB(I),I=1,NENT)
      READ(COLUNIT,20,ERR=90,END=91) (IR(I),I=1,NENT)
      READ(COLUNIT,20,ERR=90,END=91) (LEV(I),I=1,NENT)
C
C     Unpack the GB array, convert to reals and scale
C
      DO 30 J=1,NENT
         RED(J)   = IR(J)
         GREEN(J) = INT(IGB(J)/256)
         BLUE(J)  = IGB(J)- (256 * GREEN(J))      
         IF (GREEN(J).LT.0.0) GREEN(J) = GREEN(J) + 256.0
         IF (BLUE(J).LT.0.0) BLUE(J) = BLUE(J) + 256.0
         RED(J)=RED(J)/256.0
         GREEN(J)=GREEN(J)/256.0
         BLUE(J)=BLUE(J)/256.0
   30 CONTINUE
C
C     Set up the colours
C
      CALL UIS$SET_COLORS(ID_VD,0,256,RED,GREEN,BLUE)
      CALL UIS$SET_WRITING_MODE(ID_VD,0,1,UIS$C_MODE_COPY)
C
      RETURN
C
   90 STOP 'Error in colour file'
   91 STOP 'Unexpected end of colour file'
      END
