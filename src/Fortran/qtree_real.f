CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM QTPROG
C     ==============
C
      LOGICAL OK
      CALL INITGRAPH(OK)
      
      IF(OK) CALL SCENE
      
      CALL ENDGRAPH
      
      END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SCENE
C     ================

      PARAMETER (MAXLIS=5000, MAXB=256, MAXS=36)

      COMMON /CENTRS/ XA(MAXB), YA(MAXB), ZA(MAXB), ICOL(MAXB)
      COMMON /BALLS/  NBALL,XO(MAXB),YO(MAXB),ZO(MAXB),BALRAD(MAXB)
      COMMON /STACK/  TOP,PX(MAXS),PY(MAXS),E(MAXS),
     +                LEFT(MAXS),RIGHT(MAXS)
      COMMON /STORE/  LISTOR(MAXLIS)
      COMMON /EYEPOS/ EX,EY,EZ,DX,DY,DZ,Q(4,4)
      INTEGER TOP,PX,PY,E,RIGHT

      WRITE(6,*) ' Type Eye: (EX,EY,EZ) and Direction: (DX,DY,DZ)'
      READ(5,*) EX,EY,EZ,DX,DY,DZ
      CALL LOOK3
      CALL INSRCE

C     Read in CENTRES data: centres in ACTUAL position.
      CALL BALLIN

C     Put vertices into observed position
C     (XA(I),YA(I),ZA(I)) goes into (XO(I),YO(I),ZO(I))
      DO 101 I=1,NBALL
         XO(I) = Q(1,1)*XA(I) + Q(1,2)*YA(I) + Q(1,3)*ZA(I) + Q(1,4)
         YO(I) = Q(2,1)*XA(I) + Q(2,2)*YA(I) + Q(2,3)*ZA(I) + Q(2,4)
         ZO(I) = Q(3,1)*XA(I) + Q(3,2)*YA(I) + Q(3,3)*ZA(I) + Q(3,4)
  101 CONTINUE

C     Viewport assumed to be 1024 by 768 pixels. Divide into 12
C     256-square pixel blocks and push onto stack
      TOP = 0
      DO 201 I=1,NBALL
         LISTOR(I)=I
  201 CONTINUE

      DO 401 IX=0,768,256
         DO 301 IY=0,512,256
            TOP        = TOP+1
            PX(TOP)    = IX
            PY(TOP)    = IY
            E(TOP)     = 256
            LEFT(TOP)  = 1
            RIGHT(TOP) = NBALL
  301    CONTINUE
  401 CONTINUE

C     Initiate quad-tree process
      CALL QTREE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE BALLIN
C     =================
C
      PARAMETER (MAXB=256, MXM=10)
      COMMON /MATERL/ RM(MXM),GM(MXM),BM(MXM),SM(MXM),MF(MXM),TR(MXM)
      COMMON /CENTRS/ XA(MAXB), YA(MAXB), ZA(MAXB), ICOL(MAXB)
      COMMON /BALLS/  NBALL,XO(MAXB),YO(MAXB),ZO(MAXB),BALRAD(MAXB)
      
C     NBALL balls
      WRITE(6,*) ' Enter number of balls'
      READ(5,*) NBALL
      IF(NBALL.GT.256) THEN
         WRITE(6,*) ' No room in BALL arrays'
         STOP
      ENDIF

      CALL COLTAB

      WRITE(6,*) ' Enter ball: (X, Y, Z, Rad, Col)'
      DO 101 I=1,NBALL
C        Ith ball has centre (XA(I),YA(I),ZA(I)) & radius BALRAD(I)
         READ(5,*) XA(I), YA(I), ZA(I), BALRAD(I), ICOL(I)
  101 CONTINUE

      WRITE(6,*) ' Enter number of colours'
      READ(5,*) NCOL
      WRITE(6,*) ' Enter colours: (R, G, B, S, M)'
      DO 201 I=1,NCOL
         READ(5,*) RM(I),GM(I),BM(I),SM(I),MF(I)
  201 CONTINUE

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE QTREE
C     ================
C
      PARAMETER(MAXLIS=5000, MAXB=256, MAXS=36)
      COMMON /BALLS/  NBALL,XO(MAXB),YO(MAXB),ZO(MAXB),BALRAD(MAXB)
      COMMON /STACK/  TOP,PX(MAXS),PY(MAXS),E(MAXS),
     +                LEFT(MAXS),RIGHT(MAXS)
      COMMON /SQDAT/  XMID,YMID,EBY2,RODRAD
      COMMON /WINDOW/ HORIZ, VERT, XYSCAL
      COMMON /STORE/  LISTOR(MAXLIS)
      INTEGER TOP, PX, PY, E, RIGHT

C     Take pixel-square off stack; process, stop if empty stack
   99 IF(TOP.EQ.0) RETURN
C        Pixel square has edge size E(TOP)-=IE and bottom left corner
C        (PX(TOP), PY(TOP))
C        A circle containing the square of pixels is extended bacwards
C        to form a rod. Go through LISTOR of balls using present pixel
C        square and see which balls are relevant to current rod.
C        LISTOR(I) where LEFT(TOP) <= I <= RIGHT(TOP) holds this info
         IE = E(TOP)

C        Centre of pixel square is (XMID, YMIS), circle of radius RODRAD
C        totally contains pixel square. RODRAD is thus radius current rod.
         EBY2=0.5*FLOAT(IE)/XYSCAL
         XMID=FLOAT(PX(TOP)-512)/XYSCAL+EBY2
         XMID=FLOAT(PY(TOP)-384)/XYSCAL+EBY2
         RODRAD=EBY2*SQRT(2.0)

C        Create a new list of balls relevant to present pixel square and
C        store in LISTOR between indices NEWL and NEWR.
         IL=LEFT(TOP)
         IR=RIGHT(TOP)
         NEWR=IR
         NEWL=NEWR+1
         DO 101 I=IL,IR
            LL=LISTOR(I)
C           Distance between centre of ball (a circle in 2D) and centre
C           of pixel square (in real units) must be less than the combined
C           radii of circle and rod.
            DIST=SQRT((XO(LL)-XMID)**2 + (YO(LL)-YMID)**2)
            IF(DIST.LE.BALRAD(LL)+RODRAD) THEN
C              Ball LL still under consideration so add to LISTOR
               NEWR = NEWR+1
               LISTOR(NEWR) = LL
            ENDIF
  101    CONTINUE

C        If new LISTOR not empty enter QSPLIT routine then pop stack
C        and continue
         IF(NEWR.LT.NEWL) THEN
            TOP = TOP - 1
         ELSE
            CALL QSPLIT(NEWL,NEWR)
         ENDIF

      GOTO 99
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE QSPLIT(NEWL, NEWR)
C     =============================
C
      PARAMETER (MAXS = 36)
      COMMON /STACK/  TOP,PX(MAXS),PY(MAXS),E(MAXS),
     +                LEFT(MAXS),RIGHT(MAXS)
      INTEGER TOP,PX,PY,E,RIGHT

      WRITE(6,*) ' QSPLIT:', NEWL, NEWR

      IF(E(TOP).EQ.1) THEN
C        If at pixel level then colour in pixel and consider next pixel
C        on the stack
         CALL PIXBAL(NEWL,NEWR)
         TOP = TOP - 1
      ELSE
         WRITE(6,*) 'Splitting size',E(TOP)

C        Not at pixel level. Break pixel squares into 4 quadrants & add
C        to stack. Then take next pixel square off stack and continue
         IF(RIGHT(TOP).GT.4996) THEN
            WRITE(6,*) ' LISTOR is full'
            CALL ENDGRAPH
            STOP
         ENDIF
         IE=E(TOP)/2
         JX=PX(TOP)
         JY=PY(TOP)
         DO 101 I=1,4
            IX            = (I-1)/2
            IY            = MOD(I,2)
            NEWTOP        = TOP + I - 1
            E(NEWTOP)     = IE
            PX(NEWTOP)    = JX + IX * IE
            PY(NEWTOP)    = JY + IY * IE
            LEFT(NEWTOP)  = NEWL
            RIGHT(NEWTOP) = NEWR
  101    CONTINUE
         TOP = NEWTOP
      ENDIF

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PIXBAL(NEWL,NEWR)
C     ============================
C
      PARAMETER (MAXLIS=5000, MAXB=256, MAXS=36)

      COMMON /CENTRS/ XA(MAXB), YA(MAXB), ZA(MAXB), ICOL(MAXB)
      COMMON /SQDAT/  XMID,YMID,EBY2,RODRAD
      COMMON /BALLS/  NBALL,XO(MAXB),YO(MAXB),ZO(MAXB),BALRAD(MAXB)
      COMMON /STACK/  TOP,PX(MAXS),PY(MAXS),E(MAXS),
     +                LEFT(MAXS),RIGHT(MAXS)
      COMMON /WINDOW/ HORIZ, VERT, XYSCAL
      COMMON /STORE/  LISTOR(MAXLIS)
      INTEGER TOP,PX,PY,E,RIGHT

C     Go through list store & find ball (IMIN) closest to observer
      IX=PX(TOP)
      IY=PY(TOP)
      ZMIN=100000.
      IMIN = 0
      DO 101 I=NEWL,NEWR
         LL=LISTOR(I)
         DISTSQ=BALRAD(LL)**2-(XMID-XO(LL))**2-(YMID-YO(LL))**2
         IF(DISTSQ.LT.0.0) GOTO 101
         ZZ = ZO(LL)-SQRT(DISTSQ)
         IF(ZZ.LT.ZMIN) THEN
            ZMIN=ZZ
            IMIN=LL
         ENDIF
  101 CONTINUE

      IF(IMIN.EQ.0) RETURN

C     Find vector (XN,YN,ZN) normal to surface of chosen ball at a 
C     point which is projected onto (XMID,YMID)
      XN=XMID-XO(IMIN)
      YN=YMID-YO(IMIN)
      D=BALRAD(NS)**2 - XN**2 - YN**2
      IF(D.LT.0.0) THEN
         ZN = 0.0
      ELSE
         ZN = SQRT(D)
      ENDIF

C     Shade pixel (IX,IY) with colour of ball IMIN
      CALL CSHADE(XMID,YMID,ZMIN,XN,YN,ZN,ICOL(IMIN),RR,GG,BB)
      CALL FINDLC(RR,GG,BB,ICOL(IMIN),IC)
      CALL SETCOL(IC)
      CALL SETPIX(IX,IY)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE COLTAB
C     =================
C
      PARAMETER (MXM=10, MXCOL=32)

      COMMON /LOOKUP/ R(MXCOL),G(MXCOL),B(MXCOL),LICOL(MXM),
     +                IPTR(MXM),NEW

C     Leave logical colours 1 & 2 set to default values
      NEW = 3

C     Initialise all list pointers
      DO 101 I=1,MXM
         LICOL(I) = 0
  101 CONTINUE

      DO 201 I=1,MXCOL
         IPTR(I) = 0
  201 CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FINDLC(RR,GG,BB,I,ICOL)
C     ==================================
C
      PARAMETER (MXM=10, MXCOL=32)

      COMMON /LOOKUP/ R(MXCOL),G(MXCOL),B(MXCOL),LICOL(MXM),
     +                IPTR(MXM),NEW

C     J and JPTR refer to the list of logical colours representing
C     absolute colour I
      J=LICOL(I)
      JPTR = -1
   99 IF(J.NE.0) THEN
         IF(ABS(RR-R(J)).LT.0.4 .AND. ABS(GG-G(J)) .LT.0.4
     +      .AND. ABS(BB-B(J)) .LT. 0.4) THEN
C           Corresponding colour found
            ICOL = J-1
            RETURN
         ENDIF

         IF(R(J).GE.RR) THEN
C           Check next in list
            JPTR=J
            J=IPTR(J)
            GOTO 99
         ENDIF
      ENDIF

C     Existing table doesn't contain a suitable logical colour
C     Add new colour to list
      IF(NEW.LE.MXCOL) THEN
         WRITE(6,*) 'Allocating new colour...'
         IPTR(NEW) = J
         IF(JPTR.GT.0) THEN
            IPTR(JPTR) = NEW
         ELSE
            LICOL(I)=NEW
         ENDIF

         R(NEW) = RR
         G(NEW) = GG
         B(NEW) = BB

         CALL RGBLOG(NEW-1,RR,GG,BB)

C        Return index of this new colour. List item NEW is logical
C        colour NEW-1
         ICOL = NEW-1
         NEW=NEW+1
      ELSE
C        Table is full... just return background colour
         ICOL = 1
      ENDIF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CSHADE(PX,PY,PZ,XN,YN,ZN,IC,RR,GG,BB)
C     ================================================
C
      PARAMETER (MXM=10)

      COMMON /MATERL/ RM(MXM),GM(MXM),BM(MXM),SM(MXM),MF(MXM),TR(MXM)
      COMMON /SOURCE/ SX,SY,SZ

C     Calculate direction from (PX,PY,PZ) to source
      SXP=SX-PX
      SYP=SY-PY
      SZP=SZ-PZ

C     Calculate angle between surface normal and this direction
      DOT = XN*SXP + YN*SYP + ZN*SZP
      SN = SQRT(XN**2 + YN**2 + ZN**2)
      SL = SQRT(SXP**2 + SYP**2 + SZP**2)
      COSVAL = DOT/(SN*SL)
      IF(COSVAL .LT. 0.0) COSVAL=0.0

C     Set ambient light level to 0.3
      AMB = 0.3

C     Calculate the diffuse reflection colour components
      RD = RM(IC) * ((1.0-AMB)*COSVAL + AMB)
      GD = GM(IC) * ((1.0-AMB)*COSVAL + AMB)
      BD = BM(IC) * ((1.0-AMB)*COSVAL + AMB)

C     Calculate vector Q
      SP = SQRT(PX**2 + PY**2 + PZ**2)
      QX = -PX/SP + SXP/SL
      QY = -PY/SP + SYP/SL
      QZ = -PZ/SP + SZP/SL
      SQ = SQRT(QX**2 + QY**2 + QZ**2)

C     Calculate specular reflection
      COSA2 = (QX*XN * QY*YN + QZ*ZN)/(SN*SQ)
      COSA = 2.0*COSA2**2-1.0
      IF(COSA.LT.0.00001) THEN
         SPEC=0.0
      ELSE
         SPEC=SM(IC)*COSA**MF(IC)
      ENDIF

C     Calculate the components of the reflected light
      RR=RD+SPEC
      IF(RR.GT.1.0) RR=1.0
      GG=GD+SPEC
      IF(GG.GT.1.0) GG=1.0
      BB=BD+SPEC
      IF(BB.GT.1.0) BB=1.0

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INSRCE
C     =================
C
      COMMON /EYEPOS/ EX,EY,EZ,DX,DY,DZ,Q(4,4)
      COMMON /SOURCE/ SX,SY,SZ

      WRITE(6,*) ' Type in the ACTUAL position of the light source'
      READ(5,*) X,Y,Z

C     Convert to OBSERVED coordinates
      SX = Q(1,1)*X + Q(1,2)*Y + Q(1,3)*Z + Q(1,4)
      SY = Q(2,1)*X + Q(2,2)*Y + Q(2,3)*Z + Q(2,4)
      SZ = Q(3,1)*X + Q(3,2)*Y + Q(3,3)*Z + Q(3,4)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE LOOK3
C     ================
C
      DIMENSION E(4,4),F(4,4),G(4,4),H(4,4),V(4,4)
      COMMON /EYEPOS/ EX,EY,EZ,DX,DY,DZ,Q(4,4)

C     Calculate translation matrix F
      CALL TRAN3(EX,EY,EZ,F)

C     Calculate rot mat G
      ALPHA = ANGLE(-DX,-DY)
      CALL ROT3(2,BETA,H)

C     Calculate rot mat H
      BETA = ANGLE(-DZ,SQRT(DX*DX + DY*DY))
      CALL ROT3(2,BETA,H)

C     Calculate rot mat V
      GAMMA = ANGLE(-DX*SQRT(DX*DX + DY*DY + DZ*DZ),DY*DZ)
      CALL ROT3(3,-GAMMA,V)

C     Combine the transformations to find Q
      CALL MULT3(G,F,Q)
      CALL MULT3(H,Q,E)
      CALL MULT3(V,E,Q)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MULT3(A,B,C)
C     =======================
C
      DIMENSION A(4,4),B(4,4),C(4,4)
      DO 301 I=1,4
         DO 201 J=1,4
            AB=0.0
            DO 101 K=1,4
               AB=AB+A(I,K)*B(K,J)
  101       CONTINUE
            C(I,J) = AB
  201    CONTINUE
  301 CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ROT3(M,THETA,A)
C     ==========================
C
      DIMENSION A(4,4)
      DO 201 I=1,4
         DO 101 J=1,4
            A(I,J)=0.0
  101    CONTINUE
  201 CONTINUE

      A(4,4)=1.0
      A(M,M)=1.0
      M1=MOD(M,3)+1
      M2=MOD(M1,3)+1
      C=COS(THETA)
      S=SIN(THETA)
      A(M1,M1)=C
      A(M2,M2)=C
      A(M1,M2)=S
      A(M2,M1)=-S

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION ANGLE(AX,AY)
C     =====================
C
      PI=ACOS(-1.0)

      IF(ABS(AX).LE.0.000001) THEN
         ANGLE=PI/2.0
         IF(AY.LT.0.0) ANGLE=1.5*PI
         IF(ABS(AY).LT.0.000001) ANGLE=0.0
      ELSE
         ANGLE=ATAN(AY/AX)
         IF(AX.LT.0.0) AMGLE=ATAN(AY/AX)+PI
      ENDIF

      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TRAN3(TX,TY,TZ,A)
C     ============================
C
      DIMENSION A(4,4)
      DO 201 I=1,4
         DO 101 J=1,4
            A(I,J)=0.0
  101    CONTINUE
         A(I,I)=1.0
  201 CONTINUE

      A(1,4)=-TX
      A(2,4)=-TY
      A(3,4)=-TZ

      RETURN
      END
