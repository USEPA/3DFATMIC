c
c
c
C ***** TRACK3D2.F    6/19/96
      SUBROUTINE KGLOC
     I    (XPFG,CPFG,ISE,KIT,MXNPFG,MXKGL,NO2,MXNCC,XP,YP,ZP,M1,M2,
     >     KINIT,KEND,EPSX, KGL,NODE,CQ,ICODE)
C  11/17/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG,MXNCC),ISE(MXKGL,NO2)
      DIMENSION XXX(8),YYY(8),ZZZ(8),XWW(8),YWW(8),ZWW(8),CQ(7)
      DIMENSION KIT(4,6,3)
C
      ICODE=-1
      XPP=XP
      YPP=YP
      ZPP=ZP
      DO 105 MP=M1,M2
        XP=XPP
        YP=YPP
        ZP=ZPP
        CALL ELENOD(ISE(MP,5),ISE(MP,7), NODE,LS,IK)
        DO II=1,NODE
          IEM=ISE(MP,II)
          XXX(II)=XPFG(IEM,1)
          YYY(II)=XPFG(IEM,2)
          ZZZ(II)=XPFG(IEM,3)
          IF(DABS(XP-XXX(II)).LE.EPSX .AND.
     >       DABS(YP-YYY(II)).LE.EPSX .AND.
     >       DABS(ZP-ZZZ(II)).LE.EPSX)THEN
            DO KK=KINIT,KEND
              CQ(KK)=CPFG(IEM,KK)
            ENDDO
            XP=XXX(II)
            YP=YYY(II)
            ZP=ZZZ(II)
            KGL=MP
            ICODE=0
            RETURN
          ENDIF
        ENDDO
C
        DO 100 III=1,LS
          DO II1=1,3
            JJ=KIT(II1,III,IK)
            XWW(II1)=XXX(JJ)
            YWW(II1)=YYY(JJ)
            ZWW(II1)=ZZZ(JJ)
          ENDDO
          CANG=FCOS(XWW,YWW,ZWW,XP,YP,ZP,0)
          IF(LS.EQ.5 .AND. III.EQ.5)CANG=-CANG
          IF(LS.EQ.6 .AND. (III.EQ.2 .OR. III.EQ.3 .OR. III.EQ.6))
     >      CANG=-CANG
          IF(CANG.LT.0.0D0 .AND. DABS(CANG).GT.EPSX)GOTO 105
          IF(DABS(CANG).LE.EPSX)THEN
            CALL ONPLAN
     I         (XWW,YWW,ZWW,1,2,3,8,
     M          XP,YP,ZP)
          ENDIF
  100   CONTINUE
        KGL=MP
        RETURN
  105 CONTINUE
C
      RETURN
      END
c
c
c
      SUBROUTINE WRKARY(IE,X,VX,MAXEL,NO2,MAXNP,M,NODE,MXNODE,
     O                  XX,YY,ZZ,VXX,VYY,VZZ)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,NO2),X(MAXNP,3),VX(MAXNP,3)
      DIMENSION XX(MXNODE),YY(MXNODE),ZZ(MXNODE)
      DIMENSION VXX(MXNODE),VYY(MXNODE),VZZ(MXNODE)
C
          DO I1=1,NODE
            II=IE(M,I1)
            XX(I1)=X(II,1)
            YY(I1)=X(II,2)
            ZZ(I1)=X(II,3)
            VXX(I1)=VX(II,1)
            VYY(I1)=VX(II,2)
            VZZ(I1)=VX(II,3)
          ENDDO
C
      RETURN
      END
C
c
c
      SUBROUTINE BDYPLN
     I  (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,ID,
     I   XW,YW,ZW,
     M   N1,N2,N3,N4)
C  10/19/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IBW(MXNPW),IE(MAXEL,11),KIT(4,LS),XX(8),YY(8),ZZ(8)
      DIMENSION XXX(8),YYY(8),ZZZ(8),NO(4)
      DIMENSION XW(MXNPW),YW(MXNPW),ZW(MXNPW),KK(3),KNC(3)
C
      IF(N4.EQ.0)THEN
        IBWN4=0
      ELSE
        IBWN4=IBW(N4)
      ENDIF
C
      DO 100 I=1,LS
        IF(ID.EQ.1)THEN
          IF(IBW(N1).EQ.I .OR. IBW(N2).EQ.I .OR. IBW(N3).EQ.I .OR.
     >       IBWN4.EQ.I)THEN
            GOTO 20
          ELSE
            GOTO 100
          ENDIF
        ELSE
          IF(IBW(N1).EQ.I .OR. IBW(N2).EQ.I)THEN
            GOTO 20
          ELSE
            GOTO 100
          ENDIF
        ENDIF
C
  20    CONTINUE
        DO II=1,4
          KKK=KIT(II,I)
          IF(KKK.EQ.0)THEN
            NO(II)=0
          ELSE
            NO(II)=IE(M,KKK)
          ENDIF
        ENDDO
        N1=NO(1)
        N2=NO(2)
        N3=NO(3)
        N4=NO(4)
        RETURN
  100 CONTINUE
C
      KN=0
      DO 200 I=1,LS
        DO II=1,4
          KKK=KIT(II,I)
          IF(KKK.EQ.0)THEN
            NO(II)=0
          ELSE
            NO(II)=IE(M,KKK)
            XXX(II)=XX(KKK)
            YYY(II)=YY(KKK)
            ZZZ(II)=ZZ(KKK)
          ENDIF
        ENDDO
C
C *** CHECK IF POINT P IS ON THE PLANE COMPOSED WITH KK1,KK2,KK3,KK4
C
        CANG=FCOS(XXX,YYY,ZZZ,XP,YP,ZP,0)
        IF(DABS(CANG).LE.EPSX)THEN
          KN=KN+1
          KNC(KN)=I
        ENDIF
  200 CONTINUE
C
      IF(KN.EQ.1 .OR. (KN.GT.1 .AND. ID.EQ.2))THEN
        I=KNC(1)
        GOTO 340
      ENDIF
      IF(KN.GT.1)THEN
        KK(1)=N1
        KK(2)=N2
        KK(3)=N3
        DO 300 J1=1,KN
          KS=KNC(J1)
          DO J=1,3
            K=KIT(J,KS)
            XXX(J)=XX(K)
            YYY(J)=YY(K)
            ZZZ(J)=ZZ(K)
          ENDDO
          DO K=1,3
            K1=KK(K)
            CANG=FCOS(XXX,YYY,ZZZ,XW(K1),YW(K1),ZW(K1),0)
            IF(DABS(CANG).GT.EPSX)GOTO 300
          ENDDO
          I=KS
          GOTO 340
  300   CONTINUE
      ENDIF
C
      WRITE(16,*)'NOTE ||| SOMETHING WRONG WITH YOUR DATA. STOP,',
     >           ' AND CHECK IT PLEASE --- IN BDYPLN.'
      STOP
C
  340 CONTINUE
      DO II=1,4
        KKK=KIT(II,I)
        IF(KKK.EQ.0)THEN
          NO(II)=0
        ELSE
          NO(II)=IE(M,KKK)
          XXX(II)=XX(KKK)
          YYY(II)=YY(KKK)
          ZZZ(II)=ZZ(KKK)
        ENDIF
      ENDDO
      N1=NO(1)
      N2=NO(2)
      N3=NO(3)
      N4=NO(4)
      CALL ONPLAN
     I     (XXX,YYY,ZZZ,1,2,3,8,
     M      XP,YP,ZP)
C
      RETURN
      END
C       10/4/93      6:30 pm
C
c
      SUBROUTINE CKCNEL
     I                 (LRL,NLRL,M,N1,N2,N3,MAXNP,MXKBD,
     O                  ME1,KOUNT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION LRL(MXKBD,MAXNP),NLRL(MAXNP)
C
      DO 405 M1=1,NLRL(N1)
        ME1=LRL(M1,N1)
        DO 404 M2=1,NLRL(N2)
          ME2=LRL(M2,N2)
          DO 403 M3=1,NLRL(N3)
            ME3=LRL(M3,N3)
            IF(ME3.EQ.ME2 .AND. ME2.EQ.ME1 .AND.
     1         ME1.NE.M)THEN
              KOUNT=2
              GOTO 500
            ENDIF
  403     CONTINUE
  404   CONTINUE
  405 CONTINUE
      ME1=0
      KOUNT=1
  500 RETURN
      END
c
c
c
      SUBROUTINE CKCOIN
     I            (X,Y,Z,XP,YP,ZP,NLRL,NODE,N1,N2,N3,N4,MAXNP,EPSX,
     O             NN,NLRLN)
C       10/22/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP),Y(MAXNP),Z(MAXNP),NLRL(MAXNP)
C
      EPSXX=10.0D0*EPSX
C
      IF(DABS(XP-X(N1)).LE.EPSXX .AND. DABS(YP-Y(N1)).LE.EPSXX .AND.
     >   DABS(ZP-Z(N1)).LE.EPSXX)THEN
        NN=N1
        NLRLN=NLRL(N1)
      ELSEIF(DABS(XP-X(N2)).LE.EPSXX .AND. DABS(YP-Y(N2)).LE.EPSXX .AND.
     >       DABS(ZP-Z(N2)).LE.EPSXX)THEN
        NN=N2
        NLRLN=NLRL(N2)
      ELSEIF(DABS(XP-X(N3)).LE.EPSXX .AND. DABS(YP-Y(N3)).LE.EPSXX .AND.
     >       DABS(ZP-Z(N3)).LE.EPSXX)THEN
        NN=N3
        NLRLN=NLRL(N3)
      ELSEIF(NODE.EQ.8 .OR. (NODE.EQ.6 .AND. N4.NE.0))THEN
        IF(DABS(XP-X(N4)).LE.EPSXX .AND. DABS(YP-Y(N4)).LE.EPSXX .AND.
     >     DABS(ZP-Z(N4)).LE.EPSXX)THEN
          NN=N4
          NLRLN=NLRL(N4)
        ELSE
          NN=0
        ENDIF
      ELSE
        NN=0
      ENDIF
      RETURN
      END
c
c
c
      SUBROUTINE CKSIDE
     I                 (X,Y,Z,XP,YP,ZP,NODE,N1,N2,N3,N4,MAXNP,EPSX,
     O                  KON,J1,J2)
C        10/22/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP),Y(MAXNP),Z(MAXNP)
C
      EPSXX=10.0D0*EPSX
C
      IF(NODE.EQ.8)II=4
      IF(NODE.EQ.6)THEN
        IF(N4.EQ.0)THEN
          II=3
        ELSE
          II=4
        ENDIF
      ENDIF
      IF(NODE.EQ.4)II=3
      DO 401 I=1,II
        IF(I.EQ.1)THEN
          J1=N1
          J2=N2
        ELSEIF(I.EQ.2)THEN
          J1=N2
          J2=N3
        ELSEIF(I.EQ.3 .AND. NODE.EQ.8)THEN
          J1=N3
          J2=N4
        ELSEIF(I.EQ.3 .AND. NODE.EQ.6 .AND. II.EQ.3)THEN
          J1=N3
          J2=N1
        ELSEIF(I.EQ.3 .AND. NODE.EQ.6 .AND. II.EQ.4)THEN
          J1=N3
          J2=N4
        ELSEIF(I.EQ.3 .AND. NODE.EQ.4)THEN
          J1=N3
          J2=N1
        ELSE
          J1=N4
          J2=N1
        ENDIF
C
        A1=X(J2)-X(J1)
        A2=Y(J2)-Y(J1)
        A3=Z(J2)-Z(J1)
        B1=XP-X(J1)
        B2=YP-Y(J1)
        B3=ZP-Z(J1)
        C1=A2*B3-A3*B2
        C2=A3*B1-A1*B3
        C3=A1*B2-A2*B1
        D1=SQRT(C1**2+C2**2+C3**2)
        D2=SQRT(A1**2+A2**2+A3**2)
        IF(D1.LE.EPSXX*D2)THEN
          KON=1
          GOTO 500
        ENDIF
 401  CONTINUE
C
      KON=0
C
 500  CONTINUE
C
      RETURN
      END
c
c
c
      SUBROUTINE INTERP
     I       (MAXNP,MAXEL,MXNPFG,MXKGL,NO2,MXNCC,NODE,M,XQ,YQ,ZQ,X,CP,
     I        IE,XPFG,CPFG,NEL,NEFG,IDZOOM,ID,ISE,NFGMB,KINIT,KEND,EPSX,
     O        CQ)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP,3),CP(MAXNP,MXNCC),IE(MAXEL,11)
      DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG,MXNCC)
      DIMENSION NFGMB(MAXEL),ISE(MXKGL,NO2),KIT(4,6,3)
      DIMENSION XX(8),YY(8),ZZ(8),CC(8,7),DL(8),DNX(8),CQ(7)
C
      DATA KIT/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
      DO KK=1,NODE
        K=IE(M,KK)
        IF(DABS(XQ-X(K,1)).LE.EPSX .AND. DABS(YQ-X(K,2)).LE.EPSX .AND.
     >     DABS(ZQ-X(K,3)).LE.EPSX)THEN
          DO IK=KINIT,KEND
            CQ(IK)=CP(K,IK)
          ENDDO
          GOTO 900
        ENDIF
      ENDDO
C
      DO J=1,NODE
        IEM=IE(M,J)
        XX(J)=X(IEM,1)
        YY(J)=X(IEM,2)
        ZZ(J)=X(IEM,3)
        DO IK=KINIT,KEND
          CC(J,IK)=CP(IEM,IK)
        ENDDO
      ENDDO
      CALL BASE
c     I         (XX,YY,ZZ,XQ,YQ,ZQ,M,NODE,1,
     I         (XX,YY,ZZ,XQ,YQ,ZQ,NODE,1,
     O          DL,DNX,DNX,DNX)
C
C ***** FOR THE CASES OF NON-SF ELEMENTS
C
      IF(IE(M,ID).LE.0 .OR. IDZOOM.EQ.0)THEN
        DO IK=KINIT,KEND
          CQ(IK)=0.0D0
          DO J=1,NODE
            CQ(IK)=CQ(IK)+DL(J)*CC(J,IK)
          ENDDO
        ENDDO
        GOTO 900
      ENDIF
C
C ***** FOR THE CASES OF SF ELEMENTS
C NOTE: IN THE FOLLOWING, WE CHECK ALL THE POSSIBLE SUBELEMENTS FOR
C       NODE=4, CHECK SUBELEMENTS WITH ISE(MP1,6).NE.0 FOR NODE=6
C       (MP1=0 IF CANNOT FIND OUT THE GOOD ONE UNDER THIS CHECKING),
C       AND CALCULATE MP1 FOR NODE=8 (THIS CALCULATED MP1 MIGHT NOT BE
C       GOOD IF ISE(MP1,8).EQ.0).  THEREFORE, WE MIGHT NEED FURTHER
C       CHECKING FOR BOTH NODE=6 AND NODE=8.
C
      M1=NFGMB(M)+1
      IF(M.EQ.NEL)M2=NEFG
      IF(M.NE.NEL)M2=NFGMB(M+1)
      CALL KGLOC
     I    (XPFG,CPFG,ISE,KIT,MXNPFG,MXKGL,NO2,MXNCC,XQ,YQ,ZQ,M1,M2,
     >     KINIT,KEND,EPSX, KGL,NODE1,CQ,ICODE)
      IF(ICODE.EQ.0)GOTO 900
C
      DO J=1,NODE1
        ISEM=ISE(KGL,J)
        XX(J)=XPFG(ISEM,1)
        YY(J)=XPFG(ISEM,2)
        ZZ(J)=XPFG(ISEM,3)
        DO IK=KINIT,KEND
          CC(J,IK)=CPFG(ISEM,IK)
        ENDDO
      ENDDO
      CALL BASE
c     I         (XX,YY,ZZ,XQ,YQ,ZQ,M,NODE1,1,
     I         (XX,YY,ZZ,XQ,YQ,ZQ,NODE1,1,
     O          DL,DNX,DNX,DNX)
      DO IK=KINIT,KEND
        CQ(IK)=0.0D0
        DO J=1,NODE1
          CQ(IK)=CQ(IK)+DL(J)*CC(J,IK)
        ENDDO
      ENDDO
C
  900 CONTINUE
      RETURN
      END
      FUNCTION FCOS(XX,YY,ZZ,XQ,YQ,ZQ,idi)
C
C -----10/ 8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8)
C
C
      A1=XX(2)-XX(1)
      A2=XX(3)-XX(1)
      AP=XQ-XX(1)
      B1=YY(2)-YY(1)
      B2=YY(3)-YY(1)
      BP=YQ-YY(1)
      C1=ZZ(2)-ZZ(1)
      C2=ZZ(3)-ZZ(1)
      CP=ZQ-ZZ(1)
      AA=B1*C2-B2*C1
      BB=C1*A2-C2*A1
      CC=A1*B2-A2*B1
      FCOS=AA*AP+BB*BP+CC*CP
c
      if(idi.eq.0)then
        dd=dsqrt(aa*aa+bb*bb+cc*cc)
      elseif(idi.eq.1)then
        dd=dsqrt((aa*aa+bb*bb+cc*cc)*(ap*ap+bp*bp+cp*cp))
      endif
      IF(DD.EQ.0.0D0)THEN
        WRITE(16,*)'ERROR OCCURRED AT FCOS, IDI, DD=',IDI,DD,' --- STOP'
        STOP
      ENDIF
      fcos=fcos/dd
c
      RETURN
      END
c
c
c
      SUBROUTINE LOCPLN
     I    (XQ,YQ,ZQ,X,Y,Z,N1,N2,N3,N4,NODE,
     O     XSI,ETA,DL1,DL2,DL3,DL)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(8),Y(8),Z(8)
      DIMENSION XL(4),YL(4)
      DIMENSION DL(4)
C
C ----- TRANSFORM COORDINATES FROM X,Y,Z TO XL,YL
C
      XL(1)=0.0D0
      YL(1)=0.0D0
C
      XL(2)=DSQRT((X(N2)-X(N1))**2+(Y(N2)-Y(N1))**2+(Z(N2)-Z(N1))**2)
      YL(2)=0.0D0
C
      XL(3)=((X(N2)-X(N1))*(X(N3)-X(N1))+(Y(N2)-Y(N1))*(Y(N3)-Y(N1))+
     >       (Z(N2)-Z(N1))*(Z(N3)-Z(N1)))/XL(2)
      Y31=(Y(N2)-Y(N1))*(Z(N3)-Z(N1))-(Y(N3)-Y(N1))*(Z(N2)-Z(N1))
      Y32=(Z(N2)-Z(N1))*(X(N3)-X(N1))-(Z(N3)-Z(N1))*(X(N2)-X(N1))
      Y33=(X(N2)-X(N1))*(Y(N3)-Y(N1))-(X(N3)-X(N1))*(Y(N2)-Y(N1))
      YL(3)=DSQRT(Y31**2+Y32**2+Y33**2)/XL(2)
C
      IF(NODE.EQ.4)THEN
        XL(4)=((X(N2)-X(N1))*(X(N4)-X(N1))+(Y(N2)-Y(N1))*(Y(N4)-Y(N1))+
     >         (Z(N2)-Z(N1))*(Z(N4)-Z(N1)))/XL(2)
        Y41=(Y(N2)-Y(N1))*(Z(N4)-Z(N1))-(Y(N4)-Y(N1))*(Z(N2)-Z(N1))
        Y42=(Z(N2)-Z(N1))*(X(N4)-X(N1))-(Z(N4)-Z(N1))*(X(N2)-X(N1))
        Y43=(X(N2)-X(N1))*(Y(N4)-Y(N1))-(X(N4)-X(N1))*(Y(N2)-Y(N1))
        YL(4)=DSQRT(Y41**2+Y42**2+Y43**2)/XL(2)
      ENDIF
C
C ----- TRANSFORM XQ,YQ,ZQ TO XLQ,YLQ
C
      XLQ=((X(N2)-X(N1))*(XQ-X(N1))+(Y(N2)-Y(N1))*(YQ-Y(N1))+
     >     (Z(N2)-Z(N1))*(ZQ-Z(N1)))/XL(2)
      YQ1=(Y(N2)-Y(N1))*(ZQ-Z(N1))-(YQ-Y(N1))*(Z(N2)-Z(N1))
      YQ2=(Z(N2)-Z(N1))*(XQ-X(N1))-(ZQ-Z(N1))*(X(N2)-X(N1))
      YQ3=(X(N2)-X(N1))*(YQ-Y(N1))-(XQ-X(N1))*(Y(N2)-Y(N1))
      YLQ=DSQRT(YQ1**2+YQ2**2+YQ3**2)/XL(2)
C
C ----- DETERMINE BASE FUNCTION DL AND/OR LOCAL COORDINATES XSI, ETA
C       WITH BASE2D
C
c      M=0
      CALL BASE2D
c     I    (XL,YL,XLQ,YLQ,M,NODE,
     I    (XL,YL,XLQ,YLQ,NODE,
     O     DL,XSI,ETA)
      IF(NODE.EQ.3)THEN
        DL1=DL(1)
        DL2=DL(2)
        DL3=DL(3)
      ENDIF
C
      RETURN
      END
C
c
c
      SUBROUTINE LOCQ2N
     I                (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I                 N1,N2,SDT,IBF,IJUDGE,EPSX,
     O                 XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
C -----  2/19/96
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C NOTE: IF IBF=1, BACKWARD TRACKING.                                   %
C       IF IBF=2, FORWARD TRACKING.                                    %
C       IF IJUDGE=1, USING THE AVERAGE-VELOCITY APPROACH.              %
C       IF IJUDGE=2, USING THE SINGLE-VELOCITY APPROACH.               %
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8)
C
      DATA NITER/150/, EPSR/1.0D-6/, OMEGA/1.0D0/
C
      SDT1=SDT
      MXNOD=8
      NSTART=N1
      NEND=N2
      NDIM=3
      NJUMP=N2-N1
C
      VXQ=0.0D0
      VYQ=0.0D0
      VZQ=0.0D0
C
C NOTE: THE COMPUTATION IN THIS SUBROUTINE IS ACHIEVED BASED ON
C       A FORWARD PARTICLE TRACKING.  THUS, VELOCITY SIGN NEEDS
C       TO BE CHANGED IF IBF=1.
C
      IF(IBF.EQ.1)THEN
        CALL CHNGSN
     I       (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M        VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
      ENDIF
C
C ***** WHEN THE SINGLE VELOCITY APPROACH IS USED, CALL SUR2D3
C
      IF(IJUDGE.EQ.2)THEN
        CALL SUR2D3
     I     (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,XP,YP,ZP,VXP,VYP,VZP,SDT,EPSX,
     O      XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I        (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M         VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
        RETURN
      ENDIF
C
C $$$$$ DETERMINE THE LOCATION OF Q WITH DIFFERENT APPROACHES.
C
C ***** CALCULATE THE COEFFICIENTS OF XO-XP
C
      A1=XX(N2)-XP
      A2=XX(N1)-XX(N2)
C
C ***** CALCULATE THE COEFFICIENTS OF YQ-YP
C
      B1=YY(N2)-YP
      B2=YY(N1)-YY(N2)
C
C ***** CALCULATE THE COEFFICIENTS OF ZQ-ZP
C
      C1=ZZ(N2)-ZP
      C2=ZZ(N1)-ZZ(N2)
C
      IF(IJUDGE.EQ.1)THEN
C
C ***** CALCULATE THE COEFFICIENTS OF VXP+VXQ
C
        D1=VXX(N2)+VXP
        D2=VXX(N1)-VXX(N2)
C
C ***** CALCULATE THE COEFFICIENTS OF VYQ+VYP
C
        E1=VYY(N2)+VYP
        E2=VYY(N1)-VYY(N2)
C
C ***** CALCULATE THE COEFFICIENTS OF VZQ+VZP
C
        F1=VZZ(N2)+VZP
        F2=VZZ(N1)-VZZ(N2)
C
      ELSE
C
C ***** CALCULATE THE COEFFICIENTS OF VXP
C
        D1=2.0D0*VXP
        D2=0.0D0
C
C ***** CALCULATE THE COEFFICIENTS OF VYP
C
        E1=2.0D0*VYP
        E2=0.0D0
C
C ***** CALCULATE THE COEFFICIENTS OF VZP
C
        F1=2.0D0*VZP
        F2=0.0D0
      ENDIF
C
C ***** CONSTRUCT FUNCTION FOF AND FUNCTION GOF FOR ITERATION
C
C ----- CALCULATE THE COEFFICIENTS OF DX
C
C     R1=A1**2+B1**2+C1**2
      R2=2.0D0*(A1*A2+B1*B2+C1*C2)
      R3=A2**2+B2**2+C2**2
C
C ----- CALCULATE THE COEFFICIENTS OF V
C
C     S1=D1**2+E1**2+F1**2
      S2=2.0D0*(D1*D2+E1*E2+F1*F2)
      S3=D2**2+E2**2+F2**2
C
C ***** SOLVE THE ABOVE TWO EQUATIONS BY NEWTON-RAPHSON
C
C ***** START NEWTON-RAPHSON
C
      ITEST=1
  100 CONTINUE
      IF(ITEST.EQ.1)THEN
        CALL REPLAS(A1,A2,D1,D2,A1,A2,U1,U2,U3,U4,U1,U2)
      ELSEIF(ITEST.EQ.2)THEN
        CALL REPLAS(B1,B2,E1,E2,B1,B2,U1,U2,U3,U4,U1,U2)
      ELSE
        CALL REPLAS(C1,C2,F1,F2,C1,C2,U1,U2,U3,U4,U1,U2)
      ENDIF
C
      INIGES=0
  105 CONTINUE
      IF(INIGES.EQ.0)THEN
        XSIO=0.5
      ELSEIF(INIGES.EQ.1)THEN
        XSIO=0.0D0
      ELSEIF(INIGES.EQ.2)THEN
        XSIO=1.0D0
      ENDIF
C
      XSIW=XSIO
      XSI=XSIO
      XID=1.0D0
C
      DO 200 ITER=1,NITER
C
        DN1=XSIW
        DN2=1.0D0-XSIW
C
        XQ=DN1*XX(N1)+DN2*XX(N2)
        YQ=DN1*YY(N1)+DN2*YY(N2)
        ZQ=DN1*ZZ(N1)+DN2*ZZ(N2)
        VXQ=DN1*VXX(N1)+DN2*VXX(N2)
        VYQ=DN1*VYY(N1)+DN2*VYY(N2)
        VZQ=DN1*VZZ(N1)+DN2*VZZ(N2)
C
        RR=(XQ-XP)*(XQ-XP)+(YQ-YP)*(YQ-YP)+(ZQ-ZP)*(ZQ-ZP)
        IF(IJUDGE.EQ.1)THEN
          SS=(VXQ+VXP)*(VXQ+VXP)+(VYQ+VYP)*(VYQ+VYP)+(VZQ+VZP)*(VZQ+VZP)
          IF(SS.EQ.0.0D0)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
            IF(IBF.EQ.1)THEN
              CALL CHNGSN
     I             (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M              VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
            ENDIF
            SDT1=SDT
            RETURN
          ENDIF
        ELSE
          SS=(VXP+VXP)*(VXP+VXP)+(VYP+VYP)*(VYP+VYP)+(VZP+VZP)*(VZP+VZP)
        ENDIF
        RR=DSQRT(RR)
        SS=DSQRT(SS)
        T1=U1+U2*XSIW
        T2=U3+U4*XSIW
C
        FOF=SS*T1-RR*T2
        IF(FOF.EQ.0.0D0)GOTO 500
C
        DFX=SS*U2+0.5D0*T1*(S2+2.0D0*S3*XSIW)/SS
     >      -RR*U4-0.5D0*T2*(R2+2.0D0*R3*XSIW)/RR
C
C ***** COMPUTE FOR NEW XSI
C
        IF(DFX.EQ.0.0D0)THEN
          XSI=XSIO+XID*0.05D0
        ELSE
          XSI=XSIO-FOF/DFX
        ENDIF
C
C ***** TEST CONVERGENCE
C
        DIFMAX=0.0D0
        DIFMAX1=0.0D0
        IF(XSIW.NE.0.)THEN
          DIF=DABS((XSI-XSIW)/XSIW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(XSIW.EQ.0.)THEN
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DIFMAX.LE.EPSR .AND. ITER.GT.1)THEN
          IF(DIFMAX1.LE.EPSR .AND. DABS(FOF).LE.EPSX)GOTO 500
        ENDIF
C
C ***** TAKE A NEW GUESS FOR XSI AND ETA
C
        IF(XSI.GT.1.0D0)THEN
          XSI=1.0D0
          XID=-1.0D0
        ENDIF
        IF(XSI.LT.0.0D0)THEN
          XSI=0.0D0
          XID=-1.0D0
        ENDIF
        XSIW=OMEGA*XSI+(1.0D0-OMEGA)*XSIO
C
C ***** UPDATE XSIO,ETAO
C
        XSIO=XSI
  200 CONTINUE
C
      IF(ITEST.NE.3)THEN
        ITEST=ITEST+1
        GOTO 100
      ENDIF
C
      IF(DIFMAX1.GT.EPSR .OR. DABS(FOF).GT.EPSX)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I         (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M          VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
        SDT1=SDT
        RETURN
      ENDIF
C
C ***** GET THE POSITION AND VELOCITY OF POINT Q
C
  500 CONTINUE
      IF(XSI.GT.1.0D0)XSI=1.0D0
      IF(XSI.LT.0.0D0)XSI=0.0D0
C
      DN1=XSI
      DN2=1.0D0-XSI
C
      XQ=DN1*XX(N1)+DN2*XX(N2)
      YQ=DN1*YY(N1)+DN2*YY(N2)
      ZQ=DN1*ZZ(N1)+DN2*ZZ(N2)
      VXQ=DN1*VXX(N1)+DN2*VXX(N2)
      VYQ=DN1*VYY(N1)+DN2*VYY(N2)
      VZQ=DN1*VZZ(N1)+DN2*VZZ(N2)
C
C ***** CALCULATE DT
C
      RR=(XQ-XP)*(XQ-XP)+(YQ-YP)*(YQ-YP)+(ZQ-ZP)*(ZQ-ZP)
      IF(IJUDGE.EQ.1)THEN
        SS=(VXQ+VXP)*(VXQ+VXP)+(VYQ+VYP)*(VYQ+VYP)+(VZQ+VZP)*(VZQ+VZP)
        IF(SS.EQ.0.0D0)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
          IF(IBF.EQ.1)THEN
            CALL CHNGSN
     I           (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M            VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
          ENDIF
          SDT1=SDT
          RETURN
        ENDIF
      ELSE
        SS=(VXP+VXP)*(VXP+VXP)+(VYP+VYP)*(VYP+VYP)+(VZP+VZP)*(VZP+VZP)
      ENDIF
      RR=DSQRT(RR)
      SS=DSQRT(SS)
      DT=2.0D0*RR/SS
C
C -----CHECK IF A CORRECT ANSWER HAS BEEN OBTAINED
C
      IF(RR.LE.1.0D-4)GOTO 510
      IF(IJUDGE.EQ.1)THEN
        VVX=0.5D0*(VXP+VXQ)
        VVY=0.5D0*(VYP+VYQ)
        VVZ=0.5D0*(VZP+VZQ)
      ELSE
        VVX=VXP
        VVY=VYP
        VVZ=VZP
      ENDIF
      DIFFX=(XQ-XP)-VVX*DT
      DIFFY=(YQ-YP)-VVY*DT
      DIFFZ=(ZQ-ZP)-VVZ*DT
      DIFF=DSQRT(DIFFX*DIFFX+DIFFY*DIFFY+DIFFZ*DIFFZ)
      IF(DIFF.GT.2.0D0*RR*1.0D-4)THEN
        INIGES=INIGES+1
        IF(INIGES.EQ.3)THEN
          IF(ITEST.NE.3)THEN
            ITEST=ITEST+1
            GOTO 100
          ELSE
            IF(IJUDGE.EQ.2)THEN
              CALL SUR2D3
     I     (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,XP,YP,ZP,VXP,VYP,VZP,SDT,EPSX,
     O      XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
            ELSE
              SDT1=SDT
            ENDIF
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
            IF(IBF.EQ.1)THEN
              CALL CHNGSN
     I             (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M              VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
            ENDIF
            RETURN
          ENDIF
        ELSE
          GOTO 105
        ENDIF
      ENDIF
C
  510 CONTINUE
C
C
      IF(DT.GE.SDT)THEN
        SDT1=0.0D0
        XQ=XP+SDT*(XQ-XP)/DT
        YQ=YP+SDT*(YQ-YP)/DT
        ZQ=ZP+SDT*(ZQ-ZP)/DT
        DO I=N1,N2,N2-N1
          VXX(I)=-VXX(I)
          VYY(I)=-VYY(I)
          VZZ(I)=-VZZ(I)
        ENDDO
      ELSE
        SDT1=SDT-DT
        DX=DABS(VXQ*SDT1)
        DY=DABS(VYQ*SDT1)
        DZ=DABS(VZQ*SDT1)
        IF(DX.LE.EPSX .AND. DY.LE.EPSX .AND. DZ.LE.EPSX) SDT1=0.0D0
C
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I         (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M          VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
      ENDIF
C
      RETURN
      END
 
 
 
      SUBROUTINE LOCQ3N
     I                (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I                 N1,N2,N3,SDT,IBF,IJUDGE,NODE,EPSX,
     O                 XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
C -----  2/19/96
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C NOTE: IF IBF=1, BACKWARD TRACKING.                                   %
C       IF IBF=2, FORWARD TRACKING.                                    %
C       IF IJUDGE=1, USING THE AVERAGE-VELOCITY APPROACH.              %
C       IF IJUDGE=2, USING THE SINGLE-VELOCITY APPROACH.               %
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8)
C
      DATA NITER/150/, EPSR/1.0D-6/, OMEGA/1.0D0/
C
      SDT1=SDT
      MXNOD=8
      NSTART=1
      NEND=NODE
      NDIM=3
      NJUMP=1
C
      VXQ=0.0D0
      VYQ=0.0D0
      VZQ=0.0D0
C
C NOTE: THE COMPUTATION IN THIS SUBROUTINE IS ACHIEVED BASED ON
C       A FORWARD PARTICLE TRACKING.  THUS, VELOCITY SIGN NEEDS
C       TO BE CHANGED IF IBF=1.
C
      IF(IBF.EQ.1)THEN
        CALL CHNGSN
     I       (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M        VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
      ENDIF
C
C ***** WHEN THE SINGLE VELOCITY APPROACH IS USED, CALL SUR2D3
C
      IF(IJUDGE.EQ.2)THEN
        N4=0
        CALL SURE3D
     I     (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,N3,N4,XP,YP,ZP,VXP,VYP,VZP,SDT,
     O      XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I        (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M         VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
        RETURN
      ENDIF
C
C $$$$$ DETERMINE THE LOCATION OF Q WITH DIFFERENT APPROACHES.
C
C ***** CALCULATE THE COEFFICIENTS OF XQ-XP
C
      A1=XX(N3)-XP
      A2=XX(N1)-XX(N3)
      A3=XX(N2)-XX(N3)
C
C ***** CALCULATE THE COEFFICIENTS OF YQ-YP
C
      B1=YY(N3)-YP
      B2=YY(N1)-YY(N3)
      B3=YY(N2)-YY(N3)
C
C ***** CALCULATE THE COEFFICIENTS OF ZQ-ZP
C
      C1=ZZ(N3)-ZP
      C2=ZZ(N1)-ZZ(N3)
      C3=ZZ(N2)-ZZ(N3)
C
      IF(IJUDGE.EQ.1)THEN
C
C ***** CALCULATE THE COEFFICIENTS OF VXP+VXQ
C
        D1=VXX(N3)+VXP
        D2=VXX(N1)-VXX(N3)
        D3=VXX(N2)-VXX(N3)
C
C ***** CALCULATE THE COEFFICIENTS OF VYQ+VYP
C
        E1=VYY(N3)+VYP
        E2=VYY(N1)-VYY(N3)
        E3=VYY(N2)-VYY(N3)
C
C ***** CALCULATE THE COEFFICIENTS OF VZQ+VZP
C
        F1=VZZ(N3)+VZP
        F2=VZZ(N1)-VZZ(N3)
        F3=VZZ(N2)-VZZ(N3)
C
      ELSE
C
C ***** CALCULATE THE COEFFICIENTS OF VXP
C
        D1=2.0D0*VXP
        D2=0.0D0
        D3=0.0D0
C
C ***** CALCULATE THE COEFFICIENTS OF VYP
C
        E1=2.0D0*VYP
        E2=0.0D0
        E3=0.0D0
C
C ***** CALCULATE THE COEFFICIENTS OF VZP
C
        F1=2.0D0*VZP
        F2=0.0D0
        F3=0.0D0
      ENDIF
C
C ***** CONSTRUCT FUNCTION FOF AND FUNCTION GOF FOR ITERATION
C
C ----- CALCULATE THE COEFFICIENTS OF DX
C
C     R1=A1**2+B1**2+C1**2
      R2=2.0D0*(A1*A2+B1*B2+C1*C2)
      R3=2.0D0*(A1*A3+B1*B3+C1*C3)
      R4=A2**2+B2**2+C2**2
      R5=2.0D0*(A2*A3+B2*B3+C2*C3)
      R6=A3**2+B3**2+C3**2
C
C ----- CALCULATE THE COEFFICIENTS OF V
C
C     S1=D1**2+E1**2+F1**2
      S2=2.0D0*(D1*D2+E1*E2+F1*F2)
      S3=2.0D0*(D1*D3+E1*E3+F1*F3)
      S4=D2**2+E2**2+F2**2
      S5=2.0D0*(D2*D3+E2*E3+F2*F3)
      S6=D3**2+E3**2+F3**2
C
C ***** SOLVE THE ABOVE TWO EQUATIONS BY NEWTON-RAPHSON
C
C ***** START NEWTON-RAPHSON
C
      ITEST=1
  100 CONTINUE
      IF(ITEST.EQ.1)THEN
        CALL REPLAS(A1,A2,A3,D1,D2,D3,U1,U2,U3,U4,U5,U6)
        CALL REPLAS(B1,B2,B3,E1,E2,E3,V1,V2,V3,V4,V5,V6)
      ELSEIF(ITEST.EQ.2)THEN
        CALL REPLAS(B1,B2,B3,E1,E2,E3,U1,U2,U3,U4,U5,U6)
        CALL REPLAS(C1,C2,C3,F1,F2,F3,V1,V2,V3,V4,V5,V6)
      ELSE
        CALL REPLAS(C1,C2,C3,F1,F2,F3,U1,U2,U3,U4,U5,U6)
        CALL REPLAS(A1,A2,A3,D1,D2,D3,V1,V2,V3,V4,V5,V6)
      ENDIF
C
      INIGES=0
  105 CONTINUE
      IF(INIGES.EQ.0)THEN
        DN1O=0.33333333333333D0
        DN2O=0.33333333333333D0
      ELSEIF(INIGES.EQ.1)THEN
        DN1O=1.0
        DN2O=0.0
      ELSEIF(INIGES.EQ.2)THEN
        DN1O=0.0
        DN2O=1.0
      ELSEIF(INIGES.EQ.3)THEN
        DN1O=0.0
        DN2O=0.0
      ENDIF
C
      DN1W=DN1O
      DN2W=DN2O
      DN1=DN1O
      DN2=DN2O
      DN3W=1.0D0-DN1W-DN2W
C
      XID=1.0D0
      EID=1.0D0
C
      DO 200 ITER=1,NITER
C
        XQ=DN1W*XX(N1)+DN2W*XX(N2)+DN3W*XX(N3)
        YQ=DN1W*YY(N1)+DN2W*YY(N2)+DN3W*YY(N3)
        ZQ=DN1W*ZZ(N1)+DN2W*ZZ(N2)+DN3W*ZZ(N3)
        VXQ=DN1W*VXX(N1)+DN2W*VXX(N2)+DN3W*VXX(N3)
        VYQ=DN1W*VYY(N1)+DN2W*VYY(N2)+DN3W*VYY(N3)
        VZQ=DN1W*VZZ(N1)+DN2W*VZZ(N2)+DN3W*VZZ(N3)
C
        RR=(XQ-XP)*(XQ-XP)+(YQ-YP)*(YQ-YP)+(ZQ-ZP)*(ZQ-ZP)
        IF(IJUDGE.EQ.1)THEN
          SS=(VXQ+VXP)*(VXQ+VXP)+(VYQ+VYP)*(VYQ+VYP)+(VZQ+VZP)*(VZQ+VZP)
          IF(SS.EQ.0.0D0)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
            IF(IBF.EQ.1)THEN
              CALL CHNGSN
     I             (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M              VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
            ENDIF
            SDT1=SDT
            RETURN
          ENDIF
        ELSE
          SS=(VXP+VXP)*(VXP+VXP)+(VYP+VYP)*(VYP+VYP)+(VZP+VZP)*(VZP+VZP)
        ENDIF
        RR=DSQRT(RR)
        SS=DSQRT(SS)
        T1=U1+U2*DN1W+U3*DN2W
        T2=U4+U5*DN1W+U6*DN2W
        T3=V1+V2*DN1W+V3*DN2W
        T4=V4+V5*DN1W+V6*DN2W
C
        FOF=SS*T1-RR*T2
        GOF=SS*T3-RR*T4
        IF(FOF.EQ.0.0D0 .AND. GOF.EQ.0.0D0)GOTO 500
C
        DFX=SS*U2+0.5D0*T1*(S2+2.0D0*S4*DN1W+S5*DN2W)/SS
     >     -RR*U5-0.5D0*T2*(R2+2.0D0*R4*DN1W+R5*DN2W)/RR
        DFE=SS*U3+0.5D0*T1*(S3+S5*DN1W+2.0D0*S6*DN2W)/SS
     >     -RR*U6-0.5D0*T2*(R3+R5*DN1W+2.0D0*R6*DN2W)/RR
        DGX=SS*V2+0.5D0*T3*(S2+2.0D0*S4*DN1W+S5*DN2W)/SS
     >     -RR*V5-0.5D0*T4*(R2+2.0D0*R4*DN1W+R5*DN2W)/RR
        DGE=SS*V3+0.5D0*T3*(S3+S5*DN1W+2.0D0*S6*DN2W)/SS
     >     -RR*V6-0.5D0*T4*(R3+R5*DN1W+2.0D0*R6*DN2W)/RR
        DETJ=DFX*DGE-DFE*DGX
C
C ***** COMPUTE FOR NEW DN1 AND DN2
C
        CALL NEWXE(DETJ,DN1O,DN2O,DFX,DFE,DGX,DGE,XID,EID,FOF,GOF,
     O             0.05D0,  DN1,DN2)
C
C ***** TEST CONVERGENCE
C
        DIFMAX=0.0D0
        DIFMAX1=0.0D0
        IF(DN1W.NE.0.)THEN
          DIF=DABS((DN1-DN1W)/DN1W)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(DN1-DN1W)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DN1W.EQ.0.)THEN
          DIF1=DABS(DN1-DN1W)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DN2W.NE.0.)THEN
          DIF=DABS((DN2-DN2W)/DN2W)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(DN2-DN2W)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DN2W.EQ.0.)THEN
          DIF1=DABS(DN2-DN2W)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DIFMAX.LE.EPSR .AND. ITER.GT.1)THEN
          IF(DIFMAX1.LE.EPSR .AND. DABS(FOF)+DABS(GOF).LE.EPSX)GOTO 500
        ENDIF
C
C ***** TAKE A NEW GUESS FOR XSI AND ETA
C
        IF(DN1.GT.1.0D0)THEN
          DN1=1.0D0
          XID=-1.0D0
        ENDIF
        IF(DN1.LT.0.0D0)THEN
          DN1= 0.0D0
          XID=1.0D0
        ENDIF
        IF(DN2.GT.1.0D0)THEN
          DN2=1.0D0
          EID=-1.0D0
        ENDIF
        IF(DN2.LT.0.0D0)THEN
          DN2= 0.0D0
          EID=1.0D0
        ENDIF
C
        IF(DN1+DN2.GT.1.0D0)DN2=1.0D0-DN1
        DN1W=OMEGA*DN1+(1.0D0-OMEGA)*DN1O
        DN2W=OMEGA*DN2+(1.0D0-OMEGA)*DN2O
        DN3W=1.0D0-DN1W-DN2W
C
C ***** UPDATE XSIO,ETAO
C
        DN1O=DN1
        DN2O=DN2
  200 CONTINUE
C
      IF(ITEST.NE.3)THEN
        ITEST=ITEST+1
        GOTO 100
      ENDIF
C
      IF(DIFMAX1.GT.EPSR .OR. DABS(FOF)+DABS(GOF).GT.EPSX)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I         (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M          VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
        SDT1=SDT
        RETURN
      ENDIF
C
C ***** GET THE POSITION AND VELOCITY OF POINT Q
C
  500 CONTINUE
      IF(DN1.GT.1.0D0)DN1=1.0D0
      IF(DN1.LT.0.0D0)DN1=0.0D0
      IF(DN2.GT.1.0D0)DN2=1.0D0
      IF(DN2.LT.0.0D0)DN2=0.0D0
      DN3=1.0-DN1-DN2
C
      XQ=DN1*XX(N1)+DN2*XX(N2)+DN3*XX(N3)
      YQ=DN1*YY(N1)+DN2*YY(N2)+DN3*YY(N3)
      ZQ=DN1*ZZ(N1)+DN2*ZZ(N2)+DN3*ZZ(N3)
      VXQ=DN1*VXX(N1)+DN2*VXX(N2)+DN3*VXX(N3)
      VYQ=DN1*VYY(N1)+DN2*VYY(N2)+DN3*VYY(N3)
      VZQ=DN1*VZZ(N1)+DN2*VZZ(N2)+DN3*VZZ(N3)
C
C ***** CALCULATE DT
C
      RR=(XQ-XP)*(XQ-XP)+(YQ-YP)*(YQ-YP)+(ZQ-ZP)*(ZQ-ZP)
      IF(IJUDGE.EQ.1)THEN
        SS=(VXQ+VXP)*(VXQ+VXP)+(VYQ+VYP)*(VYQ+VYP)+(VZQ+VZP)*(VZQ+VZP)
        IF(SS.EQ.0.0D0)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
          IF(IBF.EQ.1)THEN
            CALL CHNGSN
     I           (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M            VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
          ENDIF
          SDT1=SDT
          RETURN
        ENDIF
      ELSE
        SS=(VXP+VXP)*(VXP+VXP)+(VYP+VYP)*(VYP+VYP)+(VZP+VZP)*(VZP+VZP)
      ENDIF
      RR=DSQRT(RR)
      SS=DSQRT(SS)
      DT=2.0D0*RR/SS
C
C -----CHECK IF A CORRECT ANSWER HAS BEEN OBTAINED
C
      IF(RR.LE.1.0D-4)GOTO 510
      IF(IJUDGE.EQ.1)THEN
        VVX=0.5D0*(VXP+VXQ)
        VVY=0.5D0*(VYP+VYQ)
        VVZ=0.5D0*(VZP+VZQ)
      ELSE
        VVX=VXP
        VVY=VYP
        VVZ=VZP
      ENDIF
      DIFFX=(XQ-XP)-VVX*DT
      DIFFY=(YQ-YP)-VVY*DT
      DIFFZ=(ZQ-ZP)-VVZ*DT
      DIFF=DSQRT(DIFFX*DIFFX+DIFFY*DIFFY+DIFFZ*DIFFZ)
      IF(DIFF.GT.2.0D0*RR*1.0D-4)THEN
        INIGES=INIGES+1
        IF(INIGES.EQ.4)THEN
          IF(ITEST.NE.4)THEN
            ITEST=ITEST+1
            GOTO 100
          ELSE
            IF(IJUDGE.EQ.2)THEN
              N4=0
              CALL SURE3D
     I       (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,N3,N4,XP,YP,ZP,VXP,VYP,VZP,SDT,
     O        XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
            ELSE
              SDT1=SDT
            ENDIF
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
            IF(IBF.EQ.1)THEN
              CALL CHNGSN
     I             (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M              VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
            ENDIF
            RETURN
          ENDIF
        ELSE
          GOTO 105
        ENDIF
      ENDIF
C
  510 CONTINUE
C
      IF(DT.GE.SDT)THEN
        SDT1=0.0D0
        XQ=XP+SDT*(XQ-XP)/DT
        YQ=YP+SDT*(YQ-YP)/DT
        ZQ=ZP+SDT*(ZQ-ZP)/DT
      ELSE
        SDT1=SDT-DT
        DX=DABS(VXQ*SDT1)
        DY=DABS(VYQ*SDT1)
        DZ=DABS(VZQ*SDT1)
        IF(DX.LE.EPSX .AND. DY.LE.EPSX .AND. DZ.LE.EPSX) SDT1=0.0D0
C
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I         (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M          VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
      ENDIF
C
      RETURN
      END
 
 
 
      SUBROUTINE LOCQ4N
     I                (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I                 N1,N2,N3,N4,SDT,IBF,IJUDGE,NODE,EPSX,
     O                 XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
C -----  2/19/96
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C NOTE: IF IBF=1, BACKWARD TRACKING.                                   %
C       IF IBF=2, FORWARD TRACKING.                                    %
C       IF IJUDGE=1, USING THE AVERAGE-VELOCITY APPROACH.              %
C       IF IJUDGE=2, USING THE SINGLE-VELOCITY APPROACH.               %
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8)
C
      DATA NITER/150/, EPSR/1.0D-6/, OMEGA/1.0D0/
C
      SDT1=SDT
      MXNOD=8
      NSTART=1
      NEND=NODE
      NDIM=3
      NJUMP=1
C
      VXQ=0.0D0
      VYQ=0.0D0
      VZQ=0.0D0
C
C NOTE: THE COMPUTATION IN THIS SUBROUTINE IS ACHIEVED BASED ON
C       A FORWARD PARTICLE TRACKING.  THUS, VELOCITY SIGN NEEDS
C       TO BE CHANGED IF IBF=1.
C
      IF(IBF.EQ.1)THEN
        CALL CHNGSN
     I       (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M        VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
      ENDIF
C
C ***** WHEN THE SINGLE VELOCITY APPROACH IS USED, CALL SUR2D3
C
      IF(IJUDGE.EQ.2)THEN
        CALL SURE3D
     I     (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,N3,N4,XP,YP,ZP,VXP,VYP,VZP,SDT,
     O      XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I        (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M         VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
        RETURN
      ENDIF
C
C ***** CALCULATE THE COEFFICIENTS OF XQ-XP
C
      A1=XX(N1)+XX(N2)+XX(N3)+XX(N4)-4.0D0*XP
      A2=XX(N2)+XX(N3)-XX(N1)-XX(N4)
      A3=XX(N3)+XX(N4)-XX(N1)-XX(N2)
      A4=XX(N1)+XX(N3)-XX(N2)-XX(N4)
C
C ***** CALCULATE THE COEFFICIENTS OF YQ-YP
C
      B1=YY(N1)+YY(N2)+YY(N3)+YY(N4)-4.0D0*YP
      B2=YY(N2)+YY(N3)-YY(N1)-YY(N4)
      B3=YY(N3)+YY(N4)-YY(N1)-YY(N2)
      B4=YY(N1)+YY(N3)-YY(N2)-YY(N4)
C
C ***** CALCULATE THE COEFFICIENTS OF ZQ-ZP
C
      C1=ZZ(N1)+ZZ(N2)+ZZ(N3)+ZZ(N4)-4.0D0*ZP
      C2=ZZ(N2)+ZZ(N3)-ZZ(N1)-ZZ(N4)
      C3=ZZ(N3)+ZZ(N4)-ZZ(N1)-ZZ(N2)
      C4=ZZ(N1)+ZZ(N3)-ZZ(N2)-ZZ(N4)
C
      IF(IJUDGE.EQ.1)THEN
C
C ***** CALCULATE THE COEFFICIENTS OF VXP+VXQ
C
        D1=VXX(N1)+VXX(N2)+VXX(N3)+VXX(N4)+4.0D0*VXP
        D2=VXX(N2)+VXX(N3)-VXX(N1)-VXX(N4)
        D3=VXX(N3)+VXX(N4)-VXX(N1)-VXX(N2)
        D4=VXX(N1)+VXX(N3)-VXX(N2)-VXX(N4)
C
C ***** CALCULATE THE COEFFICIENTS OF VXP+VXQ
C
        E1=VYY(N1)+VYY(N2)+VYY(N3)+VYY(N4)+4.0D0*VYP
        E2=VYY(N2)+VYY(N3)-VYY(N1)-VYY(N4)
        E3=VYY(N3)+VYY(N4)-VYY(N1)-VYY(N2)
        E4=VYY(N1)+VYY(N3)-VYY(N2)-VYY(N4)
C
C ***** CALCULATE THE COEFFICIENTS OF VXP+VXQ
C
        F1=VZZ(N1)+VZZ(N2)+VZZ(N3)+VZZ(N4)+4.0D0*VZP
        F2=VZZ(N2)+VZZ(N3)-VZZ(N1)-VZZ(N4)
        F3=VZZ(N3)+VZZ(N4)-VZZ(N1)-VZZ(N2)
        F4=VZZ(N1)+VZZ(N3)-VZZ(N2)-VZZ(N4)
      ELSE
C
C ***** CALCULATE THE COEFFICIENTS OF VXP
C
        D1=8.0D0*VXP
        D2=0.0D0
        D3=0.0D0
        D4=0.0D0
C
C ***** CALCULATE THE COEFFICIENTS OF VYP
C
        E1=8.0D0*VYP
        E2=0.0D0
        E3=0.0D0
        E4=0.0D0
C
C ***** CALCULATE THE COEFFICIENTS OF VZP
C
        F1=8.0D0*VZP
        F2=0.0D0
        F3=0.0D0
        F4=0.0D0
      ENDIF
C
C ***** CONSTRUCT FUNCTION FOF AND FUNCTION GOF FOR ITERATION
C
C ----- CALCULATE THE COEFFICIENTS OF DX
C
C     R1=A1**2+B1**2+C1**2
      R2=2.0D0*(A1*A2+B1*B2+C1*C2)
      R3=2.0D0*(A1*A3+B1*B3+C1*C3)
      R4=A2**2+B2**2+C2**2
      R5=2.0D0*(A1*A4+B1*B4+C1*C4+A2*A3+B2*B3+C2*C3)
      R6=A3**2+B3**2+C3**2
      R7=2.0D0*(A2*A4+B2*B4+C2*C4)
      R8=2.0D0*(A3*A4+B3*B4+C3*C4)
      R9=A4**2+B4**2+C4**2
C
C ----- CALCULATE THE COEFFICIENTS OF V
C
C     S1=D1**2+E1**2+F1**2
      S2=2.0D0*(D1*D2+E1*E2+F1*F2)
      S3=2.0D0*(D1*D3+E1*E3+F1*F3)
      S4=D2**2+E2**2+F2**2
      S5=2.0D0*(D1*D4+E1*E4+F1*F4+D2*D3+E2*E3+F2*F3)
      S6=D3**2+E3**2+F3**2
      S7=2.0D0*(D2*D4+E2*E4+F2*F4)
      S8=2.0D0*(D3*D4+E3*E4+F3*F4)
      S9=D4**2+E4**2+F4**2
C
C ***** SOLVE THE ABOVE TWO EQUATIONS BY NEWTON-RAPHSON
C
C ***** START NEWTON-RAPHSON
C
      ITEST=1
  100 CONTINUE
      IF(ITEST.EQ.1)THEN
        CALL REPLAS(A1,A2,A3,A4,D1,D2,U1,U2,U3,U4,U5,U6)
        CALL REPLAS(D3,D4,B1,B2,B3,B4,U7,U8,V1,V2,V3,V4)
        CALL REPLAS(E1,E2,E3,E4,E1,E2,V5,V6,V7,V8,V5,V6)
      ELSEIF(ITEST.EQ.2)THEN
        CALL REPLAS(B1,B2,B3,B4,E1,E2,U1,U2,U3,U4,U5,U6)
        CALL REPLAS(E3,E4,C1,C2,C3,C4,U7,U8,V1,V2,V3,V4)
        CALL REPLAS(F1,F2,F3,F4,F1,F2,V5,V6,V7,V8,V5,V6)
      ELSE
        CALL REPLAS(C1,C2,C3,C4,F1,F2,U1,U2,U3,U4,U5,U6)
        CALL REPLAS(F3,F4,A1,A2,A3,A4,U7,U8,V1,V2,V3,V4)
        CALL REPLAS(D1,D2,D3,D4,D1,D2,V5,V6,V7,V8,V5,V6)
      ENDIF
C
      INIGES=0
  105 CONTINUE
      IF(INIGES.EQ.0)THEN
        XSIO=0.0
        ETAO=0.0
      ELSEIF(INIGES.EQ.1)THEN
        XSIO=-1.0D0
        ETAO=-1.0D0
      ELSEIF(INIGES.EQ.2)THEN
        XSIO= 1.0D0
        ETAO=-1.0D0
      ELSEIF(INIGES.EQ.3)THEN
        XSIO= 1.0D0
        ETAO= 1.0D0
      ELSEIF(INIGES.EQ.4)THEN
        XSIO=-1.0D0
        ETAO= 1.0D0
      ENDIF
C
      XSIW=XSIO
      ETAW=ETAO
      XSI=XSIO
      ETA=ETAO
      XID=1.0D0
      EID=1.0D0
C
      DO 200 ITER=1,NITER
C
        DL1=0.25D0*(1.0D0-XSIW)*(1.0D0-ETAW)
        DL2=0.25D0*(1.0D0+XSIW)*(1.0D0-ETAW)
        DL3=0.25D0*(1.0D0+XSIW)*(1.0D0+ETAW)
        DL4=0.25D0*(1.0D0-XSIW)*(1.0D0+ETAW)
C
        XQ=DL1*XX(N1)+DL2*XX(N2)+DL3*XX(N3)+DL4*XX(N4)
        YQ=DL1*YY(N1)+DL2*YY(N2)+DL3*YY(N3)+DL4*YY(N4)
        ZQ=DL1*ZZ(N1)+DL2*ZZ(N2)+DL3*ZZ(N3)+DL4*ZZ(N4)
        VXQ=DL1*VXX(N1)+DL2*VXX(N2)+DL3*VXX(N3)+DL4*VXX(N4)
        VYQ=DL1*VYY(N1)+DL2*VYY(N2)+DL3*VYY(N3)+DL4*VYY(N4)
        VZQ=DL1*VZZ(N1)+DL2*VZZ(N2)+DL3*VZZ(N3)+DL4*VZZ(N4)
C
        RR=(XQ-XP)*(XQ-XP)+(YQ-YP)*(YQ-YP)+(ZQ-ZP)*(ZQ-ZP)
        IF(IJUDGE.EQ.1)THEN
          SS=(VXQ+VXP)*(VXQ+VXP)+(VYQ+VYP)*(VYQ+VYP)+(VZQ+VZP)*(VZQ+VZP)
          IF(SS.EQ.0.0D0)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
            IF(IBF.EQ.1)THEN
              CALL CHNGSN
     I             (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M              VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
            ENDIF
            SDT1=SDT
            RETURN
          ENDIF
        ELSE
          SS=(VXP+VXP)*(VXP+VXP)+(VYP+VYP)*(VYP+VYP)+(VZP+VZP)*(VZP+VZP)
        ENDIF
        RR=16.0D0*RR
        SS=16.0D0*SS
        RR=DSQRT(RR)
        SS=DSQRT(SS)
        T1=U1+U2*XSIW+U3*ETAW+U4*XSIW*ETAW
        T2=U5+U6*XSIW+U7*ETAW+U8*XSIW*ETAW
        T3=V1+V2*XSIW+V3*ETAW+V4*XSIW*ETAW
        T4=V5+V6*XSIW+V7*ETAW+V8*XSIW*ETAW
C
        FOF=SS*T1-RR*T2
        GOF=SS*T3-RR*T4
        IF(FOF.EQ.0.0D0 .AND. GOF.EQ.0.0D0)GOTO 500
C
        DFX=SS*(U2+U4*ETAW)+0.5D0*T1*(S2+2.0D0*S4*XSIW+S5*ETAW+
     >      2.0D0*S7*XSIW*ETAW+S8*ETAW**2+2.0D0*S9*XSIW*ETAW**2)/SS
     >     -RR*(U6+U8*ETAW)-0.5D0*T2*(R2+2.0D0*R4*XSIW+R5*ETAW+
     >      2.0D0*R7*XSIW*ETAW+R8*ETAW**2+2.0D0*R9*XSIW*ETAW**2)/RR
        DFE=SS*(U3+U4*XSIW)+0.5D0*T1*(S3+S5*XSIW+2.0D0*S6*ETAW+
     >      S7*XSIW**2+2.0D0*S8*XSIW*ETAW+2.0D0*S9*XSIW**2*ETAW)/SS
     >     -RR*(U7+U8*XSIW)-0.5D0*T2*(R3+R5*XSIW+2.0D0*R6*ETAW+
     >      R7*XSIW**2+2.0D0*R8*XSIW*ETAW+2.0D0*R9*XSIW**2*ETAW)/RR
        DGX=SS*(V2+V4*ETAW)+0.5D0*T3*(S2+2.0D0*S4*XSIW+S5*ETAW+
     >      2.0D0*S7*XSIW*ETAW+S8*ETAW**2+2.0D0*S9*XSIW*ETAW**2)/SS
     >     -RR*(V6+V8*ETAW)-0.5D0*T4*(R2+2.0D0*R4*XSIW+R5*ETAW+
     >      2.0D0*R7*XSIW*ETAW+R8*ETAW**2+2.0D0*R9*XSIW*ETAW**2)/RR
        DGE=SS*(V3+V4*XSIW)+0.5D0*T3*(S3+S5*XSIW+2.0D0*S6*ETAW+
     >      S7*XSIW**2+2.0D0*S8*XSIW*ETAW+2.0D0*S9*XSIW**2*ETAW)/SS
     >     -RR*(V7+V8*XSIW)-0.5D0*T4*(R3+R5*XSIW+2.0D0*R6*ETAW+
     >      R7*XSIW**2+2.0D0*R8*XSIW*ETAW+2.0D0*R9*XSIW**2*ETAW)/RR
        DETJ=DFX*DGE-DFE*DGX
C
C ***** COMPUTE FOR NEW XSI AND ETA
C
        CALL NEWXE(DETJ,XSIO,ETAO,DFX,DFE,DGX,DGE,XID,EID,FOF,GOF,
     O             0.1D0,  XSI,ETA)
C
C ***** TEST CONVERGENCE
C
        DIFMAX=0.0D0
        DIFMAX1=0.0D0
        IF(XSIW.NE.0.0D0)THEN
          DIF=DABS((XSI-XSIW)/XSIW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(XSIW.EQ.0.0D0)THEN
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(ETAW.NE.0.0D0)THEN
          DIF=DABS((ETA-ETAW)/ETAW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(ETA-ETAW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(ETAW.EQ.0.0D0)THEN
          DIF1=DABS(ETA-ETAW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DIFMAX.LE.EPSR .AND. ITER.GT.1)THEN
          IF(DIFMAX1.LE.EPSR .AND. DABS(FOF)+DABS(GOF).LE.EPSX)GOTO 500
        ENDIF
C
C ***** TAKE A NEW GUESS FOR XSI AND ETA
C
        IF(XSI.GT.1.0D0)THEN
          XSI=1.0D0
          XID=-1.0D0
        ENDIF
        IF(XSI.LT.-1.0D0)THEN
          XSI=-1.0D0
          XID=1.0D0
        ENDIF
        IF(ETA.GT.1.0D0)THEN
          ETA=1.0D0
          EID=-1.0D0
        ENDIF
        IF(ETA.LT.-1.0D0)THEN
          ETA=-1.0D0
          EID=1.0D0
        ENDIF
C
        XSIW=OMEGA*XSI+(1.0D0-OMEGA)*XSIO
        ETAW=OMEGA*ETA+(1.0D0-OMEGA)*ETAO
C
C ***** UPDATE XSIO,ETAO
C
        XSIO=XSI
        ETAO=ETA
  200 CONTINUE
C
      IF(ITEST.NE.3)THEN
        ITEST=ITEST+1
        GOTO 100
      ENDIF
C
      IF(DIFMAX1.GT.EPSR .AND. DABS(FOF)+DABS(GOF).GT.EPSX)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I         (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M          VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
        SDT1=SDT
        RETURN
      ENDIF
C
C ***** GET THE POSITION AND VELOCITY OF POINT Q
C
  500 CONTINUE
      IF(XSI.GT.1.0D0)XSI=1.0D0
      IF(XSI.LT.-1.0D0)XSI=-1.0D0
      IF(ETA.GT.1.0D0)ETA=1.0D0
      IF(ETA.LT.-1.0D0)ETA=-1.0D0
C
      DL1=0.25D0*(1.0D0-XSI)*(1.0D0-ETA)
      DL2=0.25D0*(1.0D0+XSI)*(1.0D0-ETA)
      DL3=0.25D0*(1.0D0+XSI)*(1.0D0+ETA)
      DL4=0.25D0*(1.0D0-XSI)*(1.0D0+ETA)
C
      XQ=DL1*XX(N1)+DL2*XX(N2)+DL3*XX(N3)+DL4*XX(N4)
      YQ=DL1*YY(N1)+DL2*YY(N2)+DL3*YY(N3)+DL4*YY(N4)
      ZQ=DL1*ZZ(N1)+DL2*ZZ(N2)+DL3*ZZ(N3)+DL4*ZZ(N4)
      VXQ=DL1*VXX(N1)+DL2*VXX(N2)+DL3*VXX(N3)+DL4*VXX(N4)
      VYQ=DL1*VYY(N1)+DL2*VYY(N2)+DL3*VYY(N3)+DL4*VYY(N4)
      VZQ=DL1*VZZ(N1)+DL2*VZZ(N2)+DL3*VZZ(N3)+DL4*VZZ(N4)
C
C ***** CALCULATE DT
C
      RR=(XQ-XP)*(XQ-XP)+(YQ-YP)*(YQ-YP)+(ZQ-ZP)*(ZQ-ZP)
      IF(IJUDGE.EQ.1)THEN
        SS=(VXQ+VXP)*(VXQ+VXP)+(VYQ+VYP)*(VYQ+VYP)+(VZQ+VZP)*(VZQ+VZP)
        IF(SS.EQ.0.0D0)THEN
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
          IF(IBF.EQ.1)THEN
            CALL CHNGSN
     I           (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M            VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
          ENDIF
          SDT1=SDT
          RETURN
        ENDIF
      ELSE
        SS=(VXP+VXP)*(VXP+VXP)+(VYP+VYP)*(VYP+VYP)+(VZP+VZP)*(VZP+VZP)
      ENDIF
      RR=DSQRT(RR)
      SS=DSQRT(SS)
      DT=2.0D0*RR/SS
C
C -----CHECK IF A CORRECT ANSWER HAS BEEN OBTAINED
C
      IF(RR.LE.1.0D-4)GOTO 510
      IF(IJUDGE.EQ.1)THEN
        VVX=0.5D0*(VXP+VXQ)
        VVY=0.5D0*(VYP+VYQ)
        VVZ=0.5D0*(VZP+VZQ)
      ELSE
        VVX=VXP
        VVY=VYP
        VVZ=VZP
      ENDIF
      DIFFX=(XQ-XP)-VVX*DT
      DIFFY=(YQ-YP)-VVY*DT
      DIFFZ=(ZQ-ZP)-VVZ*DT
      DIFF=DSQRT(DIFFX*DIFFX+DIFFY*DIFFY+DIFFZ*DIFFZ)
      IF(DIFF.GT.2.0D0*RR*1.0D-4)THEN
        INIGES=INIGES+1
        IF(INIGES.EQ.5)THEN
          IF(ITEST.NE.4)THEN
            ITEST=ITEST+1
            GOTO 100
          ELSE
            IF(IJUDGE.EQ.2)THEN
              CALL SURE3D
     I       (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,N3,N4,XP,YP,ZP,VXP,VYP,VZP,SDT,
     O        XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
            ELSE
              SDT1=SDT
            ENDIF
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
            IF(IBF.EQ.1)THEN
              CALL CHNGSN
     I             (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M              VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
            ENDIF
            RETURN
          ENDIF
        ELSE
          GOTO 105
        ENDIF
      ENDIF
C
  510 CONTINUE
C
      IF(DT.GE.SDT)THEN
        SDT1=0.0D0
        XQ=XP+SDT*(XQ-XP)/DT
        YQ=YP+SDT*(YQ-YP)/DT
        ZQ=ZP+SDT*(ZQ-ZP)/DT
      ELSE
        SDT1=SDT-DT
        DX=DABS(VXQ*SDT1)
        DY=DABS(VYQ*SDT1)
        DZ=DABS(VZQ*SDT1)
        IF(DX.LE.EPSX .AND. DY.LE.EPSX .AND. DZ.LE.EPSX) SDT1=0.0D0
C
C === RESET VELOCITY BY CHANGING SIGN IF IBF=1
        IF(IBF.EQ.1)THEN
          CALL CHNGSN
     I         (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M          VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
        ENDIF
      ENDIF
C
      RETURN
      END
c
c
c
      SUBROUTINE SURE3D
     I     (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2,N3,N4, XP,YP,ZP,VXP,VYP,VZP, SDT,
     O      XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 I,J
C
      DIMENSION XX(8),YY(8),ZZ(8)
      DIMENSION VXX(8),VYY(8),VZZ(8),DL(4)
C
C ===== CONSTRUCT AX+BY+CZ+D=0
C
      A=(YY(N2)-YY(N1))*(ZZ(N3)-ZZ(N1))-(YY(N3)-YY(N1))*(ZZ(N2)-ZZ(N1))
      B=(ZZ(N2)-ZZ(N1))*(XX(N3)-XX(N1))-(ZZ(N3)-ZZ(N1))*(XX(N2)-XX(N1))
      C=(XX(N2)-XX(N1))*(YY(N3)-YY(N1))-(XX(N3)-XX(N1))*(YY(N2)-YY(N1))
      D=-(A*XX(N1)+B*YY(N1)+C*ZZ(N1))
C
C ===== CONSTRUCT X+ET+F=0
C
      E=-VXP
      F=-XP
C
C ===== CONSTRUCT Y+GT+H=0
C
      G=-VYP
      H=-YP
C
C ===== CONSTRUCT Z+IT+J=0
C
      I=-VZP
      J=-ZP
C
C ===== SOLVE THE ABOVE FOUR EQUATIONS
C
      DENO=A*E+B*G+C*I
      IF(DENO.EQ.0.0D0)THEN
        WRITE(16,*)'DENO=0 IN 3DSURE'
        WRITE(16,*)'SOMETHING WRONG WITH THE DATA:'
        WRITE(16,*)'XP,YP,ZP,VXP,VYP,VZP='
        WRITE(16,100)XP,YP,ZP,VXP,VYP,VZP
        WRITE(16,*)'XX(N1),YY(N1),ZZ(N1)='
        WRITE(16,200)XX(N1),YY(N1),ZZ(N1)
        WRITE(16,*)'XX(N2),YY(N2),ZZ(N2)='
        WRITE(16,200)XX(N2),YY(N2),ZZ(N2)
        WRITE(16,*)'XX(N3),YY(N3),ZZ(N3)='
        WRITE(16,200)XX(N3),YY(N3),ZZ(N3)
  100   FORMAT(6F12.4)
  200   FORMAT(3F12.4)
        WRITE(16,*)'STOP SIMULATION'
        STOP
      ENDIF
C
      T=(D-A*F-B*H-C*J)/DENO
C
      IF(T.GE.SDT)THEN
        SDT1=0.0D0
        XQ=XP+VXP*SDT
        YQ=YP+VYP*SDT
        ZQ=ZP+VZP*SDT
        VXQ=VXP
        VYQ=VYP
        VZQ=VZP
        RETURN
      ELSE
        SDT1=SDT-T
        XQ=-F-T*E
        YQ=-H-T*G
        ZQ=-J-T*I
      ENDIF
C
C ===== COMPUTE VXQ AND VYQ BY INTERPOLATION
C
      NODE=4
      IF(N4.EQ.0)NODE=3
      CALL LOCPLN
     I    (XQ,YQ,ZQ,XX,YY,ZZ,N1,N2,N3,N4,NODE,
     O     XSI,ETA,DL1,DL2,DL3,DL)
      VXQ=0.0D0
      VYQ=0.0D0
      VZQ=0.0D0
      DO K=1,NODE
        IF(K.EQ.1)II=N1
        IF(K.EQ.2)II=N2
        IF(K.EQ.3)II=N3
        IF(K.EQ.4)II=N4
        VXQ=VXQ+VXX(II)*DL(K)
        VYQ=VYQ+VYY(II)*DL(K)
        VZQ=VZQ+VZZ(II)*DL(K)
      ENDDO
      RETURN
      END
C
c
c
      SUBROUTINE SUR2D3

     I     (XX,YY,ZZ,VXX,VYY,VZZ,N1,N2, XP,YP,ZP,VXP,VYP,VZP, SDT,EPSX,
     O      XQ,YQ,ZQ,VXQ,VYQ,VZQ, SDT1)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 I
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8)
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        XQ=XP
        YQ=YP
        ZQ=ZP
        VXQ=VXP
        VYQ=VYP
        VZQ=VZP
        SDT1=0.0D0
        RETURN
      ENDIF
C
C ===== COMPUTE THE DIRECTION OF THE TARGET LINE
C
      AA=XX(N2)-XX(N1)
      BB=YY(N2)-YY(N1)
      CC=ZZ(N2)-ZZ(N1)
      DD=XX(N1)-XP
      EE=YY(N1)-YP
      FF=ZZ(N1)-ZP
      GG=VXP*VXP+VYP*VYP+VZP*VZP
C
C ===== CONSTRUCT AX**2+2*BX+C=0
C
      A=AA*AA+BB*BB+CC*CC
      B=AA*DD+BB*EE+CC*FF
      C=DD*DD+EE*EE+FF*FF
      IF(DABS(VXP).GT. DABS(VYP) .AND. DABS(VXP).GT. DABS(VZP))THEN
        A=A-GG*AA*AA/VXP/VXP
        B=B-GG*AA*DD/VXP/VXP
        C=C-GG*DD*DD/VXP/VXP
        JMAX=1
      ELSEIF(DABS(VYP) .GT. DABS(VXP) .AND. DABS(VYP).GT. DABS(VZP))THEN
        A=A-GG*BB*BB/VYP/VYP
        B=B-GG*BB*EE/VYP/VYP
        C=C-GG*EE*EE/VYP/VYP
        JMAX=2
      ELSE
        A=A-GG*CC*CC/VZP/VZP
        B=B-GG*CC*FF/VZP/VZP
        C=C-GG*FF*FF/VZP/VZP
        JMAX=3
      ENDIF
      IF(A.EQ.0.0D0)THEN
        XSI=-C/2.0/B
      ELSE
        D=B*B-A*C
        IF(D.LT.0.0D0) THEN
          D=0.0D0
        ELSE
          D=DSQRT(D)
        ENDIF
        XSI=(-B+D)/A
      ENDIF
      IF(XSI.LT.0.0D0) XSI=0.0D0
      IF(XSI.GT.1.0D0) XSI=1.0D0
      IF(JMAX.EQ.1)THEN
        T=(AA*XSI+DD)/VXP
      ELSEIF(JMAX.EQ.2)THEN
        T=(BB*XSI+EE)/VYP
      ELSE
        T=(CC*XSI+FF)/VZP
      ENDIF
      RES1=DD+AA*XSI-VXP*T
      RES2=EE+BB*XSI-VYP*T
      RES3=FF+CC*XSI-VZP*T
      RES=RES1+RES2+RES3
      IF(T.LT.0.0D0 .OR. DABS(RES).GT.EPSX)THEN
         XSI=(-B-D)/A
         IF(JMAX.EQ.1)THEN
           T=(AA*XSI+DD)/VXP
         ELSEIF(JMAX.EQ.2)THEN
           T=(BB*XSI+EE)/VYP
         ELSE
           T=(CC*XSI+FF)/VZP
         ENDIF
      ENDIF
C
  500 CONTINUE
      IF(T.GE.SDT)THEN
        SDT1=0.0D0
        XQ=XP+VXP*SDT
        YQ=YP+VYP*SDT
        ZQ=ZP+VZP*SDT
        VXQ=VXP
        VYQ=VYP
        VZQ=VZP
        RETURN
      ELSE
        SDT1=SDT-T
        XQ=XP+VXP*T
        YQ=YP+VYP*T
        ZQ=ZP+VZP*T
      ENDIF
C
C ===== COMPUTE VXQ AND VYQ BY INTERPOLATION
C
      XSI=DSQRT((XQ-XX(N1))**2+(YQ-YY(N1))**2+(ZQ-ZZ(N1))**2)/
     >    DSQRT((XX(N2)-XX(N1))**2+(YY(N2)-YY(N1))**2+
     >          (ZZ(N2)-ZZ(N1))**2)      
      VXQ=VXX(N1)*(1.0D0-XSI)+VXX(N2)*XSI
      VYQ=VYY(N1)*(1.0D0-XSI)+VYY(N2)*XSI
      VZQ=VZZ(N1)*(1.0D0-XSI)+VZZ(N2)*XSI
      RETURN
      END
C
c
c
      SUBROUTINE REPLAS(XW,YW,ZW,VXW,VYW,VZW,XP,YP,ZP,VPX,VPY,VPZ)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      XP=XW
      YP=YW
      ZP=ZW
      VPX=VXW
      VPY=VYW
      VPZ=VZW
      RETURN
      END
c
c
c
      SUBROUTINE VALBDL(VXX,VYY,VZZ,DL,MXNODE,NODE,VPX,VPY,VPZ)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION VXX(MXNODE),VYY(MXNODE),VZZ(MXNODE),DL(MXNODE)
C
        VPX=0.0
        VPY=0.0
        VPZ=0.0
        DO I1=1,NODE
          VPX=VPX+DL(I1)*VXX(I1)
          VPY=VPY+DL(I1)*VYY(I1)
          VPZ=VPZ+DL(I1)*VZZ(I1)
        ENDDO
      RETURN
      END
C
c
c
      SUBROUTINE ONPLAN
     I    (XW,YW,ZW,N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C       10/14/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XW(MXNPW),YW(MXNPW),ZW(MXNPW)
C
      A=(YW(N2)-YW(N1))*(ZW(N3)-ZW(N1))-(YW(N3)-YW(N1))*(ZW(N2)-ZW(N1))
      B=(ZW(N2)-ZW(N1))*(XW(N3)-XW(N1))-(ZW(N3)-ZW(N1))*(XW(N2)-XW(N1))
      C=(XW(N2)-XW(N1))*(YW(N3)-YW(N1))-(XW(N3)-XW(N1))*(YW(N2)-YW(N1))
      D=A*XW(N1)+B*YW(N1)+C*ZW(N1)
      DD=A*XP+B*YP+C*ZP
      DD=DD-D
      RMAX=DMAX1(DABS(A),DABS(B))
      RMAX=DMAX1(RMAX,DABS(C))
      IF(RMAX.EQ.DABS(A))THEN
        XP=XP-DD/A
      ELSEIF(RMAX.EQ.DABS(B))THEN
        YP=YP-DD/B
      ELSE
        ZP=ZP-DD/C
      ENDIF
C
      RETURN
      END
c
c
c
      SUBROUTINE ONLINE
     I      (XW,YW,ZW,J1,J2,MXNPW,
     M       XP,YP,ZP)
C       10/14/93
       IMPLICIT REAL*8(A-H,O-Z)
C
       DIMENSION XW(MXNPW),YW(MXNPW),ZW(MXNPW)
C
       D1=DSQRT((XW(J2)-XW(J1))**2+(YW(J2)-YW(J1))**2+
     >          (ZW(J2)-ZW(J1))**2)
       D2=DSQRT((XP-XW(J1))**2+(YP-YW(J1))**2+
     >          (ZP-ZW(J1))**2)
       XSI=D2/D1
       IF(XSI.GT.1.0D0)XSI=1.0D0
       IF(XSI.LT.0.0D0)XSI=0.0D0
       XP=XW(J1)*(1.0D0-XSI)+XW(J2)*XSI
       YP=YW(J1)*(1.0D0-XSI)+YW(J2)*XSI
       ZP=ZW(J1)*(1.0D0-XSI)+ZW(J2)*XSI
C
       RETURN
       END
C
C
C
      SUBROUTINE PLANEW
     I                (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,I,J,
     I                 IBF,CIA,IJUDGE, idi,
     O                 DOT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C ----- 12/30/94
c note: if idi=1, then epsr is used to determine "DOT"
c
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),XW(8),YW(8),ZW(8)
c      DIMENSION DX(3),DY(3),DZ(3)
C
      LOGICAL DOT
C
      data epsr/3.046d-5/
C
      DOT=.FALSE.
C
      XW(1)=XP
      YW(1)=YP
      ZW(1)=ZP
      XW(2)=XX(I)
      YW(2)=YY(I)
      ZW(2)=ZZ(I)
      XW(3)=XX(J)
      YW(3)=YY(J)
      ZW(3)=ZZ(J)
C
      DD=0.0D0
      DO I2=1,3
        I1=I2+1
        IF(I1.GT.3)I1=1
        DD=DD+DABS(XW(I1)-XW(I2))+DABS(YW(I1)-YW(I2))+
     >        DABS(ZW(I1)-ZW(I2))
      ENDDO
C
C     WRITE(16,*)'DOTP=',DOTP
      IF(IJUDGE.EQ.2)THEN
        VV=DSQRT(VXP*VXP+VYP*VYP+VZP*VZP)
        IF(VV.EQ.0.0D0)VV=1.0D0
C       XQ=XP+VXP*1.0D20/VV
C       YQ=YP+VYP*1.0D20/VV
C       ZQ=ZP+VZP*1.0D20/VV
        XQ=XP+VXP*DD/VV
        YQ=YP+VYP*DD/VV
        ZQ=ZP+VZP*DD/VV
        CANG=CIA*FCOS(XW,YW,ZW,XQ,YQ,ZQ,idi)
c
        if(idi.eq.1)then
          if(dabs(cang).le.epsr)cang=0.0d0
        endif
c
        IF(IBF.EQ.1 .AND. CANG.LE.0.0D0)DOT=.TRUE.
        IF(IBF.EQ.2 .AND. CANG.GE.0.0D0)DOT=.TRUE.
      ENDIF
C
      IF(IJUDGE.EQ.1)THEN
        VP1X=0.5*(VXX(I)+VXP)
        VP1Y=0.5*(VYY(I)+VYP)
        VP1Z=0.5*(VZZ(I)+VZP)
        VP2X=0.5*(VXX(J)+VXP)
        VP2Y=0.5*(VYY(J)+VYP)
        VP2Z=0.5*(VZZ(J)+VZP)
        VV1=DSQRT(VP1X*VP1X+VP1Y*VP1Y+VP1Z*VP1Z)
        VV2=DSQRT(VP2X*VP2X+VP2Y*VP2Y+VP2Z*VP2Z)
c
        if(vv1.eq.0.0d0 .and. vv2.eq.0.0d0)return
c
        IF(VV1.EQ.0.0D0)VV1=1.0D0
        IF(VV2.EQ.0.0D0)VV2=1.0D0
C       XQ1=XP+VP1X*1.0D20/VV1
C       YQ1=YP+VP1Y*1.0D20/VV1
C       ZQ1=ZP+VP1Z*1.0D20/VV1
C       XQ2=XP+VP2X*1.0D20/VV2
C       YQ2=YP+VP2Y*1.0D20/VV2
C       ZQ2=ZP+VP2Z*1.0D20/VV2
        XQ1=XP+VP1X*DD/VV1
        YQ1=YP+VP1Y*DD/VV1
        ZQ1=ZP+VP1Z*DD/VV1
        XQ2=XP+VP2X*DD/VV2
        YQ2=YP+VP2Y*DD/VV2
        ZQ2=ZP+VP2Z*DD/VV2
        CANG1=CIA*FCOS(XW,YW,ZW,XQ1,YQ1,ZQ1,idi)
        CANG2=CIA*FCOS(XW,YW,ZW,XQ2,YQ2,ZQ2,idi)
c
        if(idi.eq.1)then
          if(dabs(cang1).le.epsr)cang1=0.0d0
          if(dabs(cang2).le.epsr)cang2=0.0d0
        endif
c
        IF(IBF.EQ.1 .AND. CANG1.LE.0.0D0 .AND. CANG2.LE.0.0D0)DOT=.TRUE.
        IF(IBF.EQ.2 .AND. CANG1.GE.0.0D0 .AND. CANG2.GE.0.0D0)DOT=.TRUE.
      ENDIF
C
C     WRITE(16,*)'DOT=',DOT
      RETURN
      END
C
C
C
      SUBROUTINE TRAK1H
     I                 (XP,YP,ZP,VXP,VYP,VZP,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,
     I                  IBF,IE,IJUDGE,J,MAXEL, EPSX, idi,
     O                  XQ,YQ,ZQ,VXQ,VYQ,VZQ,N1,N2,N3,N4,SDT1)
C
C ----- 10/8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),IE(MAXEL,8)
      DIMENSION KGB(4,6),KGB1(6,8),CIA1(3,8),KGB2(8,8)
C
      LOGICAL DOT1,DOT2,DOT(3,2)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8/
      DATA KGB1/5,1,2,3,4,6, 5,2,3,4,1,6, 5,3,4,1,2,6, 5,4,1,2,3,6,
     >          6,2,1,4,3,5, 6,3,2,1,4,5, 6,4,3,2,1,5, 6,1,4,3,2,5/
      DATA CIA1/1.0,1.0,-1.0, 1.0,-1.0,-1.0, 1.0,-1.0,1.0, 1.0,1.0,1.0,
     >       -1.0,-1.0,1.0, -1.0,-1.0,-1.0, -1.0,1.0,-1.0, -1.0,1.0,1.0/
      DATA KGB2/1,2,3,4,5,6,7,8, 2,3,4,1,6,7,8,5, 3,4,1,2,7,8,5,6,
     >          4,1,2,3,8,5,6,7, 5,8,7,6,1,4,3,2, 6,5,8,7,2,1,4,3,
     >          7,6,5,8,3,2,1,4, 8,7,6,5,4,3,2,1/
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        SDT1=0.0D0
        XQ=XP
        YQ=YP
        ZQ=ZP
        RETURN
      ENDIF
C
  10  CONTINUE
C
C ***** DETERMINE POINT Q WILL FALL IN THIS ELEMENT OR NOT
C
      DO 100 K=1,3
        MM=KGB1(K,J)
        CI=CIA1(K,J)
        DO K1=1,4
          IF(J.EQ.KGB(K1,MM))GOTO 20
        ENDDO
  20    CONTINUE
        N1=K1+1
        N2=K1+2
        N3=K1+3
        IF(N1.GT.4)N1=N1-4
        IF(N2.GT.4)N2=N2-4
        IF(N3.GT.4)N3=N3-4
        N1=KGB(N1,MM)
        N2=KGB(N2,MM)
        N3=KGB(N3,MM)
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,N1,N2,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,N2,N3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT2)
        IF(DOT1.OR.DOT2)THEN
          IF(CI.GT.0.0)THEN
            DOT(K,1)=DOT1
            DOT(K,2)=DOT2
          ELSE
            DOT(K,1)=DOT2
            DOT(K,2)=DOT1
          ENDIF
          GOTO 100
        ENDIF
        SDT1=SDT
        RETURN
 100  CONTINUE
C
C ***** DETERMINE WHICH PLANE THIS POINT COULD COME FROM
C
      IF(DOT(1,1).AND.DOT(3,2))THEN
C
C +++++ (1) PLANE KGB1(4,J)
C
        NN1=KGB2(2,J)
        NN2=KGB2(3,J)
        NN3=KGB2(7,J)
        NN4=KGB2(6,J)
        CI=1.0D0
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(DOT1)THEN
          CI=1.0D0
          CALL PLANEW
     I              (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN4,
     I               IBF,CI,IJUDGE, idi,
     O               DOT2)
          IF(DOT2)THEN
            GOTO 200
          ENDIF
        ENDIF
      ENDIF
C
      IF(DOT(1,2).AND.DOT(2,1))THEN
C
C +++++ (2) PLANE KGB1(5,J)
C
        NN1=KGB2(4,J)
        NN2=KGB2(3,J)
        NN3=KGB2(7,J)
        NN4=KGB2(8,J)
        CI=-1.0D0
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(DOT1)THEN
          CI=-1.0D0
          CALL PLANEW
     I              (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN4,
     I               IBF,CI,IJUDGE, idi,
     O               DOT2)
          IF(DOT2)THEN
            GOTO 200
          ENDIF
        ENDIF
      ENDIF
C
      IF(DOT(3,1).AND.DOT(2,2))THEN
C
C +++++ (3) PLANE KGB1(6,J)
C
        NN1=KGB2(5,J)
        NN2=KGB2(6,J)
        NN3=KGB2(7,J)
        NN4=KGB2(8,J)
        CI=1.0D0
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN4,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(DOT1)THEN
          CI=-1.0D0
          CALL PLANEW
     I              (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN2,
     I               IBF,CI,IJUDGE, idi,
     O               DOT2)
          IF(DOT2)THEN
            GOTO 200
          ENDIF
        ENDIF
      ENDIF
C
      IF(IJUDGE.EQ.1 .AND. DOT(1,1) .AND. DOT(1,2) .AND. DOT(2,1) .AND.
     >   DOT(2,2) .AND. DOT(3,1) .AND. DOT(3,2))THEN
        IJUDGE=2
        GOTO 10
      ENDIF
      SDT1=SDT
      GOTO 500
C
C ***** DETERMINE LOCATION OF Q
C
  200 CONTINUE
C
      CALL LOCQ4N
     I           (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I            NN1,NN2,NN3,NN4,SDT,IBF,IJUDGE,8,EPSX,
     O            XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
      if(sdt1.eq.sdt)then
        ijudge=2
        goto 10
      endif
C
  500 CONTINUE
C
      N1=IE(M,NN1)
      N2=IE(M,NN2)
      N3=IE(M,NN3)
      N4=IE(M,NN4)
      RETURN
      END
C
C
C
      SUBROUTINE TRAK1P
     I                 (XP,YP,ZP,VXP,VYP,VZP,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,
     I                  IBF,IE,IJUDGE,J,MAXEL, EPSX, idi,
     O                  XQ,YQ,ZQ,VXQ,VYQ,VZQ,N1,N2,N3,N4,SDT1)
C
C ----- 10/8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),IE(MAXEL,8)
      DIMENSION KGB(4,5),KGB1(5,6),CIA1(3,6),KGB2(6,6)
C
      LOGICAL DOT1,DOT2,DOT(3,2)
C
      DATA KGB/1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0/
      DATA KGB1/4,1,2,3,5,   4,2,3,1,5,   4,3,1,2,5,
     >          5,2,1,3,4,   5,3,2,1,4,   5,1,3,2,4/
      DATA CIA1/1.0,1.0,1.0,   1.0,1.0,1.0,   1.0,1.0,1.0,
     >         -1.0,1.0,1.0,  -1.0,1.0,1.0,  -1.0,1.0,1.0/
      DATA KGB2/1,2,3,4,5,6,  2,3,1,5,6,4,  3,1,2,6,4,5,
     >          4,6,5,1,3,2,  5,4,6,2,1,3,  6,5,4,3,2,1/
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        SDT1=0.0D0
        XQ=XP
        YQ=YP
        ZQ=ZP
        RETURN
      ENDIF
C
  10  CONTINUE
C
C ***** DETERMINE POINT Q WILL FALL IN THIS ELEMENT OR NOT
C
      DO 100 K=1,3
        MM=KGB1(K,J)
        CI=CIA1(K,J)
        DO K1=1,4
          IF(J.EQ.KGB(K1,MM))GOTO 20
        ENDDO
  20    CONTINUE
        N1=K1+1
        N2=K1+2
        IF(MM.EQ.4 .OR. MM.EQ.5)THEN
          IF(N1.GT.3)N1=N1-3
          IF(N2.GT.3)N2=N2-3
          N1=KGB(N1,MM)
          N2=KGB(N2,MM)
        ELSE
          N3=K1+3
          IF(N1.GT.4)N1=N1-4
          IF(N2.GT.4)N2=N2-4
          IF(N3.GT.4)N3=N3-4
          N1=KGB(N1,MM)
          N2=KGB(N2,MM)
          N3=KGB(N3,MM)
        ENDIF
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,N1,N2,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
C
        IF(MM.EQ.4 .OR. MM.EQ.5)THEN
          DOT2=.FALSE.
        ELSE
          CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,N2,N3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT2)
        ENDIF
C
        IF(DOT1.OR.DOT2)THEN
          IF(CI.GT.0.0)THEN
            DOT(K,1)=DOT1
            DOT(K,2)=DOT2
          ELSE
            DOT(K,1)=DOT2
            DOT(K,2)=DOT1
          ENDIF
          GOTO 100
        ENDIF
        SDT1=SDT
        RETURN
 100  CONTINUE
C
C ***** DETERMINE WHICH PLANE THIS POINT COULD COME FROM
C
      IF(DOT(2,2).AND.DOT(3,1))THEN
C
C +++++ (1) PLANE KGB1(5,J)
C
        NN1=KGB2(4,J)
        NN2=KGB2(5,J)
        NN3=KGB2(6,J)
        NN4=KGB2(1,J)
        CI=1.0D0
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(DOT1)THEN
          CALL LOCQ3N
     I                (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I                 NN1,NN2,NN3,SDT,IBF,IJUDGE,6,EPSX,
     O                 XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
          if(sdt1.eq.sdt)then
            ijudge=2
            goto 10
          endif
          N1=IE(M,NN1)
          N2=IE(M,NN2)
          N3=IE(M,NN3)
          N4=0
          RETURN
        ENDIF
      ENDIF
C
      IF(DOT(2,1).AND.DOT(3,2))THEN
C
C +++++ (2) PLANE KGB1(3,J)
C
        NN1=KGB2(2,J)
        NN2=KGB2(5,J)
        NN3=KGB2(6,J)
        NN4=KGB2(3,J)
        CI=-1.0D0
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(DOT1)THEN
c         CI=-1.0D0
c         CALL PLANEW
c    I              (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN4,
c    I               IBF,CI,IJUDGE, idi,
c    O               DOT2)
c         IF(DOT2)THEN
            CALL LOCQ4N
     I                (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I                 NN1,NN2,NN3,NN4,SDT,IBF,IJUDGE,6,EPSX,
     O                 XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
            if(sdt1.eq.sdt)then
              ijudge=2
              goto 10
            endif
            N1=IE(M,NN1)
            N2=IE(M,NN2)
            N3=IE(M,NN3)
            N4=IE(M,NN4)
            RETURN
c         ENDIF
        ENDIF
      ENDIF
C
      IF(IJUDGE.EQ.1 .AND. DOT(1,1) .AND. DOT(2,1) .AND.
     >   DOT(2,2) .AND. DOT(3,1) .AND. DOT(3,2))THEN
        IJUDGE=2
        GOTO 10
      ENDIF
      SDT1=SDT
C
      N1=IE(M,NN1)
      N2=IE(M,NN2)
      N3=IE(M,NN3)
      N4=IE(M,NN4)
C
      RETURN
      END
C
C
C
      SUBROUTINE TRAK1T
     I                 (XP,YP,ZP,VXP,VYP,VZP,XX,YY,ZZ,VXX,VYY,VZZ,SDT,
     I                  IBF,IE,IJUDGE,J,M,MAXEL, EPSX, idi,
     O                  XQ,YQ,ZQ,VXQ,VYQ,VZQ,N1,N2,N3,N4,SDT1)
C
C----- 10/8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),IE(MAXEL,8)
      DIMENSION KGBB(4,4)
C
      LOGICAL DOT1,DOT2,DOT3
C
      DATA KGBB/4,3,2,1, 4,1,3,2, 4,2,1,3, 1,2,3,4/
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        SDT1=0.0D0
        XQ=XP
        YQ=YP
        ZQ=ZP
        RETURN
      ENDIF
   10 CONTINUE
C
C
C ***** DETERMINE POINT Q WILL FALL IN THIS ELEMENT OR NOT
C
      NN1=KGBB(1,J)
      NN2=KGBB(2,J)
      NN3=KGBB(3,J)
      NN4=KGBB(4,J)
      CI=1.0D0
      CALL PLANEW
     I          (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN1,
     I           IBF,CI,IJUDGE, idi,
     O           DOT1)
      IF(DOT1)THEN
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN2,
     I             IBF,CI,IJUDGE, idi,
     O             DOT2)
        IF(DOT2)THEN
          CALL PLANEW
     I              (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN1,NN3,
     I               IBF,CI,IJUDGE, idi,
     O               DOT3)
          IF(DOT3)GOTO 100
        ENDIF
      ENDIF
      SDT1=SDT
      GOTO 200
C
  100 CONTINUE
C
C ***** DETERMINE LOCATION OF Q
C
      CALL LOCQ3N
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I             NN1,NN2,NN3,SDT,IBF,IJUDGE,4,EPSX,
     O             XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
      if(sdt1.eq.sdt)then
        ijudge=2
        goto 10
      endif
C
  200 CONTINUE
C
      N1=IE(M,NN1)
      N2=IE(M,NN2)
      N3=IE(M,NN3)
      N4=IE(M,NN4)
      RETURN
      END
C
C
C
      SUBROUTINE TRAK2H
     I                 (XP,YP,ZP,VXP,VYP,VZP,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,
     I                  IBF,IE,IJUDGE,ID,MAXEL, EPSX,J1,J2, idi,
     M                  N1,N2,N3,N4,
     O                  XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
C ----- 10/8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),IE(MAXEL,8)
      DIMENSION KGB(4,6),CIA(6)
C
      LOGICAL DOT1,DOT2,DOT3,DOT4
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8/
      DATA CIA/1.0, -1.0, -1.0, 1.0, 1.0, -1.0/
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        SDT1=0.0D0
        XQ=XP
        YQ=YP
        ZQ=ZP
        RETURN
      ENDIF
C
C ***** DETERMINE WHICH PLANE N1,N2,N3,N4 ARE ON
C
      KM1=0
      KM2=0
      IF(ID.EQ.1)THEN
C ----- ON THE SIDE PLANE
        DO 10 K=1,6
          DO K1=1,4
            KK=KGB(K1,K)
            KK=IE(M,KK)
            IF(N1.NE.KK .AND. N2.NE.KK .AND. N3.NE.KK .AND. N4.NE.KK)
     >         GOTO 10
          ENDDO
          KM1=K
          KM2=0
          GOTO 100
  10    CONTINUE
      ENDIF
C
      IF(ID.EQ.2)THEN
C ----- ON THE SIDE LINE
        DO 20 K=1,6
          KOUN=0
          DO K1=1,4
            KK=KGB(K1,K)
            KK=IE(M,KK)
            IF(J1.EQ.KK .OR. J2.EQ.KK)
     >         KOUN=KOUN+1
          ENDDO
          IF(KOUN.EQ.2 .AND. KM1.EQ.0)THEN
            KM1=K
            GOTO 20
          ENDIF
          IF(KOUN.EQ.2 .AND. KM1.NE.0)THEN
            KM2=K
            GOTO 100
          ENDIF
  20    CONTINUE
      ENDIF
C
C ***** START CHECKING EVERY PLANE
C
 100  CONTINUE
      DO 200 K=1,6
        IF(K.EQ.KM1 .OR. K.EQ.KM2)GOTO 200
        NN1=KGB(1,K)
        NN2=KGB(2,K)
        NN3=KGB(3,K)
        NN4=KGB(4,K)
        CI=CIA(K)
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN1,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(.NOT.DOT1)GOTO 200
C
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN2,
     I             IBF,CI,IJUDGE, idi,
     O             DOT2)
          IF(.NOT.DOT2)GOTO 200
C
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN4,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT3)
        IF(.NOT.DOT3)GOTO 200
C
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN1,NN4,
     I             IBF,CI,IJUDGE, idi,
     O             DOT4)
        IF(DOT4)GOTO 300
C
 200  CONTINUE
      SDT1=SDT
      RETURN
C
C ***** DETERMING Q
C
 300  CONTINUE
C
      CALL LOCQ4N
     I           (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I            NN1,NN2,NN3,NN4,SDT,IBF,IJUDGE,8,EPSX,
     O            XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
      if(sdt1.eq.sdt)return
C
      N1=IE(M,NN1)
      N2=IE(M,NN2)
      N3=IE(M,NN3)
      N4=IE(M,NN4)
      RETURN
      END
C
C
C
      SUBROUTINE TRAK2P
     I          (XP,YP,ZP,VXP,VYP,VZP,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,
     I           IBF,IE,IJUDGE,ID,MAXEL, EPSX, J1,J2, idi,
     M           N1,N2,N3,N4,XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
C ----- 10/8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),IE(MAXEL,8)
      DIMENSION KGB(4,5),CIA(5)
C
      LOGICAL DOT1,DOT2,DOT3,DOT4
C
      DATA KGB/1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0/
      DATA CIA/1.0, 1.0, 1.0, 1.0, -1.0/
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        SDT1=0.0D0
        XQ=XP
        YQ=YP
        ZQ=ZP
        RETURN
      ENDIF
C
C ***** DETERMINE WHICH PLANE N1,N2,N3,N4 ARE ON
C
      KM1=0
      KM2=0
      IF(ID.EQ.1)THEN
C ----- ON THE SIDE PLANE
        DO 10 K=1,5
          IF(K.EQ.4 .OR. K.EQ.5)THEN
            NODE=3
            N4=0
          ELSE
            NODE=4
          ENDIF
          DO K1=1,NODE
            KK=KGB(K1,K)
            KK=IE(M,KK)
            IF(N1.NE.KK .AND. N2.NE.KK .AND. N3.NE.KK .AND. N4.NE.KK)
     >         GOTO 10
          ENDDO
          KM1=K
          KM2=0
          GOTO 100
  10    CONTINUE
      ENDIF
C
      IF(ID.EQ.2)THEN
C ----- ON THE SIDE LINE
        DO 20 K=1,5
          IF(K.EQ.4 .OR. K.EQ.5)THEN
            NODE=3
          ELSE
            NODE=4
          ENDIF
          KOUN=0
          DO K1=1,NODE
            KK=KGB(K1,K)
            KK=IE(M,KK)
            IF(J1.EQ.KK .OR. J2.EQ.KK)
     >         KOUN=KOUN+1
          ENDDO
          IF(KOUN.EQ.2 .AND. KM1.EQ.0)THEN
            KM1=K
            GOTO 20
          ENDIF
          IF(KOUN.EQ.2 .AND. KM1.NE.0)THEN
            KM2=K
            GOTO 100
          ENDIF
  20    CONTINUE
      ENDIF
C
C ***** START CHECKING EVERY PLANE
C
 100  CONTINUE
      DO 200 K=1,5
        IF(K.EQ.KM1 .OR. K.EQ.KM2)GOTO 200
        NN1=KGB(1,K)
        NN2=KGB(2,K)
        NN3=KGB(3,K)
        NN4=KGB(4,K)
        CI=CIA(K)
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN1,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(.NOT.DOT1)GOTO 200
C
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN2,
     I             IBF,CI,IJUDGE, idi,
     O             DOT2)
        IF(.NOT.DOT2)GOTO 200
C
C ----- CONSIDER DIFFERENCE CASES
C
        IF(K.EQ.4 .OR. K.EQ.5)THEN
          CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN1,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT3)
          IF(.NOT.DOT3)GOTO 200
C
          CALL LOCQ3N
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I             NN1,NN2,NN3,SDT,IBF,IJUDGE,6,EPSX,
     O             XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
          if(sdt1.eq.sdt)return
c
          N1=IE(M,NN1)
          N2=IE(M,NN2)
          N3=IE(M,NN3)
          N4=0
          RETURN
        ELSE
          CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN4,NN3,
     I             IBF,CI,IJUDGE, idi,
     O             DOT3)
          IF(.NOT.DOT3)GOTO 200
          CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN1,NN4,
     I             IBF,CI,IJUDGE, idi,
     O             DOT4)
          IF(.NOT.DOT4)GOTO 200
C
          CALL LOCQ4N
     I               (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I                NN1,NN2,NN3,NN4,SDT,IBF,IJUDGE,6,EPSX,
     O                XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
          if(sdt1.eq.sdt)return
C
          N1=IE(M,NN1)
          N2=IE(M,NN2)
          N3=IE(M,NN3)
          N4=IE(M,NN4)
          RETURN
        ENDIF
C
 200  CONTINUE
      SDT1=SDT
      RETURN
      END
C
C
C
      SUBROUTINE TRAK2T
     I          (XP,YP,ZP,VXP,VYP,VZP,XX,YY,ZZ,VXX,VYY,VZZ,SDT,
     I           IBF,IE,IJUDGE,MAXEL,ID,M, EPSX,J1,J2, idi,
     M           N1,N2,N3,N4,XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
C----- 10/8/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),IE(MAXEL,8)
      DIMENSION KGBB(4,4)
C
      LOGICAL DOT1,DOT2,DOT3
C
      DATA KGBB/4,3,2,1, 4,1,3,2, 4,2,1,3, 1,2,3,4/
C
      IF(VXP.EQ.0.0D0 .AND. VYP.EQ.0.0D0 .AND. VZP.EQ.0.0D0)THEN
        SDT1=0.0D0
        XQ=XP
        YQ=YP
        ZQ=ZP
        RETURN
      ENDIF
C
C ***** DETERMINE WHICH PLANE N1,N2,N3 ARE ON
C
      KM1=0
      KM2=0
      IF(ID.EQ.1)THEN
C ----- ON THE SIDE PLANE
        DO 10 K=1,4
          DO K1=1,3
            KK=KGBB(K1,K)
            KK=IE(M,KK)
            IF(N1.NE.KK .AND. N2.NE.KK .AND. N3.NE.KK)GOTO 10
          ENDDO
          KM1=K
          KM2=0
          GOTO 100
  10    CONTINUE
      ENDIF
C
      IF(ID.EQ.2)THEN
C ----- ON THE SIDE LINE
        DO 20 K=1,4
          KOUN=0
          DO K1=1,3
            KK=KGBB(K1,K)
            KK=IE(M,KK)
            IF(J1.EQ.KK .OR. J2.EQ.KK)KOUN=KOUN+1
          ENDDO
          IF(KOUN.EQ.2 .AND. KM1.EQ.0)THEN
            KM1=K
            GOTO 20
          ENDIF
          IF(KOUN.EQ.2 .AND. KM1.NE.0)THEN
            KM2=K
            GOTO 100
          ENDIF
  20    CONTINUE
      ENDIF
C
C ***** START CHECKING EVERY PLANE
C
 100  CONTINUE
      DO 200 K=1,4
        IF(K.EQ.KM1 .OR. K.EQ.KM2)GOTO 200
        NN1=KGBB(1,K)
        NN2=KGBB(2,K)
        NN3=KGBB(3,K)
        NN4=KGBB(4,K)
        CI=1.0D0
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN2,NN1,
     I             IBF,CI,IJUDGE, idi,
     O             DOT1)
        IF(.NOT.DOT1)GOTO 200
C
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN3,NN2,
     I             IBF,CI,IJUDGE, idi,
     O             DOT2)
          IF(.NOT.DOT2)GOTO 200
C
        CALL PLANEW
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,NN1,NN3,
     I             IBF,CI,IJUDGE,idi,
     O             DOT3)
        IF(.NOT.DOT3)GOTO 200
        GOTO 300
 200  CONTINUE
      SDT1=SDT
      RETURN
C
C ***** DETERMINE LOCATION OF Q
C
 300  CONTINUE
      CALL LOCQ3N
     I            (XX,YY,ZZ,VXX,VYY,VZZ,XP,YP,ZP,VXP,VYP,VZP,
     I             NN1,NN2,NN3,SDT,IBF,IJUDGE,4,EPSX,
     O             XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
      if(sdt1.eq.sdt)return
C
      N1=IE(M,NN1)
      N2=IE(M,NN2)
      N3=IE(M,NN3)
      N4=IE(M,NN4)
      RETURN
      END
c
c
c
      SUBROUTINE XSI2D
c     I                (X,Y,XP,YP,M,
     I                (X,Y,XP,YP,
     O                 XSI,ETA)
C
C $$$$$ TO COMPUTE LOCAL COORDINATES BASED ON THE GIVEN ORIGINAL
C       CARTESIAN COORDINATES.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(4),Y(4)
C
      DATA NITER/100/,OMEGA/1.0D0/,EPSR/1.0D-6/
C
      X2=-1.0D38
      X1=1.0D38
      Y2=-1.0D38
      Y1=1.0D38
C
C ***** START CALCULATIONS
C
      DO 190 IQ=1,4
        X2=DMAX1(X2,X(IQ))
        Y2=DMAX1(Y2,Y(IQ))
        X1=DMIN1(X1,X(IQ))
        Y1=DMIN1(Y1,Y(IQ))
  190 CONTINUE
      DELX=X2-X1
      DELY=Y2-Y1
C
      A1=(X(1)+X(2)+X(3)+X(4))*0.25D0
      A2=(-X(1)+X(2)+X(3)-X(4))*0.25D0
      A3=(-X(1)+X(4)-X(2)+X(3))*0.25D0
      A4=(X(1)-X(4)-X(2)+X(3))*0.25D0
      B1=(Y(1)+Y(2)+Y(3)+Y(4))*0.25D0
      B2=(-Y(1)+Y(2)+Y(3)-Y(4))*0.25D0
      B3=(-Y(1)+Y(4)-Y(2)+Y(3))*0.25D0
      B4=(Y(1)-Y(4)-Y(2)+Y(3))*0.25D0
C
C *** MAKE INITIAL GUESS OF XSI AND ETA
C
      XSIO=(2.0D0*XP-(X1+X2))/DELX
      IF(XSIO.LT.-1.0D0)XSIO=-1.0D0
      IF(XSIO.GT.1.0D0)XSIO=1.0D0
      ETAO=(2.0D0*YP-(Y1+Y2))/DELY
      IF(ETAO.LT.-1.0D0)ETAO=-1.0D0
      IF(ETAO.GT.1.0D0)ETAO=1.0D0
      XSIW=XSIO
      XSI=XSIO
      ETAW=ETAO
      ETA=ETAO
C
C *** START NEWTON RALSON ITERATION LOOP
C
      DO 590 ITER=1,NITER
C
C *** COMPUTE RIGHT HAND SIDE AND THE JACOBIAN
C
        F1=XP-A1-A2*XSIW-A3*ETAW-A4*XSIW*ETAW
        F2=YP-B1-B2*XSIW-B3*ETAW-B4*XSIW*ETAW
        Z11=-A2-A4*ETAW
        Z12=-A3-A4*XSIW
        Z21=-B2-B4*ETAW
        Z22=-B3-B4*XSIW
C
C *** SOLVE FOR DELXSI AND DELETA
C
        DJAC=Z11*Z22 - Z21*Z12
        DELXSI=(F1*Z22-F2*Z12)/DJAC
        DELETA=(-F1*Z21+F2*Z11)/DJAC
C
C ***  COMPUTE FOR NEW XSI AND ETA
C
        XSI=XSIO-DELXSI
        ETA=ETAO-DELETA
C
C *** TEST CONVERGENCE
C
        DIFMAX=0.0D0
        DIFMAX1=0.0D0
        IF(XSIW.NE.0.)THEN
          DIF=DABS((XSI-XSIW)/XSIW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(ETAW.NE.0.)THEN
          DIF=DABS((ETA-ETAW)/ETAW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(ETA-ETAW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DIFMAX.LE.EPSR .AND. ITER.GT.1)GOTO 600
C
C *** TAKE A NEW GUESS FOR XSI AND ETA
C
        IF(XSI.LT.-1.0D0)XSI=-1.0D0
        IF(XSI.GT.1.0D0)XSI=1.0D0
        IF(ETA.LT.-1.0D0)ETA=-1.0D0
        IF(ETA.GT.1.0D0)ETA=1.0D0
        XSIW=OMEGA*XSI+(1.0D0-OMEGA)*XSIO
        ETAW=OMEGA*ETA+(1.0D0-OMEGA)*ETAO
C
C *** UPDATE XSIO AND ETAO
C
        XSIO=XSI
        ETAO=ETA
  590 CONTINUE
C
C *** CONVERGENCE FAILS, STOP EXECUTION
C
c     WRITE(16,2000)ITER,NITER,M,DIFMAX,DIFMAX1,EPSR
c2000 FORMAT('1','* FAIL TO CONVERGE IN COMPUTING XSI, ETA, AND ZTA:',
c    > /1X,'ITER =',I3,',   NITER =',I3,',   M =',I4/1X,
c    > 'DIFMAX =',D15.6,',  DIFMAX1=',D15.6,',   EPS =',D15.6)
c     WRITE(16,*)'IT CAN NOT REACH CONVERGENCE ----- WARNING AT XSI2D'
      IF(DIFMAX1.GT.EPSR)THEN
        WRITE(16,997)XSI,ETA
        WRITE(16,998)XP,YP
        WRITE(16,*)'XX(I),YY(I) ARE AS THE FOLLOWING:'
        WRITE(16,999)(X(I),Y(I),I=1,4)
  997 FORMAT('XSI=',E15.8,3X,'ETA=',E15.8)
  998 FORMAT('XP=',E15.8,3X,'YP=',E15.8)
  999 FORMAT(2F15.8)
        STOP
      ENDIF
C
C ------- CONVERGENT SOLUTION FOR XSI AND ETA HAS BEEN ACHIEVED.
C
  600 CONTINUE
C
      IF(XSI.LT.-1.0D0)XSI=-1.0D0
      IF(XSI.GT.1.0D0)XSI=1.0D0
      IF(ETA.LT.-1.0D0)ETA=-1.0D0
      IF(ETA.GT.1.0D0)ETA=1.0D0
C
C ***** THIS IS THE END OF CALCULATIONS
C
      RETURN
      END
c
c
c
      SUBROUTINE XSI3D
c     I                (X,Y,Z,XP,YP,ZP,M,
     I                (X,Y,Z,XP,YP,ZP,
     O                 XSI,ETA,ZTA)
C
C $$$$$ TO SOLVE FOR LOCAL COORDINATE (XSI,ETA,ZTA) OF THE PARTICLE IN
C       IN THE ELEMENT GIVEN THE GLOBAL COORDINATE (XP,YP,ZP).
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(8),Y(8),Z(8)
C
      DATA NITER/100/,EPSR/1.0D-6/,OMEGA/1.0D0/
C
      X2=-1.0D38
      X1=1.0D38
      Y2=-1.0D38
      Y1=1.0D38
      Z2=-1.0D38
      Z1=1.0D38
C
C ***** START CALCULATIONS
C
      DO 190 IQ=1,8
      X2=DMAX1(X2,X(IQ))
      Y2=DMAX1(Y2,Y(IQ))
      Z2=DMAX1(Z2,Z(IQ))
      X1=DMIN1(X1,X(IQ))
      Y1=DMIN1(Y1,Y(IQ))
      Z1=DMIN1(Z1,Z(IQ))
  190 CONTINUE
      DELX=X2-X1
      DELY=Y2-Y1
      DELZ=Z2-Z1
C
      A1=(X(1)+X(2)+X(3)+X(4)+X(5)+X(6)+X(7)+X(8))*0.125D0
      A2=(-X(1)+X(2)+X(3)-X(4)-X(5)+X(6)+X(7)-X(8))*0.125D0
      A3=(-X(1)+X(4)-X(2)+X(3)-X(5)-X(6)+X(7)+X(8))*0.125D0
      A4=(X(1)-X(4)-X(2)+X(3)+X(5)-X(6)+X(7)-X(8))*0.125D0
      A5=(-X(1)-X(2)-X(3)-X(4)+X(5)+X(6)+X(7)+X(8))*0.125D0
      A6=(X(1)+X(2)-X(3)-X(4)-X(5)-X(6)+X(7)+X(8))*0.125D0
      A7=(X(1)-X(2)-X(3)+X(4)-X(5)+X(6)+X(7)-X(8))*0.125D0
      A8=(-X(1)+X(2)-X(3)+X(4)+X(5)-X(6)+X(7)-X(8))*0.125D0
      B1=(Y(1)+Y(2)+Y(3)+Y(4)+Y(5)+Y(6)+Y(7)+Y(8))*0.125D0
      B2=(-Y(1)+Y(2)+Y(3)-Y(4)-Y(5)+Y(6)+Y(7)-Y(8))*0.125D0
      B3=(-Y(1)+Y(4)-Y(2)+Y(3)-Y(5)-Y(6)+Y(7)+Y(8))*0.125D0
      B4=(Y(1)-Y(4)-Y(2)+Y(3)+Y(5)-Y(6)+Y(7)-Y(8))*0.125D0
      B5=(-Y(1)-Y(2)-Y(3)-Y(4)+Y(5)+Y(6)+Y(7)+Y(8))*0.125D0
      B6=(Y(1)+Y(2)-Y(3)-Y(4)-Y(5)-Y(6)+Y(7)+Y(8))*0.125D0
      B7=(Y(1)-Y(2)-Y(3)+Y(4)-Y(5)+Y(6)+Y(7)-Y(8))*0.125D0
      B8=(-Y(1)+Y(2)-Y(3)+Y(4)+Y(5)-Y(6)+Y(7)-Y(8))*0.125D0
      C1=(Z(1)+Z(2)+Z(3)+Z(4)+Z(5)+Z(6)+Z(7)+Z(8))*0.125D0
      C2=(-Z(1)+Z(2)+Z(3)-Z(4)-Z(5)+Z(6)+Z(7)-Z(8))*0.125D0
      C3=(-Z(1)+Z(4)-Z(2)+Z(3)-Z(5)-Z(6)+Z(7)+Z(8))*0.125D0
      C4=(Z(1)-Z(4)-Z(2)+Z(3)+Z(5)-Z(6)+Z(7)-Z(8))*0.125D0
      C5=(-Z(1)-Z(2)-Z(3)-Z(4)+Z(5)+Z(6)+Z(7)+Z(8))*0.125D0
      C6=(Z(1)+Z(2)-Z(3)-Z(4)-Z(5)-Z(6)+Z(7)+Z(8))*0.125D0
      C7=(Z(1)-Z(2)-Z(3)+Z(4)-Z(5)+Z(6)+Z(7)-Z(8))*0.125D0
      C8=(-Z(1)+Z(2)-Z(3)+Z(4)+Z(5)-Z(6)+Z(7)-Z(8))*0.125D0
C
C *** MAKE INITIAL GUESS OF XSI, ETA, AND ZTA
C
      XSIO=(2.0D0*XP-(X1+X2))/DELX
      IF(XSIO.LT.-1.0D0) XSIO=-1.0D0
      IF(XSIO.GT.1.0D0) XSIO=1.0D0
      ETAO=(2.0D0*YP-(Y1+Y2))/DELY
      IF(ETAO.LT.-1.0D0) ETAO=-1.0D0
      IF(ETAO.GT.1.0D0) ETAO=1.0D0
      ZTAO=(2.0D0*ZP-(Z1+Z2))/DELZ
      IF(ZTAO.LT.-1.0D0) ZTAO=-1.0D0
      IF(ZTAO.GT.1.0D0) ZTAO=1.0D0
      XSIW=XSIO
      XSI=XSIO
      ETAW=ETAO
      ETA=ETAO
      ZTAW=ZTAO
      ZTA=ZTAO
C
C *** START NEWTON RALSON ITERATION LOOP
C
      DO 590 ITER=1,NITER
C
C *** COMPUTE RIGHT HAND SIDE AND THE JACOBIAN
C
        F1=XP-A1-A2*XSIW-A3*ETAW-A4*XSIW*ETAW-A5*ZTAW-A6*ETAW*ZTAW-
     1     A7*ZTAW*XSIW-A8*XSIW*ETAW*ZTAW
        F2=YP-B1-B2*XSIW-B3*ETAW-B4*XSIW*ETAW-B5*ZTAW-B6*ETAW*ZTAW-
     1     B7*ZTAW*XSIW-B8*XSIW*ETAW*ZTAW
        F3=ZP-C1-C2*XSIW-C3*ETAW-C4*XSIW*ETAW-C5*ZTAW-C6*ETAW*ZTAW-
     1     C7*ZTAW*XSIW-C8*XSIW*ETAW*ZTAW
        Z11=-A2-A4*ETAW-A7*ZTAW-A8*ETAW*ZTAW
        Z12=-A3-A4*XSIW-A6*ZTAW-A8*ZTAW*XSIW
        Z13=-A5-A6*ETAW-A7*XSIW-A8*XSIW*ETAW
        Z21=-B2-B4*ETAW-B7*ZTAW-B8*ETAW*ZTAW
        Z22=-B3-B4*XSIW-B6*ZTAW-B8*ZTAW*XSIW
        Z23=-B5-B6*ETAW-B7*XSIW-B8*XSIW*ETAW
        Z31=-C2-C4*ETAW-C7*ZTAW-C8*ETAW*ZTAW
        Z32=-C3-C4*XSIW-C6*ZTAW-C8*ZTAW*XSIW
        Z33=-C5-C6*ETAW-C7*XSIW-C8*XSIW*ETAW
C
C *** SOLVE FOR DELXSI, DELETA, AND DELZT
C
        DJAC=Z11*(Z22*Z33-Z23*Z32) - Z21*(Z12*Z33-Z32*Z13) +
     1       Z31*(Z12*Z23-Z22*Z13)
        DELXSI=(F1*(Z22*Z33-Z32*Z23)-F2*(Z12*Z33-Z32*Z13)+
     1         F3*(Z12*Z23-Z22*Z13))/DJAC
        DELETA=(-F1*(Z21*Z33-Z31*Z23)+F2*(Z11*Z33-Z31*Z13)-
     1         F3*(Z11*Z23-Z21*Z13))/DJAC
        DELZTA=(F1*(Z21*Z32-Z31*Z22)-F2*(Z11*Z32-Z31*Z12)+
     1         F3*(Z11*Z22-Z21*Z12))/DJAC
C
C *** COMPUTE FOR NEW XSI, ETA, AND ZTA
C
        XSI=XSIO-DELXSI
        ETA=ETAO-DELETA
        ZTA=ZTAO-DELZTA
C
C *** TEST CONVERGENCE
C
        DIFMAX=0.0D0
        DIFMAX1=0.0D0
        IF(XSIW.NE.0.)THEN
          DIF=DABS((XSI-XSIW)/XSIW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(ETAW.NE.0.)THEN
          DIF=DABS((ETA-ETAW)/ETAW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(ETA-ETAW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(ZTAW.NE.0.)THEN
          DIF=DABS((ZTA-ZTAW)/ZTAW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(ZTA-ZTAW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DIFMAX.LE.EPSR .AND. ITER.GT.1)GOTO 600
C
C *** TAKE A NEW GUESS FOR XSI, ETA, AND ZTA
C
        IF(XSI.LT.-1.0D0) XSI=-1.0D0
        IF(XSI.GT.1.0D0) XSI=1.0D0
        IF(ETA.LT.-1.0D0) ETA=-1.0D0
        IF(ETA.GT.1.0D0) ETA=1.0D0
        IF(ZTA.LT.-1.0D0) ZTA=-1.0D0
        IF(ZTA.GT.1.0D0) ZTA=1.0D0
        XSIW=OMEGA*XSI+(1.0D0-OMEGA)*XSIO
        ETAW=OMEGA*ETA+(1.0D0-OMEGA)*ETAO
        ZTAW=OMEGA*ZTA+(1.0D0-OMEGA)*ZTAO
C
C *** UPDATE XSIO, ETAO, AND ZTAO
C
c       IF(ITER.GT.NITER-3)THEN
c         WRITE(16,*)'ITER=',ITER
c         WRITE(16,997)XSI,ETA,ZTA
c       ENDIF
        XSIO=XSI
        ETAO=ETA
        ZTAO=ZTA
C
  590 CONTINUE
C
C *** CONVERGENCE FAILS, STOP EXECUTION
C
c     WRITE(16,2000)ITER,NITER,M,DIFMAX,DIFMAX1,EPSR
c2000 FORMAT('1','* FAIL TO CONVERGE IN COMPUTING XSI, ETA, AND ZTA:',
c    > /1X,'ITER =',I3,',   NITER =',I3,',   M =',I4/1X,
c    > 'DIFMAX =',D15.6,',  DIFMAX1=',D15.6,',   EPS =',D15.6)
c     WRITE(16,*)'IT CAN NOT REACH CONVERGENCE ----- WARNING AT XSI3D'
      IF(DIFMAX1.GT.EPSR)THEN
        WRITE(16,997)XSI,ETA,ZTA
        WRITE(16,998)XP,YP,ZP
        WRITE(16,*)'XX(I),YY(I),ZZ(I) ARE AS THE FOLLOWING:'
        WRITE(16,999)(X(I),Y(I),Z(I),I=1,8)
  997 FORMAT('XSI=',E15.8,3X,'ETA=',E15.8,3X,'ZTA=',E15.8)
  998 FORMAT('XP=',E15.8,3X,'YP=',E15.8,3X,'ZP=',E15.8)
  999 FORMAT(3F15.8)
        STOP
      ENDIF
C
C *** CONVERGENT SOLUTION FOR XSI, ETA, AND ZTA HAS BEEN ACHIEVED.
C
  600 CONTINUE
C
      IF(XSI.LT.-1.0D0) XSI=-1.0D0
      IF(XSI.GT.1.0D0) XSI=1.0D0
      IF(ETA.LT.-1.0D0) ETA=-1.0D0
      IF(ETA.GT.1.0D0) ETA=1.0D0
      IF(ZTA.LT.-1.0D0) ZTA=-1.0D0
      IF(ZTA.GT.1.0D0) ZTA=1.0D0
C
C ***** THIS IS THE END OF CALCULATIONS
C
      RETURN
      END
c
c
c
      SUBROUTINE XSI3DP
c     I        (X,Y,Z,XP,YP,ZP,M,
     I        (X,Y,Z,XP,YP,ZP,
     O         XSI,DL1,DL2,DL3)
C
C $$$$$ TO SOLVE FOR L1, L2, L3, AND XSI OF THE PARTICLE IN
C       IN THE PENTAHEDRAL ELEMENT GIVEN THE GLOBAL COORDINATE
C       (XP,YP,ZP)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(8),Y(8),Z(8)
C
      DATA NITER/100/, EPSR/1.0D-8/, OMEGA/1.0D0/
C
C ----- INITIAL GUESS OF XSIO, L1O, L2O
C
      Z1AVG=(Z(1)+Z(2)+Z(3))/3.0D0
      Z2AVG=(Z(4)+Z(5)+Z(6))/3.0D0
      ZL=DABS(Z1AVG-Z2AVG)
      IF(ZL.GT.EPSR)THEN
        XSIO=(2.0D0*ZP-Z1AVG-Z2AVG)/ZL
      ELSE
        XSIO=0.0D0
      ENDIF
      DL1O=0.33333333D0
      DL2O=0.33333333D0
C
      XSIW=XSIO
      DL1W=DL1O
      DL2W=DL2O
C
C ----- COMPUTE COOEFFICIENTS TO BE USED IN NEWTON-RALPHSON METHOD
C
      A1=X(4)+X(1)
      A2=X(5)+X(2)
      A3=X(6)+X(3)
      A4=X(4)-X(1)
      A5=X(5)-X(2)
      A6=X(6)-X(3)
      B1=Y(4)+Y(1)
      B2=Y(5)+Y(2)
      B3=Y(6)+Y(3)
      B4=Y(4)-Y(1)
      B5=Y(5)-Y(2)
      B6=Y(6)-Y(3)
      C1=Z(4)+Z(1)
      C2=Z(5)+Z(2)
      C3=Z(6)+Z(3)
      C4=Z(4)-Z(1)
      C5=Z(5)-Z(2)
      C6=Z(6)-Z(3)
C
C ------- START NEWTON RALSON ITERATION LOOP
C
      DO 590 ITER=1,NITER
C
C ------ COMPUTE RIGHT HAND SIDE AND THE JACOBIAN
C
        F1=0.5D0*((A1-A3)*DL1W+(A2-A3)*DL2W)+
     >     0.5D0*((A4-A6)*DL1W+(A5-A6)*DL2W)*XSIW+
     >     0.5D0*A6*XSIW-(XP-0.5D0*A3)
        F2=0.5D0*((B1-B3)*DL1W+(B2-B3)*DL2W)+
     >     0.5D0*((B4-B6)*DL1W+(B5-B6)*DL2W)*XSIW+
     >     0.5D0*B6*XSIW-(YP-0.5D0*B3)
        F3=0.5D0*((C1-C3)*DL1W+(C2-C3)*DL2W)+
     >     0.5D0*((C4-C6)*DL1W+(C5-C6)*DL2W)*XSIW+
     >     0.5D0*C6*XSIW-(ZP-0.5D0*C3)
C
        Z11=0.5D0*(A1-A3)+0.5D0*(A4-A6)*XSIW
        Z12=0.5D0*(A2-A3)+0.5D0*(A5-A6)*XSIW
        Z13=0.5D0*((A4-A6)*DL1W+(A5-A6)*DL2W)+0.5D0*A6
        Z21=0.5D0*(B1-B3)+0.5D0*(B4-B6)*XSIW
        Z22=0.5D0*(B2-B3)+0.5D0*(B5-B6)*XSIW
        Z23=0.5D0*((B4-B6)*DL1W+(B5-B6)*DL2W)+0.5D0*B6
        Z31=0.5D0*(C1-C3)+0.5D0*(C4-C6)*XSIW
        Z32=0.5D0*(C2-C3)+0.5D0*(C5-C6)*XSIW
        Z33=0.5D0*((C4-C6)*DL1W+(C5-C6)*DL2W)+0.5D0*C6
C
C ------- SOLVE FOR DELXSI, DELETA, AND DELZT
C
        DJAC=Z11*(Z22*Z33-Z23*Z32) - Z21*(Z12*Z33-Z32*Z13) +
     1       Z31*(Z12*Z23-Z22*Z13)
        DELDL1=(F1*(Z22*Z33-Z32*Z23)-F2*(Z12*Z33-Z32*Z13)+
     1          F3*(Z12*Z23-Z22*Z13))/DJAC
        DELDL2=(-F1*(Z21*Z33-Z31*Z23)+F2*(Z11*Z33-Z31*Z13)-
     1          F3*(Z11*Z23-Z21*Z13))/DJAC
        DELXSI=(F1*(Z21*Z32-Z31*Z22)-F2*(Z11*Z32-Z31*Z12)+
     1          F3*(Z11*Z22-Z21*Z12))/DJAC
C
C ------- COMPUTE FOR NEW XSI, ETA, AND ZTA
C
        DL1=DL1O-DELDL1
        DL2=DL2O-DELDL2
        XSI=XSIO-DELXSI
C
C ------- TEST CONVERGENCE
C
        DIFMAX=0.0D0
        DIFMAX1=0.0D0
        IF(DL1W.NE.0.0)THEN
          DIF=DABS((DL1-DL1W)/DL1W)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(DL1-DL1W)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DL2W.NE.0.0)THEN
          DIF=DABS((DL2-DL2W)/DL2W)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(DL2-DL2W)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(XSIW.NE.0.0)THEN
          DIF=DABS((XSI-XSIW)/XSIW)
          DIFMAX=DMAX1(DIF,DIFMAX)
          DIF1=DABS(XSI-XSIW)
          DIFMAX1=DMAX1(DIF1,DIFMAX1)
        ENDIF
        IF(DIFMAX.LT.EPSR .AND. ITER.GT.1) GOTO 600
C
C ------- TAKE A NEW GUESS FOR XSI, ETA, AND ZTA
C
        IF(DL1.LT.0.0D0) DL1=0.0D0
        IF(DL1.GT.1.0D0) DL1=1.0D0
        IF(DL2.LT.0.0D0) DL2=0.0D0
        IF(DL2.GT.1.0D0) DL2=1.0D0
        IF(XSI.LT.-1.0D0) XSI=-1.0D0
        IF(XSI.GT.1.0D0) XSI=1.0D0
        DL1W=OMEGA*DL1+(1.0D0-OMEGA)*DL1O
        DL2W=OMEGA*DL2+(1.0D0-OMEGA)*DL2O
        XSIW=OMEGA*XSI+(1.0D0-OMEGA)*XSIO
C
C ------- UPDATE DL1O, DL2O, AND XSIO
C
        DL1O=DL1
        DL2O=DL2
        XSIO=XSI
C
  590 CONTINUE
C
C *** CONVERGENCE FAILS, STOP EXECUTION
C
c     WRITE(16,2000)ITER,NITER,M,DIFMAX,DIFMAX1,EPSR
c2000 FORMAT('1','* FAIL TO CONVERGE IN COMPUTING DL1, DL2, AND XSI:',
c    > /1X,'ITER =',I3,',   NITER =',I3,',   M =',I4/1X,
c    > 'DIFMAX =',D15.6,',  DIFMAX1=',D15.6,',   EPS =',D15.6)
c     WRITE(16,*)'IT CAN NOT REACH CONVERGENCE ----- WARNING AT XSI3DP'
      IF(DIFMAX1.GT.EPSR)THEN
        WRITE(16,997)DL1,DL2,XSI
        WRITE(16,998)XP,YP,ZP
        WRITE(16,*)'XX(I),YY(I),ZZ(I) ARE AS THE FOLLOWING:'
        WRITE(16,999)(X(I),Y(I),Z(I),I=1,8)
  997 FORMAT('DL1=',E15.8,3X,'DL2=',E15.8,3X,'XSI=',E15.8)
  998 FORMAT('XP=',E15.8,3X,'YP=',E15.8,3X,'ZP=',E15.8)
  999 FORMAT(3F15.8)
        STOP
      ENDIF
C
C ------- CONVERGENT SOLUTION FOR XSI, ETA, AND ZTA HAS BEEN ACHIEVED.
C
  600 CONTINUE
C
      IF(DL1.LT.0.0D0) DL1=0.0D0
      IF(DL1.GT.1.0D0) DL1=1.0D0
      IF(DL2.LT.0.0D0) DL2=0.0D0
      IF(DL2.GT.1.0D0) DL2=1.0D0
      IF(XSI.LT.-1.0D0) XSI=-1.0D0
      IF(XSI.GT.1.0D0) XSI=1.0D0
C
      DL3=1.0D0-DL1-DL2
C
      RETURN
      END
C
c
c
      SUBROUTINE MOVCHK
     I  (IBF,XP,YP,ZP,VPX,VPY,VPZ,N,MM,DELT,TMAX,RAMADA,IDTI,IZOOM,
     I   X,CP,IE,XPFG,CPFG,NFGMB,ISE, MAXNP,MAXEL,MXNPFG,MXKGL,MXNEP,
     M   MXNODS,EPSX, NEL,NEFG,NPFGS,NPFGEP,NODE,XSFG,CSFG,MPLOCS,
     O   XEFG,CEFG,MPLOCE,DTIFG,  ICODE)
C
C --------------------------------------------------------------------
C **** THIS ROUTINE DETERMINES THE CONCENTRATIONS AS THE POINT
C      IS NOT GOING TO MOVE FOREVER
C --------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG),DTIFG(MXNPFG)
      DIMENSION XSFG(MXNODS,3),CSFG(MXNODS),MPLOCS(MXNODS)
      DIMENSION XEFG(MXNEP,3),CEFG(MXNEP),MPLOCE(MXNPFG)
      DIMENSION X(MAXNP,3),CP(MAXNP),IE(MAXEL,11),NFGMB(MAXEL)
      DIMENSION ISE(MXKGL,8),CQ(7)
C
      ICODE=-1
        D1=SQRT(VPX**2+VPY**2+VPZ**2)*TMAX
        IF(DABS(D1).LE.EPSX)THEN
          IF(IBF.EQ.1)THEN
            IF(IDTI.NE.1)THEN
              CALL INTERP
     I            (MAXNP,MAXEL,MXNPFG,MXKGL,8,1,NODE,MM,XP,YP,ZP,X,CP,
     I             IE,XPFG,CPFG,NEL,NEFG,IZOOM,10,ISE,NFGMB,1,1,EPSX,
     O             CQ)
              CSFG(N)=CQ(1)*DEXP(-RAMADA*DELT)
              MPLOCS(N)=MM
              IF(IDTI.EQ.2)DTIFG(N)=1.0D30
            ELSE
              DTIFG(N)=1.0D30
            ENDIF
          ELSE
            IF(MPLOCE(N).EQ.0)THEN
              NPFGS=NPFGS+1
              CALL WARMSG(NPFGS,MXNPFG,'MOVCHK','MXNPFG',1)
              XSFG(NPFGS,1)=XPFG(N,1)
              XSFG(NPFGS,2)=XPFG(N,2)
              XSFG(NPFGS,3)=XPFG(N,3)
              CSFG(NPFGS)=CPFG(N)*DEXP(-RAMADA*DELT)
              MPLOCS(NPFGS)=MM
            ELSE
              NPFGEP=NPFGEP+1
              CALL WARMSG(NPFGEP,MXNEP,'MOVCHK','MXNEP ',2)
              XEFG(NPFGEP,1)=XPFG(N,1)
              XEFG(NPFGEP,2)=XPFG(N,2)
              XEFG(NPFGEP,3)=XPFG(N,3)
              CEFG(NPFGEP)=CPFG(N)*DEXP(-RAMADA*DELT)
              MPLOCE(NPFGEP)=MM
            ENDIF
          ENDIF
          ICODE=0
        ENDIF
      RETURN
      END
c
c
c
      SUBROUTINE FIXCHK
     I    (IFIX,IDTI,N,NN,M,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I     EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I     X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I     CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O     NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C -------------------------------------------------------------------
C ***** THIS ROUTINE IS USED TO TAKE CARE THE POINT FIXED, I.E. IT
C       WON'T MOVE INTO THE REGION OF INTEREST
C $$$ IFIX: THE INDEX OF THE LOCATION IN GNTRAC OR HPTRAC CALLING
C            FIXCHK
C $$$ IDTI: THE INDEX DETERMINING DTI CALCUATIONS
C          -1 : DTI AND CS HAVE TO BE CALCULATED AND IS BACK TO GNTRAC
C           0 : ONLY CONC. BUT NO DTI AND BACK TO HPTRAC
C           1 : ONLY DTIFG AND BACK TO HPTRAC
C           2 : BOTH DTIFG AND CSFG ARE BACK TO HPTRAC WHICH IS CALLED
C               BY FGDET
C --------------------------------------------------------------------
C
      LOGICAL BOUND,SPECFY
      DIMENSION CS(MAXNP),CP(MAXNP),X(MAXNP,3),VX(MAXNP,3),DTI(MAXNP)
      DIMENSION IE(MAXEL,11),IB(MAXNP),ISB(6,MAXBES),DCOSB(3,MAXBES)
      DIMENSION NBDYB(MAXNP),IBDY(MXTUBS),XPFG(MXNPFG,3),CPFG(MXNPFG)
      DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG),ISE(MXKGL,8),NFGMB(MAXEL)
      DIMENSION MPLOCS(MXNPFG),DTIFG(MXNPFG)
      DIMENSION XEFG(MXNEP,3),CEFG(MXNEP),MPLOCE(MXNPFG)
      DIMENSION CQ(7)
C
      BOUND=.FALSE.
      SPECFY=.FALSE.
      IF(IFIX.EQ.1) THEN
        IF(IB(N).NE.0)THEN
          BOUND=.TRUE.
          IGLOBL=1
          N1=NN
        ENDIF
        IF(IB(N).LE.3)SPECFY=.TRUE.
      ELSEIF(IFIX.EQ.2) THEN
        IF(IB(NN).NE.0)THEN
          BOUND=.TRUE.
          IGLOBL=1
          N1=NN
        ENDIF
        IF(IB(NN).LE.3)SPECFY=.TRUE.
      ELSEIF(IFIX.EQ.3)THEN
        IF(NODE.EQ.8 .OR. (NODE.EQ.6 .AND. N4.NE.0))THEN
          BOUND=.TRUE.
          IGLOBL=3
          IF(IB(N1).LE.3 .AND. IB(N2).LE.3 .AND. IB(N3).LE.3 .AND.
     >       IB(N4).LE.3)SPECFY=.TRUE.
        ELSE
          BOUND=.TRUE.
          IGLOBL=3
          N4=0
          IF(IB(N1).LE.3 .AND. IB(N2).LE.3 .AND. IB(N3).LE.3)
     >       SPECFY=.TRUE.
        ENDIF
      ELSEIF(IFIX.EQ.4 .AND. IB(N1)*IB(N2).NE.0)THEN
        BOUND=.TRUE.
        IGLOBL=2
        IF(IB(N1).LE.3 .AND. IB(N2).LE.3)SPECFY=.TRUE.
      ENDIF
C
      ICODE=-1
      IF(BOUND)THEN
        IF(SPECFY)THEN
          IF(IBF.EQ.1)THEN
            IF(IFIX.EQ.1)THEN
              DTI(N)=1.0D30
            ELSEIF(IFIX.EQ.2)THEN
              DTREAL=DELT-SDT1
              IF(IDTI.LT.0)THEN
                DTI(N)=1.0D0/DTREAL
                CS(N)=CP(NN)+(CS(NN)-CP(NN))*SDT1/DELT
                CS(N)=CS(N)*DEXP(-RAMADA*DTREAL)
              ELSEIF(IDTI.EQ.1)THEN
                DTIFG(N)=1.0D0/DTREAL
                ICODE=0
              ELSE
                CSFG(N)=CP(NN)+(CS(NN)-CP(NN))*SDT1/DELT
                CSFG(N)=CSFG(N)*DEXP(-RAMADA*DTREAL)
                IF(IDTI.EQ.2)DTIFG(N)=1.0D0/DTREAL
              ENDIF
            ELSEIF(IFIX.EQ.3 .OR. IFIX.EQ.4)THEN
              IF(IDTI.NE.0)THEN
                DTREAL=DELT-SDT1
              ELSE
                DTREAL=0.0D0
              ENDIF
              IF(DTREAL.EQ.0.0D0)THEN
                DTIVS=1.0D30
              ELSE
                DTIVS=1.0D0/DTREAL
              ENDIF
              IF(IDTI.EQ.1)THEN
                DTIFG(N)=DTIVS
                ICODE=0
                RETURN
              ELSEIF(IDTI.EQ.2)THEN
                DTIFG(N)=DTIVS
              ELSEIF(IDTI.LT.0)THEN
                DTI(N)=DTIVS
              ENDIF
C
              CALL ELENOD
     I            (IE(M,5),IE(M,7),
     O             NODE,KK,KK)
              CALL INTERP
     I            (MAXNP,MAXEL,MXNPFG,MXKGL,8,1,NODE,M,XP,YP,ZP,X,CP,
     I             IE,XPFG,CPFG,NEL,NEFG,IZOOM,10,ISE,NFGMB,1,1,EPSX,
     O             CQ)
              IF(IDTI.LT.0)THEN
                CS(N)=CQ(1)*DEXP(-RAMADA*DTREAL)
              ELSE
                CSFG(N)=CQ(1)*DEXP(-RAMADA*DTREAL)
              ENDIF
            ENDIF
          ENDIF
          ICODE=0
        ELSE
          CALL ALGBDY
     I       (XP,YP,ZP,X,IE,IB,IGLOBL,EPSX,MAXNP,MAXEL,MAXBES,
     I        MXTUBS,IDETQ, IBF, NNP,VPX,VPY,VPZ,VX,SDT,
     I        NTUBS,ISB,DCOSB,NBDYB,IBDY,
c     M        N1,N2,N3,N4,
     M        N1,N2,N3,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,M)
          IF(IBF.EQ.1)THEN
            IF(SDT1.EQ.0.0D0)THEN
              DTREAL=DELT
            ELSE
              DTREAL=DELT-SDT1
            ENDIF
            IF(DTREAL.EQ.0.0D0)THEN
              DTIVS=1.0D30
            ELSE
              DTIVS=1.0D0/DTREAL
            ENDIF
            IF(IDTI.EQ.1)THEN
              DTIFG(N)=DTIVS
              ICODE=0
              RETURN
            ENDIF
            IF(IDTI.EQ.2)DTIFG(N)=DTIVS
            IF(IDTI.LT.0)DTI(N)=DTIVS
C
            CALL ELENOD
     I          (IE(M,5),IE(M,7),
     O           NODE,KK,KK)
            CALL INTERP
     I          (MAXNP,MAXEL,MXNPFG,MXKGL,8,1,NODE,M,XQ,YQ,ZQ,X,CP,
     I           IE,XPFG,CPFG,NEL,NEFG,IZOOM,10,ISE,NFGMB,1,1,EPSX,
     O           CQ)
            IF(IDTI.LT.0)THEN
              CS(N)=CQ(1)*DEXP(-RAMADA*DTREAL)
            ELSE
              CSFG(N)=CQ(1)*DEXP(-RAMADA*DTREAL)
            ENDIF
          ENDIF
          ICODE=0
        ENDIF
C
        IF(IBF.EQ.2)THEN
          IF(SDT1.EQ.0.0D0)THEN
            IF(XQ.NE.ZP .OR. YQ.NE.YP .OR. ZQ.NE.ZP)THEN
              IE(M,11)=M
C
c             IF(IDTI.GE.0)THEN
              CALL ELENOD(IE(M,5),IE(M,7),NODE1,KK,KK)
              DO KK=1,NODE1
                IEMKK=IE(M,KK)
                XR=X(IEMKK,1)
                YR=X(IEMKK,2)
                ZR=X(IEMKK,3)
                DIST=DSQRT((XQ-XR)**2+(YQ-YR)**2+(ZQ-ZR)**2)
                IF(DIST.LE.EPSX)THEN
                  IF(IDTI.LT.0)CPFGN=CPFG(N)
                  IF(IDTI.EQ.0)CPFGN=CP(N)
                  CS(IEMKK)=CPFGN*DEXP(-RAMADA*DELT)
                  GOTO 146
                ENDIF
              ENDDO
  146         CONTINUE
c           ENDIF
C
              IDO=1
              IF(IDTI.LT.0)THEN
                IF(MPLOCE(N).NE.0)IDO=0
                CPFGN=CPFG(N)
              ELSE
                CPFGN=CP(N)
              ENDIF
              IF(IDO.EQ.1)THEN
                NPFGS=NPFGS+1
                CALL WARMSG(NPFGS,MXNPFG,'FIXCHK','MXNPFG',1)
                XSFG(NPFGS,1)=XQ
                XSFG(NPFGS,2)=YQ
                XSFG(NPFGS,3)=ZQ
                CSFG(NPFGS)=CPFGN*DEXP(-RAMADA*DELT)
                MPLOCS(NPFGS)=M
              ELSE
                NPFGEP=NPFGEP+1
                CALL WARMSG(NPFGEP,MXNEP,'FIXCHK','MXNEP ',2)
                XEFG(NPFGEP,1)=XQ
                XEFG(NPFGEP,2)=YQ
                XEFG(NPFGEP,3)=ZQ
                CEFG(NPFGEP)=CPFGN*DEXP(-RAMADA*DELT)
                MPLOCE(NPFGEP)=M
              ENDIF
            ENDIF
          ENDIF
          ICODE=0
        ENDIF
C
      ENDIF
C
C
C ----- ERROR OCCURRED
C
      IF(ICODE.EQ.-1 .AND. IFIX.NE.3)THEN
        WRITE(16,*)'ERROR 1 AT FIXCHK, NODE',N,' CAN NOT BE',
     >            ' TRACKED'
        WRITE(16,1006)X(N,1),X(N,2),X(N,3),VX(N,1),VX(N,2),VX(N,3)
        WRITE(16,1007)XP,YP,ZP,VPX,VPY,VPZ
 1006   FORMAT('X(N)=',F12.6,2X,'Y(N)=',F12.6,2X,'Z(N)=',F12.6,2X,
     >         'VX(N)=',F12.6,2X,'VY(N)=',F12.6,2X,'VZ(N)=',F12.6,2X)
 1007   FORMAT('XP=',F12.6,2X,'YP=',F12.6,2X,'ZP=',F12.6,2X,'VPX=',
     >         F12.6,2X,'VPY=',F12.6,2X,'VPZ=',F12.6,2X)
        WRITE(16,*)'SDT=',SDT
        STOP
      ENDIF
C
      RETURN
      END
C
c
c
      SUBROUTINE ALGBDY
     I     (XXP,YYP,ZZP,X,IE,IB,IGLOBL,EPSX,
     I      MAXNP,MAXEL,MAXBES,MXTUBS,IDETQ,IBF,NNP,VXP,VYP,VZP,
     I      VX,SDT,NTUBS,ISB,DCOSB,NBDYB,IBDY,
c     M      N1,N2,N3,N4,
     M      N1,N2,N3,
     O      XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,M)
C  11/ 9/93
C $$$$$ TO COMPUTE PARTICAL TRACKING ALONG THOSE UNSPECIFIED BOUNDARIES
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP,3),IB(MAXNP),ISB(6,MAXBES),DCOSB(3,MAXBES)
      DIMENSION IE(MAXEL,11),NBDYB(MAXNP),IBDY(MXTUBS),VX(MAXNP,3)
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8)
c     DIMENSION KGB(4,6,3)
C
c     DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
c    >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
c    >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
      DTRES=SDT
      CALL REPLAS(XXP,YYP,ZZP,VXP,VYP,VZP,XP,YP,ZP,VPX,VPY,VPZ)
C
    5 CONTINUE
C
      IF(IGLOBL.EQ.1)THEN
C
C ***** IF THIS IS A GLOBAL NODE, ALL THE UNSPECIFIED BOUNDARY SIEDS
C       CONNECTING TO THIS NODE ARE TO BE CHECKED.
C
        N=N1
        NS1=NBDYB(N)+1
        IF(N.NE.NNP)THEN
          NS2=NBDYB(N+1)
        ELSE
          NS2=NTUBS
        ENDIF
        IF(NS1.GT.NS2)THEN
          SDT1=DTRES
          RETURN
        ENDIF
C
        IF(IDETQ.EQ.2)THEN
          IJUDGE=2
        ELSE
          IJUDGE=1
        ENDIF
C
   6    CONTINUE
C
        DO 10 NS=NS1,NS2
          KS=IBDY(NS)
          M=ISB(6,KS)
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NODE,I,IK)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >                  XX,YY,ZZ,VXX,VYY,VZZ)
C
          CALL BNDRY
     I      (MAXEL,MAXBES,M,NODE,IJUDGE,IBF,EPSX,KS,XP,YP,ZP,VPX,VPY,
     I       VPZ,DTRES,ISB,IE,DCOSB,XX,YY,ZZ,VXX,VYY,VZZ,
     O       XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,
c     M       IGLOBL,N1,N2,N3,N4)
     M       IGLOBL,N1,N2)
C
          IF(SDT1.NE.DTRES)GOTO 100
   10   CONTINUE
C
        IF(IJUDGE.EQ.1)THEN
          IJUDGE=2
          GOTO 6
        ENDIF
C
C ----- NO MORE PARTICAL TRACKING ALONG UNSPECIFIED BOUNDARY, THEN
C       RECORD THE BOUNDARY SIDE THE PARTICAL IS RIGHT ON.
C       N1, N2, N3, AND N4 ARE THE GLOBAL NODES OF THE BOUNDARY SIDE.
C       THE REASON WHY N1-N4 NEED TO BE RECORDED IS THAT WE NEED TO DO
C       INTERPOLATION AFTER GOING BACK TO  GNTRAC.
C
c       J1=KGB(1,LS,IK)
c       J2=KGB(2,LS,IK)
c       J3=KGB(3,LS,IK)
c       J4=KGB(4,LS,IK)
c       J1=IE(M,J1)
c       J2=IE(M,J2)
c       J3=IE(M,J3)
c       IF(J4.NE.0)THEN
c         J4=IE(M,J4)
c       ELSE
c         J4=0
c       ENDIF
C
c       IF(J1.EQ.N1)THEN
c         N2=J2
c         N3=J3
c         N4=J4
c       ELSEIF(J2.EQ.N1)THEN
c         N2=J3
c         IF(J4.NE.0)THEN
c           N3=J4
c           N4=J1
c         ELSE
c           N3=J1
c           N4=0
c         ENDIF
c       ELSEIF(J3.EQ.N1)THEN
c         IF(J4.NE.0)THEN
c           N2=J4
c           N3=J1
c           N4=J2
c         ELSE
c           N2=J1
c           N3=J2
c           N4=0
c         ENDIF
c       ELSE
c         N2=J1
c         N3=J2
c         N4=J3
c       ENDIF
      ENDIF
C
      IF(IGLOBL.EQ.2)THEN
C
C ***** ALL THE UNSPECIFIED BOUNDARY SIEDS CONNECTING TO BOTH GLOBAL
C       NODES N1 AND N2 ARE TO BE CHECKED.
C
        NS1=NBDYB(N1)+1
        IF(N1.NE.NNP)THEN
          NS2=NBDYB(N1+1)
        ELSE
          NS2=NTUBS
        ENDIF
        IF(NS1.GT.NS2)THEN
          SDT1=DTRES
          RETURN
        ENDIF
        NS3=NBDYB(N2)+1
        IF(N2.NE.NNP)THEN
          NS4=NBDYB(N2+1)
        ELSE
          NS4=NTUBS
        ENDIF
        IF(NS3.GT.NS4)THEN
          SDT1=DTRES
          RETURN
        ENDIF
C
        IF(IDETQ.EQ.2)THEN
          IJUDGE=2
        ELSE
          IJUDGE=1
        ENDIF
C
  16    CONTINUE
C
        DO 20 NS=NS1,NS2
          KS1=IBDY(NS)
          DO NSS=NS3,NS4
            KS2=IBDY(NSS)
            IF(KS1.EQ.KS2)THEN
              KS=KS1
              M=ISB(6,KS)
              CALL ELENOD
     I           (IE(M,5),IE(M,7),
     O            NODE,I,IK)
              CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >                    XX,YY,ZZ,VXX,VYY,VZZ)
C
              CALL BNDRY
     I          (MAXEL,MAXBES,M,NODE,IJUDGE,IBF,EPSX,
     I           KS,XP,YP,ZP,VPX,VPY,VPZ,DTRES,
     I           ISB,IE,DCOSB,XX,YY,ZZ,VXX,VYY,VZZ,
     O           XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,
c     M           IGLOBL,N1,N2,N3,N4)
     M           IGLOBL,N1,N2)
C
              IF(SDT1.NE.DTRES) GOTO 100
            ENDIF
          ENDDO
   20   CONTINUE
C
        IF(IJUDGE.EQ.1)THEN
          IJUDGE=2
          GOTO 16
        ENDIF
C
C ----- NO MORE PARTICAL TRACKING ALONG UNSPECIFIED BOUNDARY, THEN
C       RECORD THE BOUNDARY SIDE THE PARTICAL IS RIGHT ON.
C       N1, N2, N3, AND N4 ARE THE GLOBAL NODES OF THE BOUNDARY SIDE.
C       THE REASON WHY N1-N4 NEED TO BE RECORDED IS THAT WE NEED TO DO
C       INTERPOLATION AFTER GOING BACK TO GNTRAC.
C
c       J1=KGB(1,LS,IK)
c       J2=KGB(2,LS,IK)
c       J3=KGB(3,LS,IK)
c       J4=KGB(4,LS,IK)
c       J1=IE(M,J1)
c       J2=IE(M,J2)
c       J3=IE(M,J3)
c       IF(J4.NE.0)THEN
c         J4=IE(M,J4)
c       ELSE
c         J4=0
c       ENDIF
c
c       IF(J1.EQ.N1)THEN
c         IF(J2.EQ.N2)THEN
c           N3=J3
c           N4=J4
c         ENDIF
c         IF(J4.EQ.N2)THEN
c           N3=J3
c           N4=J2
c         ENDIF
c         IF(J3.EQ.N2)THEN
c           N3=J2
c           N4=0
c         ENDIF
c       ELSEIF(J2.EQ.N1)THEN
c         IF(J3.EQ.N2)THEN
c           IF(J4.NE.0)THEN
c             N3=J4
c             N4=J1
c           ELSE
c             N3=J1
c             N4=0
c           ENDIF
c         ENDIF
c         IF(J1.EQ.N2)THEN
c           IF(J4.NE.0)THEN
c             N3=J4
c             N4=J3
c           ELSE
c             N3=J3
c             N4=0
c           ENDIF
c         ENDIF
c       ELSEIF(J3.EQ.N1)THEN
c         IF(J4.EQ.N2)THEN
c           N3=J1
c           N4=J2
c         ENDIF
c         IF(J1.EQ.N2)THEN
c           N3=J2
c           N4=0
c         ENDIF
c         IF(J2.EQ.N2)THEN
c           N3=J1
c           N4=J4
c         ENDIF
c       ELSE
c         IF(J1.EQ.N2)THEN
c           N3=J2
c           N4=J3
c         ENDIF
c         IF(J3.EQ.N2)THEN
c           N3=J2
c           N4=J1
c         ENDIF
c       ENDIF
      ENDIF
c
      IF(IGLOBL.EQ.3)THEN
C
C ***** THE UNSPECIFIED BOUNDARY SIED CONNECTING TO GLOBAL NODES N1, N2,
C       N3, AND N4 ARE TO BE CHECKED.
C
        NS1=NBDYB(N1)+1
        IF(N1.NE.NNP)THEN
          NS2=NBDYB(N1+1)
        ELSE
          NS2=NTUBS
        ENDIF
        IF(NS1.GT.NS2)THEN
          SDT1=DTRES
          RETURN
        ENDIF
        NS3=NBDYB(N2)+1
        IF(N2.NE.NNP)THEN
          NS4=NBDYB(N2+1)
        ELSE
          NS4=NTUBS
        ENDIF
        IF(NS3.GT.NS4)THEN
          SDT1=DTRES
          RETURN
        ENDIF
        NS5=NBDYB(N3)+1
        IF(N3.NE.NNP)THEN
          NS6=NBDYB(N3+1)
        ELSE
          NS6=NTUBS
        ENDIF
        IF(NS5.GT.NS6)THEN
          SDT1=DTRES
          RETURN
        ENDIF
C
        IF(IDETQ.EQ.2)THEN
          IJUDGE=2
        ELSE
          IJUDGE=1
        ENDIF
C
  26    CONTINUE
C
        DO 30 NS=NS1,NS2
          KS1=IBDY(NS)
          DO 29 NSS=NS3,NS4
            KS2=IBDY(NSS)
            IF(KS1.NE.KS2)GOTO 29
            DO 28 NSSS=NS5,NS6
              KS3=IBDY(NSSS)
              IF(KS3.NE.KS1)GOTO 28
              KS=KS1
              M=ISB(6,KS)
              CALL ELENOD
     I            (IE(M,5),IE(M,7),
     O             NODE,I,I)
              CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >                    XX,YY,ZZ,VXX,VYY,VZZ)
C
              CALL BNDRY
     I          (MAXEL,MAXBES,M,NODE,IJUDGE,IBF,EPSX,
     I           KS,XP,YP,ZP,VPX,VPY,VPZ,DTRES,
     I           ISB,IE,DCOSB,XX,YY,ZZ,VXX,VYY,VZZ,
     O           XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,
c     M           IGLOBL,N1,N2,N3,N4)
     M           IGLOBL,N1,N2)
C
              IF(SDT1.NE.DTRES) GOTO 100
   28       CONTINUE
   29     CONTINUE
   30   CONTINUE
C
        IF(IJUDGE.EQ.1)THEN
          IJUDGE=2
          GOTO 26
        ENDIF
      ENDIF
C
C ----- CANNOT BE TRACKED ALONG THE BOUNDARY ANY MORE, THEN SDT1 IS
C       FORCE TO BE ZERO NO MATTER HOW MUCH DTRES IS.
C
      SDT1=0.0D0
      XQ=XP
      YQ=YP
      ZQ=ZP
      RETURN
C
  100 CONTINUE
C
      IF(SDT1.EQ.0.0D0)RETURN
C
      IF(IGLOBL.EQ.1)THEN
        IF(IB(N1).LE.3)RETURN
      ELSE
        IF(IB(N1).LE.3 .AND. IB(N2).LE.3)THEN
          NS1=NBDYB(N1)+1
          IF(N1.NE.NNP)THEN
            NS2=NBDYB(N1+1)
          ELSE
            NS2=NTUBS
          ENDIF
          IF(NS1.GT.NS2)THEN
            SDT1=DTRES
            RETURN
          ENDIF
          NS3=NBDYB(N2)+1
          IF(N2.NE.NNP)THEN
            NS4=NBDYB(N2+1)
          ELSE
            NS4=NTUBS
          ENDIF
          IF(NS3.GT.NS4)THEN
            SDT1=DTRES
            RETURN
          ENDIF
          DO 120 NS=NS1,NS2
            KS1=IBDY(NS)
            DO NSS=NS3,NS4
              KS2=IBDY(NSS)
              IF(KS1.EQ.KS2 .AND. KS1.NE.KS)GOTO 150
            ENDDO
  120     CONTINUE
          RETURN
        ENDIF
      ENDIF
C
C ----- TRY NEXT TRACKING
C
  150 CONTINUE
      DTRES=SDT1
      CALL REPLAS(XQ,YQ,ZQ,VQX,VQY,VQZ,XP,YP,ZP,VPX,VPY,VPZ)
      GOTO 5
C
c      RETURN
      END
C
c
c
      SUBROUTINE BNDRY
     I    (MAXEL,MAXBES,M,NODE,IJUDGE,IBF,EPSX,KS,XP,YP,ZP,VPX,VPY,VPZ,
     I     SDT,ISB,IE,DCOSB,XX,YY,ZZ,VXX,VYY,VZZ,
     O     XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,
c     M     ID,J1,J2,J3,J4)
     M     ID,J1,J2)
C  10/14/93
C $$$$$ TO COMPUTE PARTICAL TRACKING ALONG AN UNSPECIFIED BOUNDARY SIDE
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION KGB(4,6,3)
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES),IE(MAXEL,11)
      DIMENSION XXX(8),YYY(8),ZZZ(8),VXXX(8),VYYY(8),VZZZ(8)
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
      IF(IJUDGE.EQ.1)IJUDGE1=1
      IF(IJUDGE.EQ.2)IJUDGE1=2
C
      LS=ISB(5,KS)
      DN1=DCOSB(1,KS)
      DN2=DCOSB(2,KS)
      DN3=DCOSB(3,KS)
C
      IF(NODE.EQ.6 .AND. LS.EQ.5)GOTO 100
      IF(NODE.EQ.8 .AND. (LS.EQ.2 .OR. LS.EQ.3 .OR. LS.EQ.6))GOTO 100
      DN1=-DN1
      DN2=-DN2
      DN3=-DN3
C
  100 CONTINUE
      KNODE=4
      IF(NODE.EQ.6 .AND. (LS.EQ.4 .OR. LS.EQ.5))KNODE=3
      IF(NODE.EQ.4)KNODE=3
C
C ----- COMPUTE THE VELOCITIES ALONG THE BOUNDARY SIDE
C
      IF(NODE.EQ.4)THEN
        IK=3
      ELSEIF(NODE.EQ.6)THEN
        IK=2
      ELSEIF(NODE.EQ.8)THEN
        IK=1
      ENDIF
      K1=KGB(1,LS,IK)
      N1=IE(M,K1)
      XXX(1)=XX(K1)
      YYY(1)=YY(K1)
      ZZZ(1)=ZZ(K1)
      VXV=VYY(K1)*DN3-VZZ(K1)*DN2
      VYV=VZZ(K1)*DN1-VXX(K1)*DN3
      VZV=VXX(K1)*DN2-VYY(K1)*DN1
      VXXX(1)=DN2*VZV-DN3*VYV
      VYYY(1)=DN3*VXV-DN1*VZV
      VZZZ(1)=DN1*VYV-DN2*VXV
C
      K2=KGB(2,LS,IK)
      N2=IE(M,K2)
      XXX(2)=XX(K2)
      YYY(2)=YY(K2)
      ZZZ(2)=ZZ(K2)
      VXV=VYY(K2)*DN3-VZZ(K2)*DN2
      VYV=VZZ(K2)*DN1-VXX(K2)*DN3
      VZV=VXX(K2)*DN2-VYY(K2)*DN1
      VXXX(2)=DN2*VZV-DN3*VYV
      VYYY(2)=DN3*VXV-DN1*VZV
      VZZZ(2)=DN1*VYV-DN2*VXV
C
      K3=KGB(3,LS,IK)
c      N3=IE(M,K3)
      XXX(3)=XX(K3)
      YYY(3)=YY(K3)
      ZZZ(3)=ZZ(K3)
      VXV=VYY(K3)*DN3-VZZ(K3)*DN2
      VYV=VZZ(K3)*DN1-VXX(K3)*DN3
      VZV=VXX(K3)*DN2-VYY(K3)*DN1
      VXXX(3)=DN2*VZV-DN3*VYV
      VYYY(3)=DN3*VXV-DN1*VZV
      VZZZ(3)=DN1*VYV-DN2*VXV
C
      IF(KNODE.EQ.4)THEN
        K4=KGB(4,LS,IK)
c        N4=IE(M,K4)
        XXX(4)=XX(K4)
        YYY(4)=YY(K4)
        ZZZ(4)=ZZ(K4)
        VXV=VYY(K4)*DN3-VZZ(K4)*DN2
        VYV=VZZ(K4)*DN1-VXX(K4)*DN3
        VZV=VXX(K4)*DN2-VYY(K4)*DN1
        VXXX(4)=DN2*VZV-DN3*VYV
        VYYY(4)=DN3*VXV-DN1*VZV
        VZZZ(4)=DN1*VYV-DN2*VXV
      ELSE
        K4=0
c        N4=0
      ENDIF
C
      VXV=VPY*DN3-VPZ*DN2
      VYV=VPZ*DN1-VPX*DN3
      VZV=VPX*DN2-VPY*DN1
      VPXX=DN2*VZV-DN3*VYV
      VPYY=DN3*VXV-DN1*VZV
      VPZZ=DN1*VYV-DN2*VXV
C
C ----- CHECK THE VELOCITY,ALONG THE SURFACE, OF POINT P
C
      IF(VPXX.EQ.0.0D0 .AND. VPYY.EQ.0.0D0 .AND. VPZZ.EQ.0.0D0)THEN
        SDT1=0.0D0
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
c       J1=N1
c       J2=N2
c       J3=N3
c       J4=N4
        RETURN
      ENDIF
C
C ----- START CHECKING
C
      IF(ID.EQ.1)THEN
        DO K=1,NODE
          IF(J1.EQ.IE(M,K))THEN
            JJ=K
            GOTO 110
          ENDIF
        ENDDO
  110   CONTINUE
      ENDIF
C
      IF(ID.EQ.2)THEN
        DO K=1,NODE
          IF(J1.EQ.IE(M,K))JJ1=K
          IF(J2.EQ.IE(M,K))JJ2=K
        ENDDO
      ENDIF
C
  120 CONTINUE
      DO 200 K=1,KNODE
        K1=K+1
        IF(K1.GT.KNODE)K1=1
C
C ----- CONSIDER DIFFERENT CASES
C
        KK=KGB(K,LS,IK)
        KK1=KGB(K1,LS,IK)
        IF(ID.EQ.1)THEN
          IF(JJ.EQ.KK .OR. JJ.EQ.KK1)GOTO 200
        ENDIF
        IF(ID.EQ.2)THEN
          IF(JJ1.EQ.KK .AND. JJ2.EQ.KK1)GOTO 200
          IF(JJ1.EQ.KK1 .AND. JJ2.EQ.KK)GOTO 200
        ENDIF
C
C ----- CHECK IF POINT P CAME FROM THIS SIDE
C
  140   CONTINUE
        CALL REPLAS(XXX(K),YYY(K),ZZZ(K),VXXX(K),VYYY(K),VZZZ(K),
     >              X1,Y1,Z1,VX1,VY1,VZ1)
        IF(IJUDGE1.EQ.1)THEN
          DOTX=(YP-Y1)*(VPZZ+VZ1)-(ZP-Z1)*(VPYY+VY1)
          DOTY=(ZP-Z1)*(VPXX+VX1)-(XP-X1)*(VPZZ+VZ1)
          DOTZ=(XP-X1)*(VPYY+VY1)-(YP-Y1)*(VPXX+VX1)
          DOTL=DSQRT(DOTX*DOTX+DOTY*DOTY+DOTZ*DOTZ)
          IF(DOTL.EQ.0.0D0)DOTL=1.0D0
          DOT1=(DOTX*DN1+DOTY*DN2+DOTZ*DN3)/DOTL
        ELSE
          DOTX=(YP-Y1)*VPZZ-(ZP-Z1)*VPYY
          DOTY=(ZP-Z1)*VPXX-(XP-X1)*VPZZ
          DOTZ=(XP-X1)*VPYY-(YP-Y1)*VPXX
          DOTL=DSQRT(DOTX*DOTX+DOTY*DOTY+DOTZ*DOTZ)
          IF(DOTL.EQ.0.0D0)DOTL=1.0D0
          DOT1=(DOTX*DN1+DOTY*DN2+DOTZ*DN3)/DOTL
        ENDIF
        IF(DOT1.LT.0.0D0)GOTO 200
C
  145   CONTINUE
        CALL REPLAS(XXX(K1),YYY(K1),ZZZ(K1),VXXX(K1),VYYY(K1),VZZZ(K1),
     >              X2,Y2,Z2,VX2,VY2,VZ2)
        IF(IJUDGE1.EQ.1)THEN
          DOTX=(YP-Y2)*(VPZZ+VZ2)-(ZP-Z2)*(VPYY+VY2)
          DOTY=(ZP-Z2)*(VPXX+VX2)-(XP-X2)*(VPZZ+VZ2)
          DOTZ=(XP-X2)*(VPYY+VY2)-(YP-Y2)*(VPXX+VX2)
          DOTL=DSQRT(DOTX*DOTX+DOTY*DOTY+DOTZ*DOTZ)
          IF(DOTL.EQ.0.0D0)DOTL=1.0D0
          DOT2=(DOTX*DN1+DOTY*DN2+DOTZ*DN3)/DOTL
        ELSE
          DOTX=(YP-Y2)*VPZZ-(ZP-Z2)*VPYY
          DOTY=(ZP-Z2)*VPXX-(XP-X2)*VPZZ
          DOTZ=(XP-X2)*VPYY-(YP-Y2)*VPXX
          DOTL=DSQRT(DOTX*DOTX+DOTY*DOTY+DOTZ*DOTZ)
          IF(DOTL.EQ.0.0D0)DOTL=1.0D0
          DOT2=(DOTX*DN1+DOTY*DN2+DOTZ*DN3)/DOTL
        ENDIF
        IF(DOT2.GT.0.0D0)GOTO 200
C
C ----- HERE WE GO, SOLVE FOR POINT Q
C NOTE: (XQ-X2):(X1-X2) = XSI:1
C
        N1=KK
        N2=KK1
        CALL REPLAS(VX1,VY1,VZ1,VX2,VY2,VZ2,VXXX(N1),VYYY(N1),VZZZ(N1),
     >              VXXX(N2),VYYY(N2),VZZZ(N2))
        CALL LOCQ2N
     I           (XX,YY,ZZ,VXXX,VYYY,VZZZ,XP,YP,ZP,VPXX,VPYY,VPZZ,
     I            N1,N2,SDT,IBF,IJUDGE1,EPSX,
     O            XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
        if(sdt1.eq.sdt)goto 200
c
        IF(SDT1.GT.0.0D0)THEN
          IF(DABS(XQ-XX(N1)).LE.EPSX .AND. DABS(YQ-YY(N1)).LE.EPSX .AND.
     >      DABS(ZQ-ZZ(N1)).LE.EPSX)THEN
            ID=1
            CALL REPLAS(XX(N1),YY(N1),ZZ(N1),VXX(N1),VYY(N1),VZZ(N1),
     >                   XQ,YQ,ZQ,VQX,VQY,VQZ)
            J1=IE(M,N1)
          ELSEIF(DABS(XQ-XX(N2)).LE.EPSX .AND. DABS(YQ-YY(N2)).LE.EPSX
     >           .AND. DABS(ZQ-ZZ(N2)).LE.EPSX)THEN
            ID=1
            CALL REPLAS(XX(N2),YY(N2),ZZ(N2),VXX(N2),VYY(N2),VZZ(N2),
     >                   XQ,YQ,ZQ,VQX,VQY,VQZ)
            J1=IE(M,N2)
          ELSE
            ID=2
            J1=IE(M,N1)
            J2=IE(M,N2)
          ENDIF
        ENDIF
        GOTO 300
  200 CONTINUE
C
      IF(IJUDGE1.EQ.1)THEN
        IJUDGE1=2
        GOTO 120
      ENDIF
      SDT1=SDT
      RETURN
C
  300 CONTINUE
C
      RETURN
      END
C
c
c
      SUBROUTINE BASE1
     I                (SS,TT,UU,
     O                 DL)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION DL(8)
C
      DO 10 I=1,8
        DL(I)=0.0D0
 10   CONTINUE
C
C ***** START CALCULATIONS
C
      SM=1.0D0-SS
      SP=1.0D0+SS
      TM=1.0D0-TT
      TP=1.0D0+TT
      UM=1.0D0-UU
      UP=1.0D0+UU
      DL(1)=0.125D0*SM*TM*UM
      DL(2)=0.125D0*SP*TM*UM
      DL(3)=0.125D0*SP*TP*UM
      DL(4)=0.125D0*SM*TP*UM
      DL(5)=0.125D0*SM*TM*UP
      DL(6)=0.125D0*SP*TM*UP
      DL(7)=0.125D0*SP*TP*UP
      DL(8)=0.125D0*SM*TP*UP
C
C ***** THIS IS THE END OF CALCULATIONS
C
      RETURN
      END
C
C
c
      SUBROUTINE BASE2D
c     I    (XX,YY,XQ,YQ,M,NODE,
     I    (XX,YY,XQ,YQ,NODE,
     O     DL,XSI,ETA)
C
C $$$$$ TO COMPUTE BASE FUNCTIONS BASED ON THE GIVEN ORIGINAL CARTESIAN
C       COORDINATES.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(4),YY(4),DL(4)
C
      IF(NODE.EQ.4)THEN
        CALL XSI2D
c     I    (XX,YY,XQ,YQ,M,
     I    (XX,YY,XQ,YQ,
     O     XSI,ETA)
        DL(1)=0.25*(1.-XSI)*(1.-ETA)
        DL(2)=0.25*(1.+XSI)*(1.-ETA)
        DL(3)=0.25*(1.+XSI)*(1.+ETA)
        DL(4)=0.25*(1.-XSI)*(1.+ETA)
      ENDIF
      IF(NODE.EQ.3)THEN
        A1=XX(2)*YY(3)-YY(2)*XX(3)
        A2=XX(3)*YY(1)-YY(3)*XX(1)
        A3=XX(1)*YY(2)-YY(1)*XX(2)
        B1=YY(2)-YY(3)
        B2=YY(3)-YY(1)
        B3=YY(1)-YY(2)
        C1=XX(3)-XX(2)
        C2=XX(1)-XX(3)
        C3=XX(2)-XX(1)
        A=A1+A2+A3
        DL(1)=(A1+B1*XQ+C1*YQ)/A
        DL(2)=(A2+B2*XQ+C2*YQ)/A
        DL(3)=(A3+B3*XQ+C3*YQ)/A
      ENDIF
      RETURN
      END
C
C
c
      SUBROUTINE SFDET
     I          (MAXNP,MAXEL,MXNODE,MXNPFG,NPFGS,ADPARM,ADPEPS,CMX,
     I           X,CS,XSFG,CSFG,MPLOCS,
     M           IE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP,3),CS(MAXNP)
      DIMENSION XSFG(MXNODE,3),CSFG(MXNODE)
      DIMENSION MPLOCS(MXNPFG),IE(MAXEL,11)
      DIMENSION XX(8),YY(8),ZZ(8),CC(8),DL(8),DNX(8)
C
C ***** DETERMINE SHARP FRONT
C
      DO 225 NP=1,NPFGS
        M=MPLOCS(NP)
        IF(IE(M,11).NE.0)GOTO 225
        XP=XSFG(NP,1)
        YP=XSFG(NP,2)
        ZP=XSFG(NP,3)
        CP=CSFG(NP)
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,I,I)
        DO I=1,NODE
          IEM=IE(M,I)
          XX(I)=X(IEM,1)
          YY(I)=X(IEM,2)
          ZZ(I)=X(IEM,3)
          CC(I)=CS(IEM)
        ENDDO
C
        CALL BASE
c     I           (XX,YY,ZZ,XP,YP,ZP,M,NODE,1,
     I           (XX,YY,ZZ,XP,YP,ZP,NODE,1,
     O            DL,DNX,DNX,DNX)
C
        CIJFE=0.0D0
        DO I=1,NODE
          CIJFE=CIJFE+CC(I)*DL(I)
        ENDDO
        CDIF=CP-CIJFE
        IF(DABS(CDIF).GT.ADPARM*CP .AND. DABS(CDIF).GT.ADPEPS*
     >     CMX)IE(M,11)=M
  225 CONTINUE
C
      RETURN
      END
C
C
c
        SUBROUTINE FGDET
     I      (IE,IBCHK,X,IB,CP,CS,VX,LRL,NLRL,ISE,NFGMB,XPFG,CPFG,MPLOC,
     I       IEW,IBW,LRLW,NLRLW,XW,CW,MWLOC,VXW,NEL,NEFG,NPFG,NNP,NTUBS,
     I       NCC,MICONF,IEPC,NXG,NYG,NZG,NPFGEP,NXW,NYW,NZW,IDETQ,TMAX,
     I       DELT,RAMADA,ADPARM,ADPEPS,EPSX,IZOOM,CMX,DTI,DTIFG,XEFG,
     I       CEFG,MPLOCE,IBE,ISB,DCOSB,NBDYB,IBDY,MAXBES,MXTUBS,
     I       MAXNP,MAXEL,MXNPFG,MXKGL,MXNPW,MXELW,MXKBD,MXNEP,
     i       MXNCC,DL468,nwnp,npw,iwtyp,wss,mxwnp,mxwpr,
     o       NPFGS,XSFG,CSFG,MPLOCS)
C
        IMPLICIT REAL*8(A-H,O-Z)
C
        DIMENSION IE(MAXEL,11),IBCHK(MAXEL),X(MAXNP,3),IBE(MAXEL)
        DIMENSION CS(MAXNP,MXNCC),DTI(MAXNP,MXNCC),DTIFG(MXNPFG,MXNCC)
        DIMENSION CP(MAXNP,MXNCC),IB(MAXNP),CMX(MXNCC),IBDY(MXTUBS)
        DIMENSION VX(MAXNP,3,MXNCC),LRL(MXKBD,MAXNP),NLRL(MAXNP)
        DIMENSION ISB(6,MAXBES),DCOSB(3,MAXBES),NBDYB(MAXNP)
        DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC)
        DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG,MXNCC)
        DIMENSION XEFG(MXNEP,3),CEFG(MXNEP,MXNCC)
        DIMENSION MPLOC(MXNPFG),MPLOCS(MXNPFG),MPLOCE(MXNPFG)
        DIMENSION ISE(MXKGL,8),NFGMB(MAXEL)
        DIMENSION XW(MXNPW,3),MWLOC(MXNPW),VXW(MXNPW,3),CW(MXNPW)
        DIMENSION LRLW(24,MXNPW,3),NLRLW(MXNPW,3)
        DIMENSION IEW(MXELW,8,3),IBW(MXNPW,3),DL468(8,MXNPW,3)
        dimension npw(mxwnp),iwtyp(mxwnp,mxncc),wss(mxwpr)
        DIMENSION XX(8),YY(8),ZZ(8),DL(8),IM6(6),ips(8)
C
        LOGICAL BINDEX
C
        NPFGS=0
        NPTET=(NXG+1)*(NXG+2)*(NXG+3)/6
c       NETET=NXG**3
        NPPRI=(NXG+1)*(NXG+2)*(NZG+1)/2
c       NEPRI=NXG**2*NZG
        NPHEX=(NXG+1)*(NYG+1)*(NZG+1)
c       NEHEX=NXG*NYG*NZG
C
        DO 300 M=1,NEL
C
          IF(IEPC.EQ.1)THEN
            IF(IE(M,11).EQ.0 .AND. IBCHK(M).EQ.0)GOTO 300
          ELSE
            IF(IE(M,11).EQ.0)GOTO 300
          ENDIF
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NODE,I,I)
          CALL WRKARY(IE,X,X,MAXEL,11,MAXNP,M,NODE,8,
     >                XX,YY,ZZ,XX,YY,ZZ)
C
          do i=1,8
            ips(i)=0
          enddo
          do 10 i=1,nwnp
            np=npw(i)
            do ii=1,node
              iem=ie(m,ii)
              if(iem.eq.np)then
                ips(ii)=ii
                goto 10
              endif
            enddo
  10      continue
c
          NPS=NPFGS
          NEPP=NPFGEP
C
          DO I=6,1,-1
            IM6(I)=IBE(M)/10**(I-1)
            DO J=6,I+1,-1
              IM6(I)=IM6(I)-IM6(J)*10**(J-I)
            ENDDO
          ENDDO
C
C ----- REFINE THIS SHARP-FRONT ELEMENT ===> XSFG(N,1..3)
C ----- FOR THE CASES OF HEXAHEDRAL ELEMENTS
C
          IF(NODE.EQ.8)THEN
            DO 240 K=1,NZG+1
              ZI=2.0*DBLE(K-1)/DBLE(NZG)-1.0
              DO 230 I=1,NYG+1
                YI=2.0*DBLE(I-1)/DBLE(NYG)-1.0
                DO 220 J=1,NXG+1
                  XI=2.0*DBLE(J-1)/DBLE(NXG)-1.0
                  N=J+(I-1)*(NXG+1)+(K-1)*(NXG+1)*(NYG+1)+NPFGS
                  CALL WARMSG(N,MXNPFG,'FGDET ','MXNPFG',1)
C
C ----- DETERMINE WORKING ARRAYS XW, YW, ZW, VXW, VYW, AND VZW
C
                  CALL BASE1
     I                      (XI,YI,ZI,
     O                       DL)
                  CALL VALBDL(XX,YY,ZZ,DL,8,8,XSFG(N,1),XSFG(N,2),
     >                        XSFG(N,3))
                  MPLOCS(N)=M
C
C ----- CHECK NEW EPCOF POINTS DUE TO INFLOW BOUNDARIES
C
                  NEP=NPFGEP
                  IF(K.EQ.NZG .OR. K.EQ.NZG+1)ISF=6
                  IF(K.EQ.1 .OR. K.EQ.2)ISF=5
                  IF(I.EQ.NYG .OR. I.EQ.NYG+1)ISF=4
                  IF(J.EQ.NXG .OR. J.EQ.NXG+1)ISF=3
                  IF(I.EQ.1 .OR. I.EQ.2)ISF=2
                  IF(J.EQ.1 .OR. J.EQ.2)ISF=1
                  BINDEX=.FALSE.
                  IJK=0
  215             CONTINUE
                    IJK=IJK+1
                    IF(IM6(IJK).EQ.ISF)BINDEX=.TRUE.
                  IF(.NOT. BINDEX .AND. IJK.LT.6)GOTO 215
                  if(bindex)goto 217
C
                  if(nwnp.ne.0)then
                    do 216 ijk=1,8
                      if(ips(ijk).eq.0)goto 216
                      if(ips(ijk).eq.1)then
                        if(j.eq.1 .or. j.eq.2)then
                          if(i.eq.1 .or. i.eq.2)then
                            if(k.eq.1 .or. k.eq.2)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.2)then
                        if(j.eq.nxg .or. j.eq.nxg+1)then
                          if(i.eq.1 .or. i.eq.2)then
                            if(k.eq.1 .or. k.eq.2)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.3)then
                        if(j.eq.nxg .or. j.eq.nxg+1)then
                          if(i.eq.nyg .or. i.eq.nyg+1)then
                            if(k.eq.1 .or. k.eq.2)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.4)then
                        if(j.eq.1 .or. j.eq.2)then
                          if(i.eq.nyg .or. i.eq.nyg+1)then
                            if(k.eq.1 .or. k.eq.2)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.5)then
                        if(j.eq.1 .or. j.eq.2)then
                          if(i.eq.1 .or. i.eq.2)then
                            if(k.eq.nzg .or. k.eq.nzg+1)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.6)then
                        if(j.eq.nxg .or. j.eq.nxg+1)then
                          if(i.eq.1 .or. i.eq.2)then
                            if(k.eq.nzg .or. k.eq.nzg+1)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.7)then
                        if(j.eq.nxg .or. j.eq.nxg+1)then
                          if(i.eq.nyg .or. i.eq.nyg+1)then
                            if(k.eq.nzg .or. k.eq.nzg+1)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      elseif(ips(ijk).eq.8)then
                        if(j.eq.1 .or. j.eq.2)then
                          if(i.eq.nyg .or. i.eq.nyg+1)then
                            if(k.eq.nzg .or. k.eq.nzg+1)then
                              bindex=.true.
                            endif
                          endif
                        endif
                      endif
                      if(bindex)goto 217
  216               continue
                  endif
  217             continue
C
                  IF(BINDEX)NPFGEP=NPFGEP+1
C
                  IF(NPFGEP.GT.NEP .AND. IEPC.NE.0)THEN
                    CALL WARMSG(NPFGEP,MXNEP,'FGDET ','MXNEP ',2)
                    XEFG(NPFGEP,1)=XSFG(N,1)
                    XEFG(NPFGEP,2)=XSFG(N,2)
                    XEFG(NPFGEP,3)=XSFG(N,3)
                    MPLOCE(NPFGEP)=M
                    MPLOC(NPFGEP)=N
                  ENDIF
C
  220           CONTINUE
  230         CONTINUE
  240       CONTINUE
            NPFGS=NPFGS+NPHEX
            NPWRK=NPHEX
C
C ----- DETERMINE WORKING ARRAYS XW, YW, ZW, VXW, VYW, AND VZW
C
          ELSEIF(NODE.EQ.6)THEN
c           NWRK=NPPRI
            NN=0
            DO 248 K=1,NZG+1
              ZI=2.0D0*DBLE(K-1)/DBLE(NZG)-1.0D0
              DO 247 I=NXG+1,1,-1
                DL1=DBLE(I-1)/DBLE(NXG)
                DO 246 J=NXG+1,1,-1
                  DL2=DBLE(J-1)/DBLE(NXG)
                  IF(DL1+DL2.GT.1.0D0 .AND. DABS(DL1+DL2-1.0D0)
     >               .GT.EPSX)GOTO 246
                  DO 245 L=NXG+1,1,-1
                    DL3=DBLE(L-1)/DBLE(NXG)
                    IF(DABS(DL1+DL2+DL3-1.0D0).GT.EPSX)GOTO 245
                    NN=NN+1
                    N=NN+NPFGS
                    CALL WARMSG(N,MXNPFG,'FGDET ','MXNPFG',3)
C
                    DL(1)=DL1*(1.0D0-ZI)*0.5D0
                    DL(2)=DL2*(1.0D0-ZI)*0.5D0
                    DL(3)=DL3*(1.0D0-ZI)*0.5D0
                    DL(4)=DL1*(1.0D0+ZI)*0.5D0
                    DL(5)=DL2*(1.0D0+ZI)*0.5D0
                    DL(6)=DL3*(1.0D0+ZI)*0.5D0
                    CALL VALBDL(XX,YY,ZZ,DL,8,6,XSFG(N,1),XSFG(N,2),
     >                          XSFG(N,3))
                    MPLOCS(N)=M
C
C ----- CHECK NEW EPCOF POINTS DUE TO INFLOW BOUNDARIES
C
                    NEP=NPFGEP
                    IF(L.EQ.1 .OR. L.EQ.2)ISF=1
                    IF(J.EQ.1 .OR. J.EQ.2)ISF=3
                    IF(I.EQ.1 .OR. I.EQ.2)ISF=2
                    IF(K.EQ.NZG .OR. K.EQ.NZG+1)ISF=5
                    IF(K.EQ.1 .OR. K.EQ.2)ISF=4
                    BINDEX=.FALSE.
                    IJK=0
  241               CONTINUE
                      IJK=IJK+1
                      IF(IM6(IJK).EQ.ISF)BINDEX=.TRUE.
                    IF(.NOT. BINDEX .AND. IJK.LT.6)GOTO 241
c
                    if(bindex)goto 243
C
                    if(nwnp.ne.0)then
                      do 242 ijk=1,6
                        if(ips(ijk).eq.0)goto 242
                        if(ips(ijk).eq.1)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(j.eq.1 .or. j.eq.2)then
                              if(k.eq.1 .or. k.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.2)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(i.eq.1 .or. i.eq.2)then
                              if(k.eq.1 .or. k.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.3)then
                          if(j.eq.1 .or. j.eq.2)then
                            if(i.eq.1 .or. i.eq.2)then
                              if(k.eq.1 .or. k.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.4)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(j.eq.1 .or. j.eq.2)then
                              if(k.eq.nzg .or. k.eq.nzg+1)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.5)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(i.eq.1 .or. i.eq.2)then
                              if(k.eq.nzg .or. k.eq.nzg+1)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.6)then
                          if(j.eq.1 .or. j.eq.2)then
                            if(i.eq.1 .or. i.eq.2)then
                              if(k.eq.nzg .or. k.eq.nzg+1)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        endif
                        if(bindex)goto 243
  242                 continue
                    endif
  243               continue
C
                    IF(BINDEX)NPFGEP=NPFGEP+1
C
                    IF(NPFGEP.GT.NEP .AND. IEPC.NE.0)THEN
                      CALL WARMSG(NPFGEP,MXNEP,'FGDET ','MXNEP ',4)
                      XEFG(NPFGEP,1)=XSFG(N,1)
                      XEFG(NPFGEP,2)=XSFG(N,2)
                      XEFG(NPFGEP,3)=XSFG(N,3)
                      MPLOCE(NPFGEP)=M
                      MPLOC(NPFGEP)=N
                    ENDIF
C
  245             CONTINUE
  246           CONTINUE
  247         CONTINUE
  248       CONTINUE
            NPFGS=NPFGS+NPPRI
            NPWRK=NPPRI
C
C----- FOR THE CASES OF TETRAHEDRAL ELEMENTS
C
          ELSE
            NN=0
            DO 290 K=NXG+1,1,-1
              DL1=DBLE(K-1)/DBLE(NXG)
              DO 280 I=NXG+1,1,-1
                DL2=DBLE(I-1)/DBLE(NXG)
                IF(DL1+DL2.GT.1.0D0 .AND. DABS(DL1+DL2-1.0D0).GT.EPSX)
     >             GOTO 280
                DO 270 J=NXG+1,1,-1
                  DL3=DBLE(J-1)/DBLE(NXG)
                  IF(DL1+DL2+DL3.GT.1.0D0 .AND. DABS(DL1+DL2+DL3-1.0D0)
     >               .GT.EPSX)GOTO 270
                  DO 260 L=NXG+1,1,-1
                    DL4=DBLE(L-1)/DBLE(NXG)
                    IF(DABS(DL1+DL2+DL3+DL4-1.0D0).GT.EPSX)GOTO 260
                    NN=NN+1
                    N=NN+NPFGS
                    CALL WARMSG(N,MXNPFG,'FGDET ','MXNPFG',5)
C
                    DL(1)=DL1
                    DL(2)=DL2
                    DL(3)=DL3
                    DL(4)=DL4
                    CALL VALBDL(XX,YY,ZZ,DL,8,4,XSFG(N,1),XSFG(N,2),
     >                          XSFG(N,3))
                    MPLOCS(N)=M
C
C ----- CHECK NEW EPCOF POINTS DUE TO INFLOW BOUNDARIES
C
                    NEP=NPFGEP
                    IF(L.EQ.1 .OR. L.EQ.2)ISF=4
                    IF(J.EQ.1 .OR. J.EQ.2)ISF=3
                    IF(I.EQ.1 .OR. I.EQ.2)ISF=2
                    IF(K.EQ.1 .OR. K.EQ.2)ISF=1
                    BINDEX=.FALSE.
                    IJK=0
  255               CONTINUE
                      IJK=IJK+1
                      IF(IM6(IJK).EQ.ISF)BINDEX=.TRUE.
                    IF(.NOT. BINDEX .AND. IJK.LT.6)GOTO 255
c
                    if(bindex)goto 257
C
                    if(nwnp.ne.0)then
                      do 256 ijk=1,4
                        if(ips(ijk).eq.0)goto 256
                        if(ips(ijk).eq.1)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(j.eq.1 .or. j.eq.2)then
                              if(i.eq.1 .or. i.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.2)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(j.eq.1 .or. j.eq.2)then
                              if(k.eq.1 .or. k.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.3)then
                          if(l.eq.1 .or. l.eq.2)then
                            if(i.eq.1 .or. i.eq.2)then
                              if(k.eq.1 .or. k.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        elseif(ips(ijk).eq.4)then
                          if(i.eq.1 .or. i.eq.2)then
                            if(j.eq.1 .or. j.eq.2)then
                              if(k.eq.1 .or. k.eq.2)then
                                bindex=.true.
                              endif
                            endif
                          endif
                        endif
                        if(bindex)goto 257
  256                 continue
                    endif
  257               continue
C
                    IF(BINDEX)NPFGEP=NPFGEP+1
C
                    IF(NPFGEP.GT.NEP .AND. IEPC.NE.0)THEN
                      CALL WARMSG(NPFGEP,MXNEP,'FGDET ','MXNEP ',6)
                      IF(NPFGEP.GT.MXNEP)THEN
                        WRITE(16,*)'NPFGEP=',NPFGEP,'--> INCREASE MXNEP'
                        STOP
                      ENDIF
                      XEFG(NPFGEP,1)=XSFG(N,1)
                      XEFG(NPFGEP,2)=XSFG(N,2)
                      XEFG(NPFGEP,3)=XSFG(N,3)
                      MPLOCE(NPFGEP)=M
                      MPLOC(NPFGEP)=N
                    ENDIF
C
  260             CONTINUE
  270           CONTINUE
  280         CONTINUE
  290       CONTINUE
            NPFGS=NPFGS+NPTET
            NPWRK=NPTET
          ENDIF
C
          NINIT=NPS+1
          IBF=1
          DO K=1,NCC
            DO NP=NINIT,NPFGS
              DTIFG(NP,K)=1.0D0/DELT
              IF(MICONF.EQ.1 .AND. K.LE.3)THEN
                IDTI=0
              ELSE
                IDTI=2
              ENDIF
            ENDDO
            CALL HPTRAK
     I       (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,MXTUBS,
     I        MXNEP,MXNEP,MXNPFG,
     I        NNP,NEL,NEFG,NTUBS,NXW,NYW,NZW, NPFG,NINIT,
     I        IBF,IDETQ, DELT,TMAX,RAMADA, IDTI,IZOOM,EPSX,
     I        CP(1,K),X,VX(1,1,K),CPFG(1,K),XPFG,MPLOC,ISE,
     I        NFGMB,IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,
     I        ISB,DCOSB,NBDYB,IBDY,DL468,
     I        nwnp,npw,iwtyp(1,k),wss,mxwnp,mxwpr,
     O        CS(1,K),DTI(1,K),DTIFG(1,K),NPFGS,XSFG,CSFG(1,K),MPLOCS,
     O        NPFGEP,XEFG,CEFG(1,K),MPLOCE,
     M        XW,VXW)
          ENDDO
C
C ----- CHECK IF ELEMENT M IS REALLY A SF ELEMENT DUE TO INFLOW
C       BOUNDARIES
C
c         IF(IEPC.EQ.1)THEN
c           IF(IBCHK(M).NE.0 .AND. IE(M,11).EQ.0 .AND.
c    >         NPFGEP.GT.NEPP)THEN
            IF(IBCHK(M).NE.0 .AND. IE(M,11).EQ.0) then
              DO I=1,NPWRK
                XW(I,1)=XSFG(NPS+I,1)
                XW(I,2)=XSFG(NPS+I,2)
                XW(I,3)=XSFG(NPS+I,3)
                MWLOC(I)=MPLOCS(NPS+I)
              ENDDO
              DO K=1,NCC
                DO I=1,NPWRK
                  CW(I)=CSFG(NPS+I,K)
                ENDDO
                CALL SFDET
     I             (MAXNP,MAXEL,MXNPW,MXNPW,NPWRK,ADPARM,ADPEPS,CMX(K),
     I              X,CS(1,K),XW,CW,MWLOC,
     M              IE)
                IF(IE(M,11).NE.0)THEN
                  IF(NPFGEP.GT.NEPP)THEN
                    if(iepc.ne.0)then
                      DO N=NEPP+1,NPFGEP
                        N1=MPLOC(N)
                        DO IK=1,NCC
                          CEFG(N,IK)=CSFG(N1,IK)
                        ENDDO
                      ENDDO
                    else
                      npfgep=nepp
                    endif
                  ENDIF
                ELSE
                  NPFGS=NPS
                  NPFGEP=NEPP
                ENDIF
              ENDDO
            ELSE
              NPFGEP=NEPP
            ENDIF
c         ELSE
c           NPFGEP=NEPP
c         ENDIF
  300   CONTINUE
      RETURN
      END
C
C
c
      SUBROUTINE ISEHIL
     I  (MAXEL,MXKGL,MXNPFG,MXELW,MXNPW,MXEPW,MXNCC,NEL,NCC,
     I   NPFGW,EPSX,NFGM,NFGMB,MPLOCS,IE,XWFG,CWFG,MPLOC,
     I   XSFG,CSFG, NXG,NYG,NZG,
     M   NEPWN,NEPW,
     O   NPFGS,NEFGS,NLGELM,ISE,MAXFGW,MINFGW,CMAXFG,CMINFG,MPLOCE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION NFGMB(MAXEL),IE(MAXEL,11)
      DIMENSION XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC),MPLOC(MXNPFG)
      DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC),ISE(MXKGL,8)
      DIMENSION MAXFGW(MXELW,MXNCC),MINFGW(MXELW,MXNCC)
      DIMENSION CMAXFG(MXELW,MXNCC),CMINFG(MXELW,MXNCC)
      DIMENSION NFGM(MAXEL),MPLOCS(MXNPFG),MPLOCE(MXNPFG)
      DIMENSION NEPWN(MXELW),NEPW(MXELW,MXEPW),CQ(7),KIT(4,6,3)
C
      DATA KIT/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
C ----- DETERMINE ISE( , )
C
        NP=0
        NE=0
        NLGELM=0
        NELWK6=(NXG+1)*(NXG+2)/2
        NELWK8=(NXG+1)*(NYG+1)
        DO 400 M=1,NEL
C
C NOTE: THE FOLLOWING ARRAY NFGMB IS A WORKING ARRAY
C
          NB=NFGMB(M)
C
C NOTE: THE FOLLOWING ARRAY NFGMB IS NOT A WORKING ARRAY
C
          NFGMB(M)=NE
          IF(IE(M,11).EQ.0)GOTO 400
          NLGELM=NLGELM+1
          CALL ELENOD(IE(M,5),IE(M,7),NODE,I,I)
          IF(NODE.EQ.4)THEN
            NPWRK=(NXG+1)*(NXG+2)*(NXG+3)/6
            NEWRK=NXG**3
          ELSEIF(NODE.EQ.6)THEN
            NPWRK=(NXG+1)*(NXG+2)*(NZG+1)/2
            NEWRK=NXG**2*NZG
          ELSE
            NPWRK=(NXG+1)*(NYG+1)*(NZG+1)
            NEWRK=NXG*NYG*NZG
          ENDIF
          CALL WARMSG(NEWRK,MXELW,'ISEHIL','MXELW ',1)
          CALL WARMSG(NPWRK,MXNPW,'ISEHIL','MXNPW ',2)
          DO I=1,NEWRK
            NEPWN(I)=0
            DO J=1,MXEPW
              NEPW(I,J)=0
            ENDDO
          ENDDO
C
C ----- FOR THE CASES OF HEXAHEDRAL ELEMENTS
C
          IF(NODE.EQ.8)THEN
            DO 87 K=1,NZG
              DO 86 I=1,NYG
                DO 85 J=1,NXG
                  MW=(I-1)*NXG+J+(K-1)*NXG*NYG+NE
                  ISE(MW,1)=(I-1)*(NXG+1)+J+(K-1)*NELWK8+NP
                  ISE(MW,2)=(I-1)*(NXG+1)+J+1+(K-1)*NELWK8+NP
                  ISE(MW,3)=I*(NXG+1)+J+1+(K-1)*NELWK8+NP
                  ISE(MW,4)=I*(NXG+1)+J+(K-1)*NELWK8+NP
                  ISE(MW,5)=ISE(MW,1)+NELWK8
                  ISE(MW,6)=ISE(MW,2)+NELWK8
                  ISE(MW,7)=ISE(MW,3)+NELWK8
                  ISE(MW,8)=ISE(MW,4)+NELWK8
  85            CONTINUE
  86          CONTINUE
  87        CONTINUE
C
C ----- FOR THE CASES OF PRIZM ELEMENTS
C
          ELSEIF(NODE.EQ.6)THEN
            DO 90 K=1,NZG
              NC=0
              NI=0
              NK=(K-1)*NXG**2
              NKI=(K-1)*NELWK6
              DO 89 I=1,NXG
                NT=2*I-1
                NII=NI
                DO 88 J=NC+1,NC+NT,2
                  MW=J+NE+NK
                  NI=NI+1
                  NJ=NI+NKI
                  ISE(MW,1)=NJ+NP
                  ISE(MW,2)=NJ+I+NP
                  ISE(MW,3)=NJ+I+1+NP
                  ISE(MW,4)=ISE(MW,1)+NELWK6
                  ISE(MW,5)=ISE(MW,2)+NELWK6
                  ISE(MW,6)=ISE(MW,3)+NELWK6
                  ISE(MW,7)=0
                  ISE(MW,8)=0
  88            CONTINUE
                IF(I.GT.1)THEN
                  DO J=NC+2,NC+NT-1,2
                    NII=NII+1
                    NIJ=NII+NKI
                    MW=J+NE+NK
                    ISE(MW,1)=NIJ+NP
                    ISE(MW,2)=NIJ+I+1+NP
                    ISE(MW,3)=NIJ+1+NP
                    ISE(MW,4)=ISE(MW,1)+NELWK6
                    ISE(MW,5)=ISE(MW,2)+NELWK6
                    ISE(MW,6)=ISE(MW,3)+NELWK6
                    ISE(MW,7)=0
                    ISE(MW,8)=0
                  ENDDO
                ENDIF
                NC=I**2
  89          CONTINUE
  90        CONTINUE
C
C----- FOR THE CASES OF TETRAHEDRAL ELEMENTS
C
          ELSE
            J2=1
            J1=0
            DO 92 L=1,NXG
              JUMP1=L
              JUMP2=(L-1)*L/2
              J1=J1+JUMP1
              J2=J2+JUMP2
c             NCOUNT=L**3-(L-1)**3
              M1=(L-1)**3+1+NE
              ISE(M1,1)=J2+NP
              ISE(M1,2)=J2+J1+NP
              ISE(M1,3)=J2+J1+1+NP
              ISE(M1,4)=J2+J1+2+NP
C
              IF(L.EQ.1)GOTO 92
              J3=J2
              DO 91 L1=2,L
                M2=(L1-2)*(L1-1)*3+2+(L-1)**3+NE
                J3=J3+L1-1
                ISE(M2,1)=J3+NP
                ISE(M2,2)=J3+J1+NP
                ISE(M2,3)=J3+J1+L1+NP
                ISE(M2,4)=J3+J1+L1+1+NP
                ISE(M2+1,1)=J3+NP
                ISE(M2+1,2)=J3+J1+1+NP
                ISE(M2+1,3)=J3-L1+1+NP
                ISE(M2+1,4)=J3+J1+NP
                ISE(M2+2,1)=J3+NP
                ISE(M2+2,2)=J3+J1+1+NP
                ISE(M2+2,3)=J3+J1+NP
                ISE(M2+2,4)=J3+J1+L1+1+NP
                ISE(M2+3,1)=J3+NP
                ISE(M2+3,2)=J3+J1+1+NP
                ISE(M2+3,3)=J3+J1+L1+1+NP
                ISE(M2+3,4)=J3+1+NP
                ISE(M2+4,1)=J3+NP
                ISE(M2+4,2)=J3+J1+1+NP
                ISE(M2+4,3)=J3+1+NP
                ISE(M2+4,4)=J3-L1+1+NP
                ISE(M2+5,1)=J3+1+NP
                ISE(M2+5,2)=J3+J1+1+NP
                ISE(M2+5,3)=J3+J1+L1+1+NP
                ISE(M2+5,4)=J3+J1+L1+2+NP
                IF(L.GE.3)THEN
                  M5=M2+5
                  J4=J3
                  DO L2=1,L-2
                    J4=J4+1
                    M5=M5+1
                    ISE(M5,1)=J4+NP
                    ISE(M5,2)=J4-L1+1+NP
                    ISE(M5,3)=J4-L1+NP
                    ISE(M5,4)=J4+J1+NP
                    M5=M5+1
                    ISE(M5,1)=J4+NP
                    ISE(M5,2)=J4+J1+1+NP
                    ISE(M5,3)=J4-L1+1+MP
                    ISE(M5,4)=J4+J1+NP
                    M5=M5+1
                    ISE(M5,1)=J4+NP
                    ISE(M5,2)=J4+J1+1+NP
                    ISE(M5,3)=J4+J1+NP
                    ISE(M5,4)=J4+J1+L1+1+NP
                    M5=M5+1
                    ISE(M5,1)=J4+NP
                    ISE(M5,2)=J4+J1+1+NP
                    ISE(M5,3)=J4+J1+L1+1+NP
                    ISE(M5,4)=J4+1+NP
                    M5=M5+1
                    ISE(M5,1)=J4+NP
                    ISE(M5,2)=J4+J1+1+NP
                    ISE(M5,3)=J4+1+NP
                    ISE(M5,4)=J4-L1+1+NP
                    M5=M5+1
                    ISE(M5,1)=J4+1+NP
                    ISE(M5,2)=J4+J1+1+NP
                    ISE(M5,3)=J4+J1+L1+1+NP
                    ISE(M5,4)=J4+J1+L1+2+NP
                  ENDDO
                ENDIF
  91          CONTINUE
  92        CONTINUE
            DO 93 I1=1,NXG**3
              II1=I1+NE
              DO I2=5,8
                ISE(II1,I2)=0
              ENDDO
  93        CONTINUE
          ENDIF
C
          NP=NP+NPWRK
          CALL WARMSG(NP,MXNPFG,'ISEHIL','MXNPFG',3)
          NE=NE+NEWRK
          CALL WARMSG(NE,MXKGL,'ISEHIL','MXKGL ',4)
C
C ***** START CHECKING THE HIGH AND LOW OF SUBELEMENTS
C
          DO I=1,MXELW
            DO K=1,NCC
              MAXFGW(I,K)=0
              MINFGW(I,K)=0
              CMAXFG(I,K)=-1.0D38
              CMINFG(I,K)=1.0D38
            ENDDO
          ENDDO
c         DO I=1,NODE
c           IEM=IE(M,I)
c           XX(I)=X(IEM,1)
c           YY(I)=X(IEM,2)
c           ZZ(I)=X(IEM,3)
c           DO K=1,NCC
c             CC(I,K)=CS(IEM,K)
c           ENDDO
c         ENDDO
C
C CHECK EVERY FORWARD-TRACKED POINT, WHICK FELL IN THIS ELEMENT, ONE
C BY ONE
C
          DO 398 N=1,NFGM(M)
            NM=N+NB
C
C NOTE: ARRAY MPLOC IN THE NEXT LINE IS THE WORKING ARRAY GENERATED JUST
C       NOW
C
            NN=MPLOC(NM)
            XP=XWFG(NN,1)
            YP=XWFG(NN,2)
            ZP=XWFG(NN,3)
C
            M1=NFGMB(M)+1
            M2=NE
            CALL KGLOC
     I         (XSFG,CSFG,ISE,KIT,MXNPFG,MXKGL,8,MXNCC,
     >          XP,YP,ZP,M1,M2,1,1,EPSX, MP1,NODE,CQ,ICODE)
            XWFG(NN,1)=XP
            XWFG(NN,2)=YP
            XWFG(NN,3)=ZP
            MP=MP1-NFGMB(M)
C
            IF(NN.GT.NPFGW)THEN
              DO I=1,NODE
                IEM=ISE(MP1,I)
                IF(DABS(XP-XSFG(IEM,1)).LE.EPSX .AND.
     >             DABS(YP-XSFG(IEM,2)).LE.EPSX .AND.
     >             DABS(ZP-XSFG(IEM,3)).LE.EPSX)THEN
                  MPLOCE(IEM)=IEM
                  DO K=1,NCC
                    CSFG(IEM,K)=CWFG(NN,K)
                  ENDDO
                  GOTO 398
                ENDIF
              ENDDO
              NEPWN(MP)=NEPWN(MP)+1
              CALL WARMSG(NEPWN(MP),MXEPW,'ISEHIL','MXEPW ',5)
              NBN=NEPWN(MP)
              NEPW(MP,NBN)=NN
              GOTO 398
            ENDIF
C
            DO K=1,NCC
              CCP=CWFG(NN,K)
              MXK=0
              MNK=0
              DO I=1,NODE
                IEM=ISE(MP1,I)
                IF(DABS(XP-XSFG(IEM,1)).LE.EPSX .AND.
     >             DABS(YP-XSFG(IEM,2)).LE.EPSX .AND.
     >             DABS(YP-XSFG(IEM,3)).LE.EPSX)GOTO 398
                IF(CCP.GE.CSFG(IEM,K))MXK=MXK+1
                IF(CCP.LE.CSFG(IEM,K))MNK=MNK+1
              ENDDO
              IF(MXK.EQ.NODE)THEN
                IF(CCP.GT.CMAXFG(MP,K))THEN
                  DO KK=1,NCC
                    MAXFGW(MP,KK)=NN
                    CMAXFG(MP,KK)=CWFG(NN,KK)
                  ENDDO
                  GOTO 398
                ENDIF
              ENDIF
              IF(MNK.EQ.NODE)THEN
                IF(CCP.LT.CMINFG(MP,K))THEN
                  DO KK=1,NCC
                    MINFGW(MP,KK)=NN
                    CMINFG(MP,KK)=CWFG(NN,KK)
                  ENDDO
                  GOTO 398
                ENDIF
              ENDIF
            ENDDO
  398     CONTINUE
C
C ----- CHECK IF THOSE HIGHS AND LOWS ARE REALLY NEEDED
C
          DO 399 I=1,NEWRK
            IDO=1
            DO K=1,NCC
              IF(MAXFGW(I,K).NE.0 .AND. CMAXFG(I,K).GT.1.0D-20)IDO=0
            ENDDO
            IF(IDO.EQ.1)THEN
              DO K=1,NCC
                MAXFGW(I,K)=0
                CMAXFG(I,K)=-1.0D38
              ENDDO
            ENDIF
            IDO=1
            DO K=1,NCC
              IF(MINFGW(I,K).NE.0)THEN
                DO J=1,NODE
                  IEM=ISE(NFGMB(M)+I,J)
                  IF(CSFG(IEM,K).GT.1.0D-20)IDO=0
                ENDDO
              ENDIF
            ENDDO
            IF(IDO.EQ.1)THEN
              DO K=1,NCC
                MINFGW(I,K)=0
                CMINFG(I,K)=1.0D38
              ENDDO
            ENDIF
  399     CONTINUE
C
          CALL TRIANG
     I        (NEWRK,MAXFGW,MINFGW,XWFG,CWFG,NFGMB,NEPWN,NEPW,
     I         NODE,MXELW,MXNPFG,MAXEL,MXKGL,MXEPW,MXNCC,M,NCC,EPSX,
     M         XSFG,CSFG,ISE,NE,NPFGS,MPLOCS,MPLOCE)
C
          CALL WARMSG(NE,MXKGL,'ISEHIL','MXKGL ',6)
  400   CONTINUE
        NEFGS=NE
        RETURN
        END
c
c
c
      SUBROUTINE TRIANG
     I    (NEWRK,MAXFGW,MINFGW,XWFG,CWFG,NFGMB,NEPWN,NEPW,
     I     NODE,MXELW,MXNPFG,MAXEL,MXKGL,MXEPW,MXNCC,M,NCC,EPSX,
     M     XSFG,CSFG,ISE,NE,NPFGS,MPLOCS,MPLOCE)
C      11/ 5/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (MXNPK=50,MXNEK=200)
      DIMENSION XWRK(MXNPK,3),NWRK(MXNPK),IZE(MXNEK,8),IAP(MXNPK)
      DIMENSION MAXFGW(MXELW,MXNCC),MINFGW(MXELW,MXNCC),NFGMB(MAXEL)
      DIMENSION XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC)
      DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC)
      DIMENSION ISE(MXKGL,8),MPLOCS(MXNPFG),MPLOCE(MXNPFG)
      DIMENSION NEPWN(MXELW),NEPW(MXELW,MXEPW),XV(MXNEK,3)
      DIMENSION XX(8),YY(8),ZZ(8),IBP(MXNEK,4),N(3),KID(4,6,3)
C
      DATA KID /1,2,3,4, 3,4,5,6, 2,3,4,5, 0,0,0,0, 0,0,0,0, 0,0,0,0,
     >          4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0,
     >          1,2,4,5, 5,6,4,8, 5,2,4,6, 2,3,4,6, 6,3,4,8, 6,3,8,7/
C
      DO 400 MP=1,NEWRK
        CALL ELENOD
     I      (ISE(NFGMB(M)+MP,5),ISE(NFGMB(M)+MP,7),
     O       NODE,I,I)
C
C ----- CHECK EPCOF POINTS AND THE HIGH AND LOW OF EVERY SUBELEMENT
C
        NK=0
        IF(NEPWN(MP).NE.0)THEN
          DO 10 I=1,NEPWN(MP)
            J=NEPW(MP,I)
            DO II=1,NK
              IF(DABS(XWRK(II,1)-XWFG(J,1)).LE.EPSX .AND.
     >           DABS(XWRK(II,2)-XWFG(J,2)).LE.EPSX .AND.
     >           DABS(XWRK(II,3)-XWFG(J,3)).LE.EPSX)GOTO 10
            ENDDO
            NK=NK+1
            CALL WARMSG(NK,MXNPK,'TRIANG','MXNPK ',1)
            NPFGS=NPFGS+1
            CALL WARMSG(NPFGS,MXNPFG,'TRIANG','MXNPFG',2)
            XSFG(NPFGS,1)=XWFG(J,1)
            XSFG(NPFGS,2)=XWFG(J,2)
            XSFG(NPFGS,3)=XWFG(J,3)
            DO K=1,NCC
              CSFG(NPFGS,K)=CWFG(J,K)
            ENDDO
            MPLOCS(NPFGS)=M
            MPLOCE(NPFGS)=NPFGS
            XWRK(NK,1)=XSFG(NPFGS,1)
            XWRK(NK,2)=XSFG(NPFGS,2)
            XWRK(NK,3)=XSFG(NPFGS,3)
            NWRK(NK)=NPFGS
  10      CONTINUE
        ENDIF
C
        IDO=1
        K=0
  15    CONTINUE
          K=K+1
          IF(MAXFGW(MP,K).NE.0)IDO=0
        IF(IDO.EQ.1 .AND. K.LT.NCC)GOTO 15
        IF(IDO.EQ.0)THEN
          DO K=1,NCC
            J=MAXFGW(MP,K)
            DO I=1,NK
              IF(DABS(XWRK(I,1)-XWFG(J,1)).LE.EPSX .AND.
     >           DABS(XWRK(I,2)-XWFG(J,2)).LE.EPSX .AND.
     >           DABS(XWRK(I,3)-XWFG(J,3)).LE.EPSX)GOTO 20
            ENDDO
            NK=NK+1
            NPFGS=NPFGS+1
            CALL WARMSG(NPFGS,MXNPFG,'TRIANG','MXNPFG',3)
            XSFG(NPFGS,1)=XWFG(J,1)
            XSFG(NPFGS,2)=XWFG(J,2)
            XSFG(NPFGS,3)=XWFG(J,3)
            DO IK=1,NCC
              CSFG(NPFGS,IK)=CWFG(J,IK)
            ENDDO
            MPLOCS(NPFGS)=M
            XWRK(NK,1)=XSFG(NPFGS,1)
            XWRK(NK,2)=XSFG(NPFGS,2)
            XWRK(NK,3)=XSFG(NPFGS,3)
            NWRK(NK)=NPFGS
          ENDDO
        ENDIF
C
  20    CONTINUE
        IDO=1
        K=0
  25    CONTINUE
          K=K+1
          IF(MINFGW(MP,K).NE.0)IDO=0
        IF(IDO.EQ.1 .AND. K.LT.NCC)GOTO 25
        IF(IDO.EQ.0)THEN
          DO K=1,NCC
            J=MINFGW(MP,K)
            DO I=1,NK
              IF(DABS(XWRK(I,1)-XWFG(J,1)).LE.EPSX .AND.
     >           DABS(XWRK(I,2)-XWFG(J,2)).LE.EPSX .AND.
     >           DABS(XWRK(I,3)-XWFG(J,3)).LE.EPSX)GOTO 30
            ENDDO
            NK=NK+1
            NPFGS=NPFGS+1
            CALL WARMSG(NPFGS,MXNPFG,'TRIANG','MXNPFG',4)
            XSFG(NPFGS,1)=XWFG(J,1)
            XSFG(NPFGS,2)=XWFG(J,2)
            XSFG(NPFGS,3)=XWFG(J,3)
            DO IK=1,NCC
              CSFG(NPFGS,IK)=CWFG(J,IK)
            ENDDO
            MPLOCS(NPFGS)=M
            XWRK(NK,1)=XSFG(NPFGS,1)
            XWRK(NK,2)=XSFG(NPFGS,2)
            XWRK(NK,3)=XSFG(NPFGS,3)
            NWRK(NK)=NPFGS
          ENDDO
        ENDIF
C
  30    CONTINUE
        IF(NK.EQ.0)GOTO 400
C
C ----- WHEN THE HIGH OR THE LOW EXISTS, START TRIANGULATION IN THIS
C       SUBELEMENT
C
        NKK=NK
C
        IF(NODE.EQ.4)THEN
          NV=1
          DO J1=1,4
            NK=NK+1
            CALL WARMSG(NK,MXNPK,'TRIANG','MXNPK ',5)
            IEM=ISE(MP+NFGMB(M),J1)
            XWRK(NK,1)=XSFG(IEM,1)
            XWRK(NK,2)=XSFG(IEM,2)
            XWRK(NK,3)=XSFG(IEM,3)
            NWRK(NK)=IEM
            IZE(1,J1)=NK
            XX(J1)=XWRK(NK,1)
            YY(J1)=XWRK(NK,2)
            ZZ(J1)=XWRK(NK,3)
          ENDDO
          CALL CENTER
     I                 (XX,YY,ZZ,
     O                  XCN,YCN,ZCN)
          XV(1,1)=XCN
          XV(1,2)=YCN
          XV(1,3)=ZCN
          DO J1=5,8
            IZE(1,J1)=0
          ENDDO
        ELSEIF(NODE.EQ.6)THEN
          DO J1=1,6
            NK=NK+1
            CALL WARMSG(NK,MXNPK,'TRIANG','MXNPK ',6)
            IEM=ISE(MP+NFGMB(M),J1)
            XWRK(NK,1)=XSFG(IEM,1)
            XWRK(NK,2)=XSFG(IEM,2)
            XWRK(NK,3)=XSFG(IEM,3)
            NWRK(NK)=IEM
          ENDDO
          NV=0
          DO 175 J1=1,3
            NV=NV+1
            DO 174 J2=1,4
              JJ=KID(J2,J1,1)
              DO 179 J4=NKK+1,NK
                IEM=ISE(MP+NFGMB(M),JJ)
                NJ=NWRK(J4)
                IF(IEM.EQ.NJ)THEN
                  IZE(NV,J2)=J4
                  XX(J2)=XWRK(J4,1)
                  YY(J2)=XWRK(J4,2)
                  ZZ(J2)=XWRK(J4,3)
                  GOTO 174
                ENDIF
  179         CONTINUE
  174       CONTINUE
C
            CALL CENTER
     I                   (XX,YY,ZZ,
     O                    XCN,YCN,ZCN)
            XV(NV,1)=XCN
            XV(NV,2)=YCN
            XV(NV,3)=ZCN
  175     CONTINUE
C
C ***** CHECK IZE(M,5...8)
C
          DO 185 J1=1,NV
            DO 181 J2=5,8
              DO I=1,3
                I1=KID(I,J2-4,2)
                N(I)=IZE(J1,I1)
              ENDDO
              DO 201 JJ1=1,NV
                IF(J1.EQ.JJ1)GOTO 201
                DO 202 JJ2=1,4
                  IF(N(1).EQ.IZE(JJ1,JJ2))THEN
                    DO 203 JJ3=1,4
                      IF(N(2).EQ.IZE(JJ1,JJ3))THEN
                        DO 204 JJ4=1,4
                          IF(N(3).EQ.IZE(JJ1,JJ4))THEN
                            IZE(J1,J2)=JJ1
                            GOTO 181
                          ENDIF
  204                   CONTINUE
                      ENDIF
  203               CONTINUE
                  ENDIF
  202           CONTINUE
  201         CONTINUE
              IZE(J1,J2)=0
  181       CONTINUE
  185     CONTINUE
        ELSE
          DO J1=1,8
            NK=NK+1
            CALL WARMSG(NK,MXNPK,'TRIANG','MXNPK ',7)
            IEM=ISE(MP+NFGMB(M),J1)
            XWRK(NK,1)=XSFG(IEM,1)
            XWRK(NK,2)=XSFG(IEM,2)
            XWRK(NK,3)=XSFG(IEM,3)
            NWRK(NK)=IEM
          ENDDO
          NV=0
          DO 275 J1=1,6
            NV=NV+1
            DO 274 J2=1,4
              JJ=KID(J2,J1,3)
              DO 279 J4=NKK+1,NK
                IEM=ISE(MP+NFGMB(M),JJ)
                NJ=NWRK(J4)
                IF(IEM.EQ.NJ)THEN
                  IZE(NV,J2)=J4
                  XX(J2)=XWRK(J4,1)
                  YY(J2)=XWRK(J4,2)
                  ZZ(J2)=XWRK(J4,3)
                  GOTO 274
                ENDIF
  279         CONTINUE
  274       CONTINUE
C
            CALL CENTER
     I                     (XX,YY,ZZ,
     O                      XCN,YCN,ZCN)
            XV(NV,1)=XCN
            XV(NV,2)=YCN
            XV(NV,3)=ZCN
  275     CONTINUE
C
C ***** CHECK IZE(M,5...8)
C
          DO 285 J1=1,NV
            DO 281 J2=5,8
              DO I=1,3
                I1=KID(I,J2-4,2)
                N(I)=IZE(J1,I1)
              ENDDO
              DO 301 JJ1=1,NV
                IF(J1.EQ.JJ1)GOTO 301
                DO 302 JJ2=1,4
                  IF(N(1).EQ.IZE(JJ1,JJ2))THEN
                    DO 303 JJ3=1,4
                      IF(N(2).EQ.IZE(JJ1,JJ3))THEN
                        DO 304 JJ4=1,4
                          IF(N(3).EQ.IZE(JJ1,JJ4))THEN
                            IZE(J1,J2)=JJ1
                            GOTO 281
                          ENDIF
  304                   CONTINUE
                      ENDIF
  303               CONTINUE
                  ENDIF
  302           CONTINUE
  301         CONTINUE
              IZE(J1,J2)=0
  281       CONTINUE
  285     CONTINUE
        ENDIF
C
        DO I=1,NKK
          IAP(I)=I
        ENDDO
C
        CALL GRID3DN
     I            (NKK,XWRK,MXNPK,MXNEK,IAP,EPSX,
     M             IZE,NV,XV,IBP)
        NLNK=NKK+1
        DO NL=NPFGS,NPFGS-NKK+1,-1
          NLNK=NLNK-1
          XSFG(NL,1)=XWRK(NLNK,1)
          XSFG(NL,2)=XWRK(NLNK,2)
          XSFG(NL,3)=XWRK(NLNK,3)
        ENDDO
C
C ***** RECORD INFORMATION OF ISE
C
        DO 289 NL=1,NV
          IF(NL.EQ.1)THEN
            DO I1=1,4
              IZEM=IZE(NL,I1)
              ISE(MP+NFGMB(M),I1)=NWRK(IZEM)
            ENDDO
            ISE(MP+NFGMB(M),5)=0
            ISE(MP+NFGMB(M),6)=0
            ISE(MP+NFGMB(M),7)=0
            ISE(MP+NFGMB(M),8)=0
          ELSE
            NE=NE+1
            CALL WARMSG(NE,MXKGL,'TRIANG','MXKGL ',8)
            DO I1=1,4
              IZEM=IZE(NL,I1)
              ISE(NE,I1)=NWRK(IZEM)
            ENDDO
            ISE(NE,5)=0
            ISE(NE,6)=0
            ISE(NE,7)=0
            ISE(NE,8)=0
          ENDIF
  289   CONTINUE
  400 CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE CENTER
     I                   (XX,YY,ZZ,
     O                    XCN,YCN,ZCN)
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8)
      DIMENSION A(3),B(3),C(3),D(3)
C
C ***** FROM 2,3 ----- DETERMIN PLANE A : A1*X+B1*Y+C1*Z=D1
C
      A(1)=XX(3)-XX(2)
      B(1)=YY(3)-YY(2)
      C(1)=ZZ(3)-ZZ(2)
      D(1)=(A(1)*(XX(3)+XX(2))+B(1)*(YY(3)+YY(2))+C(1)*(ZZ(3)+ZZ(2)))/2.
C
C ***** FROM 3,4 ----- DETERMINE PLANE B : A2*X+B2*Y+C2*Z=D2
C
      A(2)=XX(4)-XX(3)
      B(2)=YY(4)-YY(3)
      C(2)=ZZ(4)-ZZ(3)
      D(2)=(A(2)*(XX(4)+XX(3))+B(2)*(YY(4)+YY(3))+C(2)*(ZZ(4)+ZZ(3)))/2.
C
C ***** FROM 2,3,4 ----- DETERMINE PLANE CONTAINING 2,3,4 : A3*X+B3*Y+
C                        C3*Z=D3
C
      A(3)=(YY(3)-YY(2))*(ZZ(4)-ZZ(2))-(YY(4)-YY(2))*(ZZ(3)-ZZ(2))
      B(3)=(ZZ(3)-ZZ(2))*(XX(4)-XX(2))-(ZZ(4)-ZZ(2))*(XX(3)-XX(2))
      C(3)=(XX(3)-XX(2))*(YY(4)-YY(2))-(XX(4)-XX(2))*(YY(3)-YY(2))
      D(3)=A(3)*XX(2)+B(3)*YY(2)+C(3)*ZZ(2)
C
C ***** DETERMINE LOCATION OF P ----- P IS THE COMMON POINT OF THOSE
C                                     THREE PLANES
      DJAB=0.0D0
      DJABX=0.0D0
      DJABY=0.0D0
      DJABZ=0.0D0
      DO 10 I1=1,3
        I2=I1+1
        I3=I2+1
        IF(I2.GT.3)I2=I2-3
        IF(I3.GT.3)I3=I3-3
        DJAB=DJAB+A(I1)*B(I2)*C(I3)-C(I1)*B(I2)*A(I3)
        DJABX=DJABX+D(I1)*B(I2)*C(I3)-C(I1)*B(I2)*D(I3)
        DJABY=DJABY+A(I1)*D(I2)*C(I3)-C(I1)*D(I2)*A(I3)
        DJABZ=DJABZ+A(I1)*B(I2)*D(I3)-D(I1)*B(I2)*A(I3)
   10 CONTINUE
      XP=DJABX/DJAB
      YP=DJABY/DJAB
      ZP=DJABZ/DJAB
C
C ***** DETERMINE THE VECTOR OF THE LINE THAT CONTAINED BY BOTH PLANE
C       A AND PLANE B
C
      AA=A(3)
      BB=B(3)
      CC=C(3)
C
C ***** DETERMINE THE COORDINATE OF SPHERE CENTROID
C
      E1=XP-XX(1)
      E2=YP-YY(1)
      E3=ZP-ZZ(1)
      F1=XP-XX(4)
      F2=YP-YY(4)
      F3=ZP-ZZ(4)
      DDEN=2.0D0*(AA*(E1-F1)+BB*(E2-F2)+CC*(E3-F3))
      DNUM=F1**2+F2**2+F3**2-(E1**2+E2**2+E3**2)
      TT=DNUM/DDEN
      XCN=XP+TT*AA
      YCN=YP+TT*BB
      ZCN=ZP+TT*CC
c     DO 500 I=1,4
c       DIST=(XX(I)-XCN)**2+(YY(I)-YCN)**2+(ZZ(I)-ZCN)**2
c500  CONTINUE
C
      RETURN
      END
c
c
c
      SUBROUTINE GRID3DN
     I                  (KFN,XWRK,MXNPFG,MXKGL,IAP,EPSX,
     O                   IZE,NV,
     M                   XV,IBP)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NOWRK=300,LAYER=300)
      DIMENSION XWRK(MXNPFG,3),XV(MXKGL,3)
      DIMENSION IZE(MXKGL,8),IZEW(NOWRK,8),IAP(MXNPFG)
      DIMENSION XX(8),YY(8),ZZ(8),KD(NOWRK),JJ(LAYER),KK(LAYER)
      DIMENSION XVW(NOWRK),YVW(NOWRK),ZVW(NOWRK)
      DIMENSION IBP(MXKGL,4)
C
C ***** CHECK BOUNDARY PLANES OF SF REGION
C
      NBP=0
      DO 30 M=1,NV
        IF(IZE(M,5).EQ.0)THEN
          NBP=NBP+1
          IBP(NBP,1)=IZE(M,2)
          IBP(NBP,2)=IZE(M,3)
          IBP(NBP,3)=IZE(M,4)
          IBP(NBP,4)=IZE(M,1)
        ENDIF
        IF(IZE(M,6).EQ.0)THEN
          NBP=NBP+1
          IBP(NBP,1)=IZE(M,3)
          IBP(NBP,2)=IZE(M,1)
          IBP(NBP,3)=IZE(M,4)
          IBP(NBP,4)=IZE(M,2)
        ENDIF
        IF(IZE(M,7).EQ.0)THEN
          NBP=NBP+1
          IBP(NBP,1)=IZE(M,4)
          IBP(NBP,2)=IZE(M,1)
          IBP(NBP,3)=IZE(M,2)
          IBP(NBP,4)=IZE(M,3)
        ENDIF
        IF(IZE(M,8).EQ.0)THEN
          NBP=NBP+1
          IBP(NBP,1)=IZE(M,1)
          IBP(NBP,2)=IZE(M,3)
          IBP(NBP,3)=IZE(M,2)
          IBP(NBP,4)=IZE(M,4)
        ENDIF
  30  CONTINUE
C
C ##### STEP 1
C
C     WRITE(16,*)'NV=',NV,'-----BEFORE GRID3D'
      NCOUNT=0
      DO 900 NNN=KFN,1,-1
        N=IAP(NNN)
        ICHECK=0
C       PRINT *,'NNN=',NNN,' N=',N
        XPP=XWRK(N,1)
        YPP=XWRK(N,2)
        ZPP=XWRK(N,3)
C
C ##### STEP2 & STEP3
C
  40    CONTINUE
        L=0
        ND=0
C
C
        IJK1=1
        IJK2=NV
        DO 60 I=IJK1,IJK2
          XP=XPP
          YP=YPP
          ZP=ZPP
          J=NV-I+1
          DO 55 IJ=1,4
            IF(IJ.EQ.1)THEN
              K1=2
              K2=4
              K3=3
            ELSEIF(IJ.EQ.2)THEN
              K1=3
              K2=4
              K3=1
            ELSEIF(IJ.EQ.3)THEN
              K1=4
              K2=2
              K3=1
            ELSE
              K1=1
              K2=2
              K3=3
            ENDIF
            DO 51 J1=1,3
              IF(J1.EQ.1)JI=K1
              IF(J1.EQ.2)JI=K2
              IF(J1.EQ.3)JI=K3
              NN=IZE(J,JI)
              XX(J1)=XWRK(NN,1)
              YY(J1)=XWRK(NN,2)
              ZZ(J1)=XWRK(NN,3)
   51       CONTINUE
            CANG=FCOS(XX,YY,ZZ,XP,YP,ZP,0)
            IF(CANG.LT.0.0D0)GOTO 60
            IF(DABS(CANG).LE.EPSX)THEN
              CALL ONPLAN
     I           (XX,YY,ZZ,1,2,3,8,
     M            XP,YP,ZP)
            ENDIF
   55     CONTINUE
          J2=IZE(J,1)
          D1=DSQRT((XP-XV(J,1))**2+(YP-XV(J,2))**2+(ZP-XV(J,3))**2)
          D2=DSQRT((XWRK(J2,1)-XV(J,1))**2+(XWRK(J2,2)-XV(J,2))**2+
     >             (XWRK(J2,3)-XV(J,3))**2)
          IF(D1.LT.D2 .OR. DABS(D1-D2).LE.EPSX)GOTO 70
  60    CONTINUE
C
C ***** THIS POINT COULD NOT BE FOUND OUT IN ANY TETRAHEDRAL
C       SUBELEMENTS, BECAUSE THIS POINT IS ON (VERY CLOSE TO) THE
C       BOUNDARY OF THE WORKING SF REGION.  WE ARE GOING TO SKIP THIS
C       POINT
C
        NCOUNT=NCOUNT+1
        GOTO 900
  70    CONTINUE
        XWRK(N,1)=XP
        XWRK(N,2)=YP
        XWRK(N,3)=ZP
        ND=ND+1
        CALL WARMSG(ND,NOWRK,'GRID3D','NOWRK ',1)
        KD(ND)=J
C
C ***** COMPOSE NEW TETRAHEDRAL SUBELMENTS AFTER POINT P ADDED INTO
C       THIS SF REGION
        CALL NEWTETRA
     I      (ND,XP,YP,ZP,KD,IZE,XWRK,NOWRK,MXKGL,MXNPFG,N,EPSX,
     O       NW,IZEW,XVW,YVW,ZVW,INEW)
        IF(INEW.EQ.1)THEN
          NCOUNT=NCOUNT+1
          GOTO 900
        ENDIF
C
C ***** CHECK IF ANY GEOMETRIC VIOLATIONS EXIST WHEN NEW TETRAHEDRAL
C       SUBELEMENTS ARE GENERATED
        CALL CKTETRA
     I                (NW,IZEW,XWRK,NOWRK,MXNPFG,
     O                 KCK)
C
        IF(KCK.EQ.1)THEN
          NCOUNT=NCOUNT+1
          GOTO 900
        ENDIF
C
        J2=IZE(J,1)
        L=L+1
        CALL WARMSG(L,LAYER,'GRID3D','LAYER ',2)
        JJ(L)=J
c       JL=0
        JM=JJ(L)
        K=0
  80    CONTINUE
        IF(K.EQ.4)GOTO 100
        K=K+1
        KK(L)=K
        J1=IZE(JM,K+4)
        DO 85 JD=1,ND
          IF(J1.EQ.KD(JD))GOTO 100
  85    CONTINUE
        IF(J1.EQ.0)GOTO 100
        J2=IZE(J1,1)
        J3=IZE(J1,2)
c       J4=IZE(J1,3)
c       J5=IZE(J1,4)
        D1=DSQRT((XP-XV(J1,1))**2+(YP-XV(J1,2))**2+(ZP-XV(J1,3))**2)
        D2=DSQRT((XWRK(J2,1)-XV(J1,1))**2+(XWRK(J2,2)-XV(J1,2))**2+
     >     (XWRK(J2,3)-XV(J1,3))**2)
c       D3=DSQRT((XWRK(J3,1)-XV(J1,1))**2+(XWRK(J3,2)-XV(J1,2))**2+
c    >     (XWRK(J3,3)-XV(J1,3))**2)
c       D4=DSQRT((XWRK(J4,1)-XV(J1,1))**2+(XWRK(J4,2)-XV(J1,2))**2+
c    >     (XWRK(J4,3)-XV(J1,3))**2)
c       D5=DSQRT((XWRK(J5,1)-XV(J1,1))**2+(XWRK(J5,2)-XV(J1,2))**2+
c    >     (XWRK(J5,3)-XV(J1,3))**2)
        IF(D1.GT.D2 .OR. DABS(D1-D2).LE.EPSX)GOTO 100
        DO IID=1,ND
          IF(J1.EQ.KD(IID))GOTO 100
        ENDDO
C
        IF(ICHECK.EQ.1)THEN
          IF(K.EQ.1)THEN
            JM1=4
            JM2=3
            JM3=2
          ELSEIF(K.EQ.2)THEN
            JM1=4
            JM2=1
            JM3=3
          ELSEIF(K.EQ.3)THEN
            JM1=4
            JM2=2
            JM3=1
          ELSE
            JM1=1
            JM2=2
            JM3=3
          ENDIF
          JJ1=IZE(JM,JM1)
          JJ2=IZE(JM,JM2)
          JJ3=IZE(JM,JM3)
          DO IID=1,4
            IF(IZE(J1,IID).NE.JJ1 .AND. IZE(J1,IID).NE.JJ2 .AND.
     >        IZE(J1,IID).NE.JJ3)GOTO 90
          ENDDO
  90      CONTINUE
          JCK=IZE(J1,IID)
          XQ=XWRK(JCK,1)
          YQ=XWRK(JCK,2)
          ZQ=XWRK(JCK,3)
          DO 95 IID=1,NBP
            DO IJD=1,4
              JID=IBP(IID,IJD)
              XX(IJD)=XWRK(JID,1)
              YY(IJD)=XWRK(JID,2)
              ZZ(IJD)=XWRK(JID,3)
            ENDDO
C
C ***** CHECK IF ANY GEOMETRIC VIOLATIONS EXIST IN THE PROCESS OF
C       DETERMINING DELETED TETRAHEDRAL SUBELEMENTS.
C
            CALL CKSFBP
     I                 (XX,YY,ZZ,XP,YP,ZP,XQ,YQ,ZQ,IID,MXKGL,IBP,NBP,
     I                  MXNPFG,XWRK,
     O                  KCK)
C           PRINT *,'------------------'
            IF(KCK.EQ.1)GOTO 100
  95      CONTINUE
        ENDIF
C
        ND=ND+1
        CALL WARMSG(ND,NOWRK,'GRID3D','NOWRK ',3)
        KD(ND)=J1
C
        IF(ICHECK.EQ.1)THEN
C ***** COMPOSE NEW TETRAHEDRAL SUBELMENTS AFTER POINT P ADDED INTO
C       THIS SF REGION
          CALL NEWTETRA
     I        (ND,XP,YP,ZP,KD,IZE,XWRK,NOWRK,MXKGL,MXNPFG,N,EPSX,
     O         NW,IZEW,XVW,YVW,ZVW,INEW)
          IF(INEW.EQ.1)THEN
            ND=ND-1
            GOTO 100
          ENDIF
C
C ***** CHECK IF ANY GEOMETRIC VIOLATIONS EXIST WHEN NEW TETRAHEDRAL
C       SUBELEMENTS ARE GENERATED
          CALL CKTETRA
     I                  (NW,IZEW,XWRK,NOWRK,MXNPFG,
     O                   KCK)
          IF(KCK.EQ.1)THEN
            ND=ND-1
            GOTO 100
          ENDIF
        ENDIF
C
        L=L+1
        CALL WARMSG(L,LAYER,'GRID3D','LAYER ',4)
        JJ(L)=J1
        JM=JJ(L)
c       IF(L.NE.1)THEN
c         JL=JJ(L-1)
c       ELSE
c         JL=0
c       ENDIF
        K=0
        GOTO 80
  100   CONTINUE
        IF(K.LT.4)THEN
          GOTO 80
        ELSEIF(L.EQ.1)THEN
          GOTO 200
        ELSE
          L=L-1
          K=KK(L)
          JM=JJ(L)
c         IF(L.NE.1)THEN
c           JL=JJ(L-1)
c         ELSE
c           JL=0
c         ENDIF
          GOTO 80
        ENDIF
  200   CONTINUE
C
C ##### STEP 4,5,6
C
C ***** COMPOSE NEW TETRAHEDRAL SUBELMENTS AFTER POINT P ADDED INTO
C       THIS SF REGION
C
        CALL NEWTETRA
     I      (ND,XP,YP,ZP,KD,IZE,XWRK,NOWRK,MXKGL,MXNPFG,N,EPSX,
     O       NW,IZEW,XVW,YVW,ZVW,INEW)
        IF(INEW.EQ.1 .AND. ICHECK.EQ.1)THEN
          NCOUNT=NCOUNT+1
          GOTO 900
c         WRITE(16,*)'ERROR OCCURRED AFTER NEWTETRA, NNN,N=',NNN,N
c         STOP
        ENDIF
        IF(INEW.EQ.1 .AND. ICHECK.EQ.0)THEN
          NCOUNT=NCOUNT+1
          GOTO 900
c         ICHECK=1
c         GOTO 40
        ENDIF
C
C ***** CHECK IF ANY GEOMETRIC VIOLATIONS EXIST WHEN NEW TETRAHEDRAL
C       SUBELEMENTS ARE GENERATED
        CALL CKTETRA
     I                  (NW,IZEW,XWRK,NOWRK,MXNPFG,
     O                   KCK)
        IF(ICHECK.EQ.1 .AND. KCK.EQ.1)THEN
          NCOUNT=NCOUNT+1
          GOTO 900
c         WRITE(16,*)'ERROR OCCURRED AFTER CKTETRA, NNN,N=', NNN,N
c         STOP
        ENDIF
        IF(ICHECK.EQ.0 .AND. KCK.EQ.1)THEN
          NCOUNT=NCOUNT+1
          GOTO 900
c         ICHECK=1
c         GOTO 40
        ENDIF
C
C ##### STEP 7
C
C ***** CHECK IF THERE ARE INAPPROPRIATE CONNECTIONS BETWEEN NEWLY
C       COMPOSED TETRAHEDRONS ----- CHECK IF ANY TWO NEW TETRAHEDRONS
C       HAVE UNUSUAL CONNECTIONS (I.E. THEY CONNECT TO EACH OTHER WITH
C       MORE THAN ONE COMMON TRIANGLE)
C
C ***** DETERMINE IZEW(NW,6),IZEW(NW,7), AND IZEW(NW,8)
C
        DO 400 I=1,NW
          DO 390 I1=6,8
            IF(I1.EQ.6)THEN
              K1=3
              K2=4
              K3=1
            ELSEIF(I1.EQ.7)THEN
              K1=4
              K2=2
              K3=1
            ELSE
              K1=1
              K2=2
              K3=3
            ENDIF
            N1=IZEW(I,K1)
            N2=IZEW(I,K2)
            N3=IZEW(I,K3)
            DO 350 J=1,NW
              IF(J.EQ.I)GOTO 350
              DO 310 J1=1,4
                IF(N1.EQ.IZEW(J,J1))THEN
                  DO 320 J2=1,4
                    IF(N2.EQ.IZEW(J,J2))THEN
                      DO 330 J3=1,4
                        IF(N3.EQ.IZEW(J,J3))THEN
                          IZEW(I,I1)=J
                          GOTO 390
                        ENDIF
  330                 CONTINUE
                    ENDIF
  320             CONTINUE
                  GOTO 350
                ENDIF
  310         CONTINUE
  350       CONTINUE
            IZEW(I,I1)=0
  390     CONTINUE
  400   CONTINUE
C
C ##### STEP 8 OVERWRITE
C
        NVV=NV
        IF(NW.LT.ND)THEN
          NDW=NW
        ELSE
          NDW=ND
        ENDIF
        DO 500 I=1,NDW
          M1=KD(I)
          XV(M1,1)=XVW(I)
          XV(M1,2)=YVW(I)
          XV(M1,3)=ZVW(I)
          DO 410 I1=1,8
            IZE(M1,I1)=IZEW(I,I1)
            IF(I1.GE.6)THEN
              IF(IZE(M1,I1).EQ.0)GOTO 410
              IF(IZE(M1,I1).LE.ND)THEN
                MM=IZE(M1,I1)
                IZE(M1,I1)=KD(MM)
              ELSE
                MM=IZE(M1,I1)
                IZE(M1,I1)=NVV+MM-ND
              ENDIF
            ENDIF
  410     CONTINUE
C         IF(M1.EQ.430 .OR. M1.EQ.1468)PRINT *,'NNN=',NNN,'-----'
          M2=IZE(M1,5)
          IF(M2.EQ.0)GOTO 500
          N1=IZE(M1,2)
          N2=IZE(M1,3)
          N3=IZE(M1,4)
          DO 420 I1=5,8
            IF(I1.EQ.5)THEN
              NN1=IZE(M2,2)
              NN2=IZE(M2,3)
              NN3=IZE(M2,4)
            ELSEIF(I1.EQ.6)THEN
              NN1=IZE(M2,3)
              NN2=IZE(M2,1)
              NN3=IZE(M2,4)
            ELSEIF(I1.EQ.7)THEN
              NN1=IZE(M2,4)
              NN2=IZE(M2,1)
              NN3=IZE(M2,2)
            ELSE
              NN1=IZE(M2,1)
              NN2=IZE(M2,3)
              NN3=IZE(M2,2)
            ENDIF
            IF(N1.EQ.NN1 .OR. N1.EQ.NN2 .OR. N1.EQ.NN3)THEN
              IF(N2.EQ.NN1 .OR. N2.EQ.NN2 .OR. N2.EQ.NN3)THEN
                IF(N3.EQ.NN1 .OR. N3.EQ.NN2 .OR. N3.EQ.NN3)THEN
                  IZE(M2,I1)=M1
                  GOTO 500
                ENDIF
              ENDIF
            ENDIF
  420     CONTINUE
  500   CONTINUE
C
        IF(NW.EQ.ND)GOTO 890
        IF(NW.LT.ND)THEN
          NDIF=ND-NW
          M1=NV+1
          DO 480 I=1,NDIF
            M=KD(ND-I+1)
            IF(M.GT.NV-NDIF)GOTO 480
  450       CONTINUE
            M1=M1-1
            DO IJ=1,NDIF
              MM=KD(ND-IJ+1)
              IF(M1.EQ.MM)GOTO 450
            ENDDO
            XV(M,1)=XV(M1,1)
            XV(M,2)=XV(M1,2)
            XV(M,3)=XV(M1,3)
            DO IJ=1,8
              IZE(M,IJ)=IZE(M1,IJ)
            ENDDO
            DO 460 II=1,NV
              DO IJ=5,8
                IF(IZE(II,IJ).EQ.M1)IZE(II,IJ)=M
              ENDDO
  460       CONTINUE
  480     CONTINUE
          NV=NV-NDIF
          GOTO 890
        ENDIF
C
        DO 600 I=ND+1,NW
          NV=NV+1
          CALL WARMSG(NV,MXKGL,'GRID3D','MXKGL ',5)
          XV(NV,1)=XVW(I)
          XV(NV,2)=YVW(I)
          XV(NV,3)=ZVW(I)
          DO 510 I1=1,8
            IZE(NV,I1)=IZEW(I,I1)
            IF(I1.GE.6)THEN
              IF(IZE(NV,I1).EQ.0)GOTO 510
              IF(IZE(NV,I1).LE.ND)THEN
                MM=IZE(NV,I1)
                IZE(NV,I1)=KD(MM)
              ELSE
                MM=IZE(NV,I1)
                IZE(NV,I1)=NVV+MM-ND
              ENDIF
            ENDIF
  510     CONTINUE
          M2=IZE(NV,5)
          IF(M2.EQ.0)GOTO 600
          N1=IZE(NV,2)
          N2=IZE(NV,3)
          N3=IZE(NV,4)
          DO 520 I1=5,8
            IF(I1.EQ.5)THEN
              NN1=IZE(M2,2)
              NN2=IZE(M2,3)
              NN3=IZE(M2,4)
            ELSEIF(I1.EQ.6)THEN
              NN1=IZE(M2,3)
              NN2=IZE(M2,1)
              NN3=IZE(M2,4)
            ELSEIF(I1.EQ.7)THEN
              NN1=IZE(M2,4)
              NN2=IZE(M2,1)
              NN3=IZE(M2,2)
            ELSE
              NN1=IZE(M2,1)
              NN2=IZE(M2,3)
              NN3=IZE(M2,2)
            ENDIF
            IF(N1.EQ.NN1 .OR. N1.EQ.NN2 .OR. N1.EQ.NN3)THEN
              IF(N2.EQ.NN1 .OR. N2.EQ.NN2 .OR. N2.EQ.NN3)THEN
                IF(N3.EQ.NN1 .OR. N3.EQ.NN2 .OR. N3.EQ.NN3)THEN
                  IZE(M2,I1)=NV
                  GOTO 600
                ENDIF
              ENDIF
            ENDIF
  520     CONTINUE
  600   CONTINUE
C
C ##### STEP 9
C
  890 CONTINUE
C     PRINT *,'NNN,N,ND,NW,NV=',NNN,N,ND,NW,NV
  900 CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE CKTETRA
     I                  (NW,IZEW,XWRK,NOWRK,MXNPFG,
     O                   KCK)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IZEW(NOWRK,8),XWRK(MXNPFG,3)
      DIMENSION XX(8),YY(8),ZZ(8)
C
      KCK=0
        DO 888 I=1,NW
          N1=IZEW(I,4)
          N2=IZEW(I,3)
          N3=IZEW(I,2)
          X1=XWRK(N1,1)
          Y1=XWRK(N1,2)
          Z1=XWRK(N1,3)
          X2=XWRK(N2,1)
          Y2=XWRK(N2,2)
          Z2=XWRK(N2,3)
          X3=XWRK(N3,1)
          Y3=XWRK(N3,2)
          Z3=XWRK(N3,3)
          DO 777 IK=1,NW
            IF(IK.EQ.I)GOTO 777
            DO 667 IKK=1,3
              IF(IKK.EQ.1)THEN
                NPN=N1
                XPP=X1
                YPP=Y1
                ZPP=Z1
              ELSEIF(IKK.EQ.2)THEN
                NPN=N2
                XPP=X2
                YPP=Y2
                ZPP=Z2
              ELSE
                NPN=N3
                XPP=X3
                YPP=Y3
                ZPP=Z3
              ENDIF
              DO 666 KI=1,3
                IF(KI.EQ.1)THEN
                  NN1=IZEW(IK,4)
                  NN2=IZEW(IK,2)
                  NN3=IZEW(IK,1)
                ELSEIF(KI.EQ.2)THEN
                  NN1=IZEW(IK,4)
                  NN2=IZEW(IK,1)
                  NN3=IZEW(IK,3)
                ELSE
                  NN1=IZEW(IK,1)
                  NN2=IZEW(IK,2)
                  NN3=IZEW(IK,3)
                ENDIF
                IF(NPN.EQ.NN1 .OR. NPN.EQ.NN2 .OR. NPN.EQ.NN3)GOTO 667
                XX(1)=XWRK(NN1,1)
                YY(1)=XWRK(NN1,2)
                ZZ(1)=XWRK(NN1,3)
                XX(2)=XWRK(NN2,1)
                YY(2)=XWRK(NN2,2)
                ZZ(2)=XWRK(NN2,3)
                XX(3)=XWRK(NN3,1)
                YY(3)=XWRK(NN3,2)
                ZZ(3)=XWRK(NN3,3)
                CANG=FCOS(XX,YY,ZZ,XPP,YPP,ZPP,0)
                IF(CANG.LT.0.0D0)GOTO 777
 666          CONTINUE
              KCK=1
 667        CONTINUE
 777      CONTINUE
 888    CONTINUE
      RETURN
      END
c
c
c
      SUBROUTINE CKSFBP
     I                 (XX,YY,ZZ,XP,YP,ZP,XQ,YQ,ZQ,IID,MXKGL,IBP,NBP,
     I                  MXNPFG,XWRK,
     O                  KCK)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),DL(8),DNX(8),IBP(MXKGL,4)
      DIMENSION XWRK(MXNPFG,3)
C
      DATA EPS/1.0D-6/,EPSR/1.0D-10/
C
      KCK=0
      A=(YY(2)-YY(1))*(ZZ(3)-ZZ(1))-(YY(3)-YY(1))*(ZZ(2)-ZZ(1))
      B=(ZZ(2)-ZZ(1))*(XX(3)-XX(1))-(ZZ(3)-ZZ(1))*(XX(2)-XX(1))
      C=(XX(2)-XX(1))*(YY(3)-YY(1))-(XX(3)-XX(1))*(YY(2)-YY(1))
      D=-(A*XX(2)+B*YY(2)+C*ZZ(2))
      DEN=A*(XQ-XP)+B*(YQ-YP)+C*(ZQ-ZP)
      IF(DABS(DEN).LE.EPSR)RETURN
      T=-(D+A*XP+B*YP+C*ZP)/DEN
C     PRINT *,T
      IF((T.GT.1.0D0 .AND. DABS(T-1.0D0).GT.EPSR) .OR.
     >   (T.LT.0.0D0 .AND. DABS(T).GT.EPSR))RETURN
      XNN=XP+T*(XQ-XP)
      YNN=YP+T*(YQ-YP)
      ZNN=ZP+T*(ZQ-ZP)
C     PRINT *,XNN,YNN,ZNN
      CALL BASE
c     I         (XX,YY,ZZ,XNN,YNN,ZNN,989,3,1,
     I         (XX,YY,ZZ,XNN,YNN,ZNN,3,1,
     O          DL,DNX,DNX,DNX)
C     PRINT *,'DL',(DL(I),I=1,4)
      DO I=1,3
        IF((DL(I).LT.0.0D0 .AND. DABS(DL(I)).GT.EPS) .OR.
     >     (DL(I).GT.1.0D0 .AND. DABS(DL(I)-1.0D0).GT.EPS))RETURN
      ENDDO
      IF(DABS(T).LE.EPSR)THEN
        T=T+EPSR*1.0D1
        XN=XP+T*(XQ-XP)
        YN=YP+T*(YQ-YP)
        ZN=ZP+T*(ZQ-ZP)
        CANG=FCOS(XX,YY,ZZ,XN,YN,ZN,0)
        IF(CANG.LE.0.0D0) RETURN
      ELSEIF(DABS(T-1.0D0).LE.EPSR)THEN
        T=T-EPSR*1.0D1
        XN=XP+T*(XQ-XP)
        YN=YP+T*(YQ-YP)
        ZN=ZP+T*(ZQ-ZP)
        CANG=FCOS(XX,YY,ZZ,XN,YN,ZN,0)
        IF(CANG.LE.0.0D0)RETURN
      ELSE
        T1=T+EPSR*1.0D1
        IF(T1.GT.1.0D0)GOTO 50
        XN=XP+T1*(XQ-XP)
        YN=YP+T1*(YQ-YP)
        ZN=ZP+T1*(ZQ-ZP)
        CANG=FCOS(XX,YY,ZZ,XN,YN,ZN,0)
        IF(CANG.GT.0.0D0)GOTO 100
C
   50   CONTINUE
        T2=T-EPSR*1.0D1
        IF(T2.LT.0.0D0)RETURN
        XN=XP+T2*(XQ-XP)
        YN=YP+T2*(YQ-YP)
        ZN=ZP+T2*(ZQ-ZP)
        CANG=FCOS(XX,YY,ZZ,XN,YN,ZN,0)
        IF(CANG.LE.0.0D0) RETURN
C
  100   CONTINUE
        DO I=1,3
          I1=I+1
          IF(I1.EQ.4)I1=1
          X1=XX(I)
          Y1=YY(I)
          Z1=ZZ(I)
          X2=XX(I1)
          Y2=YY(I1)
          Z2=ZZ(I1)
          D1=DSQRT((X2-X1)**2+(Y2-Y1)**2+(Z2-Z1)**2)
          DETX=(Y2-Y1)*(ZNN-Z1)-(YNN-Y1)*(Z2-Z1)
          DETY=(Z2-Z1)*(XNN-X1)-(ZNN-Z1)*(X2-X1)
          DETZ=(X2-X1)*(YNN-Y1)-(XNN-X1)*(Y2-Y1)
          DET=DSQRT(DETX**2+DETY**2+DETZ**2)/D1
          IF(DABS(DET).LE.EPSR)GOTO 150
        ENDDO
        KCK=1
        RETURN
C
  150   CONTINUE
        IF(I.EQ.1)THEN
          N1=IBP(IID,1)
          N2=IBP(IID,2)
        ELSEIF(I.EQ.2)THEN
          N1=IBP(IID,2)
          N2=IBP(IID,3)
        ELSE
          N1=IBP(IID,3)
          N2=IBP(IID,1)
        ENDIF
C       PRINT *,'N1,N2=',N1,N2
        DO 200 I=1,NBP
          IF(I.EQ.IID)GOTO 200
          DO 180 J=1,3
            IF(IBP(I,J).EQ.N1)THEN
              DO 170 K=1,3
                IF(IBP(I,K).EQ.N2)THEN
                  DO KK=1,3
                    IK=IBP(I,KK)
                    XX(KK)=XWRK(IK,1)
                    YY(KK)=XWRK(IK,2)
                    ZZ(KK)=XWRK(IK,3)
                  ENDDO
                  CANG=FCOS(XX,YY,ZZ,XN,YN,ZN,0)
                  IF(CANG.LE.0.0D0) RETURN
                ENDIF
 170          CONTINUE
            ENDIF
 180      CONTINUE
 200    CONTINUE
      ENDIF
      KCK=1
      RETURN
      END
c
c
c
      SUBROUTINE NEWTETRA
     I    (ND,XP,YP,ZP,KD,IZE,XWRK,NOWRK,MXKGL,MXNPFG,N,EPSX,
     O     NW,IZEW,XVW,YVW,ZVW,INEW)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION KD(NOWRK),IZE(MXKGL,8),XWRK(MXNPFG,3)
      DIMENSION IZEW(NOWRK,8),XVW(NOWRK),YVW(NOWRK)
      DIMENSION ZVW(NOWRK),XX(8),YY(8),ZZ(8),XW(8),YW(8),ZW(8)
C
        INEW=0
        NW=0
        DO 250 I=1,ND
          M1=KD(I)
          DO 245 I1=1,4
            IF(I1.EQ.1)THEN
              K1=2
              K2=3
              K3=4
            ELSEIF(I1.EQ.2)THEN
              K1=3
              K2=1
              K3=4
            ELSEIF(I1.EQ.3)THEN
              K1=4
              K2=1
              K3=2
            ELSE
              K1=1
              K2=3
              K3=2
            ENDIF
            N1=IZE(M1,K1)
            N2=IZE(M1,K2)
            N3=IZE(M1,K3)
            DO 230 J=1,ND
              M2=KD(J)
              IF(M2.EQ.M1)GOTO 230
              DO 210 J1=1,4
                IF(N1.EQ.IZE(M2,J1))THEN
                  DO 220 J2=1,4
                    IF(N2.EQ.IZE(M2,J2))THEN
                      DO 225 J3=1,4
                        IF(N3.EQ.IZE(M2,J3))GOTO 245
  225                 CONTINUE
                    ENDIF
  220             CONTINUE
                  GOTO 230
                ENDIF
  210         CONTINUE
  230       CONTINUE
            NW=NW+1
            IF(NW.GT.NOWRK)THEN
              WRITE(16,*)'NW.GT.NOWRK ----- IN GRID3DN'
              WRITE(16,*)'N=',N
              INEW=1
              RETURN
            ENDIF
            IZEW(NW,1)=N
            IZEW(NW,2)=N1
            IZEW(NW,3)=N2
            IZEW(NW,4)=N3
C
C
            DO 234 K=1,4
              NN=IZEW(NW,K)
              XX(K)=XWRK(NN,1)
              YY(K)=XWRK(NN,2)
              ZZ(K)=XWRK(NN,3)
  234       CONTINUE
            DO 235 K=2,4
              NN=IZEW(NW,K)
              KK=5-K
              XW(KK)=XWRK(NN,1)
              YW(KK)=XWRK(NN,2)
              ZW(KK)=XWRK(NN,3)
  235       CONTINUE
            CANG=FCOS(XW,YW,ZW,XP,YP,ZP,0)
            IF(DABS(CANG).LE.EPSX)THEN
              CALL ONPLAN(XW,YW,ZW,1,2,3,8,
     >                    XWRK(N,1),XWRK(N,2),XWRK(N,3))
              XP=XWRK(N,1)
              YP=XWRK(N,2)
              ZP=XWRK(N,3)
              NW=NW-1
              GOTO 245
            ENDIF
            CALL CENTER
     I                   (XX,YY,ZZ,
     O                    XCN,YCN,ZCN)
            XVW(NW)=XCN
            YVW(NW)=YCN
            ZVW(NW)=ZCN
C
C ##### CHECK IZEW(NW,5)
C
            DO 240 K=5,8
              MM=IZE(M1,K)
              IF(MM.EQ.0)GOTO 240
              DO 236 K1=1,4
                IF(N1.EQ.IZE(MM,K1))THEN
                  DO 237 K2=1,4
                    IF(N2.EQ.IZE(MM,K2))THEN
                      DO 238 K3=1,4
                        IF(N3.EQ.IZE(MM,K3))THEN
                          IZEW(NW,5)=MM
                          GOTO 245
                        ENDIF
  238                 CONTINUE
                    ENDIF
  237             CONTINUE
                  GOTO 240
                ENDIF
  236         CONTINUE
  240       CONTINUE
            IZEW(NW,5)=0
  245     CONTINUE
  250   CONTINUE
      RETURN
      END
C
c
c
      SUBROUTINE DFPREP
     I    (MAXNP,MAXEL,MXNPFG,MXKGL,MXKGLD,MXNPW,MXADNP,MXJBD,MXKBD,
     I     MXNCC,MXNDB,IDZOOM,NCC,NEL,NNP,NXD,NYD,NZD,EPSX,X,CS,IB,IE,
     I     NLRN,LRL,NLRL,ISE,NFGMB,XSFG,CSFG,NEFGS, MXMSV,MXLSV,IPNTST,
     M     NFGM,XW,MWLOC,MPLOCW,XWFG,CWFG,
     O     NLS, ILSV,IMSV,
     O     MPLOC,ISED,NDBD,MPLOCD,NDFG,NEFGD,NDB,NLRND,LRN,ND)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP,3),CS(MAXNP,MXNCC),IE(MAXEL,11),NLRL(MAXNP)
      DIMENSION IB(MAXNP),NLRN(MAXNP),XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC)
      DIMENSION LRL(MXKBD,MAXNP),ISE(MXKGL,8),NFGMB(MAXEL),NFGM(MAXEL)
      DIMENSION XW(MXNPW,3),MWLOC(MXNPW)
      DIMENSION MPLOCW(MXNPFG),XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC)
      DIMENSION MPLOC(MXNPFG),ISED(MXKGLD,9),NLRND(MXADNP)
      DIMENSION LRN(MXJBD,MXADNP),NDBD(MXNDB),MPLOCD(MXNDB)
c      DIMENSION MCON(6),XX(8),YY(8),ZZ(8),CC(8,7),MCB(3),KID(3,6,3)
      DIMENSION MCON(6),XX(8),YY(8),ZZ(8),MCB(3),KID(3,6,3)
      DIMENSION ILSV(MXLSV,5),IMSV(MXMSV,4),ND(MXADNP)
      DIMENSION KK(100),KJ(100),MCON1(100)
C
      DATA KID/4,8,5, 1,6,5, 2,7,6, 3,8,4, 1,2,3, 5,8,7,
     >         3,1,6, 1,2,4, 2,3,5, 1,2,3, 4,5,6, 0,0,0,
     >         4,3,2, 4,1,3, 4,2,1, 1,2,3, 0,0,0, 0,0,0/
C
      NDFG=NNP
      NEFGD=0
      NMS=0
      NLS=0
C
      DO I=1,NNP
        XWFG(I,1)=X(I,1)
        XWFG(I,2)=X(I,2)
        XWFG(I,3)=X(I,3)
        DO K=1,NCC
          CWFG(I,K)=CS(I,K)
        ENDDO
      ENDDO
      IF(IDZOOM.EQ.0)GOTO 501
C
      DO I=1,NNP
        MPLOC(I)=LRL(1,I)
      ENDDO
C
C ----- FIND OUT ALL THE NON SF ELEMENTS WHICH SURROUND SF REGIONS
C
      nlgelm=0
      DO 10 M=1,NEL
        if(ie(m,11).ne.0)nlgelm=nlgelm+1
        IF(IE(M,11).EQ.0 .OR. IE(M,11).EQ.-1)GOTO 10
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,I,I)
        DO I=1,NODE
          NP=IE(M,I)
          DO J=1,NLRL(NP)
            MJ=LRL(J,NP)
            IF(IE(MJ,11).EQ.0)IE(MJ,11)=-1
          ENDDO
        ENDDO
   10 CONTINUE
c
ccc      print *,'the no. of extended rough elements is ',nlgelm
C
C ----- IN THE FOLLOWING, ARRAY NFGM(M) REPRESENTS THE ACCUMULATED
C       DIFFUSION FINE GRIDS IN THE FIRST M-1 ELEMENTS
      DO M=1,NEL
        NFGM(M)=0
      ENDDO
      NDB=0
      DO 490 M=1,NEL
        IF(IE(M,11).EQ.0)THEN
          IF(M.EQ.1)NFGM(M)=NNP
          IF(M.NE.1)NFGM(M)=NFGM(M-1)
          GOTO 490
        ENDIF
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,II,ID)
        DO I=1,NODE
          IEM=IE(M,I)
          XX(I)=X(IEM,1)
          YY(I)=X(IEM,2)
          ZZ(I)=X(IEM,3)
c         DO K=1,NCC
c           CC(I,K)=CS(IEM,K)
c         ENDDO
        ENDDO
C
C ----- CHECK THE ELEMENTS CONNECTING TO ELEMENT M, AND STORE THE
C       INFORMATION INTO ARRAY MCON(I).  MCON(I) IS THE GLOBAL ELEMENT
C       NUMBER OF THE I-TH ELEMENT CONNECTED TO ELEMENT M. MCON(I) NEEDS
C       TO BE SMALLER THAN M BECAUSE THE ELEMENT IS CHECKED FROM 1 TO
C       NEL.  IT IS NECESSARY TO CHECK ONLY THE ELEMENTS WHICH HAVE BEEN
C       ZOOMED WITH DIFFUSION FINE GRIDS.
C
        DO I=1,100
          MCON1(I)=0
        ENDDO
C
        I=0
        DO 15 J=1,NODE
          IEM=IE(M,J)
          DO 14 K=1,NLRL(IEM)
            MM=LRL(K,IEM)
            IF(MM.EQ.0 .OR. MM.GE.M .OR. IE(MM,11).EQ.0)GOTO 14
            IF(I.EQ.0)THEN
              I=I+1
              MCON1(I)=MM
            ELSE
              DO L=1,I
                IF(MM.EQ.MCON1(I))GOTO 14
              ENDDO
              I=I+1
              MCON1(I)=MM
            ENDIF
  14      CONTINUE
  15    CONTINUE
        NCON1=I
C
C ----- CHECK THE ELEMENTS CONNECTING TO ELEMENT M, AND STORE THE
C       INFORMATION INTO ARRAY MCON(LS) WHICH MEANS THE NUMBER OF
C       OF GLOBAL ELEMENT CONNECTED TO THE LS-TH SURFACE OF ELEMENT
C       M
C
        DO 20 LS=1,II
          MCON(LS)=0
          I1=KID(1,LS,ID)
          N1=IE(M,I1)
          I2=KID(2,LS,ID)
          N2=IE(M,I2)
          I3=KID(3,LS,ID)
          N3=IE(M,I3)
          DO 19 J1=1,NLRL(N1)
            M1=LRL(J1,N1)
            DO 18 J2=1,NLRL(N2)
              M2=LRL(J2,N2)
              IF(M1.NE.M2)GOTO 18
              DO J3=1,NLRL(N3)
                M3=LRL(J3,N3)
                IF(M3.EQ.M1 .AND. M1.NE.M)THEN
                  MCON(LS)=M1
                  GOTO 20
                ENDIF
              ENDDO
   18       CONTINUE
   19     CONTINUE
          MCON(LS)=0
   20   CONTINUE
C
C ----- CHECK THOSE DIFFUSION FINE GRIDS WHICH LOCATE ON THE BOUNDARY OF
C       EACH ELEMENTS
C
        IF(NODE.EQ.8)THEN
          DO 35 K=1,NZD+1
            ZI=2.0D0*DBLE(K-1)/DBLE(NZD)-1.0D0
            DO 34 I=1,NYD+1
              YI=2.0D0*DBLE(I-1)/DBLE(NYD)-1.0D0
              DO 33 J=1,NXD+1
                XI=2.0D0*DBLE(J-1)/DBLE(NXD)-1.0D0
                N=(I-1)*(NXD+1)+J+(K-1)*(NXD+1)*(NYD+1)
                CALL WARMSG(N,MXNPW,'DFPREP','MXNPW ',1)
                CALL BASEXI
     I              (XI,YI,ZI,DL4,XX,YY,ZZ,NODE,
     O               XW(N,1),XW(N,2),XW(N,3))
C
                IF(J.EQ.1 .OR. J.EQ.NXD+1 .OR.  I.EQ.1 .OR. I.EQ.NYD+1
     >             .OR. K.EQ.1 .OR. K.EQ.NZD+1)THEN
C ------------    CHECK THOSE POINTS ON THE BOUNDARY WITH BOTH S-F
C                 ELEMENTS AND THE INNERMOST NON-SF ELEMENTS
                  CALL GLBCHK
     I                (NCON1,N,NNP,NEL,NEFGS,NODE,EPSX,NCC,IDZOOM,
     I                 XW(N,1),XW(N,2),XW(N,3),MCON,MCON1,M,
     I                 NXD,NYD,NZD,MAXNP,MAXEL,MXNPFG,MXKGL,MXNDB,MXNCC,
     I                 MXKBD,MXNPW,MXADNP,I,J,K,DL1,DL2,DL3,DL4,
     I                 NFGM,X,IE,IB,NLRL,LRL,CS,XSFG,CSFG,ISE,NFGMB,
     M                 NDB,NDFG,NCB,NDBD,MPLOCW,MPLOCD,
     M                 XWFG,CWFG,MPLOC,MWLOC,MCB)
                ELSE
C ------ refine the mesh for diffusion zooming
                  CALL FPLUS1
     I                (MAXNP,MAXEL,MXNPFG,MXKGL,MXADNP,MXNCC,MXNPW,
     I                 NODE,EPSX,NCC,IDZOOM,N,M,XW(N,1),XW(N,2),
     I                 XW(N,3),X,CS,IE,XSFG,CSFG,NEL,NEFGS,ISE,NFGMB,
     O                 NDFG,XWFG,CWFG,MPLOC,MWLOC)
                ENDIF
 33           CONTINUE
 34         CONTINUE
 35       CONTINUE
C
C ----- DETERMINE THE DIFFUSION SUBELEMENTS FOR THE CURRENT TIME STEP
          NFGM(M)=NDFG
          CALL GRISED
     I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
     I         MAXNP,MAXEL,MXKBD,MXLSV,MXMSV,X,IE,MCON,II,
     I         XWFG,MXNPFG, EPSX,  LRL,NLRL,
     O         NEFGD,MWLOC,ISED,NMS,NLS,ILSV,IMSV)
C         CALL GRISED
C    I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
C    O         NEFGD,MWLOC,ISED)
C
        ELSEIF(NODE.EQ.6)THEN
          N=0
          DO 46 K=1,NZD+1
            ZI=2.0D0*DBLE(K-1)/DBLE(NZD)-1.0D0
            DO 45 I=NXD+1,1,-1
              DL1=DBLE(I-1)/DBLE(NXD)
              DO 44 J=NXD+1,1,-1
                DL2=DBLE(J-1)/DBLE(NXD)
                IF(DL1+DL2.GT.1.0D0 .AND. DABS(DL1+DL2-1.0D0) .GT. EPSX)
     >            GOTO 44
                DO 43 L=NXD+1,1,-1
                  DL3=DBLE(L-1)/DBLE(NXD)
                  IF(DABS(DL1+DL2+DL3-1.0D0).GT.EPSX)GOTO 43
                  N=N+1
                  CALL WARMSG(N,MXNPW,'DFPREP','MXNPW ',1)
                  CALL BASEXI
     I                (ZI,DL1,DL2,DL3,XX,YY,ZZ,NODE,
     O                 XW(N,1),XW(N,2),XW(N,3))
                  IF(K.EQ.1 .OR. K.EQ.NZD+1 .OR. DABS(DL1).LE.1.0D-6
     >               .OR. DABS(DL2).LE.1.0D-6 .OR. DABS(DL3).LE.1.0D-6)
     >               THEN
C ------------    CHECK THOSE POINTS ON THE BOUNDARY WITH BOTH S-F
C                 ELEMENTS AND THE INNERMOST NON-SF ELEMENTS
                    CALL GLBCHK
     I                (NCON1,N,NNP,NEL,NEFGS,NODE,EPSX,NCC,IDZOOM,
     I                 XW(N,1),XW(N,2),XW(N,3),MCON,MCON1,M,
     I                 NXD,NYD,NZD,MAXNP,MAXEL,MXNPFG,MXKGL,MXNDB,MXNCC,
     I                 MXKBD,MXNPW,MXADNP,K,I,J,DL1,DL2,DL3,DL4,
     I                 NFGM,X,IE,IB,NLRL,LRL,CS,XSFG,CSFG,ISE,NFGMB,
     M                 NDB,NDFG,NCB,NDBD,MPLOCW,MPLOCD,
     M                 XWFG,CWFG,MPLOC,MWLOC,MCB)
                  ELSE
C ------ refine the mesh for diffusion zooming
                    CALL FPLUS1
     I                  (MAXNP,MAXEL,MXNPFG,MXKGL,MXADNP,MXNCC,MXNPW,
     I                   NODE,EPSX,NCC,IDZOOM,N,M,XW(N,1),XW(N,2),
     I                   XW(N,3),X,CS,IE,XSFG,CSFG,NEL,NEFGS,ISE,NFGMB,
     O                   NDFG,XWFG,CWFG,MPLOC,MWLOC)
                  ENDIF
  43            CONTINUE
  44          CONTINUE
  45        CONTINUE
  46      CONTINUE
C
          NFGM(M)=NDFG
          CALL GRISED
     I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
     I         MAXNP,MAXEL,MXKBD,MXLSV,MXMSV,X,IE,MCON,II,
     I         XWFG,MXNPFG, EPSX,  LRL,NLRL,
     O         NEFGD,MWLOC,ISED,NMS,NLS,ILSV,IMSV)
C         CALL GRISED
C    I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
C    O         NEFGD,MWLOC,ISED)
        ELSE
          N=0
          DO 56 K=NXD+1,1,-1
            DL1=DBLE(K-1)/DBLE(NXD)
            DO 55 I=NXD+1,1,-1
              DL2=DBLE(I-1)/DBLE(NXD)
            IF(DL1+DL2.GT.1.0D0 .AND. DABS(DL1+DL2-1.0D0).GT.EPSX)
     >         GOTO 55
              DO 54 J=NXD+1,1,-1
                DL3=DBLE(J-1)/DBLE(NXD)
                IF(DL1+DL2+DL3.GT.1.0D0 .AND. DABS(DL1+DL2+DL3-1.0D0)
     >             .GT.EPSX)GOTO 54
                DO 53 L=NXD+1,1,-1
                  DL4=DBLE(L-1)/DBLE(NXD)
                  IF(DABS(DL1+DL2+DL3+DL4-1.0D0).GT.EPSX)GOTO 53
                  N=N+1
                  CALL WARMSG(N,MXNPW,'DFPREP','MXNPW ',1)
                  CALL BASEXI
     I                (DL1,DL2,DL3,DL4,XX,YY,ZZ,NODE,
     O                 XW(N,1),XW(N,2),XW(N,3))
C
                  IF(DABS(DL1).LE.1.0D-6 .OR. DABS(DL2).LE.1.0D-6 .OR.
     >               DABS(DL3).LE.1.0D-6 .OR. DABS(DL4).LE.1.0D-6)THEN
C ------------    CHECK THOSE POINTS ON THE BOUNDARY WITH BOTH S-F
C                 ELEMENTS AND THE INNERMOST NON-SF ELEMENTS
                    CALL GLBCHK
     I                (NCON1,N,NNP,NEL,NEFGS,NODE,EPSX,NCC,IDZOOM,
     I                 XW(N,1),XW(N,2),XW(N,3),MCON,MCON1,M,
     I                 NXD,NYD,NZD,MAXNP,MAXEL,MXNPFG,MXKGL,MXNDB,MXNCC,
     I                 MXKBD,MXNPW,MXADNP,I,J,K,DL1,DL2,DL3,DL4,
     I                 NFGM,X,IE,IB,NLRL,LRL,CS,XSFG,CSFG,ISE,NFGMB,
     M                 NDB,NDFG,NCB,NDBD,MPLOCW,MPLOCD,
     M                 XWFG,CWFG,MPLOC,MWLOC,MCB)
                  ELSE
                    CALL FPLUS1
     I                  (MAXNP,MAXEL,MXNPFG,MXKGL,MXADNP,MXNCC,MXNPW,
     I                   NODE,EPSX,NCC,IDZOOM,N,M,XW(N,1),XW(N,2),
     I                   XW(N,3),X,CS,IE,XSFG,CSFG,NEL,NEFGS,ISE,NFGMB,
     O                   NDFG,XWFG,CWFG,MPLOC,MWLOC)
                  ENDIF
  53            CONTINUE
  54          CONTINUE
  55        CONTINUE
  56      CONTINUE
          NFGM(M)=NDFG
          CALL GRISED
     I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
     I         MAXNP,MAXEL,MXKBD,MXLSV,MXMSV,X,IE,MCON,II,
     I         XWFG,MXNPFG, EPSX,  LRL,NLRL,
     O         NEFGD,MWLOC,ISED,NMS,NLS,ILSV,IMSV)
C         CALL GRISED
C    I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
C    O         NEFGD,MWLOC,ISED)
        ENDIF
C
  490 CONTINUE
c     NELD=NEFGD
ccc      print *,'NMS=',nms
ccc      print *,'NLS=',nls
C
C ----- CREATE ARRAY NLRND(NNP) FOR DIFFUSION
C
  501 CONTINUE
C
      DO I=1,NNP
        NLRND(I)=NLRN(I)
      ENDDO
      IF(IDZOOM.NE.0)THEN
        DO I=NNP+1,NDFG
          NLRND(I)=0
        ENDDO
        mxnlrn=0
        DO 505 M=1,NEFGD
          CALL ELENOD
     I        (ISED(M,5),ISED(M,7),
     O         NODE,I,I)
          DO 504 I=1,NODE
            IEM=ISED(M,I)
            DO 503 J=1,NODE
              IEMJ=ISED(M,J)
              IF(NLRND(IEM).EQ.0)THEN
                NLRND(IEM)=NLRND(IEM)+1
                LRN(1,IEM)=IEMJ
              ELSE
                DO K=1,NLRND(IEM)
                  N=LRN(K,IEM)
                  IF(N.EQ.IEMJ)GOTO 503
                ENDDO
                NLRND(IEM)=NLRND(IEM)+1
                NLRNN=NLRND(IEM)
                CALL WARMSG(NLRNN,MXJBD,'DFPREP','MXJBD ',1)
                mxnlrn=max0(mxnlrn,nlrnn)
                LRN(NLRNN,IEM)=IEMJ
              ENDIF
 503        CONTINUE
 504      CONTINUE
 505    CONTINUE
ccc        print *,'MXJBD=',mxnlrn
C
        IF(IPNTST.EQ.3)THEN
C
C ----- REARRANGE LRN IN ACENDING ORDER AND COMPUTE ND IF IPNTST=3
C
          DO 550 NP=1,NDFG
            NCOUNT=0
            IF(NP.LE.NNP)THEN
              I1=NLRN(NP)+1
            ELSE
              I1=1
            ENDIF
            DO 540 I=I1,NLRND(NP)
              NI=LRN(I,NP)
              IF(NCOUNT.EQ.0)THEN
                KK(1)=I
                KJ(1)=NI
                NCOUNT=NCOUNT+1
              ELSE
                DO J=NCOUNT,1,-1
                  NJ=KJ(J)
                  IF(J.EQ.1)THEN
                    IF(NI.LT.NJ)GOTO 510
                  ELSE
                    NK=KJ(J-1)
                    IF(NI.LT.NJ .AND. NI.GT.NK)GOTO 510
                  ENDIF
                ENDDO
                NCOUNT=NCOUNT+1
                CALL WARMSG(NCOUNT,100,'DFPREP','  100 ',2)
                KK(NCOUNT)=I
                KJ(NCOUNT)=NI
                GOTO 540
  510           CONTINUE
                NCOUNT=NCOUNT+1
                DO JJ=NCOUNT,J+1,-1
                  KK(JJ)=KK(JJ-1)
                  KJ(JJ)=KJ(JJ-1)
                ENDDO
                KK(J)=I
                KJ(J)=NI
              ENDIF
  540       CONTINUE
C
            DO I=I1,NLRND(NP)
              II1=I-I1+1
              LRN(I,NP)=KJ(II1)
              IF(NP.GT.NNP .AND. KJ(II1).EQ.NP)ND(NP)=I
            ENDDO
C
  550     CONTINUE
        ENDIF
C
        NFGM(1)=0
        DO M=2,NEL
          IF(IE(M-1,11).NE.0)THEN
            IF(IE(M,8).EQ.0)THEN
              IF(IE(M,6).EQ.0)THEN
                NDD=NXD*NXD*NXD
              ELSE
                NDD=NZD*NXD*NXD
              ENDIF
            ELSE
              NDD=NXD*NYD*NZD
            ENDIF
            NFGM(M)=NFGM(M-1)+NDD
          ELSE
            NFGM(M)=NFGM(M-1)
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
c
c
c
      SUBROUTINE BASEXI
     I     (DL1,DL2,DL3,DL4,XX,YY,ZZ,NODE,
     O      XW1,XW2,XW3)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),DL(8)
C
      IF(NODE.EQ.8)THEN
        DL(1)=0.125D0*(1.0D0-DL1)*(1.0D0-DL2)*(1.0D0-DL3)
        DL(2)=0.125D0*(1.0D0+DL1)*(1.0D0-DL2)*(1.0D0-DL3)
        DL(3)=0.125D0*(1.0D0+DL1)*(1.0D0+DL2)*(1.0D0-DL3)
        DL(4)=0.125D0*(1.0D0-DL1)*(1.0D0+DL2)*(1.0D0-DL3)
        DL(5)=0.125D0*(1.0D0-DL1)*(1.0D0-DL2)*(1.0D0+DL3)
        DL(6)=0.125D0*(1.0D0+DL1)*(1.0D0-DL2)*(1.0D0+DL3)
        DL(7)=0.125D0*(1.0D0+DL1)*(1.0D0+DL2)*(1.0D0+DL3)
        DL(8)=0.125D0*(1.0D0-DL1)*(1.0D0+DL2)*(1.0D0+DL3)
      ELSEIF(NODE.EQ.6)THEN
        DL(1)=DL2*(1.0D0-DL1)*0.5D0
        DL(2)=DL3*(1.0D0-DL1)*0.5D0
        DL(3)=DL4*(1.0D0-DL1)*0.5D0
        DL(4)=DL2*(1.0D0+DL1)*0.5D0
        DL(5)=DL3*(1.0D0+DL1)*0.5D0
        DL(6)=DL4*(1.0D0+DL1)*0.5D0
      ELSE
        DL(1)=DL1
        DL(2)=DL2
        DL(3)=DL3
        DL(4)=DL4
      ENDIF
C
      XW1=0.0D0
      XW2=0.0D0
      XW3=0.0D0
      DO IJ=1,NODE
        XW1=XW1+DL(IJ)*XX(IJ)
        XW2=XW2+DL(IJ)*YY(IJ)
        XW3=XW3+DL(IJ)*ZZ(IJ)
      ENDDO
      RETURN
      END
c
c
c
      SUBROUTINE GLBCHK
     I    (NCON1,N,NNP,NEL,NEFGS,NODE,EPSX,NCC,IDZOOM,XW1,XW2,XW3,MCON,
     I     MCON1,M,NXD,NYD,NZD,MAXNP,MAXEL,MXNPFG,MXKGL,MXNDB,MXNCC,
     I     MXKBD,MXNPW,MXADNP,L1,L2,L3,DL1,DL2,DL3,DL4,NFGM,X,IE,IB,
     I     NLRL,LRL,CS,XSFG,CSFG,ISE,NFGMB,
     M     NDB,NDFG,NCB,NDBD,MPLOCW,MPLOCD,XWFG,CWFG,MPLOC,MWLOC,MCB)
C
C ----------------------------------------------------------------------
C ***** THIS ROUTINE CHECKS THE GENERATED FINE GRID WHICH IS LOCATED
C ***** ON GLOBAL NODES TO SEE IF THEY ARE ALSO ON THE BOUNDARY OF
C ***** ROUGH REGION
C ----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP,3),CS(MAXNP,MXNCC),IE(MAXEL,11),IB(MAXNP)
      DIMENSION NLRL(MAXNP),LRL(MXKBD,MAXNP),ISE(MXKGL,8),NFGMB(MAXEL)
      DIMENSION NFGM(MAXEL),XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC)
      DIMENSION NDBD(MXNDB),MPLOCW(MXNPFG),MPLOCD(MXNDB),MPLOC(MXNPFG)
      DIMENSION XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC),MWLOC(MXNPW)
      DIMENSION MCON(6),MCB(3),CQ(7),MCON1(100)
C
      DO 31 KI=1,NCON1
        MP=MCON1(KI)
        IF(MP.EQ.0)GOTO 31
        IF(IE(MP,11).NE.0)THEN
          IF(MP.EQ.1)THEN
            NM1=NNP
            NM2=NFGM(MP)
          ELSE
            NM1=NFGM(MP-1)
            NM2=NFGM(MP)
          ENDIF
          IF(NM2.NE.0)THEN
            DO KK=NM1+1,NM2
              IF(DABS(XW1-XWFG(KK,1)).LE.EPSX .AND.
     >           DABS(XW2-XWFG(KK,2)).LE.EPSX .AND.
     >           DABS(XW3-XWFG(KK,3)).LE.EPSX)THEN
C ----- MWLOC(N)=KK MEANS THE N-TH DIFFUSION FINE GRID IN ELEMENT M IS
C       THE KK-TH DIFFUSION POINT OF ALL DIFFUSION POINTS
                MWLOC(N)=KK
                GOTO 33
              ENDIF
            ENDDO
          ENDIF
        ENDIF
   31 CONTINUE
C
      DO JI=1,NODE
        IEMP=IE(M,JI)
        IF(IEMP.NE.0)THEN
          IF(DABS(XW1-X(IEMP,1)).LE.EPSX .AND.
     >       DABS(XW2-X(IEMP,2)).LE.EPSX .AND.
     >       DABS(XW3-X(IEMP,3)).LE.EPSX)THEN
C ------- the n-th diffusion point in the M-th element is located on
C         global node
            MWLOC(N)=IEMP
            IF(IB(IEMP).EQ.0)THEN
              DO KM=1,NLRL(IEMP)
                MK=LRL(KM,IEMP)
                IF(IE(MK,11).NE.0 .AND. IE(MK,11).NE.-1)GOTO 33
              ENDDO
            ENDIF
C
C ------- boundary point or the outest layer of rough region
            DO KM=1,NLRL(IEMP)
              MK=LRL(KM,IEMP)
C ------- this element has already been taken care before
              IF(MK.LT.M .AND. IE(MK,11).NE.0)GOTO 33
            ENDDO
C
            NDB=NDB+1
            CALL WARMSG(NDB,MXNDB,'GLBCHK','MXNDB ',1)
C ------ store the diffusion boundary point
            NDBD(NDB)=IEMP
            MPLOCW(NDB)=M
C ------- store the plane on which the diffusion boundary point fallen
            IF(NODE.EQ.8)THEN
              IF(JI.EQ.1 .OR. JI.EQ.2 .OR. JI.EQ.3 .OR.JI.EQ.4)
     >           MPLOCD(NDB)=5
              IF(JI.EQ.5 .OR. JI.EQ.6 .OR. JI.EQ.7 .OR.JI.EQ.8)
     >           MPLOCD(NDB)=6
            ELSEIF(NODE.EQ.6)THEN
              IF(JI.EQ.1 .OR. JI.EQ.2 .OR. JI.EQ.3)MPLOCD(NDB)=4
              IF(JI.EQ.4 .OR. JI.EQ.5 .OR. JI.EQ.6)MPLOCD(NDB)=5
            ELSE
              IF(JI.EQ.1)MPLOCD(NDB)=2
              IF(JI.EQ.2 .OR. JI.EQ.3 .OR. JI.EQ.4)MPLOCD(NDB)=1
            ENDIF
            GOTO 33
C
          ENDIF
        ENDIF
      ENDDO
      NDFG=NDFG+1
      CALL WARMSG(NDFG,MXADNP,'GLBCHK','MXADNP',1)
      MWLOC(N)=NDFG
C ----- DETERMINE THE CONCENTRATION OF THIS DIFFUSION FINE GRID BY
C       INTERPOLATION FROM CS AND CSFG
      CALL INTERP
     I    (MAXNP,MAXEL,MXNPFG,MXKGL,8,MXNCC,NODE,M,XW1,XW2,XW3,X,CS,
     I     IE,XSFG,CSFG,NEL,NEFGS,IDZOOM,11,ISE,NFGMB,1,NCC,EPSX,
     O     CQ)
      XWFG(NDFG,1)=XW1
      XWFG(NDFG,2)=XW2
      XWFG(NDFG,3)=XW3
      DO K=1,NCC
        CWFG(NDFG,K)=CQ(K)
      ENDDO
      MPLOC(NDFG)=M
C
C ----- CHECK THOSE POINTS ON THE BOUNDARY INTERSECTING SF & NON-SF
C       ELEMENT OR ON THE GLOBAL BOUNDARIES
C
      NCB=0
c     MP1=-10
c     MP2=-10
c     MP3=-10
      IEMP1=-100
      IEMP2=-100
      IEMP3=-100
C
      IF(NODE.EQ.8)THEN
        IF(L1.EQ.1)THEN
          NCB=NCB+1
          MCB(NCB)=2
        ENDIF
        IF(L1.EQ.NYD+1)THEN
          NCB=NCB+1
          MCB(NCB)=4
        ENDIF
        IF(L2.EQ.1)THEN
          NCB=NCB+1
          MCB(NCB)=1
        ENDIF
        IF(L2.EQ.NXD+1)THEN
          NCB=NCB+1
          MCB(NCB)=3
        ENDIF
        IF(L3.EQ.1)THEN
          NCB=NCB+1
          MCB(NCB)=5
        ENDIF
        IF(L3.EQ.NZD+1)THEN
          NCB=NCB+1
          MCB(NCB)=6
        ENDIF
      ELSEIF(NODE.EQ.6)THEN
        IF(L1.EQ.1)THEN
          NCB=NCB+1
          MCB(NCB)=4
        ENDIF
        IF(L1.EQ.NZD+1)THEN
         NCB=NCB+1
         MCB(NCB)=5
        ENDIF
        IF(DABS(DL1).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=2
        ENDIF
        IF(DABS(DL2).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=3
        ENDIF
        IF(DABS(DL3).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=1
        ENDIF
      ELSE
        IF(DABS(DL1).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=4
        ENDIF
        IF(DABS(DL2).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=2
        ENDIF
        IF(DABS(DL3).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=3
        ENDIF
        IF(DABS(DL4).LE.1.0D-6)THEN
          NCB=NCB+1
          MCB(NCB)=1
        ENDIF
      ENDIF
C
      IF(NCB.GE.1)THEN
        KK=MCB(1)
        MP=MCON(KK)
        IF(MP.EQ.0)GOTO 32
        IEMP1=IE(MP,11)
        IF(NCB.GE.2)THEN
          KK=MCB(2)
          MP=MCON(KK)
          IF(MP.EQ.0)GOTO 32
          IEMP2=IE(MP,11)
          IF(NCB.GE.3)THEN
            KK=MCB(3)
            MP=MCON(KK)
            IF(MP.EQ.0)GOTO 32
            IEMP3=IE(MP,11)
            IF(IEMP3.EQ.0)GOTO 32
          ENDIF
          IF(IEMP2.EQ.0)GOTO 32
        ENDIF
        IF(IEMP1.EQ.0)GOTO 32
      ENDIF
      GOTO 33
C
  32  CONTINUE
      NDB=NDB+1
      CALL WARMSG(NDB,MXNDB,'GLBCHK','MXNDB ',3)
      NDBD(NDB)=NDFG
      MPLOCW(NDB)=M
      MPLOCD(NDB)=KK
  33  RETURN
      END
c
c
c
      SUBROUTINE FPLUS1
     I     (MAXNP,MAXEL,MXNPFG,MXKGL,MXADNP,MXNCC,MXNPW,NODE,EPSX,NCC,
     I      IDZOOM,N,M,XQ,YQ,ZQ,X,CS,IE,XSFG,CSFG,NEL,NEFGS,ISE,NFGMB,
     O      NDFG,XWFG,CWFG,MPLOC,MWLOC)
C
C ---------------------------------------------------------------------
C ***** THIS ROUTINE STORES THE GENERATED FINE GRID FOR DIFFUSION
C ***** ZOOMING.  THE ASSOCIATED INFORMATION IS STORED IN XWFG,CWFG,
C ***** MPLOC AND MWLOC ARRAYS
C ---------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP,3),CS(MAXNP,MXNCC),IE(MAXEL,11),NFGMB(MAXEL)
      DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC),ISE(MXKGL,8)
      DIMENSION XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC),MPLOC(MXNPFG)
      DIMENSION MWLOC(MXNPW),CQ(7)
C
      NDFG=NDFG+1
      CALL WARMSG(NDFG,MXADNP,'FPLUS1','MXADNP',1)
      MWLOC(N)=NDFG
C ----- DETERMINE THE CONCENTRATION OF THIS DIFFUSION FINE GRID BY
C       INTERPOLATION FROM CS AND CSFG
      CALL INTERP
     I    (MAXNP,MAXEL,MXNPFG,MXKGL,8,MXNCC,NODE,M,XQ,YQ,ZQ,X,CS,
     I     IE,XSFG,CSFG,NEL,NEFGS,IDZOOM,11,ISE,NFGMB,1,NCC,EPSX,
     O     CQ)
      XWFG(NDFG,1)=XQ
      XWFG(NDFG,2)=YQ
      XWFG(NDFG,3)=ZQ
      DO K=1,NCC
        CWFG(NDFG,K)=CQ(K)
      ENDDO
      MPLOC(NDFG)=M
      RETURN
      END
C
c
c
      SUBROUTINE GRISED
     I        (NXD,NYD,NZD,MXKGLD,MXNPW,NODE,M,
     I         MAXNP,MAXEL,MXKBD,MXLSV,MXMSV,X,IE,MCON,II,
     I         XWFG,MXNPFG, EPSX,  LRL,NLRL,
     O         NEFGD,MWLOC,ISED,NMS,NLS,ILSV,IMSV)
C 11/23/93
C ----------------------------------------------------------------------
C ***** THIS ROUTINE GENERATES ISED ARRAY AND STORE MWLOC INFORMATION
C ----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION MWLOC(MXNPW),ISED(MXKGLD,9)
      DIMENSION IE(MAXEL,11),MCON(6),X(MAXNP,3)
      DIMENSION ILSV(MXLSV,5),IMSV(MXMSV,4)
      DIMENSION KID(4,6,3),CID(6,3),MSI(4)
      DIMENSION XWFG(MXNPFG,3),LRL(MXKBD,MAXNP),NLRL(MAXNP)
C
      DATA KID /1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >          1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >          4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
      DATA CID /1.0D0, -1.0D0, -1.0D0, 1.0D0,  1.0D0, -1.0D0,
     >          1.0D0,  1.0D0,  1.0D0, 1.0D0, -1.0D0,  0.0D0,
     >          1.0D0,  1.0D0,  1.0D0, 1.0D0,  0.0D0,  0.0D0/
C
      IF(NODE.EQ.8)THEN
        DO 38 K=1,NZD
          DO 37 I=1,NYD
            DO J=1,NXD
              MD=(K-1)*NXD*NYD+(I-1)*NXD+J+NEFGD
              N11=(K-1)*(NXD+1)*(NYD+1)+(I-1)*(NXD+1)
              N22=N11+(NXD+1)*(NYD+1)
              ISED(MD,1)=MWLOC(N11+J)
              ISED(MD,2)=MWLOC(N11+J+1)
              ISED(MD,3)=MWLOC(N11+NXD+1+J+1)
              ISED(MD,4)=MWLOC(N11+NXD+1+J)
              ISED(MD,5)=MWLOC(N22+J)
              ISED(MD,6)=MWLOC(N22+J+1)
              ISED(MD,7)=MWLOC(N22+NXD+1+J+1)
              ISED(MD,8)=MWLOC(N22+NXD+1+J)
              ISED(MD,9)=M
            ENDDO
   37     CONTINUE
   38   CONTINUE
C
C ----- PREPARE INFORMATION FOR INTRA-BOUNDARY SIDES
C
c       IF(IFLUX.EQ.0)GOTO 46
        DO 45 MS=1,II
          IF(MCON(MS).EQ.0 .OR. IE(M,11).NE.-1)GOTO 45
          NMS=NMS+1
          CALL WARMSG(NMS,MXMSV,'GRISED','MXMSV ',1)
          DO KK=1,4
            MSI(KK)=KID(KK,MS,1)
            IF(MSI(KK).NE.0)THEN
              IMSV(NMS,KK)=IE(M,MSI(KK))
            ELSE
              IMSV(NMS,KK)=0
            ENDIF
          ENDDO
C
C ----- CHECK IF THIS INTRA-BOUNDARY IS A GLOBAL BOUNDARY
C
          N1=IMSV(NMS,1)
          N2=IMSV(NMS,2)
          N3=IMSV(NMS,3)
          CALL CKCNEL
     I             (LRL,NLRL,M,N1,N2,N3,MAXNP,MXKBD,
     O              ME1,KOUNT)
          IF(KOUNT.EQ.1)THEN
            NMS=NMS-1
            GOTO 45
          ENDIF
C
          IF(IE(ME1,11).NE.0)THEN
            NMS=NMS-1
            GOTO 45
          ENDIF
C
c         NI=IMSV(NMS,1)
c         NJ=IMSV(NMS,2)
c         A1=X(NJ,1)-X(NI,1)
c         A2=X(NJ,2)-X(NI,2)
c         A3=X(NJ,3)-X(NI,3)
c         NK=IMSV(NMS,3)
c         B1=X(NK,1)-X(NI,1)
c         B2=X(NK,2)-X(NI,2)
c         B3=X(NK,3)-X(NI,3)
c         AB23=A2*B3-A3*B2
c         AB31=A3*B1-A1*B3
c         AB12=A1*B2-A2*B1
c         AREA=DSQRT(AB23*AB23+AB31*AB31+AB12*AB12)
c         DMSV(NMS,1)=AB23/AREA*CID(MS,1)
c         DMSV(NMS,2)=AB31/AREA*CID(MS,1)
c         DMSV(NMS,3)=AB12/AREA*CID(MS,1)
C
          K1=1
          K2=NZD
          I1=1
          I2=NYD
          J1=1
          J2=NXD
          IF(MS.EQ.1)J2=1
          IF(MS.EQ.2)I2=1
          IF(MS.EQ.3)J1=NXD
          IF(MS.EQ.4)I1=NYD
          IF(MS.EQ.5)K2=1
          IF(MS.EQ.6)K1=NZD
          DO 44 K=K1,K2
            DO 43 I=I1,I2
              DO 42 J=J1,J2
                MD=(K-1)*NXD*NYD+(I-1)*NXD+J+NEFGD
                NLS=NLS+1
                CALL WARMSG(NLS,MXLSV,'GRISED','MXLSV ',1)
                DO KK=1,4
                  IF(MSI(KK).NE.0)THEN
                    ILSV(NLS,KK)=ISED(MD,MSI(KK))
                  ELSE
                    ILSV(NLS,KK)=0
                  ENDIF
                ENDDO
                ILSV(NLS,5)=NMS
c               DO KK=1,3
c                 DLSV(NLS,KK)=-DMSV(NMS,KK)
c               ENDDO
  42          CONTINUE
  43        CONTINUE
  44      CONTINUE
  45    CONTINUE
C
c 46    CONTINUE
        NEFGD=NEFGD+NXD*NYD*NZD
C
      ELSEIF(NODE.EQ.6)THEN
        N11=(NXD+1)*(NXD+2)/2
        DO 49 K=1,NZD
          NC=0
          NI=0
          NK=(K-1)*NXD**2
          NKI=(K-1)*(NXD+1)*(NXD+2)/2
          DO 48 I=1,NXD
            NT=2*I-1
            NII=NI
            DO 47 J=NC+1,NC+NT,2
              MW=J+NEFGD+NK
              NI=NI+1
              NJ=NI+NKI
              ISED(MW,1)=MWLOC(NJ)
              ISED(MW,2)=MWLOC(NJ+I)
              ISED(MW,3)=MWLOC(NJ+I+1)
              ISED(MW,4)=MWLOC(NJ+N11)
              ISED(MW,5)=MWLOC(NJ+I+N11)
              ISED(MW,6)=MWLOC(NJ+I+1+N11)
              ISED(MW,7)=0
              ISED(MW,8)=0
              ISED(MW,9)=M
  47        CONTINUE
            IF(I.GT.1)THEN
              DO J=NC+2,NC+NT-1,2
                NII=NII+1
                NIJ=NII+NKI
                MW=J+NEFGD+NK
                ISED(MW,1)=MWLOC(NIJ)
                ISED(MW,2)=MWLOC(NIJ+I+1)
                ISED(MW,3)=MWLOC(NIJ+1)
                ISED(MW,4)=MWLOC(NIJ+N11)
                ISED(MW,5)=MWLOC(NIJ+I+1+N11)
                ISED(MW,6)=MWLOC(NIJ+1+N11)
                ISED(MW,7)=0
                ISED(MW,8)=0
                ISED(MW,9)=M
              ENDDO
            ENDIF
            NC=I**2
  48      CONTINUE
  49    CONTINUE
C
C ----- PREPARE INFORMATION FOR INTRA-BOUNDARY SIDES
C
c       IF(IFLUX.EQ.0)GOTO 56
        DO 55 MS=1,II
          IF(MCON(MS).NE.0 .OR. IE(M,11).NE.-1)GOTO 55
          NMS=NMS+1
          CALL WARMSG(NMS,MXMSV,'GRISED','MXMSV ',2)
          DO KK=1,4
            MSI(KK)=KID(KK,MS,2)
            IF(MSI(KK).NE.0)THEN
              IMSV(NMS,KK)=IE(M,MSI(KK))
            ELSE
              IMSV(NMS,KK)=0
            ENDIF
          ENDDO
C
C ----- CHECK IF THIS INTRA-BOUNDARY IS A GLOBAL BOUNDARY
C
          N1=IMSV(NMS,1)
          N2=IMSV(NMS,2)
          N3=IMSV(NMS,3)
          CALL CKCNEL
     I             (LRL,NLRL,M,N1,N2,N3,MAXNP,MXKBD,
     O              ME1,KOUNT)
          IF(KOUNT.EQ.1)THEN
            NMS=NMS-1
            GOTO 55
          ENDIF
C
c         NI=IMSV(NMS,1)
c         NJ=IMSV(NMS,2)
c         A1=X(NJ,1)-X(NI,1)
c         A2=X(NJ,2)-X(NI,2)
c         A3=X(NJ,3)-X(NI,3)
c         NK=IMSV(NMS,3)
c         B1=X(NK,1)-X(NI,1)
c         B2=X(NK,2)-X(NI,2)
c         B3=X(NK,3)-X(NI,3)
c         AB23=A2*B3-A3*B2
c         AB31=A3*B1-A1*B3
c         AB12=A1*B2-A2*B1
c         AREA=DSQRT(AB23*AB23+AB31*AB31+AB12*AB12)
c         DMSV(NMS,1)=AB23/AREA*CID(MS,2)
c         DMSV(NMS,2)=AB31/AREA*CID(MS,2)
c         DMSV(NMS,3)=AB12/AREA*CID(MS,2)
C
          K1=1
          K2=NZD
          I1=1
          I2=NXD
          IF(MS.EQ.4)K2=1
          IF(MS.EQ.5)K1=NZD
          IF(MS.EQ.3)I1=NXD
          DO 54 K=K1,K2
            NC=0
            NK=(K-1)*NXD**2
            NKI=(K-1)*(NXD+1)*(NXD+2)/2
            DO 53 I=I1,I2
              NT=2*I-1
              NC=(I-1)*(I-1)
              J1=NC+1
              J2=NC+NT
              IF(MS.EQ.1)THEN
                J1=NC+NT
                J2=J1
              ELSEIF(MS.EQ.2)THEN
                J1=NC+1
                J2=J1
              ENDIF
              DO 52 J=J1,J2,2
                MD=J+NEFGD+NK
                NLS=NLS+1
                CALL WARMSG(NLS,MXLSV,'GRISED','MXLSV ',2)
                DO KK=1,4
                  IF(MSI(KK).NE.0)THEN
                    ILSV(NLS,KK)=ISED(MD,MSI(KK))
                  ELSE
                    ILSV(NLS,KK)=0
                  ENDIF
                ENDDO
                ILSV(NLS,5)=NMS
c               DO KK=1,3
c                 DLSV(NLS,KK)=-DMSV(NMS,KK)
c               ENDDO
  52          CONTINUE
              IF(I.GT.1 .AND. MS.GE.4)THEN
                DO J=NC+2,NC+NT-1,2
                  MD=J+NEFGD+NK
                  NLS=NLS+1
                  DO KK=1,4
                    IF(MSI(KK).NE.0)THEN
                      ILSV(NLS,KK)=ISED(MD,MSI(KK))
                    ELSE
                      ILSV(NLS,KK)=0
                    ENDIF
                  ENDDO
                  ILSV(NLS,5)=NMS
                ENDDO
              ENDIF
  53        CONTINUE
  54      CONTINUE
  55    CONTINUE
C
c 56    CONTINUE
        NEFGD=NEFGD+NXD*NXD*NZD
C
      ELSE
        J2=1
        J1=0
        DO 58 L=1,NXD
          JUMP1=L
          JUMP2=(L-1)*L/2
          J1=J1+JUMP1
          J2=J2+JUMP2
c         NCOUNT=L**3-(L-1)**3
          M1=(L-1)**3+1+NEFGD
          ISED(M1,1)=MWLOC(J2)
          ISED(M1,2)=MWLOC(J2+J1)
          ISED(M1,3)=MWLOC(J2+J1+1)
          ISED(M1,4)=MWLOC(J2+J1+2)
          ISED(M1,9)=M
C
          IF(L.EQ.1)GOTO 58
          J3=J2
          DO 57 L1=2,L
            M2=(L1-2)*(L1-1)*3+2+(L-1)**3+NEFGD
            J3=J3+L1-1
            ISED(M2,1)=MWLOC(J3)
            ISED(M2,2)=MWLOC(J3+J1)
            ISED(M2,3)=MWLOC(J3+J1+L1)
            ISED(M2,4)=MWLOC(J3+J1+L1+1)
            ISED(M2,9)=M
            ISED(M2+1,1)=MWLOC(J3)
            ISED(M2+1,2)=MWLOC(J3+J1+1)
            ISED(M2+1,3)=MWLOC(J3-L1+1)
            ISED(M2+1,4)=MWLOC(J3+J1)
            ISED(M2+1,9)=M
            ISED(M2+2,1)=MWLOC(J3)
            ISED(M2+2,2)=MWLOC(J3+J1+1)
            ISED(M2+2,3)=MWLOC(J3+J1)
            ISED(M2+2,4)=MWLOC(J3+J1+L1+1)
            ISED(M2+2,9)=M
            ISED(M2+3,1)=MWLOC(J3)
            ISED(M2+3,2)=MWLOC(J3+J1+1)
            ISED(M2+3,3)=MWLOC(J3+J1+L1+1)
            ISED(M2+3,4)=MWLOC(J3+1)
            ISED(M2+3,9)=M
            ISED(M2+4,1)=MWLOC(J3)
            ISED(M2+4,2)=MWLOC(J3+J1+1)
            ISED(M2+4,3)=MWLOC(J3+1)
            ISED(M2+4,4)=MWLOC(J3-L1+1)
            ISED(M2+4,9)=M
            ISED(M2+5,1)=MWLOC(J3+1)
            ISED(M2+5,2)=MWLOC(J3+J1+1)
            ISED(M2+5,3)=MWLOC(J3+J1+L1+1)
            ISED(M2+5,4)=MWLOC(J3+J1+L1+2)
            ISED(M2+5,9)=M
            IF(L.GE.3)THEN
              M5=M2+5
              J4=J3
              DO L2=1,L-2
                J4=J4+1
                M5=M5+1
                ISED(M5,1)=MWLOC(J4)
                ISED(M5,2)=MWLOC(J4-L1+1)
                ISED(M5,3)=MWLOC(J4-L1)
                ISED(M5,4)=MWLOC(J4+J1)
                ISED(M5,9)=M
                M5=M5+1
                ISED(M5,1)=MWLOC(J4)
                ISED(M5,2)=MWLOC(J4+J1+1)
                ISED(M5,3)=MWLOC(J4-L1+1)
                ISED(M5,4)=MWLOC(J4+J1)
                ISED(M5,9)=M
                M5=M5+1
                ISED(M5,1)=MWLOC(J4)
                ISED(M5,2)=MWLOC(J4+J1+1)
                ISED(M5,3)=MWLOC(J4+J1)
                ISED(M5,4)=MWLOC(J4+J1+L1+1)
                ISED(M5,9)=M
                M5=M5+1
                ISED(M5,1)=MWLOC(J4)
                ISED(M5,2)=MWLOC(J4+J1+1)
                ISED(M5,3)=MWLOC(J4+J1+L1+1)
                ISED(M5,4)=MWLOC(J4+1)
                ISED(M5,9)=M
                M5=M5+1
                ISED(M5,1)=MWLOC(J4)
                ISED(M5,2)=MWLOC(J4+J1+1)
                ISED(M5,3)=MWLOC(J4+1)
                ISED(M5,4)=MWLOC(J4-L1+1)
                ISED(M5,9)=M
                M5=M5+1
                ISED(M5,1)=MWLOC(J4+1)
                ISED(M5,2)=MWLOC(J4+J1+1)
                ISED(M5,3)=MWLOC(J4+J1+L1)
                ISED(M5,4)=MWLOC(J4+J1+L1+2)
                ISED(M5,9)=M
              ENDDO
            ENDIF
  57      CONTINUE
  58    CONTINUE
        DO 59 I1=1,NXD**3
          II1=I1+NEFGD
          DO I2=5,8
            ISED(II1,I2)=0
          ENDDO
  59    CONTINUE
C
C ----- PREPARE INFORMATION FOR INTRA-BOUNDARY SIDES
C
c       IF(IFLUX.EQ.0)GOTO 66
        DO 65 MS=1,II
          IF(MCON(MS).NE.0 .OR. IE(M,11).NE.-1)GOTO 65
          NMS=NMS+1
          CALL WARMSG(NMS,MXMSV,'GRISED','MXMSV ',3)
          DO KK=1,4
            MSI(KK)=KID(KK,MS,3)
            IF(MSI(KK).NE.0)THEN
              IMSV(NMS,KK)=IE(M,MSI(KK))
            ELSE
              IMSV(NMS,KK)=0
            ENDIF
          ENDDO
C
C ----- CHECK IF THIS INTRA-BOUNDARY IS A GLOBAL BOUNDARY
C
          N1=IMSV(NMS,1)
          N2=IMSV(NMS,2)
          N3=IMSV(NMS,3)
          CALL CKCNEL
     I             (LRL,NLRL,M,N1,N2,N3,MAXNP,MXKBD,
     O              ME1,KOUNT)
          IF(KOUNT.EQ.1)THEN
            NMS=NMS-1
            GOTO 65
          ENDIF
C
          NI=IMSV(NMS,1)
          NJ=IMSV(NMS,2)
          A1=X(NJ,1)-X(NI,1)
          A2=X(NJ,2)-X(NI,2)
          A3=X(NJ,3)-X(NI,3)
          NK=IMSV(NMS,3)
          B1=X(NK,1)-X(NI,1)
          B2=X(NK,2)-X(NI,2)
          B3=X(NK,3)-X(NI,3)
          AB23=A2*B3-A3*B2
          AB31=A3*B1-A1*B3
          AB12=A1*B2-A2*B1
          AREA=DSQRT(AB23*AB23+AB31*AB31+AB12*AB12)
          DMSV1=AB23/AREA*CID(MS,3)
          DMSV2=AB31/AREA*CID(MS,3)
          DMSV3=AB12/AREA*CID(MS,3)
C
          DO 64 MK=1,NXD*NXD*NXD
            MD=MK+NEFGD
            DO 63 LS=1,4
              DO KK=1,4
                MSI(KK)=KID(KK,LS,3)
                IF(MSI(KK).NE.0)THEN
                  IKK=ISED(MD,MSI(KK))
                ELSE
                  IKK=0
                ENDIF
                IF(IKK.NE.0)THEN
                  XP=XWFG(IKK,1)-X(NI,1)
                  YP=XWFG(IKK,2)-X(NI,2)
                  ZP=XWFG(IKK,3)-X(NI,3)
                  CANG=XP*DMSV1+YP*DMSV2+ZP*DMSV3
                  IF(DABS(CANG).GT.EPSX)GOTO 63
                ENDIF
              ENDDO
              NLS=NLS+1
              CALL WARMSG(NLS,MXLSV,'GRISED','MXLSV ',3)
              DO KK=1,4
                IF(MSI(KK).NE.0)THEN
                  ILSV(NLS,KK)=ISED(MD,MSI(KK))
                ELSE
                  ILSV(NLS,KK)=0
                ENDIF
              ENDDO
              ILSV(NLS,5)=NMS
c             DO KK=1,3
c               DLSV(NLS,KK)=-DMSV(NMS,KK)
c             ENDDO
              GOTO 64
  63        CONTINUE
  64      CONTINUE
  65    CONTINUE
C
c 66    CONTINUE
        NEFGD=NEFGD+NXD*NXD*NXD
C
      ENDIF
      CALL WARMSG(NEFGD,MXKGLD,'FPLUS1','MXKGLD',2)
      RETURN
      END
c
c
c
      SUBROUTINE NEWXE(DETJ,XSIO,ETAO,DFX,DFE,DGX,DGE,XID,EID,FOF,GOF,
     O                 SCALAR,  XSI,ETA)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
        IF(DETJ.EQ.0.0D0)THEN
          IF(DFX.EQ.0.0D0 .AND. DGX.EQ.0.0D0)THEN
            IF(DFE.EQ.0.0D0 .AND. DGE.NE.0.0D0)THEN
              ETA=ETAO-GOF/DGE
            ELSEIF(DFE.NE.0.0D0 .AND. DGE.EQ.0.0D0)THEN
              ETA=ETAO-FOF/DFE
            ELSEIF(DFE.NE.0.0D0 .AND. DGE.NE.0.0D0)THEN
              ETA=ETAO-0.5D0*(FOF/DFE+GOF/DGE)
            ELSE
              ETA=ETAO+EID*SCALAR
            ENDIF
            XSI=XSIO
          ELSEIF(DFE.EQ.0.0D0 .AND. DGE.EQ.0.0D0)THEN
            IF(DFX.EQ.0.0D0 .AND. DGX.NE.0.0D0)THEN
              XSI=XSIO-GOF/DGX
            ELSEIF(DFX.NE.0.0D0 .AND. DGX.EQ.0.0D0)THEN
              XSI=XSIO-FOF/DFX
            ELSEIF(DFX.NE.0.0D0 .AND. DGX.NE.0.0D0)THEN
              XSI=XSIO-0.5D0*(FOF/DFX+GOF/DGX)
            ELSE
              XSI=XSIO+XID*SCALAR
            ENDIF
            ETA=ETAO
          ELSE
            XSI=XSIO+XID*SCALAR
            ETA=ETAO+EID*SCALAR
          ENDIF
        ELSE
          XSI=XSIO-(DGE*FOF-DFE*GOF)/DETJ
          ETA=ETAO-(-DGX*FOF+DFX*GOF)/DETJ
        ENDIF
C
        RETURN
        END
C
c
c
      SUBROUTINE IBWCK4
     I          (IBWJ1,IBWJ2,XP,YP,ZP,VPX,VPY,VPZ,NXW,M,IE,SDT,MAXEL,
     O           N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
C       10/10/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,9)
C
      IF(IBWJ1*IBWJ2.NE.0)THEN
        IF(NXW.EQ.1)THEN
          IBSUM=IBWJ1+IBWJ2
          IF(IBSUM.EQ.256)GOTO 160
          IF(IBSUM.EQ.267)GOTO 160
          IF(IBSUM.EQ.275)GOTO 160
          IF(IBSUM.EQ.358)GOTO 170
          IF(IBSUM.EQ.366)GOTO 170
          IF(IBSUM.EQ.377)GOTO 180
          WRITE(16,*)'SOMETHING WRONG WITH IBWJ1,IBWJ2'
          GOTO 240
        ELSE
          IBPROD=IBWJ1*IBWJ2
          IF(IBWJ1.EQ.1 .OR. IBWJ2.EQ.1)THEN
            GOTO 160
          ELSEIF(IBWJ1.EQ.2 .OR. IBWJ2.EQ.2)THEN
            GOTO 170
          ELSEIF(IBWJ1.EQ.3 .OR. IBWJ2.EQ.3)THEN
            GOTO 180
          ELSEIF(IBWJ1.EQ.4 .OR. IBWJ2.EQ.4)THEN
            GOTO 190
          ELSEIF(IBPROD.EQ.782 .OR. IBPROD.EQ.966 .OR. IBPROD.EQ.1428)
     >           THEN
            GOTO 160
          ELSEIF(IBPROD.EQ.182 .OR. IBPROD.EQ.442 .OR. IBPROD.EQ.476)
     >           THEN
            GOTO 170
          ELSEIF(IBPROD.EQ.168 .OR. IBPROD.EQ.504 .OR. IBPROD.EQ.588)
     >           THEN
            GOTO 180
          ELSEIF(IBPROD.EQ.156 .OR. IBPROD.EQ.276 .OR. IBPROD.EQ.299)
     >           THEN
            GOTO 190
          ELSEIF(IBWJ1.EQ.34 .OR. IBWJ2.EQ.34)THEN
            GOTO 160
          ELSEIF(IBWJ1.EQ.42 .OR. IBWJ2.EQ.42)THEN
            GOTO 160
          ELSEIF(IBWJ1.EQ.23 .OR. IBWJ2.EQ.23)THEN
            GOTO 160
          ELSEIF(IBWJ1.EQ.13 .OR. IBWJ2.EQ.13)THEN
            GOTO 170
          ELSEIF(IBWJ1.EQ.14 .OR. IBWJ2.EQ.14)THEN
            GOTO 170
          ELSE
            GOTO 180
          ENDIF
        ENDIF
C
 160    CONTINUE
        N1=IE(M,2)
        N2=IE(M,3)
        N3=IE(M,4)
        N4=0
        GOTO 220
C
 170    CONTINUE
        N1=IE(M,1)
        N2=IE(M,4)
        N3=IE(M,3)
        N4=0
        GOTO 220
C
 180    CONTINUE
        N1=IE(M,4)
        N2=IE(M,1)
        N3=IE(M,2)
        N4=0
        GOTO 220
C
 190    CONTINUE
        N1=IE(M,1)
        N2=IE(M,3)
        N3=IE(M,2)
        N4=0
C
 220    CONTINUE
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
  240 CONTINUE
C
      WRITE(16,*)'ERROR MESSAGE AT IBWCK4: The Following Point Can ',
     > 'Not be Tracked - STOP'
      WRITE(16,1002)XP,YP,ZP,VPX,VPY,VPZ
 1002 FORMAT('XP=',F12.6,2X,'YP=',F12.6,2X,'ZP=',F12.6,2X,'VPX=',
     >       F12.6,2X,'VPY=',F12.6,2X,'VPZ=',F12.6,2X)
      WRITE(16,*)'SDT=',SDT
      STOP
C
c      RETURN
      END
C
c
c
      SUBROUTINE IBWCK6
     I      (IBWJ1,IBWJ2,XP,YP,ZP,VPX,VPY,VPZ,NXW,NZW,M,IE,SDT,MAXEL,
     O       N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
C       10/10/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,9)
C
      IF(IBWJ1*IBWJ2.NE.0)THEN
        IBPROD=IBWJ1*IBWJ2
        IF(NXW.EQ.1)THEN
          IF(IBPROD.EQ.16616 .OR. IBPROD.EQ.16875 .OR. IBPROD.EQ.156)
     >       GOTO 160
          IF(IBPROD.EQ.29016 .OR. IBPROD.EQ.29375 .OR. IBPROD.EQ.276)
     >       GOTO 170
          IF(IBPROD.EQ.31356 .OR. IBPROD.EQ.31725 .OR. IBPROD.EQ.299)
     >       GOTO 180
        ENDIF
        IF(NZW.EQ.1)THEN
          IF(IBPROD.EQ.15500 .OR. IBPROD.EQ.210)GOTO 160
          IF(IBPROD.EQ.54990 .OR. IBPROD.EQ.600)GOTO 170
          IF(IBPROD.EQ.18090 .OR. IBPROD.EQ.1190)GOTO 180
        ENDIF
        IF(IBWJ1.EQ.1 .OR. IBWJ2.EQ.1)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.2 .OR. IBWJ2.EQ.2)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.3 .OR. IBWJ2.EQ.3)THEN
          GOTO 180
        ELSEIF(IBWJ1.EQ.4 .OR. IBWJ2.EQ.4)THEN
          GOTO 190
        ELSEIF(IBWJ1.EQ.5 .OR. IBWJ2.EQ.5)THEN
          GOTO 200
        ELSEIF(IBWJ1.EQ.12 .OR. IBWJ2.EQ.12)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.13 .OR. IBWJ2.EQ.13)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.14 .OR. IBWJ2.EQ.14)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.15 .OR. IBWJ2.EQ.15)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.23 .OR. IBWJ2.EQ.23)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.24 .OR. IBWJ2.EQ.24)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.25 .OR. IBWJ2.EQ.25)THEN
          GOTO 170
        ELSE
          GOTO 180
        ENDIF
C
 160    CONTINUE
        N1=IE(M,1)
        N2=IE(M,3)
        N3=IE(M,6)
        N4=IE(M,4)
        GOTO 220
C
 170    CONTINUE
        N1=IE(M,1)
        N2=IE(M,4)
        N3=IE(M,5)
        N4=IE(M,2)
        GOTO 220
C
 180    CONTINUE
        N1=IE(M,2)
        N2=IE(M,5)
        N3=IE(M,6)
        N4=IE(M,3)
        GOTO 220
C
 190    CONTINUE
        N1=IE(M,1)
        N2=IE(M,2)
        N3=IE(M,3)
        N4=0
        GOTO 220
C
 200    CONTINUE
        N1=IE(M,4)
        N2=IE(M,5)
        N3=IE(M,6)
        N4=0
C
 220    CONTINUE
        XQ=XP
        YQ=YP
        ZQ=ZP
        VQX=VPX
        VQY=VPY
        VQZ=VPZ
        SDTT1=SDT
        RETURN
      ENDIF
c
      WRITE(16,*)'ERROR MESSAGE AT IBWCK6: The Following Point Can ',
     > 'Not be Tracked - STOP'
      WRITE(16,1002)XP,YP,ZP,VPX,VPY,VPZ
 1002 FORMAT('XP=',F12.6,2X,'YP=',F12.6,2X,'ZP=',F12.6,2X,'VPX=',
     >       F12.6,2X,'VPY=',F12.6,2X,'VPZ=',F12.6,2X)
      WRITE(16,*)'SDT=',SDT
      STOP
C
c      RETURN
      END
C
c
c
      SUBROUTINE IBWCK8
     I          (IBWJ1,IBWJ2,XP,YP,ZP,VPX,VPY,VPZ,NXW,NYW,NZW,M,IE,SDT,
     O           MAXEL,  N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
C       10/10/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,9)
C
      IF(IBWJ1*IBWJ2.NE.0)THEN
        IBPROD=IBWJ1*IBWJ2
        IF(NXW.EQ.1)THEN
          IF(IBPROD.EQ.132 .OR. IBPROD.EQ.182)GOTO 200
          IF(IBPROD.EQ.240 .OR. IBPROD.EQ.306)GOTO 210
        ENDIF
        IF(NYW.EQ.1)THEN
          IF(IBPROD.EQ.154 .OR. IBPROD.EQ.270)GOTO 160
          IF(IBPROD.EQ.156 .OR. IBPROD.EQ.272)GOTO 180
        ENDIF
        IF(NZW.EQ.1)THEN
          IF(IBPROD.EQ.165 .OR. IBPROD.EQ.192)GOTO 170
          IF(IBPROD.EQ.221 .OR. IBPROD.EQ.252)GOTO 190
        ENDIF
        IF(IBWJ1.EQ.1 .OR. IBWJ2.EQ.1)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.2 .OR. IBWJ2.EQ.2)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.3 .OR. IBWJ2.EQ.3)THEN
          GOTO 180
        ELSEIF(IBWJ1.EQ.4 .OR. IBWJ2.EQ.4)THEN
          GOTO 190
        ELSEIF(IBWJ1.EQ.5 .OR. IBWJ2.EQ.5)THEN
          GOTO 200
        ELSEIF(IBWJ1.EQ.6 .OR. IBWJ2.EQ.6)THEN
          GOTO 210
        ELSEIF(IBWJ1.EQ.21 .OR. IBWJ2.EQ.21)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.22 .OR. IBWJ2.EQ.22)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.23 .OR. IBWJ2.EQ.23)THEN
          GOTO 180
        ELSEIF(IBWJ1.EQ.24 .OR. IBWJ2.EQ.24)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.25 .OR. IBWJ2.EQ.25)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.26 .OR. IBWJ2.EQ.26)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.27 .OR. IBWJ2.EQ.27)THEN
          GOTO 180
        ELSEIF(IBWJ1.EQ.28 .OR. IBWJ2.EQ.28)THEN
          GOTO 190
        ELSEIF(IBWJ1.EQ.29 .OR. IBWJ2.EQ.29)THEN
          GOTO 160
        ELSEIF(IBWJ1.EQ.30 .OR. IBWJ2.EQ.30)THEN
          GOTO 170
        ELSEIF(IBWJ1.EQ.31 .OR. IBWJ2.EQ.31)THEN
          GOTO 180
        ELSE
          GOTO 190
        ENDIF
C
 160    CONTINUE
        N1=IE(M,1)
        N2=IE(M,4)
        N3=IE(M,8)
        N4=IE(M,5)
        GOTO 220
C
 170    CONTINUE
        N1=IE(M,1)
        N2=IE(M,2)
        N3=IE(M,6)
        N4=IE(M,5)
        GOTO 220
C
 180    CONTINUE
        N1=IE(M,2)
        N2=IE(M,3)
        N3=IE(M,7)
        N4=IE(M,6)
        GOTO 220
C
 190    CONTINUE
        N1=IE(M,4)
        N2=IE(M,3)
        N3=IE(M,7)
        N4=IE(M,8)
        GOTO 220
C
 200    CONTINUE
        N1=IE(M,1)
        N2=IE(M,2)
        N3=IE(M,3)
        N4=IE(M,4)
        GOTO 220
C
 210    CONTINUE
        N1=IE(M,5)
        N2=IE(M,6)
        N3=IE(M,7)
        N4=IE(M,8)
C
 220    CONTINUE
        XQ=XP
        YQ=YP
        ZQ=ZP
        VQX=VPX
        VQY=VPY
        VQZ=VPZ
        SDTT1=SDT
        RETURN
      ENDIF
c
      WRITE(16,*)'ERROR MESSAGE AT IBWCK8: The Following Point Can ',
     > 'Not be Tracked - STOP'
      WRITE(16,1002)XP,YP,ZP,VPX,VPY,VPZ
 1002 FORMAT('XP=',F12.6,2X,'YP=',F12.6,2X,'ZP=',F12.6,2X,'VPX=',
     >       F12.6,2X,'VPY=',F12.6,2X,'VPZ=',F12.6,2X)
      WRITE(16,*)'SDT=',SDT
      STOP
C
c      RETURN
      END
C
C
C
      SUBROUTINE CHNGSN
     I   (MXNOD,NSTART,NEND,NJUMP, NDIM,
     M    VXP,VYP,VZP,VXQ,VYQ,VZQ,VXX,VYY,VZZ)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION VXX(MXNOD),VYY(MXNOD),VZZ(MXNOD)
C
      VXP=-VXP
      VYP=-VYP
      VXQ=-VXQ
      VYQ=-VYQ
      IF(NDIM.EQ.3)THEN
        VZP=-VZP
        VZQ=-VZQ
      ENDIF
      DO I=NSTART,NEND,NJUMP
        VXX(I)=-VXX(I)
        VYY(I)=-VYY(I)
        IF(NDIM.EQ.3)VZZ(I)=-VZZ(I)
      ENDDO
      RETURN
      END
