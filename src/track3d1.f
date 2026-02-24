C $$$$$$$$$$$ TRACK3D1.F     3/15/95
C     10/4/93    6:00 pm
C
      SUBROUTINE CKBDY
     I   (NPBB,ISB,ISV,ISC,IB,
     O    NBDYB,IBDY,
     M    IWRK)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /TCBC/ MXCNP,MXCES,MXCPR,MXCDP,
     .              NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,
     .              NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /TDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
C
      DIMENSION NBDYB(MAXNP),IBDY(MXTUBS),IWRK(MAXBES)
      DIMENSION NPBB(MAXBNP),IB(MAXNP)
      DIMENSION ISV(5,MXVES),ISC(5,MXCES),ISB(6,MAXBES)
C
C ----- INITIATION
C NOTE: NBDYB(NP) REPRESENTS THE ACCUMULATED NUMBER OF UNSPECIFIED
C       BOUNDARY SIDES THAT THE 1-ST THROUGH THE (NP-1)-TH GLOBAL NODE
C       CONNECT.
C       IBDY(NBS) REPRESENTS THE BOUNDARY SIDE WHICH RELATES TO THE NBS-
C       TH ACCUMULATED UNSPECIFIED BOUNDARY SIDE.
C
      DO NP=1,MAXNP
        NBDYB(NP)=0
      ENDDO
      DO NBS=1,MXTUBS
        IBDY(NBS)=0
      ENDDO
C
C ----- CREATE WORKING ARRAY
C
      DO IBS=1,NBES
        IWRK(IBS)=0
      ENDDO
      DO IV=1,NVES
        NV=ISV(5,IV)
        IWRK(NV)=1
      ENDDO
      DO IC=1,NCES
        NC=ISC(5,IC)
        IWRK(NC)=1
      ENDDO
C
C ----- CHECK ALL THE UNSPECIFIED BOUNDARY SIDES
C
      DO 100 IBS=1,NBES
        IF(IWRK(IBS).EQ.1)GOTO 100
C
C ----- THIS IS AN UNSPECIFIED BOUNDADRY SIDE
C
        N1=ISB(1,IBS)
        N2=ISB(2,IBS)
        N3=ISB(3,IBS)
        N4=ISB(4,IBS)
        NN1=NPBB(N1)
        NN2=NPBB(N2)
        NN3=NPBB(N3)
        IF(N4.NE.0)THEN
          NN4=NPBB(N4)
          IF(IB(NN1).LE.3 .AND. IB(NN2).LE.3 .AND. IB(NN3).LE.3 .AND.
     >       IB(NN4).LE.3)THEN
            IWRK(IBS)=1
            GOTO 100
          ENDIF
        ELSE
          NN4=0
          IF(IB(NN1).LE.3 .AND. IB(NN2).LE.3 .AND. IB(NN3).LE.3)THEN
            IWRK(IBS)=1
            GOTO 100
          ENDIF
        ENDIF
C
        NBDYB(NN1)=NBDYB(NN1)+1
        NBDYB(NN2)=NBDYB(NN2)+1
        NBDYB(NN3)=NBDYB(NN3)+1
        IF(NN4.NE.0)THEN
          NBDYB(NN4)=NBDYB(NN4)+1
        ENDIF
  100 CONTINUE
C
C ----- REARRANGE NBDYB
C
      DO NP=2,NNP
        IF(NBDYB(NP).EQ.0)THEN
          NBDYB(NP)=NBDYB(NP-1)
        ELSE
          NBDYB(NP)=NBDYB(NP-1)+NBDYB(NP)
        ENDIF
      ENDDO
      NTUBS=NBDYB(NNP)
      CALL WARMSG(NTUBS,MXTUBS,' CKBDY','MXTUBS',1)
c
ccc    print *,'the total number of unspecified boundaries,NTUBS, is ',
ccc     >        NTUBS
c
      DO NP=NNP,2,-1
        NBDYB(NP)=NBDYB(NP-1)
      ENDDO
      NBDYB(1)=0
C
C ----- CREATE IBDY
C
      DO 200 IBS=1,NBES
        IF(IWRK(IBS).EQ.1)GOTO 200
C
C ----- THIS IS AN UNSPECIFIED BOUNDADRY SIDE
C
        N1=ISB(1,IBS)
        N2=ISB(2,IBS)
        N3=ISB(3,IBS)
        N4=ISB(4,IBS)
        NN1=NPBB(N1)
        NN2=NPBB(N2)
        NN3=NPBB(N3)
        IF(N4.NE.0)THEN
          NN4=NPBB(N4)
        ELSE
          NN4=0
        ENDIF
C
        NBDYB(NN1)=NBDYB(NN1)+1
        NB=NBDYB(NN1)
        IBDY(NB)=IBS
C
        NBDYB(NN2)=NBDYB(NN2)+1
        NB=NBDYB(NN2)
        IBDY(NB)=IBS
C
        NBDYB(NN3)=NBDYB(NN3)+1
        NB=NBDYB(NN3)
        IBDY(NB)=IBS
C
        IF(NN4.NE.0)THEN
          NBDYB(NN4)=NBDYB(NN4)+1
          NB=NBDYB(NN4)
          IBDY(NB)=IBS
        ENDIF
  200 CONTINUE
C
C ----- READJUST NBDYB
C
      DO NP=NNP,2,-1
        NBDYB(NP)=NBDYB(NP-1)
      ENDDO
      NBDYB(1)=0
C
      RETURN
      END
C
c
c
      SUBROUTINE IBE3D
     I     (MAXEL,MAXNP,MXKBD,IE,NLRL,LRL,IB,NEL,
     O      IBE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,11),NLRL(MAXNP),LRL(MXKBD,MAXNP)
      DIMENSION IB(MAXNP),IBE(MAXEL)
      DIMENSION MS(6),LS(3,6,3)
C
      DATA LS /4,1,8, 1,2,5, 2,3,6, 3,4,7, 1,2,3, 5,6,7,
     >         1,2,4, 2,3,5, 3,1,6, 1,2,3, 4,5,6, 0,0,0,
     >         2,3,4, 1,4,3, 1,2,4, 1,3,2, 0,0,0, 0,0,0/
C
      DO 400 M=1,NEL
        NBS=0
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,II,IK)
C
C ----- CHECK BOUNDARY SIDES FOR THE ELEMENT
C
        DO 200 I=1,II
          N1=LS(1,I,IK)
          N2=LS(2,I,IK)
          N3=LS(3,I,IK)
C
          NN1=IE(M,N1)
          NN2=IE(M,N2)
          NN3=IE(M,N3)
          IF(IB(NN1).EQ.0 .OR. IB(NN2).EQ.0 .OR. IB(NN3).EQ.0)GOTO 200
C
C ----- CHECK IF THIS SIDE IS A BOUNDARY SIDE
C
          DO 100 J1=1,NLRL(NN1)
            M1=LRL(J1,NN1)
            IF(M1.EQ.M)GOTO 100
            DO 90 J2=1,NLRL(NN2)
              M2=LRL(J2,NN2)
              IF(M2.NE.M1)GOTO 90
              DO 80 J3=1,NLRL(NN3)
                M3=LRL(J3,NN3)
                IF(M3.EQ.M2)GOTO 200
   80         CONTINUE
   90       CONTINUE
  100     CONTINUE
C
C ----- HERE WE GO, THIS SIDE IS A BOUNDARY SIDE
C
          NBS=NBS+1
          MS(NBS)=I
  200   CONTINUE
C
        IBE(M)=0
        DO I=1,NBS
          IBE(M)=IBE(M)+MS(I)*10**(NBS-I)
        ENDDO
  400 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE ADVW3D
     I                 (MXNPW,MXELW,NXW,NYW,NZW,ISHAPE,EPSX,
     O                  IBW,IEW,NLRLW,LRLW,DL468)
C
C $$$$$ TO GENERATE FINE GRIDS FOR TRACKING.
C  NOTE: IF ISHAPE=4, TETRAHEDRAL ELEMENTS
C        IF ISHAPE=6, PENTAHEDRAL ELEMENTS
C        IF ISHAPE=8, HEXAHEDRAL ELEMENTS
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IBW(MXNPW),IEW(MXELW,8),NLRLW(MXNPW),LRLW(24,MXNPW)
      DIMENSION DL468(8,MXNPW)
C
      IF(ISHAPE.EQ.4 .OR. ISHAPE.EQ.0)THEN
C ----- FOR THE CASES OF TRETRAHEFRAL ELEMENTS
        N=0
        DO 120 K=NXW+1,1,-1
          DL1=DBLE(K-1)/DBLE(NXW)
          DO 115 I=NXW+1,1,-1
            DL2=DBLE(I-1)/DBLE(NXW)
            IF(DL1+DL2.GT.1.0D0 .AND. DABS(DL1+DL2-1.0D0).GT.EPSX)
     >         GOTO 115
            DO 110 J=NXW+1,1,-1
              DL3=DBLE(J-1)/DBLE(NXW)
              IF(DL1+DL2+DL3.GT.1.0D0 .AND. DABS(DL1+DL2+DL3-1.0D0)
     >           .GT.EPSX)GOTO 110
              DO 105 L=NXW+1,1,-1
                DL4=DBLE(L-1)/DBLE(NXW)
                IF(DABS(DL1+DL2+DL3+DL4-1.0D0).GT.EPSX)GOTO 105
                N=N+1
C
                DL468(1,N)=DL1
                DL468(2,N)=DL2
                DL468(3,N)=DL3
                DL468(4,N)=DL4
                DL468(5,N)=0.0D0
                DL468(6,N)=0.0D0
                DL468(7,N)=0.0D0
                DL468(8,N)=0.0D0
                DO II=1,8
                  DL468(II,N)=DINT(DL468(II,N)*1.0D6)
                ENDDO
                SUM=0.0D0
                DO II=1,8
                  SUM=SUM+DL468(II,N)
                  IF(SUM.GT.1.0D6)DL468(II,N)=DL468(II,N)-SUM+1.0D6
                ENDDO
                IF(SUM.LT.1.0D6)THEN
                  DO II=1,8
                    IF(DL468(II,N).NE.0.0D0)THEN
                      DL468(II,N)=DL468(II,N)+1.0D6-SUM
                      GOTO 100
                    ENDIF
                  ENDDO
                ENDIF
  100           CONTINUE
                DO II=1,8
                  DL468(II,N)=DINT(DL468(II,N))
                ENDDO
C
                IF(K.EQ.1)THEN
                  IF(I.EQ.1)THEN
                    IF(J.EQ.1)THEN
                      IBW(N)=132
                    ELSEIF(L.EQ.1)THEN
                      IBW(N)=124
                    ELSE
                      IBW(N)=34
                    ENDIF
                  ELSEIF(J.EQ.1)THEN
                    IF(L.EQ.1)THEN
                      IBW(N)=143
                    ELSE
                      IBW(N)=42
                    ENDIF
                  ELSEIF(L.EQ.1)THEN
                    IBW(N)=23
                  ELSE
                    IBW(N)=1
                  ENDIF
                ELSEIF(I.EQ.1)THEN
                  IF(J.EQ.1)THEN
                    IF(L.EQ.1)THEN
                      IBW(N)=234
                    ELSE
                      IBW(N)=14
                    ENDIF
                  ELSEIF(L.EQ.1)THEN
                    IBW(N)=13
                  ELSE
                    IBW(N)=2
                  ENDIF
                ELSEIF(J.EQ.1)THEN
                  IF(L.EQ.1)THEN
                    IBW(N)=12
                  ELSE
                    IBW(N)=3
                  ENDIF
                ELSEIF(L.EQ.1)THEN
                  IBW(N)=4
                ELSE
                  IBW(N)=0
                ENDIF
  105         CONTINUE
  110       CONTINUE
  115     CONTINUE
  120   CONTINUE
        J2=1
        J1=0
        DO 130 L=1,NXW
          JUMP1=L
          JUMP2=(L-1)*L/2
          J1=J1+JUMP1
          J2=J2+JUMP2
c         NCOUNT=L**3-(L-1)**3
          M1=(L-1)**3+1
          IEW(M1,1)=J2
          IEW(M1,2)=J2+J1
          IEW(M1,3)=J2+J1+1
          IEW(M1,4)=J2+J1+2
C
          IF(L.EQ.1)GOTO 130
          J3=J2
          DO 125 L1=2,L
            M2=(L1-2)*(L1-1)*3+2+(L-1)**3
            J3=J3+L1-1
            IEW(M2,1)=J3
            IEW(M2,2)=J3+J1
            IEW(M2,3)=J3+J1+L1
            IEW(M2,4)=J3+J1+L1+1
            IEW(M2+1,1)=J3
            IEW(M2+1,2)=J3+J1+1
            IEW(M2+1,3)=J3-L1+1
            IEW(M2+1,4)=J3+J1
            IEW(M2+2,1)=J3
            IEW(M2+2,2)=J3+J1+1
            IEW(M2+2,3)=J3+J1
            IEW(M2+2,4)=J3+J1+L1+1
            IEW(M2+3,1)=J3
            IEW(M2+3,2)=J3+J1+1
            IEW(M2+3,3)=J3+J1+L1+1
            IEW(M2+3,4)=J3+1
            IEW(M2+4,1)=J3
            IEW(M2+4,2)=J3+J1+1
            IEW(M2+4,3)=J3+1
            IEW(M2+4,4)=J3-L1+1
            IEW(M2+5,1)=J3+1
            IEW(M2+5,2)=J3+J1+1
            IEW(M2+5,3)=J3+J1+L1+1
            IEW(M2+5,4)=J3+J1+L1+2
            IF(L.GE.3)THEN
              M5=M2+5
              J4=J3
              DO L2=1,L-2
                J4=J4+1
                M5=M5+1
                IEW(M5,1)=J4
                IEW(M5,2)=J4-L1+1
                IEW(M5,3)=J4-L1
                IEW(M5,4)=J4+J1
                M5=M5+1
                IEW(M5,1)=J4
                IEW(M5,2)=J4+J1+1
                IEW(M5,3)=J4-L1+1
                IEW(M5,4)=J4+J1
                M5=M5+1
                IEW(M5,1)=J4
                IEW(M5,2)=J4+J1+1
                IEW(M5,3)=J4+J1
                IEW(M5,4)=J4+J1+L1+1
                M5=M5+1
                IEW(M5,1)=J4
                IEW(M5,2)=J4+J1+1
                IEW(M5,3)=J4+J1+L1+1
                IEW(M5,4)=J4+1
                M5=M5+1
                IEW(M5,1)=J4
                IEW(M5,2)=J4+J1+1
                IEW(M5,3)=J4+1
                IEW(M5,4)=J4-L1+1
                M5=M5+1
                IEW(M5,1)=J4+1
                IEW(M5,2)=J4+J1+1
                IEW(M5,3)=J4+J1+L1+1
                IEW(M5,4)=J4+J1+L1+2
              ENDDO
            ENDIF
  125     CONTINUE
  130   CONTINUE
        DO 135 I1=1,NXW**3
          DO I2=5,8
            IEW(I1,I2)=0
          ENDDO
  135   CONTINUE
C
        NELW=NXW**3
        NP1=NXW+1
        NNPW=NP1*(NP1+1)*(NP1+2)/6
        CALL LRL3D
     I            (IEW,NELW,NNPW,MXELW,MXNPW,24,8,
     O             NLRLW,LRLW)
      ENDIF
C
      IF(ISHAPE.EQ.6)THEN
C ----- FOR THE CASES OF PENTAHEDRAL ELEMENTS
        N=0
        DO 220 L=1,NZW+1
          XI1=DBLE(NZW+1-L)/DBLE(NZW)
          XI2=DBLE(L-1)/DBLE(NZW)
          DO 215 I=NXW+1,1,-1
            DL1=DBLE(I-1)/DBLE(NXW)
            DO 210 J=NXW+1,1,-1
              DL2=DBLE(J-1)/DBLE(NXW)
              IF(DL1+DL2.GT.1.0D0 .AND. DABS(DL1+DL2-1.0D0).GT.EPSX)
     >          GOTO 210
              DO 205 K=NXW+1,1,-1
                DL3=DBLE(K-1)/DBLE(NXW)
                IF(DABS(DL1+DL2+DL3-1.0D0).LE.1.0D-6)THEN
                  N=N+1
C
                  DL468(1,N)=DL1*XI1
                  DL468(2,N)=DL2*XI1
                  DL468(3,N)=DL3*XI1
                  DL468(4,N)=DL1*XI2
                  DL468(5,N)=DL2*XI2
                  DL468(6,N)=DL3*XI2
                  DL468(7,N)=0.0D0
                  DL468(8,N)=0.0D0
                  DO II=1,8
                    DL468(II,N)=DINT(DL468(II,N)*1.0D6)
                  ENDDO
                  SUM=0.0D0
                  DO II=1,8
                    SUM=SUM+DL468(II,N)
                    IF(SUM.GT.1.0D6)DL468(II,N)=DL468(II,N)-SUM+1.0D6
                  ENDDO
                  IF(SUM.LT.1.0D6)THEN
                    DO II=1,8
                      IF(DL468(II,N).NE.0.0D0)THEN
                        DL468(II,N)=DL468(II,N)+1.0D6-SUM
                        GOTO 200
                      ENDIF
                    ENDDO
                  ENDIF
  200             CONTINUE
                  DO II=1,8
                    DL468(II,N)=DINT(DL468(II,N))
                  ENDDO
C
                  IF(L.EQ.1)THEN
                    IF(I.EQ.1)THEN
                      IF(J.EQ.1)THEN
                        IBW(N)=134
                      ELSEIF(K.EQ.1)THEN
                        IBW(N)=234
                      ELSE
                        IBW(N)=34
                      ENDIF
                    ELSEIF(J.EQ.1)THEN
                      IF(K.EQ.1)THEN
                        IBW(N)=124
                      ELSE
                        IBW(N)=14
                      ENDIF
                    ELSEIF(K.EQ.1)THEN
                      IBW(N)=24
                    ELSE
                      IBW(N)=4
                    ENDIF
                  ELSEIF(L.EQ.NZW+1)THEN
                    IF(I.EQ.1)THEN
                      IF(J.EQ.1)THEN
                        IBW(N)=135
                      ELSEIF(K.EQ.1)THEN
                        IBW(N)=235
                      ELSE
                        IBW(N)=35
                      ENDIF
                    ELSEIF(J.EQ.1)THEN
                      IF(K.EQ.1)THEN
                        IBW(N)=125
                      ELSE
                        IBW(N)=15
                      ENDIF
                    ELSEIF(K.EQ.1)THEN
                      IBW(N)=25
                    ELSE
                      IBW(N)=5
                    ENDIF
                  ELSE
                    IF(I.EQ.1)THEN
                      IF(J.EQ.1)THEN
                        IBW(N)=13
                      ELSEIF(K.EQ.1)THEN
                        IBW(N)=23
                      ELSE
                        IBW(N)=3
                      ENDIF
                    ELSEIF(J.EQ.1)THEN
                      IF(K.EQ.1)THEN
                        IBW(N)=12
                      ELSE
                        IBW(N)=1
                      ENDIF
                    ELSEIF(K.EQ.1)THEN
                      IBW(N)=2
                    ELSE
                      IBW(N)=0
                    ENDIF
                  ENDIF
                ENDIF
  205         CONTINUE
  210       CONTINUE
  215     CONTINUE
  220   CONTINUE
C
        NW=(NXW**2+3*NXW+2)/2
        DO 250 L=1,NZW
          LL=(L-1)*(NXW**2)
          NN=(L-1)*(NXW**2+3*NXW+2)/2
          NC=0
          NI=0
          DO 240 I=1,NXW
            NT=2*I-1
            NII=NI
            DO J=NC+1,NC+NT,2
              NI=NI+1
              IEW(J+LL,1)=NI+NN
              IEW(J+LL,2)=NI+I+NN
              IEW(J+LL,3)=NI+I+1+NN
              IEW(J+LL,4)=NI+NN+NW
              IEW(J+LL,5)=NI+I+NN+NW
              IEW(J+LL,6)=NI+I+1+NN+NW
              IEW(J+LL,7)=0
              IEW(J+LL,8)=0
            ENDDO
            IF(I.GT.1)THEN
              DO J=NC+2,NC+NT-1,2
                NII=NII+1
                IEW(J+LL,1)=NII+NN
                IEW(J+LL,2)=NII+I+1+NN
                IEW(J+LL,3)=NII+1+NN
                IEW(J+LL,4)=NII+NN+NW
                IEW(J+LL,5)=NII+I+1+NN+NW
                IEW(J+LL,6)=NII+1+NN+NW
                IEW(J+LL,7)=0
                IEW(J+LL,8)=0
              ENDDO
            ENDIF
            NC=I**2
  240     CONTINUE
  250   CONTINUE
C
        NELW=NZW*NXW**2
        NNPW=(NZW+1)*(NXW**2+3*NXW+2)/2
        CALL LRL3D
     I            (IEW,NELW,NNPW,MXELW,MXNPW,24,8,
     O             NLRLW,LRLW)
      ENDIF
C
      IF(ISHAPE.EQ.8)THEN
C ----- FOR THE CASES OF HEXAHEDRAL ELEMENTS
        DO 460 K=1,NZW+1
          ZI1=DBLE(NZW+1-K)/DBLE(NZW)
          ZI2=DBLE(K-1)/DBLE(NZW)
          DO 455 I=1,NYW+1
            YI1=DBLE(NYW+1-I)/DBLE(NYW)
            YI2=DBLE(I-1)/DBLE(NYW)
            DO 450 J=1,NXW+1
              XI1=DBLE(NXW+1-J)/DBLE(NXW)
              XI2=DBLE(J-1)/DBLE(NXW)
              N=J+(I-1)*(NXW+1)+(K-1)*(NXW+1)*(NYW+1)
C
              DL468(1,N)=XI1*YI1*ZI1
              DL468(2,N)=XI2*YI1*ZI1
              DL468(3,N)=XI2*YI2*ZI1
              DL468(4,N)=XI1*YI2*ZI1
              DL468(5,N)=XI1*YI1*ZI2
              DL468(6,N)=XI2*YI1*ZI2
              DL468(7,N)=XI2*YI2*ZI2
              DL468(8,N)=XI1*YI2*ZI2
              DO II=1,8
                DL468(II,N)=DINT(DL468(II,N)*1.0D6)
              ENDDO
              SUM=0.0D0
              DO II=1,8
                SUM=SUM+DL468(II,N)
                IF(SUM.GT.1.0D6)DL468(II,N)=DL468(II,N)-SUM+1.0D6
              ENDDO
              IF(SUM.LT.1.0D6)THEN
                DO II=1,8
                  IF(DL468(II,N).NE.0.0D0)THEN
                    DL468(II,N)=DL468(II,N)+1.0D6-SUM
                    GOTO 445
                  ENDIF
                ENDDO
              ENDIF
  445         CONTINUE
              DO II=1,8
                DL468(II,N)=DINT(DL468(II,N))
              ENDDO
C
              IF(J.EQ.1)THEN
                IF(I.EQ.1)THEN
                  IBW(N)=21
                  IF(K.EQ.1)IBW(N)=11
                  IF(K.EQ.NZW+1)IBW(N)=15
                ELSEIF(I.EQ.NYW+1)THEN
                  IBW(N)=24
                  IF(K.EQ.1)IBW(N)=14
                  IF(K.EQ.NZW+1)IBW(N)=18
                ELSEIF(K.EQ.1)THEN
                  IBW(N)=25
                ELSEIF(K.EQ.NZW+1)THEN
                  IBW(N)=29
                ELSE
                  IBW(N)=1
                ENDIF
              ELSEIF(I.EQ.1)THEN
                IF(J.EQ.NXW+1)THEN
                  IBW(N)=22
                  IF(K.EQ.1)IBW(N)=12
                  IF(K.EQ.NZW+1)IBW(N)=16
                ELSEIF(K.EQ.1)THEN
                  IBW(N)=26
                ELSEIF(K.EQ.NZW+1)THEN
                  IBW(N)=30
                ELSE
                  IBW(N)=2
                ENDIF
              ELSEIF(J.EQ.NXW+1)THEN
                IF(I.EQ.NYW+1)THEN
                  IBW(N)=23
                  IF(K.EQ.1)IBW(N)=13
                  IF(K.EQ.NZW+1)IBW(N)=17
                ELSEIF(K.EQ.1)THEN
                  IBW(N)=27
                ELSEIF(K.EQ.NZW+1)THEN
                  IBW(N)=31
                ELSE
                  IBW(N)=3
                ENDIF
              ELSEIF(I.EQ.NYW+1)THEN
                IF(K.EQ.1)THEN
                  IBW(N)=28
                ELSEIF(K.EQ.NZW+1)THEN
                  IBW(N)=32
                ELSE
                  IBW(N)=4
                ENDIF
              ELSEIF(K.EQ.1)THEN
                IBW(N)=5
              ELSEIF(K.EQ.NZW+1)THEN
                IBW(N)=6
              ELSE
                IBW(N)=0
              ENDIF
  450       CONTINUE
  455     CONTINUE
  460   CONTINUE
C ----- DETERMINE WORKING ARRAY IEW
        DO 475 K=1,NZW
          DO 470 I=1,NYW
            DO 465 J=1,NXW
              MW=(I-1)*NXW+J+(K-1)*NXW*NYW
              IEW(MW,1)=(I-1)*(NXW+1)+J+(K-1)*(NXW+1)*(NYW+1)
              IEW(MW,2)=(I-1)*(NXW+1)+J+1+(K-1)*(NXW+1)*(NYW+1)
              IEW(MW,3)=I*(NXW+1)+J+1+(K-1)*(NXW+1)*(NYW+1)
              IEW(MW,4)=I*(NXW+1)+J+(K-1)*(NXW+1)*(NYW+1)
              IEW(MW,5)=IEW(MW,1)+(NXW+1)*(NYW+1)
              IEW(MW,6)=IEW(MW,2)+(NXW+1)*(NYW+1)
              IEW(MW,7)=IEW(MW,3)+(NXW+1)*(NYW+1)
              IEW(MW,8)=IEW(MW,4)+(NXW+1)*(NYW+1)
  465       CONTINUE
  470     CONTINUE
  475   CONTINUE
C      DETERMINE WORKING ARRAYS NLRLW, AND LRLW
        NELW=NXW*NYW*NZW
        NNPW=(NXW+1)*(NYW+1)*(NZW+1)
        CALL LRL3D
     I            (IEW,NELW,NNPW,MXELW,MXNPW,24,8,
     O             NLRLW,LRLW)
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE GNTRAK
     I     (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,MXTUBS,
     I      MXNEP,NNP,NEL,NEFG,NTUBS,NXW,NYW,NZW,IDTI,IZOOM,EPSX,IBF,
     I      IDETQ, DELT,RAMADA,CP,X,VX,CPFG,XPFG,ISE,NFGMB,IE,IB,
     I      LRL,NLRL,IEW,IBW,LRLW,NLRLW,ISB,DCOSB,NBDYB,IBDY,DL468,
     i      ieww,ibww,lrlww,nlrlww,dl468w,mxnpww,mxelww,nxg,nyg,nzg,
     i      nwnp,npw,iwtyp,wss,mxwnp,mxwpr,
     O      CS,DTI,DTIFG,NPFGS,XSFG,CSFG,MPLOCS,IBCHK,XEFG,CEFG,MPLOCE,
     M      XW,VXW,xww,vxww)
C
C $$$$$ TO COMPUTE PARTICLE TRACKING STARTING FROM GLOBAL NODES
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CS(MAXNP),CP(MAXNP)
      DIMENSION DTI(MAXNP),X(MAXNP,3),VX(MAXNP,3)
      DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG),MPLOCS(MXNPFG)
      DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG),DTIFG(MXNPFG)
      DIMENSION ISE(MXKGL,8),NFGMB(MAXEL),DL468(8,MXNPW,3)
      DIMENSION XEFG(MXNEP,3),CEFG(MXNEP),MPLOCE(MXNPFG)
C
      DIMENSION IE(MAXEL,11),IB(MAXNP),LRL(MXKBD,MAXNP),NLRL(MAXNP)
      DIMENSION IBCHK(MAXEL),IEW(MXELW,8,3),IBW(MXNPW,3)
      DIMENSION LRLW(24,MXNPW,3),NLRLW(MXNPW,3),XW(MXNPW,3),VXW(MXNPW,3)
      DIMENSION ISB(6,MAXBES),DCOSB(3,MAXBES),NBDYB(MAXNP),IBDY(MXTUBS)
c
      DIMENSION IEWw(MXELWw,8,3),IBWw(MXNPWw,3),LRLWw(24,MXNPWw,3)
      DIMENSION NLRLWw(MXNPWw,3),XWw(MXNPWw,3),VXWw(MXNPWw,3)
      DIMENSION DL468w(8,MXNPWw,3)
c
      dimension npw(mxwnp),iwtyp(mxwnp),wss(mxwpr)
c
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),MK(8),CQ(7)
C
      DATA N1/1/, N2/2/, N3/3/, N4/4/
C
C ##### START TRACKINGS
C
      DO 910 N=1,NNP
C
        idi=0
        idps=0
C
        IF(IBF.EQ.1)DTI(N)=1.0D0/DELT
C
C ----- CHECK IF THIS POINT HAS ZERO VELOCITY OR IT IS NEGLIGIBLE
C
        CALL REPLAS(X(N,1),X(N,2),X(N,3),VX(N,1),VX(N,2),VX(N,3),
     >              XP,YP,ZP,VPX,VPY,VPZ)
        IF(VPX.EQ.0.0D0 .AND. VPY.EQ.0.0D0 .AND. VPZ.EQ.0.0D0)THEN
          IF(IBF.EQ.1)THEN
            CS(N)=CP(N)*DEXP(-RAMADA*DELT)
            XQ=X(N,1)
            YQ=X(N,2)
            ZQ=X(N,3)
          ENDIF
          GOTO 900
        ENDIF
C
C ----- PREPARE INFORMATION OF SOURCE POINT P
C
        NLRLN=NLRL(N)
        SDT=DELT
c
c       if(izoom.ne.0)then                                               3/15/95
          do 200 i=1,nwnp
            if(n.ne.npw(i))goto 200
            if(ibf.eq.1)then
              if(ib(n).eq.-1)then
                dti(n)=1.0d38
                iprof=iwtyp(i)
                cs(n)=wss(iprof)
                goto 900
              endif
              goto 210
            else
              if(ib(n).eq.-2)goto 900
              idps=1
              iprof=iwtyp(i)
              cps=wss(iprof)
              cpsp=cp(n)
              goto 210
            endif
  200     continue
c       endif                                                            3/15/95
C
  210   continue
c
        DO 250 I=1,NLRLN
          M=LRL(I,N)
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NODE,I1,I1)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          if(idps.eq.1)then
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,
     I          M,IBF,IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,IBWw(1,1),
     I          NLRLWw(1,1),LRLWw(1,1,1),IEWw(1,1,1),DL468w(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,
     I          M,IBF,IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,NZg,
     I          IBWw(1,2),NLRLWw(1,2),LRLWw(1,1,2),IEWw(1,1,2),
     i          DL468w(1,1,2),idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,NYg,NZg,
     I          IBWw(1,3),NLRLWw(1,3),LRLWw(1,1,3),IEWw(1,1,3),
     i          DL468w(1,1,3),idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
          else
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,
     I          M,IBF,IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,
     I          M,IBF,IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
          endif
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
C
  250   CONTINUE
C
C ***** THIS POINT WON'T MOVE INTO THE REGION OF INTEREST
C
        IFIX=1
        CALL FIXCHK
     I     (IFIX,IDTI,N,N,M,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
        IF(ICODE.EQ.0)THEN
          GO TO 900
        ELSE
          if(idi.eq.0)then
            idi=1
            goto 210
          endif
          WRITE(16,*)'ERROR IN FIXCHK, ICODE=',ICODE
          STOP
        ENDIF
C
C ##### START THE SUBSEQUENT TRACKINGS
C
  260   CONTINUE
c
        idi=0
c
        IF(IBF.EQ.2 .AND. IB(N).NE.0 .AND. IB(N).LE.3)IBCHK(M)=M
        SDT=SDT1
        CALL REPLAS(XQ,YQ,ZQ,VQX,VQY,VQZ,XP,YP,ZP,VPX,VPY,VPZ)
c       NN1=N1
c       NN2=N2
c       NN3=N3
c       NN4=N4
C
C ----- CHECK IF POINT P COINSIDES WITH ANY GLOBAL NODES
C
        CALL CKCOIN
     I      (X(1,1),X(1,2),X(1,3),XP,YP,ZP,NLRL,NODE,N1,N2,N3,N4,MAXNP,
     >       EPSX, NN,NLRLN)
        IF(NN.EQ.0)GOTO 400
C
C ----- FOR THE CASES THAT POINT P IS A GLOBAL NODE
C
        CALL REPLAS(X(NN,1),X(NN,2),X(NN,3),VX(NN,1),VX(NN,2),VX(NN,3),
     >              XP,YP,ZP,VPX,VPY,VPZ)
c
        if(izoom.ne.0)then
          do 300 i=1,nwnp
            if(nn.ne.npw(i))goto 300
            if(ibf.eq.1)then
              dtreal=delt-sdt
              xsi=(delt-dtreal)/delt
              iprof=iwtyp(i)
              cc=cp(nn)*xsi+wss(iprof)*(1.0d0-xsi)
              dti(n)=1.0d0/dtreal
              cs(n)=cc*dexp(-ramada*dtreal)
            endif
            goto 900
  300     continue
        endif
c
  310   continue
c
        DO 350 I=1,NLRLN
          M=LRL(I,NN)
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NODE,I1,I1)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          if(idps.eq.1)then
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,IBWw(1,1),
     I          NLRLWw(1,1),LRLWw(1,1,1),IEWw(1,1,1),DL468w(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,NZg,IBWw(1,2),
     I          NLRLWw(1,2),LRLWw(1,1,2),IEWw(1,1,2),DL468w(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,NYg,NZg,
     I          IBWw(1,3),NLRLWw(1,3),LRLWw(1,1,3),IEWw(1,1,3),
     i          DL468w(1,1,3), idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
          else
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
          endif
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
C
  350   CONTINUE
C
C ----- CHECK IF THIS POINT IS A BOUNDARY NODE
C
        IFIX=2
        CALL FIXCHK
     I    (IFIX,IDTI,N,NN,M,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I     EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I     X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I     CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O     NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
        IF(ICODE.EQ.0)THEN
          GO TO 900
        ELSE
          if(idi.eq.0)then
            idi=1
            goto 310
          endif
          WRITE(16,*)'ERROR IN FIXCHK, ICODE=',ICODE
          STOP
        ENDIF
C
C ----- FOR THE CASES THAT POINT P IS NOT A GLOBAL NODE
C
  400   CONTINUE
C
C ----- CHECK IF P IS ON ANY SIDES OF TETRAGON N1,N2,N3,N4
C
        CALL CKSIDE
     I      (X(1,1),X(1,2),X(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MAXNP,EPSX,
     O       KON,J1,J2)
        IF(KON.EQ.1)GOTO 491
C
C ----- NO, IT IS NOT ON ANY SIDE
C
        CALL ONPLAN
     I    (X(1,1),X(1,2),X(1,3),N1,N2,N3,MAXNP,
     M     XP,YP,ZP)
C ----- CHECK IF N1,N2,N3,N4 ARE ON A BOUNDARY PLANE
C
        CALL CKCNEL
     I              (LRL,NLRL,M,N1,N2,N3,MAXNP,MXKBD,
     O               ME1,KOUNT)
        IF(KOUNT.EQ.1)GOTO 490
C
C ----- N1,N2,N3,N4 ARE NOT ON A BOUNDARY PLANE
C
C ----- ELEMENT M1 IS WHAT WE NEED, THEN CHECK WHICH PLANE POINT P IS
C       ONTO
C
        ME=M
        M=ME1
C
  431   CONTINUE
C
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,I1,I1)
        CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >              XX,YY,ZZ,VXX,VYY,VZZ)
C
  490   CONTINUE
C
        if(idps.eq.1)then
        IF(NODE.EQ.4)THEN
          CALL ELTRK4
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,IBWw(1,1),
     I        NLRLWw(1,1),LRLWw(1,1,1),IEWw(1,1,1),DL468w(1,1,1),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XWw,VXWw,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSEIF(NODE.EQ.6)THEN
          CALL ELTRK6
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,NZg,IBWw(1,2),
     I        NLRLWw(1,2),LRLWw(1,1,2),IEWw(1,1,2),DL468w(1,1,2),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XWw,VXWw,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSE
          CALL ELTRK8
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELWw,MXNPWw,NXg,NYg,NZg,IBWw(1,3),
     I        NLRLWw(1,3),LRLWw(1,1,3),IEWw(1,1,3),DL468w(1,1,3),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4, XWw,VXWw,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ENDIF
        else
        IF(NODE.EQ.4)THEN
          CALL ELTRK4
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,
     I        IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSEIF(NODE.EQ.6)THEN
          CALL ELTRK6
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NZW,
     I        IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSE
          CALL ELTRK8
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NYW,NZW,
     I        IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4, XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ENDIF
        endif
C
        IF(SDT1.EQ.0.)GOTO 500
        IF(SDT1.NE.SDT)GOTO 260
C
C ----- CHECK IF P IS ON THE BOUNDARY PLANE N1,N2,N3,N4
C
        IF(KOUNT.EQ.1)THEN
C
C ----- FOR THE CASES THAT THIS BOUNDARY SIDE IS COMPOSED BY 4 POINTS
C       OR BY 3 POINTS
C
          IFIX=3
          CALL FIXCHK
     I     (IFIX,IDTI,N,N,M,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
          IF(ICODE.EQ.0)GOTO 900
        ENDIF
C
        IF(M.EQ.ME1)THEN
          M=ME
          GOTO 431
        ENDIF
c
        if(idi.eq.0)then
          idi=1
          m=me1
          goto 431
        endif
C
C ----- ERROR OCCURRED
C
        WRITE(16,*)'ERROR 3 AT GNTRAC, NODE',N,' CAN NOT BE',
     >            ' TRACKED'
        WRITE(16,1006)X(N,1),X(N,2),X(N,3),VX(N,1),VX(N,2),VX(N,3)
        WRITE(16,1007)XP,YP,ZP,VPX,VPY,VPZ
        WRITE(16,*)'SDT=',SDT
 1006   FORMAT('X(N)=',F12.6,2X,'Y(N)=',F12.6,2X,'Z(N)=',F12.6,2X,
     >         'VX(N)=',F12.6,2X,'VY(N)=',F12.6,2X,'VZ(N)=',F12.6,2X)
 1007   FORMAT('XP=',F12.6,2X,'YP=',F12.6,2X,'ZP=',F12.6,2X,'VPX=',
     >         F12.6,2X,'VPY=',F12.6,2X,'VPZ=',F12.6,2X)
        STOP
C
C *** YES, IT IS ON J1,J2 SIDE
C
  491   CONTINUE
C
        CALL ONLINE
     I      (X(1,1),X(1,2),X(1,3),J1,J2,MAXNP,
     M       XP,YP,ZP)
        NLRL1=NLRL(J1)
        NLRL2=NLRL(J2)
        KOUNT=0
        DO 493 I1=1,NLRL1
          M1=LRL(I1,J1)
          DO 492 I2=1,NLRL2
            M2=LRL(I2,J2)
            IF(M1.EQ.M2)THEN
              KOUNT=KOUNT+1
              MK(KOUNT)=M1
            ENDIF
  492     CONTINUE
  493   CONTINUE
C
  494   continue
c
        DO 499 I=1,KOUNT
          M=MK(I)
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NODE,I1,I1)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,M,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          if(idps.eq.1)then
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELWw,MXNPWw,NXg,IBWw(1,1),
     I          NLRLWw(1,1),LRLWw(1,1,1),IEWw(1,1,1),DL468w(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELWw,MXNPWw,NXg,NZg,IBWw(1,2),
     I          NLRLWw(1,2),LRLWw(1,1,2),IEWw(1,1,2),DL468w(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELWw,MXNPWw,NXg,NYg,NZg,
     I          IBWw(1,3),NLRLWw(1,3),LRLWw(1,1,3),IEWw(1,1,3),
     i          DL468w(1,1,3), idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XWw,VXWw,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
          else
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,M,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
          endif
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
C
  499   CONTINUE
C
        IFIX=4
        CALL FIXCHK
     I    (IFIX,IDTI,N,N,M,J1,J2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
        IF(ICODE.EQ.0)THEN
          GO TO 900
        ELSE
          if(idi.eq.0)then
            idi=1
            goto 494
          endif
          WRITE(16,*)'ERROR IN FIXCHK, ICODE=',ICODE
          STOP
        ENDIF
C
C ***** RECORDING INFORMATION
C
  500   CONTINUE
C
C ----- FOR THE CASES OF BACKWARD TRACKING
C
        IF(IBF.EQ.1)THEN
          CALL INTERP
     I               (MAXNP,MAXEL,MXNPFG,MXKGL,8,1,NODE,M,XQ,YQ,ZQ,X,CP,
     I                IE,XPFG,CPFG,NEL,NEFG,IZOOM,10,ISE,NFGMB,1,1,EPSX,
     O                CQ)
          CS(N)=CQ(1)*DEXP(-RAMADA*DELT)
          GOTO 900
        ENDIF
C
C ----- FOR FORWARD TRACKING NODES
C
        IF(IB(N).NE.0 .AND. IB(N).LE.3)IBCHK(M)=M                        1/06/95
        IE(M,11)=M
        NPFGS=NPFGS+1
        CALL WARMSG(NPFGS,MXNPFG,'GNTRAC','MXNPFG',1)
        XSFG(NPFGS,1)=XQ
        XSFG(NPFGS,2)=YQ
        XSFG(NPFGS,3)=ZQ
        CSFG(NPFGS)=CP(N)*DEXP(-RAMADA*DELT)
        MPLOCS(NPFGS)=M
C
  900   CONTINUE
  910 CONTINUE
C
  990 FORMAT(3F15.4)
  992 FORMAT(I5,3X,2F15.4)
C
      RETURN
      END
C
C
C
      SUBROUTINE HPTRAK
     I     (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,MXTUBS,
     I      MXNEP,MXNODE,MXNODS,NNP,NEL,NEFG,NTUBS,NXW,NYW,NZW,
     I      NPFG,NINIT,IBF,IDETQ, DELT,TMAX,RAMADA, IDTI,IZOOM,EPSX,
     I      CP,X,VX,CPFG,XPFG,MPLOC,ISE,NFGMB,IE,IB,LRL,NLRL,IEW,IBW,
     I      LRLW,NLRLW,ISB,DCOSB,NBDYB,IBDY,DL468,
     i      nwnp,npw,iwtyp,wss,mxwnp,mxwpr,
     O      CS,DTI,DTIFG,NPFGS,XSFG,CSFG,MPLOCS,NPFGEP,XEFG,CEFG,MPLOCE,
     M      XW,VXW)
C
C $$$$$ TO COMPUTE PARTICLE TRACKING STARTING FROM HIDDEN POINTS
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION CS(MAXNP),CP(MAXNP)
      DIMENSION DTI(MAXNP),DTIFG(MXNPFG),X(MAXNP,3),VX(MAXNP,3)
      DIMENSION XSFG(MXNODS,3),CSFG(MXNODS),MPLOCS(MXNODS)
      DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG),MPLOC(MXNPFG)
      DIMENSION ISE(MXKGL,8),NFGMB(MAXEL)
      DIMENSION XEFG(MXNODE,3),CEFG(MXNODE),MPLOCE(MXNPFG)
C
      DIMENSION IE(MAXEL,11),IB(MAXNP),LRL(MXKBD,MAXNP),NLRL(MAXNP)
      DIMENSION IEW(MXELW,8,3),IBW(MXNPW,3)
      DIMENSION LRLW(24,MXNPW,3),NLRLW(MXNPW,3)
      DIMENSION XW(MXNPW,3),VXW(MXNPW,3)
      DIMENSION ISB(6,MAXBES),DCOSB(3,MAXBES)
      DIMENSION NBDYB(MAXNP),IBDY(MXTUBS),DL468(8,MXNPW,3)
      dimension npw(mxwnp),iwtyp(mxwnp),wss(mxwpr)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),CQ(7),IEMM(4)
      DIMENSION MK(8),DL(8),KIT(4,6,3),DNX(8)
C
      DATA N1/1/, N2/2/, N3/3/, N4/4/
      DATA KIT/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
C ##### START TRACKINGS
C
      IF(IBF.EQ.1)THEN
        NWFG=NPFGS
      ELSE
        NWFG=NPFG
      ENDIF
C
      DO 900 N=NINIT,NWFG
c
        idi=0
        idps=0
c
        IF(IBF.EQ.1)THEN
          XP=XSFG(N,1)
          YP=XSFG(N,2)
          ZP=XSFG(N,3)
          M=MPLOCS(N)
        ELSE
          XP=XPFG(N,1)
          YP=XPFG(N,2)
          ZP=XPFG(N,3)
          M=MPLOC(N)
        ENDIF
        MM=M
        CALL ELENOD
     I      (IE(MM,5),IE(MM,7),
     O       NODE,II,ID)
C
C ----- CHECK IF THIS POINT COINSIDES WITH ANY GLOBAL NODES
C
        DO 120 I=1,NODE
          IEM=IE(MM,I)
          IF(DABS(XP-X(IEM,1)).LE.EPSX .AND. DABS(YP-X(IEM,2)).LE.EPSX
     >       .AND. DABS(ZP-X(IEM,3)).LE.EPSX)THEN
            IF(IBF.EQ.1)THEN
              CSFG(N)=CS(IEM)
              MPLOCS(N)=MM
              IF(IDTI.GT.0)DTIFG(N)=DTI(IEM)
              GOTO 900
            ENDIF
            IF(IBF.EQ.2)THEN
C             MM=0
              NN=IEM
              NLRLN=NLRL(IEM)
              CALL REPLAS(X(IEM,1),X(IEM,2),X(IEM,3),VX(IEM,1),
     >                    VX(IEM,2),VX(IEM,3),XP,YP,ZP,VPX,VPY,VPZ)
              SDT=DELT
C
              CALL MOVCHK
     I         (IBF,XP,YP,ZP,VPX,VPY,VPZ,N,MM,DELT,TMAX,RAMADA,IDTI,
     I          IZOOM, X,CP,IE,XPFG,CPFG,NFGMB,ISE, MAXNP,MAXEL,MXNPFG,
     M          MXKGL,MXNODE,MXNODS,EPSX,NEL,NEFG,NPFGS,NPFGEP,NODE,
     O          XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTIFG, ICODE)
              IF(ICODE.EQ.0)GOTO 900
c
              do ii=1,nwnp
                if(iem.eq.npw(ii))goto 900
              enddo
c
              GOTO 349
            ENDIF
          ENDIF
  120   continue
C
        SDT=DELT
C ----- CHECK IF POINT P IS ON ANY SIDE PLANE OF ELEMENT MM
        DO 130 LS=1,II
          DO IQ=1,4
            I=KIT(IQ,LS,ID)
            IF(I.EQ.0)THEN
              NI=0
            ELSE
              NI=IE(MM,I)
              XX(IQ)=X(NI,1)
              YY(IQ)=X(NI,2)
              ZZ(IQ)=X(NI,3)
            ENDIF
            IEMM(IQ)=NI
          ENDDO
C ----- CHECK IF POINT P IS ON THE PLANE COMPOSED WITH KK1,KK2,KK3,KK4
          CANG=FCOS(XX,YY,ZZ,XP,YP,ZP,0)
          IF(DABS(CANG).LE.EPSX*10.0D0)GOTO 140
  130   CONTINUE
C ----- POINT P IS NOT ON ANY SIDE PLANE OF ELEMENT MM, IT IS INSIDE
C       ELEMENT MM
        CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,MM,NODE,8,
     >              XX,YY,ZZ,VXX,VYY,VZZ)
C
        CALL BASE
c     I           (XX,YY,ZZ,XP,YP,ZP,MM,NODE,1,
     I           (XX,YY,ZZ,XP,YP,ZP,NODE,1,
     O            DL,DNX,DNX,DNX)
C
        CALL VALBDL(VXX,VYY,VZZ,DL,8,NODE,VPX,VPY,VPZ)
C
        CALL MOVCHK
     I    (IBF,XP,YP,ZP,VPX,VPY,VPZ,N,MM,DELT,TMAX,RAMADA,IDTI,IZOOM,
     I     X,CP,IE,XPFG,CPFG,NFGMB,ISE, MAXNP,MAXEL,MXNPFG,MXKGL,MXNODE,
     M     MXNODS,EPSX,NEL,NEFG,NPFGS,NPFGEP,NODE,XSFG,CSFG,MPLOCS,
     O     XEFG,CEFG,MPLOCE,DTIFG,  ICODE)
        IF(ICODE.EQ.0)GOTO 900
c
  131   continue
C
        IF(NODE.EQ.4)THEN
          CALL ELTRK4
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I        IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,
     I        IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSEIF(NODE.EQ.6)THEN
          CALL ELTRK6
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I        IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NZW,
     I        IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSE
          CALL ELTRK8
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I        IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NYW,NZW,
     I        IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4, XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ENDIF
C
        IF(SDT1.EQ.0.)GOTO 500
        IF(SDT1.NE.SDT)GOTO 260
C
C ----- ERROR OCCURRED
C
        if(idi.eq.0)then
          idi=1
          goto 131
        endif
        WRITE(16,*)'ERROR 1 AT HPTRAC: POINT',N,' CAN NOT BE',
     >            ' TRACKED'
        STOP
C
C ----- YES, POINT P IS ON A SIDE PLANE N1,N2,N3,N4
  140   CONTINUE
        N1=IEMM(1)
        N2=IEMM(2)
        N3=IEMM(3)
        N4=IEMM(4)
C ----- CHECK IF THIS POINT IS ON THE SIDE OF TETRAGON N1,N2,N3,N4
        CALL CKSIDE
     I      (X(1,1),X(1,2),X(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MAXNP,EPSX,
     O       KON,J1,J2)
        IF(KON.EQ.1)GOTO 150
C ----- NO, POINT P IS NOT ON ANY SIDE
C
        CALL ONPLAN
     I    (X(1,1),X(1,2),X(1,3),N1,N2,N3,MAXNP,
     M     XP,YP,ZP)
C
        CALL CKCNEL
     I             (LRL,NLRL,MM,N1,N2,N3,MAXNP,MXKBD,
     O              ME1,KOUNT)
C
  141   continue
c
        DO 145 I=1,KOUNT
          IF(I.EQ.2)MM=ME1
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,MM,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          CALL BASE
c     I               (XX,YY,ZZ,XP,YP,ZP,MM,NODE,1,
     I               (XX,YY,ZZ,XP,YP,ZP,NODE,1,
     O                DL,DNX,DNX,DNX)
C
          CALL VALBDL(VXX,VYY,VZZ,DL,8,NODE,VPX,VPY,VPZ)
C
          CALL MOVCHK
     I     (IBF,XP,YP,ZP,VPX,VPY,VPZ,N,MM,DELT,TMAX,RAMADA,IDTI,IZOOM,
     I      X,CP,IE,XPFG,CPFG,NFGMB,ISE,MAXNP,MAXEL,MXNPFG,MXKGL,MXNODE,
     M      MXNODS,EPSX,NEL,NEFG,NPFGS,NPFGEP,NODE,XSFG,CSFG,MPLOCS,
     O      XEFG,CEFG,MPLOCE,DTIFG,  ICODE)
          IF(ICODE.EQ.0)GOTO 900
C
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
  145   CONTINUE
C
C ----- CHECK IF P IS ON THE BOUNDARY PLANE N1,N2,N3,N4
C
        IF(KOUNT.EQ.1)THEN
C
C ----- FOR THE CASES THAT THIS BOUNDARY SIDE IS COMPOSED BY 4 POINTS
C       OR 3 POINTS
C
          IFIX=3
          CALL FIXCHK
     I    (IFIX,IDTI,N,N,MM,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
          IF(ICODE.EQ.0)GOTO 900
        ENDIF
C
C ----- ERROR OCCURRED
C
        if(idi.eq.0)then
          idi=1
          goto 141
        endif
        WRITE(16,*)'ERROR OCCURRED AT HPTRAC,ICODE=',ICODE
        STOP
C
  150   CONTINUE
C
C ----- YES, POINT P IS ON THE SIDE WITH NODES J1 AND J2
C
        CALL ONLINE
     I      (X(1,1),X(1,2),X(1,3),J1,J2,MAXNP,
     M       XP,YP,ZP)
        NLRL1=NLRL(J1)
        NLRL2=NLRL(J2)
        KOUNT=0
        DO 153 I1=1,NLRL1
          M1=LRL(I1,J1)
          DO I2=1,NLRL2
            M2=LRL(I2,J2)
            IF(M1.EQ.M2)THEN
              KOUNT=KOUNT+1
              MK(KOUNT)=M1
            ENDIF
          ENDDO
  153   CONTINUE
C
  154   continue
c
        DO 155 I=1,KOUNT
          MM=MK(I)
          CALL ELENOD
     I      (IE(MM,5),IE(MM,7),
     O       NODE,I1,I1)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,MM,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          CALL BASE
c     I             (XX,YY,ZZ,XP,YP,ZP,MM,NODE,1,
     I             (XX,YY,ZZ,XP,YP,ZP,NODE,1,
     O              DL,DNX,DNX,DNX)
C
          CALL VALBDL(VXX,VYY,VZZ,DL,8,NODE,VPX,VPY,VPZ)
C
          CALL MOVCHK
     I     (IBF,XP,YP,ZP,VPX,VPY,VPZ,N,MM,DELT,TMAX,RAMADA,IDTI,IZOOM,
     I      X,CP,IE,XPFG,CPFG,NFGMB,ISE,MAXNP,MAXEL,MXNPFG,MXKGL,MXNODE,
     M      MXNODS,EPSX,NEL,NEFG,NPFGS,NPFGEP,NODE,XSFG,CSFG,MPLOCS,
     O      XEFG,CEFG,MPLOCE,DTIFG,  ICODE)
          IF(ICODE.EQ.0)GOTO 900
C
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
C
  155   CONTINUE
C
        IFIX=4
        CALL FIXCHK
     I    (IFIX,IDTI,N,N,MM,J1,J2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
        IF(ICODE.EQ.0)THEN
          GO TO 900
        ELSE
          if(idi.eq.0)then
            idi=1
            goto 154
          endif
          WRITE(16,*)'ERROR IN FIXCHK, ICODE=',ICODE
          STOP
        ENDIF
C
C
C ##### START THE SUBSEQUENT TRACKINGS
C
  260   CONTINUE
        SDT=SDT1
        CALL REPLAS(XQ,YQ,ZQ,VQX,VQY,VQZ,XP,YP,ZP,VPX,VPY,VPZ)
c       NN1=N1
c       NN2=N2
c       NN3=N3
c       NN4=N4
C
C ----- CHECK IF POINT P COINSIDES WITH ANY GLOBAL NODES
C
        CALL CKCOIN
     I      (X(1,1),X(1,2),X(1,3),XP,YP,ZP,NLRL,NODE,N1,N2,N3,N4,MAXNP,
     >       EPSX, NN,NLRLN)
        IF(NN.EQ.0)GOTO 400
C
C ----- FOR THE CASES THAT POINT P IS A GLOBAL NODE
C
        CALL REPLAS(X(NN,1),X(NN,2),X(NN,3),VX(NN,1),VX(NN,2),VX(NN,3),
     >              XP,YP,ZP,VPX,VPY,VPZ)
c
        do 300 i=1,nwnp
          if(nn.ne.npw(i))goto 300
          if(ibf.eq.1)then
            dtreal=delt-sdt
            if(idti.eq.1 .or. idti.eq.2)dtifg(n)=1.0d0/dtreal
            if(idti.ne.1)then
              iprof=iwtyp(i)
              csfg(n)=wss(iprof)*dexp(-ramada*dtreal)
            endif
          endif
          goto 900
  300   continue
c
  349   CONTINUE
C
        DO 350 I=1,NLRLN
          MM=LRL(I,NN)
          CALL ELENOD
     I        (IE(MM,5),IE(MM,7),
     O         NODE,I1,I1)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,MM,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
C
  350   CONTINUE
C
C ----- CHECK IF THIS POINT IS A BOUNDARY NODE
C
        IFIX=2
        CALL FIXCHK
     I    (IFIX,IDTI,N,NN,MM,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
        IF(ICODE.EQ.0)THEN
          GO TO 900
        ELSE
          if(idi.eq.0)then
            idi=1
            goto 349
          endif
          WRITE(16,*)'ERROR IN FIXCHK, ICODE=',ICODE
          STOP
        ENDIF
C
C ----- FOR THE CASES THAT POINT P IS NOT A GLOBAL NODE
C
  400   CONTINUE
C
C ----- CHECK IF P IS ON ANY SIDES OF TETRAGON N1,N2,N3,N4
C
        CALL CKSIDE
     I      (X(1,1),X(1,2),X(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MAXNP,EPSX,
     O       KON,J1,J2)
        IF(KON.EQ.1)GOTO 491
C
C ----- NO, IT IS NOT ON ANY SIDE
        CALL ONPLAN
     I    (X(1,1),X(1,2),X(1,3),N1,N2,N3,MAXNP,
     M     XP,YP,ZP)
C
C ----- CHECK IF N1,N2,N3,N4 ARE ON A BOUNDARY PLANE
C
        CALL CKCNEL
     I              (LRL,NLRL,MM,N1,N2,N3,MAXNP,MXKBD,
     O               ME1,KOUNT)
        IF(KOUNT.EQ.1)GOTO 480
C
C ----- N1,N2,N3,N4 ARE NOT ON A BOUNDARY PLANE
C
C ----- ELEMENT M1 IS WHAT WE NEED, THEN CHECK WHICH PLANE POINT P IS
C       ONTO
C
        ME=MM
        MM=ME1
C
  431   CONTINUE
C
        CALL ELENOD
     I      (IE(MM,5),IE(MM,7),
     O       NODE,I1,I1)
        CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,MM,NODE,8,
     >              XX,YY,ZZ,VXX,VYY,VZZ)
C
  480   CONTINUE
C
        IF(NODE.EQ.4)THEN
          CALL ELTRK4
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,
     I        IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSEIF(NODE.EQ.6)THEN
          CALL ELTRK6
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NZW,
     I        IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4,XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ELSE
          CALL ELTRK8
     I       (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I        IE,IDETQ,MAXEL,EPSX,0,MXELW,MXNPW,NXW,NYW,NZW,
     I        IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i        idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M        N1,N2,N3,N4, XW,VXW,
     O        XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
        ENDIF
C
        IF(SDT1.EQ.0.)GOTO 500
        IF(SDT1.NE.SDT)GOTO 260
C
C ----- CHECK IF P IS ON THE BOUNDARY PLANE N1,N2,N3,N4
C
        IF(KOUNT.EQ.1)THEN
          IFIX=3
          CALL FIXCHK
     I      (IFIX,IDTI,N,N,MM,N1,N2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I       EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I       X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I       CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O       NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
          IF(ICODE.EQ.0)GOTO 900
        ENDIF
C
        IF(M.EQ.ME1)THEN
          M=ME
          GOTO 431
        ENDIF
C
C ----- ERROR OCCURRED
C
        if(idi.eq.0)then
          idi=1
          m=me1
          goto 431
        endif
        WRITE(16,*)'ERROR OCCURRED AT HPTRAC'
        STOP
C
C *** YES, IT IS ON J1,J2 SIDE
C
  491   CONTINUE
        CALL ONLINE
     I      (X(1,1),X(1,2),X(1,3),J1,J2,MAXNP,
     M       XP,YP,ZP)
C
        NLRL1=NLRL(J1)
        NLRL2=NLRL(J2)
        KOUNT=0
        DO 493 I1=1,NLRL1
          M1=LRL(I1,J1)
          DO 492 I2=1,NLRL2
            M2=LRL(I2,J2)
            IF(M1.EQ.M2)THEN
              KOUNT=KOUNT+1
              MK(KOUNT)=M1
            ENDIF
  492     CONTINUE
  493   CONTINUE
C
  494   continue
c
        DO 499 I=1,KOUNT
          MM=MK(I)
          CALL ELENOD
     I        (IE(MM,5),IE(MM,7),
     O         NODE,I1,I1)
          CALL WRKARY(IE,X,VX,MAXEL,11,MAXNP,MM,NODE,8,
     >                XX,YY,ZZ,VXX,VYY,VZZ)
C
          IF(NODE.EQ.4)THEN
            CALL ELTRK4
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,
     I          IBW(1,1),NLRLW(1,1),LRLW(1,1,1),IEW(1,1,1),DL468(1,1,1),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSEIF(NODE.EQ.6)THEN
            CALL ELTRK6
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NZW,
     I          IBW(1,2),NLRLW(1,2),LRLW(1,1,2),IEW(1,1,2),DL468(1,1,2),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4,XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ELSE
            CALL ELTRK8
     I         (XP,YP,ZP,VPX,VPY,VPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDT,MM,IBF,
     I          IE,IDETQ,MAXEL,EPSX,1,MXELW,MXNPW,NXW,NYW,NZW,
     I          IBW(1,3),NLRLW(1,3),LRLW(1,1,3),IEW(1,1,3),DL468(1,1,3),
     i          idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M          N1,N2,N3,N4, XW,VXW,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1,npfgs,xsfg,csfg,mplocs)
          ENDIF
C
          IF(SDT1.EQ.0.)GOTO 500
          IF(SDT1.NE.SDT)GOTO 260
C
  499   CONTINUE
C
        IFIX=4
        CALL FIXCHK
     I    (IFIX,IDTI,N,N,MM,J1,J2,N3,N4,NODE,IBF,XP,YP,ZP,VPX,VPY,VPZ,
     I      EPSX,DELT,SDT,SDT1,NNP,NEL,NEFG,IDETQ,NTUBS,NPFGEP,
     I      X,IE,IB,VX,MAXNP,MAXEL,MAXBES,MXTUBS,MXKGL,MXNPFG,MXNEP,
     I      CP,ISB,DCOSB,NBDYB,IBDY,XPFG,CPFG,ISE,NFGMB,IZOOM,RAMADA,
     O      NPFGS,CS,XSFG,CSFG,MPLOCS,XEFG,CEFG,MPLOCE,DTI,DTIFG,ICODE)
        IF(ICODE.EQ.0)THEN
          GO TO 900
        ELSE
          if(idi.eq.0)then
            idi=1
            goto 494
          endif
          WRITE(16,*)'ERROR IN FIXCHK, ICODE=',ICODE
          STOP
        ENDIF
C
C ***** RECORDING INFORMATION
C
  500   CONTINUE
C
C ----- FOR THE CASES OF BACKWARD TRACKING
C
        IF(IBF.EQ.1)THEN
          IF(IDTI.NE.1)THEN
            CALL INTERP
     I             (MAXNP,MAXEL,MXNPFG,MXKGL,8,1,NODE,MM,XQ,YQ,ZQ,X,CP,
     I              IE,XPFG,CPFG,NEL,NEFG,IZOOM,10,ISE,NFGMB,1,1,EPSX,
     O              CQ)
            CSFG(N)=CQ(1)*DEXP(-RAMADA*DELT)
            IF(IDTI.EQ.2)DTIFG(N)=1.0D0/DELT
          ELSE
            DTIFG(N)=1.0D0/DELT
          ENDIF
          GOTO 900
        ENDIF
C
C ----- FOR FORWARD TRACKING NODES
C
        IE(MM,11)=MM
C
        CALL ELENOD
     I      (IE(MM,5),IE(MM,7),
     O       NODE1,KK,KK)
         DO KK=1,NODE1
           IEMKK=IE(MM,KK)
           XR=X(IEMKK,1)
           YR=X(IEMKK,2)
           ZR=X(IEMKK,3)
           DIST=DSQRT((XQ-XR)**2+(YQ-YR)**2+(ZQ-ZR)**2)
           IF(DIST.LE.EPSX)THEN
             CS(IEMKK)=CPFG(N)*DEXP(-RAMADA*DELT)
             GOTO 506
           ENDIF
         ENDDO
  506    CONTINUE
C
        IF(MPLOCE(N).EQ.0)THEN
          NPFGS=NPFGS+1
          CALL WARMSG(NPFGS,MXNPFG,'HPTRAC','MXNPFG',1)
          XSFG(NPFGS,1)=XQ
          XSFG(NPFGS,2)=YQ
          XSFG(NPFGS,3)=ZQ
          CSFG(NPFGS)=CPFG(N)*DEXP(-RAMADA*DELT)
          MPLOCS(NPFGS)=MM
        ELSE
          NPFGEP=NPFGEP+1
          CALL WARMSG(NPFGEP,MXNEP,'HPTRAC','MXNEP ',2)
          XEFG(NPFGEP,1)=XQ
          XEFG(NPFGEP,2)=YQ
          XEFG(NPFGEP,3)=ZQ
          CEFG(NPFGEP)=CPFG(N)*DEXP(-RAMADA*DELT)
          MPLOCE(NPFGEP)=MM
        ENDIF
C
  900 CONTINUE
C
  990 FORMAT(3F15.4)
  992 FORMAT(I5,3X,2F15.4)
C
      RETURN
      END
C
C
C
      SUBROUTINE ELTRK4
     I          (XPP,YPP,ZPP,VPPX,VPPY,VPPZ,XX,YY,ZZ,VXX,VYY,VZZ,
     I           SDTT,M,IBF,IE,IDETQ,MAXEL,EPSX,ICALL,
     I           MXELW,MXNPW,NXW,IBW,NLRLW,LRLW,IEW,DL468,
     i           idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M           N1,N2,N3,N4,XW,VXW,
     O           XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1,npfgs,xsfg,csfg,mplocs)
C     10/10/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),MK(24)
      DIMENSION XW(MXNPW,3),VXW(MXNPW,3),IEW(MXELW,8),IE(MAXEL,11)
      DIMENSION NLRLW(MXNPW),LRLW(24,MXNPW),IBW(MXNPW),KIT(4,4),NNA(3)
      DIMENSION XXX(8),YYY(8),ZZZ(8),VXXX(8),VYYY(8),VZZZ(8),KIN(4)
      DIMENSION KNC(3),DL468(8,MXNPW)
      dimension xsfg(mxnpfg,3),csfg(mxnpfg),mplocs(mxnpfg)
C
      DATA KIT/4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0/,LS/4/,NODE/4/
C
      NELW=NXW*NXW*NXW
      NPW=(NXW+1)*(NXW+2)*(NXW+3)/6
      SDT=SDTT
      CALL REPLAS(XPP,YPP,ZPP,VPPX,VPPY,VPPZ,XP,YP,ZP,VPX,VPY,VPZ)
C
      IF(VPX.EQ.0.0D0 .AND. VPY.EQ.0.0D0 .AND. VPZ.EQ.0.0D0)THEN
        XQ=XP
        YQ=YP
        ZQ=ZP
        SDTT1=SDT
        RETURN
      ENDIF
C
C***** PREPARE INFORMATION FOR ELEMENT TRACKING
C
      DO 60 N=1,NPW
C ----- DETERMINE WORKING ARRAYS XW, YW, ZW, VXW, VYW, AND VZW
        CALL VALBDL(XX,YY,ZZ,DL468(1,N),8,NODE,
     >              XW(N,1),XW(N,2),XW(N,3))
        CALL VALBDL
     >             (VXX,VYY,VZZ,DL468(1,N),8,NODE,
     >              VXW(N,1),VXW(N,2),VXW(N,3))
        DO II=1,3
          XW(N,II)=XW(N,II)*1.0D-6
          VXW(N,II)=VXW(N,II)*1.0D-6
        ENDDO
  60  CONTINUE
C
C ***** DETERMINE THE SUBELEMENT AT WHICH POINT P LOCATES
C
      CALL MMLOC
     I    (XW,IEW,VXW,NLRLW,KIT,NELW,MXNPW,MXELW,LS,EPSX,NODE,
     I     XP,YP,ZP,VPX,VPY,VPZ,
     O     KIN,MM,NN,NLRLN,ICODE)
C
C ----- START TRACKING
C
      NA1=N1
      NA2=N2
      NA3=N3
      NA4=N4
      DO I=1,NODE
        IF(NA1.EQ.IE(M,I))NNA(1)=I
        IF(NA2.EQ.IE(M,I))NNA(2)=I
        IF(NA3.EQ.IE(M,I))NNA(3)=I
      ENDDO
C
C ----- POINT P COINCIDES WITH ANY NODES OF WORKING ELEMENT MM
C
      IF(ICODE.EQ.0)GOTO 290
C
C----- CHECK IF POINT P IS ON ANY SIDE PLANE OF WORKING ELEMENT MM
C
      KN=0
      DO I=1,LS
        IF(KIN(I).EQ.1)THEN
          KN=KN+1
          KNC(KN)=I
        ENDIF
      ENDDO
      IF(KN.EQ.1 .OR. (KN.GT.1 .AND. ICALL.EQ.1))THEN
        I=KNC(1)
        GOTO 140
      ENDIF
      IF(KN.GT.1)THEN
        DO I=1,3
          II=NNA(I)
          XXX(I)=XX(II)
          YYY(I)=YY(II)
          ZZZ(I)=ZZ(II)
        ENDDO
        DO 100 J=1,KN
          KS=KNC(J)
          DO I=1,3
            II=KIT(I,KS)
            II=IEW(MM,II)
            CANG=FCOS(XXX,YYY,ZZZ,XW(II,1),XW(II,2),XW(II,3),0)
            IF(DABS(CANG).GT.EPSX)GOTO 100
          ENDDO
          I=KS
          GOTO 140
  100   CONTINUE
        I=KNC(1)
        GOTO 140
c       WRITE(16,*)'ERROR OCCURRED AT ELTRK4 --- AFTER MMLOC'
c       STOP
      ENDIF
C
C ***** POINT P IS NOT ON ANY SIDE PLANE OF ELEMENT MM, IT IS INSIDE
C ***** ELEMENT MM
C
      ID=0
      CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O            XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  131 CONTINUE
C
      CALL TRAK2T
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             IBF,IEW,IJUDGE,MXELW,ID,MM, EPSX, J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
      IF(SDT1.EQ.0.)THEN
        SDTT1=0.0D0
        RETURN
      ENDIF
      IF(SDT1.NE.SDT)GOTO 250
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 131
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 1 IN ELTRK4: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C *** YES, POINT P IS ON A SIDE PLANE N1,N2,N3,N4
C *** CHECK IF THIS POINT IS ON THE SIDE OF TETRAGON N1,N2,N3,N4
C
  140 CONTINUE
      N1=KIT(1,I)
      N2=KIT(2,I)
      N3=KIT(3,I)
      N1=IEW(MM,N1)
      N2=IEW(MM,N2)
      N3=IEW(MM,N3)
C
      CALL ONPLAN
     I    (XW(1,1),XW(1,2),XW(1,3),N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C
      CALL CKSIDE
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MXNPW,EPSX,
     O     KON,J1,J2)
      IF(KON.EQ.1)GOTO 150
C
C *** NO, POINT P IS NOT ON ANY SIDE
C
      ID=1
      CALL CKCNEL
     I            (LRLW,NLRLW,MM,N1,N2,N3,MXNPW,24,
     O             ME1,KOUNT)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  141 CONTINUE
C
      DO 145 I=1,KOUNT
        IF(I.EQ.2)MM=ME1
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2T
     I          (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I           IBF,IEW,IJUDGE,MXELW,ID,MM, EPSX, J1,J2, idi,
     M           N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  145 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 141
      ENDIF
      IF(KOUNT.EQ.1)THEN
        N1=NA1
        N2=NA2
        N3=NA3
        N4=NA4
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 2 AT ELTRK4: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
  150 CONTINUE
C
C *** YES, POINT P IS ON THE SIDE WITH NODES J1 AND J2
C
      CALL ONLINE
     I      (XW(1,1),XW(1,2),XW(1,3),J1,J2,MXNPW,
     M       XP,YP,ZP)
C
      J12=J1
  151 CONTINUE
      IF(DABS(XP-XW(J12,1)).LE.EPSX .AND. DABS(YP-XW(J12,2)).LE.EPSX
     >    .AND. DABS(ZP-XW(J12,3)).LE.EPSX)THEN
        CALL REPLAS(XW(J12,1),XW(J12,2),XW(J12,3),VXW(J12,1),VXW(J12,2),
     >              VXW(J12,3),XP,YP,ZP,VPX,VPY,VPZ)
        NLRLN=NLRLW(J12)
        NN=J12
        GOTO 290
      ELSEIF(J12.EQ.J1)THEN
        J12=J2
        GOTO 151
      ENDIF
C
      ID=2
      NLRL1=NLRLW(J1)
      NLRL2=NLRLW(J2)
      KOUNT=0
      DO 153 I1=1,NLRL1
        M1=LRLW(I1,J1)
        DO 152 I2=1,NLRL2
          M2=LRLW(I2,J2)
          IF(M1.EQ.M2)THEN
            KOUNT=KOUNT+1
            MK(KOUNT)=M1
          ENDIF
  152   CONTINUE
  153 CONTINUE
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  154 CONTINUE
C
      DO 155 I=1,KOUNT
        MM=MK(I)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2T
     I          (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I           IBF,IEW,IJUDGE,MXELW,ID,MM, EPSX, J1,J2, idi,
     M           N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  155 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 154
      ENDIF
C
      N1=NA1
      N2=NA2
      N3=NA3
      N4=NA4
      CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
      SDTT1=SDT
c     CALL IBWCK4
c    I    (IBW(J1),IBW(J2),XP,YP,ZP,VPX,VPY,VPX,NXW,M,IE,SDT,MAXEL,
c    O     N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
      RETURN
C
C $$$$$ START CONSEQUENT TRACKINGS
C
  250 CONTINUE
      if(idps.eq.1)then
        npfgs=npfgs+1
        xsfg(npfgs,1)=0.5d0*(xp+xq)
        xsfg(npfgs,2)=0.5d0*(yp+yq)
        xsfg(npfgs,3)=0.5d0*(zp+zq)
        mplocs(npfgs)=m
        dtreal=delt-0.5d0*(sdt+sdt1)
        xsi=(delt-dtreal)/delt
        cc=cpsp*xsi+cps*(1.0d0-xsi)
        csfg(npfgs)=cc*dexp(-ramada*dtreal)
      endif
      SDT=SDT1
      CALL REPLAS(XQ,YQ,ZQ,VQX,VQY,VQZ,XP,YP,ZP,VPX,VPY,VPZ)
C
C ***** CHECK IF POINT P COINSIDES WITH ANY GLOBAL NODES
C
      CALL CKCOIN
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NLRLW,NODE,
     >     N1,N2,N3,N4,MXNPW,EPSX,  NN,NLRLN)
      IF(NN.EQ.0)GOTO 400
C
C ***** YES, P COINSIDES WITH NODE NN, THEN PREPARE DATA FOR TRACK1
C
      CALL REPLAS(XW(NN,1),XW(NN,2),XW(NN,3),VXW(NN,1),VXW(NN,2),
     >            VXW(NN,3),  XP,YP,ZP,VPX,VPY,VPZ)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
  290 CONTINUE
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  300 CONTINUE
C
      DO 350 I=1,NLRLN
        MM=LRLW(I,NN)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
        DO 305 J=1,NODE
          IF(NN.EQ.IEW(MM,J))GOTO 310
  305   CONTINUE
C
  310   CONTINUE
C
        CALL TRAK1T
     I          (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I           IBF,IEW,IJUDGE,J,MM,MXELW, EPSX, idi,
     O           XQ,YQ,ZQ,VQX,VQY,VQZ,N1,N2,N3,N4,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  350 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 300
      ENDIF
C
      IF(IBW(NN).NE.0)THEN
        IF(IBW(NN).EQ.1 .OR. IBW(NN).EQ.34 .OR. IBW(NN).EQ.42 .OR.
     >     IBW(NN).EQ.23 .OR. IBW(NN).EQ.132 .OR. IBW(NN).EQ.143 .OR.
     >     IBW(NN).EQ.124)THEN
          I=1
        ELSEIF(IBW(NN).EQ.2 .OR. IBW(NN).EQ.13 .OR. IBW(NN).EQ.14 .OR.
     >     IBW(NN).EQ.234)THEN
          I=2
        ELSEIF(IBW(NN).EQ.3)THEN
          I=3
        ELSE
          I=4
        ENDIF
        N1=KIT(1,I)
        N2=KIT(2,I)
        N3=KIT(3,I)
        N1=IE(M,N1)
        N2=IE(M,N2)
        N3=IE(M,N3)
        N4=0
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 4 AT ELTRK4: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C ***** NO, P DOESN'T COINSIDE WITH ANY NODES, THEN PREPARE DATA FOR
C       TRACK2
C
  400 CONTINUE
C
C ***** CHECK IF P IS ON ANY SIDES OF TETRAGON N1,N2,N3,N4
C
      CALL ONPLAN
     I    (XW(1,1),XW(1,2),XW(1,3),N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C
      CALL CKSIDE
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MXNPW,EPSX,
     O     KON,J1,J2)
      IF(KON.EQ.1)GOTO 491
C
C *** NO, IT IS NOT ON ANY SIDE.
C
C === CHECK IF N1,N2,N3,N4 ARE ON A BOUNDARY PLANE
C
      CALL CKCNEL
     I            (LRLW,NLRLW,MM,N1,N2,N3,MXNPW,24,
     O             ME1,KOUNT)
      IF(KOUNT.EQ.1)GOTO 490
C
C --- N1,N2,N3,N4 ARE NOT ON A BOUNDARY PLANE. THEN CHECK WHICH PLANE
C --- POINT P IS ONTO
C
      ME=MM
      MM=ME1
C
  431 CONTINUE
C
      CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O            XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
  490 CONTINUE
C
        ID=1
C
C === DETERMINE HOW TO TRACK THIS POINT
C
  484 CONTINUE
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  485 CONTINUE
C
      CALL TRAK2T
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             IBF,IEW,IJUDGE,MXELW,ID,MM, EPSX, J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
      IF(SDT1.EQ.0.)THEN
        SDTT1=0.0D0
        RETURN
      ENDIF
      IF(SDT1.NE.SDT)GOTO 250
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 485
      ENDIF
      IF(KOUNT.EQ.1)THEN
        CALL BDYPLN
     I        (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,1,
     I         XW(1,1),XW(1,2),XW(1,3),
     O         N1,N2,N3,N4)
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
C
      IF(MM.EQ.ME1)THEN
        MM=ME
        GOTO 431
      ENDIF
C
      WRITE(16,*)'ERROR MESSAGE 5 AT ELTRK4 : POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C ***** YES
C
  491 CONTINUE
C
      CALL ONLINE
     I      (XW(1,1),XW(1,2),XW(1,3),J1,J2,MXNPW,
     M       XP,YP,ZP)
C
      J12=J1
  495 CONTINUE
      IF(DABS(XP-XW(J12,1)).LE.EPSX .AND. DABS(YP-XW(J12,2)).LE.EPSX
     >    .AND. DABS(ZP-XW(J12,3)).LE.EPSX)THEN
        CALL REPLAS(XW(J12,1),XW(J12,2),XW(J12,3),VXW(J12,1),VXW(J12,2),
     >              VXW(J12,3),XP,YP,ZP,VPX,VPY,VPZ)
        NLRLN=NLRLW(J12)
        NN=J12
        GOTO 290
      ELSEIF(J12.EQ.J1)THEN
        J12=J2
        GOTO 495
      ENDIF
C
      ID=2
      NLRL1=NLRLW(J1)
      NLRL2=NLRLW(J2)
      KOUNT=0
      DO 493 I1=1,NLRL1
        M1=LRLW(I1,J1)
        DO 492 I2=1,NLRL2
          M2=LRLW(I2,J2)
          IF(M1.EQ.M2)THEN
            KOUNT=KOUNT+1
            MK(KOUNT)=M1
          ENDIF
  492   CONTINUE
  493 CONTINUE
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  494 CONTINUE
C
      DO 499 I=1,KOUNT
        MM=MK(I)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2T
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             IBF,IEW,IJUDGE,MXELW,ID,MM, EPSX, J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VXQ,VYQ,VZQ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  499 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 494
      ENDIF
C
      N1=J1
      N2=J2
      CALL BDYPLN
     I         (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,2,
     I          XW(1,1),XW(1,2),XW(1,3),
     O          N1,N2,N3,N4)
      CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
      SDTT1=SDT
c     CALL IBWCK4
c    I    (IBW(J1),IBW(J2),XP,YP,ZP,VPX,VPY,VPX,NXW,M,IE,SDT,MAXEL,
c    O     N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
C
      RETURN
      END
C
c
c
      SUBROUTINE ELTRK6
     I           (XPP,YPP,ZPP,VPPX,VPPY,VPPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDTT,
     I            M,IBF,IE,IDETQ,MAXEL, EPSX,ICALL,
     I            MXELW,MXNPW,NXW,NZW,IBW,NLRLW,LRLW,IEW,DL468,
     i            idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M            N1,N2,N3,N4,XW,VXW,
     O            XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1,npfgs,xsfg,csfg,mplocs)
C
C        10/10/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),MK(10)
      DIMENSION XW(MXNPW,3),VXW(MXNPW,3),IEW(MXELW,8),IE(MAXEL,11)
      DIMENSION NLRLW(MXNPW),LRLW(24,MXNPW),IBW(MXNPW),KIN(5),NNA(3)
      DIMENSION XXX(8),YYY(8),ZZZ(8),VXXX(8),VYYY(8),VZZZ(8),KIT(4,5)
      DIMENSION KNC(3),DL468(8,MXNPW)
      dimension xsfg(mxnpfg,3),csfg(mxnpfg),mplocs(mxnpfg)
c
      DATA KIT/1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0/
      DATA LS/5/,NODE/6/
C
C
      NELW=NXW*NXW*NZW
      NPW=(NXW+1)*(NXW+2)*(NZW+1)/2
      SDT=SDTT
      CALL REPLAS(XPP,YPP,ZPP,VPPX,VPPY,VPPZ,XP,YP,ZP,VPX,VPY,VPZ)
      IF(VPX.EQ.0.0D0 .AND. VPY.EQ.0.0D0 .AND. VPZ.EQ.0.0D0)THEN
        XQ=XP
        YQ=YP
        ZQ=ZP
        SDTT1=SDT
        RETURN
      ENDIF
C
C***** PREPARE INFORMATION FOR ELEMENT TRACKING
C
      DO 60 N=1,NPW
        CALL VALBDL(XX,YY,ZZ,DL468(1,N),8,NODE,
     >              XW(N,1),XW(N,2),XW(N,3))
        CALL VALBDL
     >             (VXX,VYY,VZZ,DL468(1,N),8,NODE,
     >              VXW(N,1),VXW(N,2),VXW(N,3))
        DO II=1,3
          XW(N,II)=XW(N,II)*1.0D-6
          VXW(N,II)=VXW(N,II)*1.0D-6
        ENDDO
  60  CONTINUE
C
C ***** DETERMINE THE SUBELEMENT AT WHICH POINT P LOCATES
C
      CALL MMLOC
     I    (XW,IEW,VXW,NLRLW,KIT,NELW,MXNPW,MXELW,LS,EPSX,NODE,
     I     XP,YP,ZP,VPX,VPY,VPZ,
     O     KIN,MM,NN,NLRLN,ICODE)
C
C ----- START TRACKING
C
      NA1=N1
      NA2=N2
      NA3=N3
      NA4=N4
      DO I=1,NODE
        IF(NA1.EQ.IE(M,I))NNA(1)=I
        IF(NA2.EQ.IE(M,I))NNA(2)=I
        IF(NA3.EQ.IE(M,I))NNA(3)=I
      ENDDO
C
C ----- POINT P COINCIDES WITH ANY NODES OF WORKING ELEMENT MM
      IF(ICODE.EQ.0)GOTO 290
C
C----- CHECK IF POINT P IS ON ANY SIDE PLANE OF WORKING ELEMENT MM
C
      KN=0
      DO I=1,LS
        IF(KIN(I).EQ.1)THEN
          KN=KN+1
          KNC(KN)=I
        ENDIF
      ENDDO
      IF(KN.EQ.1 .OR. (KN.GT.1 .AND. ICALL.EQ.1))THEN
        I=KNC(1)
        GOTO 140
      ENDIF
      IF(KN.GT.1)THEN
        DO I=1,3
          II=NNA(I)
          XXX(I)=XX(II)
          YYY(I)=YY(II)
          ZZZ(I)=ZZ(II)
        ENDDO
        DO 100 J=1,KN
          KS=KNC(J)
          DO I=1,3
            II=KIT(I,KS)
            II=IEW(MM,II)
            CANG=FCOS(XXX,YYY,ZZZ,XW(II,1),XW(II,2),XW(II,3),0)
            IF(DABS(CANG).GT.EPSX)GOTO 100
          ENDDO
          I=KS
          GOTO 140
  100   CONTINUE
        I=KNC(1)
        GOTO 140
c       WRITE(16,*)'ERROR OCCURRED AT ELTRK6 --- AFTER MMLOC'
c       STOP
      ENDIF
C
C ***** POINT P IS NOT ON ANY SIDE PLANE OF ELEMENT MM, IT IS INSIDE
C ***** ELEMENT MM
C
      ID=0
      CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O            XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  131 CONTINUE
C
      CALL TRAK2P
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
      IF(SDT1.EQ.0.)THEN
        SDTT1=0.0D0
        RETURN
      ENDIF
      IF(SDT1.NE.SDT)GOTO 250
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 131
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 1 IN ELTRK6: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
 
  140 CONTINUE
      N1=KIT(1,I)
      N2=KIT(2,I)
      N3=KIT(3,I)
      N4=KIT(4,I)
      N1=IEW(MM,N1)
      N2=IEW(MM,N2)
      N3=IEW(MM,N3)
      IF(N4.NE.0)N4=IEW(MM,N4)
C
C *** YES, POINT P IS ON A SIDE PLANE N1,N2,N3,N4
C *** CHECK IF THIS POINT IS ON THE SIDE OF TETRAGON N1,N2,N3,N4
C
      CALL ONPLAN
     I    (XW(1,1),XW(1,2),XW(1,3),N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C
      CALL CKSIDE
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MXNPW,EPSX,
     O     KON,J1,J2)
      IF(KON.EQ.1)GOTO 150
C
C *** NO, POINT P IS NOT ON ANY SIDE
C
      ID=1
      CALL CKCNEL
     I            (LRLW,NLRLW,MM,N1,N2,N3,MXNPW,24,
     O             ME1,KOUNT)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  141 CONTINUE
C
      DO 145 I=1,KOUNT
        IF(I.EQ.2)MM=ME1
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2P
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M          N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  145 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 141
      ENDIF
C
      IF(KOUNT.EQ.1)THEN
        N1=NA1
        N2=NA2
        N3=NA3
        N4=NA4
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 2 AT ELTRK6: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
  150 CONTINUE
C
C *** YES, POINT P IS ON THE SIDE WITH NODES J1 AND J2
C
      CALL ONLINE
     I      (XW(1,1),XW(1,2),XW(1,3),J1,J2,MXNPW,
     M       XP,YP,ZP)
C
      J12=J1
  151 CONTINUE
      IF(DABS(XP-XW(J12,1)).LE.EPSX .AND. DABS(YP-XW(J12,2)).LE.EPSX
     >    .AND. DABS(ZP-XW(J12,3)).LE.EPSX)THEN
        CALL REPLAS(XW(J12,1),XW(J12,2),XW(J12,3),VXW(J12,1),VXW(J12,2),
     >              VXW(J12,3),XP,YP,ZP,VPX,VPY,VPZ)
        NLRLN=NLRLW(J12)
        NN=J12
        GOTO 290
      ELSEIF(J12.EQ.J1)THEN
        J12=J2
        GOTO 151
      ENDIF
C
      ID=2
      NLRL1=NLRLW(J1)
      NLRL2=NLRLW(J2)
      KOUNT=0
      DO 153 I1=1,NLRL1
        M1=LRLW(I1,J1)
        DO 152 I2=1,NLRL2
          M2=LRLW(I2,J2)
          IF(M1.EQ.M2)THEN
            KOUNT=KOUNT+1
            MK(KOUNT)=M1
          ENDIF
  152   CONTINUE
  153 CONTINUE
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  154 CONTINUE
C
      DO 155 I=1,KOUNT
        MM=MK(I)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2P
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX, J1,J2, idi,
     M          N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  155 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 154
      ENDIF
C
      N1=NA1
      N2=NA2
      N3=NA3
      N4=NA4
      CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
      SDTT1=SDT
c     CALL IBWCK6
c    I    (IBW(J1),IBW(J2),XP,YP,ZP,VPX,VPY,VPZ,NXW,NZW,M,IE,SDT,MAXEL,
c    O     N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
      RETURN
C
C $$$$$ START CONSEQUENT TRACKINGS
C
  250 CONTINUE
      if(idps.eq.1)then
        npfgs=npfgs+1
        xsfg(npfgs,1)=0.5d0*(xp+xq)
        xsfg(npfgs,2)=0.5d0*(yp+yq)
        xsfg(npfgs,3)=0.5d0*(zp+zq)
        mplocs(npfgs)=m
        dtreal=delt-0.5d0*(sdt+sdt1)
        xsi=(delt-dtreal)/delt
        cc=cpsp*xsi+cps*(1.0d0-xsi)
        csfg(npfgs)=cc*dexp(-ramada*dtreal)
      endif
      SDT=SDT1
      CALL REPLAS(XQ,YQ,ZQ,VQX,VQY,VQZ,XP,YP,ZP,VPX,VPY,VPZ)
C
C ***** CHECK IF POINT P COINSIDES WITH ANY GLOBAL NODES
C
      CALL CKCOIN
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NLRLW,NODE,
     O     N1,N2,N3,N4,MXNPW,EPSX,  NN,NLRLN)
      IF(NN.EQ.0)GOTO 400
C
C ***** YES, P COINSIDES WITH NODE NN, THEN PREPARE DATA FOR TRACK1
C
      CALL REPLAS(XW(NN,1),XW(NN,2),XW(NN,3),VXW(NN,1),VXW(NN,2),
     >            VXW(NN,3),  XP,YP,ZP,VPX,VPY,VPZ)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
  290 CONTINUE
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  300 CONTINUE
C
      DO 350 I=1,NLRLN
        MM=LRLW(I,NN)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
        DO 305 J=1,NODE
          IF(NN.EQ.IEW(MM,J))GOTO 310
  305   CONTINUE
C
  310   CONTINUE
C
        CALL TRAK1P
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,J,MXELW, EPSX, idi,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,N1,N2,N3,N4,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  350 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 300
      ENDIF
C
      IF(IBW(NN).NE.0)THEN
        IF(IBW(NN).EQ.1 .OR. IBW(NN).EQ.12 .OR. IBW(NN).EQ.13 .OR.
     >     IBW(NN).EQ.14 .OR. IBW(NN).EQ.15 .OR. IBW(NN).EQ.124 .OR.
     >     IBW(NN).EQ.125.OR. IBW(NN).EQ.134 .OR. IBW(NN).EQ.135)THEN
          I=1
        ELSEIF(IBW(NN).EQ.2 .OR. IBW(NN).EQ.23 .OR. IBW(NN).EQ.24 .OR.
     >     IBW(NN).EQ.25 .OR. IBW(NN).EQ.234 .OR. IBW(NN).EQ.235)THEN
          I=2
        ELSEIF(IBW(NN).EQ.3 .OR. IBW(NN).EQ.34 .OR. IBW(NN).EQ.35)THEN
          I=3
        ELSEIF(IBW(NN).EQ.4)THEN
          I=4
        ELSE
          I=5
        ENDIF
        N1=KIT(1,I)
        N2=KIT(2,I)
        N3=KIT(3,I)
        N4=KIT(4,I)
        N1=IE(M,N1)
        N2=IE(M,N2)
        N3=IE(M,N3)
        N4=IE(M,N4)
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 2 AT ELTRK6: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C ***** NO, P DOESN'T COINSIDE WITH ANY NODES, THEN PREPARE DATA FOR
C       TRACK2
C
  400 CONTINUE
C
C ***** CHECK IF P IS ON ANY SIDES OF TETRAGON N1,N2,N3,N4
C
      CALL ONPLAN
     I    (XW(1,1),XW(1,2),XW(1,3),N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C
      CALL CKSIDE
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MXNPW,EPSX,
     O     KON,J1,J2)
      IF(KON.EQ.1)GOTO 491
C
C *** NO, IT IS NOT ON ANY SIDE.
C
C === CHECK IF N1,N2,N3,N4 ARE ON A BOUNDARY PLANE
C
      CALL CKCNEL
     I            (LRLW,NLRLW,MM,N1,N2,N3,MXNPW,24,
     O             ME1,KOUNT)
      IF(KOUNT.EQ.1)GOTO 490
C
C --- N1,N2,N3,N4 ARE NOT ON A BOUNDARY PLANE. THEN CHECK WHICH PLANE
C --- POINT P IS ONTO
C
      ME=MM
      MM=ME1
C
  431 CONTINUE
C
      CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O            XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
  490 CONTINUE
C
        ID=1
C
C === DETERMINE HOW TO TRACK THIS POINT
C
  484 CONTINUE
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  485 CONTINUE
C
      CALL TRAK2P
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
      IF(SDT1.EQ.0.)THEN
        SDTT1=0.0D0
        RETURN
      ENDIF
      IF(SDT1.NE.SDT)GOTO 250
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 485
      ENDIF
      IF(KOUNT.EQ.1)THEN
        CALL BDYPLN
     I        (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,1,
     I         XW(1,1),XW(1,2),XW(1,3),
     O         N1,N2,N3,N4)
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
C
      IF(MM.EQ.ME1)THEN
        MM=ME
        GOTO 431
      ENDIF
C
      WRITE(16,*)'ERROR MESSAGE 3 AT ELTRK6 --- POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C ***** YES
C
  491 CONTINUE
C
      CALL ONLINE
     I      (XW(1,1),XW(1,2),XW(1,3),J1,J2,MXNPW,
     M       XP,YP,ZP)
C
      J12=J1
  495 CONTINUE
      IF(DABS(XP-XW(J12,1)).LE.EPSX .AND. DABS(YP-XW(J12,2)).LE.EPSX
     >    .AND. DABS(ZP-XW(J12,3)).LE.EPSX)THEN
        CALL REPLAS(XW(J12,1),XW(J12,2),XW(J12,3),VXW(J12,1),VXW(J12,2),
     >              VXW(J12,3),XP,YP,ZP,VPX,VPY,VPZ)
        NLRLN=NLRLW(J12)
        NN=J12
        GOTO 290
      ELSEIF(J12.EQ.J1)THEN
        J12=J2
        GOTO 495
      ENDIF
C
      ID=2
      NLRL1=NLRLW(J1)
      NLRL2=NLRLW(J2)
      KOUNT=0
      DO 493 I1=1,NLRL1
        M1=LRLW(I1,J1)
        DO 492 I2=1,NLRL2
          M2=LRLW(I2,J2)
          IF(M1.EQ.M2)THEN
            KOUNT=KOUNT+1
            MK(KOUNT)=M1
          ENDIF
  492   CONTINUE
  493 CONTINUE
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  494 CONTINUE
C
      DO 499 I=1,KOUNT
        MM=MK(I)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2P
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M          N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  499 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 494
      ENDIF
C
      N1=J1
      N2=J2
      CALL BDYPLN
     I         (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,2,
     I          XW(1,1),XW(1,2),XW(1,3),
     O          N1,N2,N3,N4)
      CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
      SDTT1=SDT
c     CALL IBWCK6
c    I    (IBW(J1),IBW(J2),XP,YP,ZP,VPX,VPY,VPZ,NXW,NZW,M,IE,SDT,MAXEL,
c    O     N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
C
      RETURN
      END
c
c
c
      SUBROUTINE ELTRK8
     I           (XPP,YPP,ZPP,VPPX,VPPY,VPPZ,XX,YY,ZZ,VXX,VYY,VZZ,SDTT,
     I            M,IBF,IE,IDETQ,MAXEL, EPSX,ICALL,
     I            MXELW,MXNPW,NXW,NYW,NZW,IBW,NLRLW,LRLW,IEW,DL468,
     i            idi,idps,cps,cpsp,ramada,delt,mxnpfg,
     M            N1,N2,N3,N4,XW,VXW,
     O            XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1,npfgs,xsfg,csfg,mplocs)
C
C        10/10/93
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XX(8),YY(8),ZZ(8),VXX(8),VYY(8),VZZ(8),MK(4)
      DIMENSION XW(MXNPW,3),VXW(MXNPW,3),IEW(MXELW,8),IE(MAXEL,11)
      DIMENSION NLRLW(MXNPW),LRLW(24,MXNPW),IBW(MXNPW),KIN(6),NNA(3)
      DIMENSION XXX(8),YYY(8),ZZZ(8),VXXX(8),VYYY(8),VZZZ(8),KIT(4,6)
      DIMENSION KNC(3),DL468(8,MXNPW)
      dimension xsfg(mxnpfg,3),csfg(mxnpfg),mplocs(mxnpfg)
C
      DATA KIT/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8/
      DATA LS/6/,NODE/8/
C
      NELW=NXW*NYW*NZW
      NPW=(NXW+1)*(NYW+1)*(NZW+1)
      SDT=SDTT
      CALL REPLAS(XPP,YPP,ZPP,VPPX,VPPY,VPPZ,XP,YP,ZP,VPX,VPY,VPZ)
      IF(VPX.EQ.0.0D0 .AND. VPY.EQ.0.0D0 .AND. VPZ.EQ.0.0D0)THEN
        XQ=XP
        YQ=YP
        ZQ=ZP
        SDTT1=SDT
        RETURN
      ENDIF
C
C***** PREPARE INFORMATION FOR ELEMENT TRACKING
C
      DO 60 N=1,NPW
C ----- DETERMINE WORKING ARRAYS XW, YW, ZW, VXW, VYW, AND VZW
        CALL VALBDL(XX,YY,ZZ,DL468(1,N),8,NODE,
     >              XW(N,1),XW(N,2),XW(N,3))
        CALL VALBDL
     >      (VXX,VYY,VZZ,DL468(1,N),8,NODE,
     >       VXW(N,1),VXW(N,2),VXW(N,3))
        DO II=1,3
          XW(N,II)=XW(N,II)*1.0D-6
          VXW(N,II)=VXW(N,II)*1.0D-6
        ENDDO
  60  CONTINUE
C
C ***** DETERMINE THE SUBELEMENT AT WHICH POINT P LOCATES
C
      CALL MMLOC
     I    (XW,IEW,VXW,NLRLW,KIT,NELW,MXNPW,MXELW,LS,EPSX,
     I     NODE,  XP,YP,ZP,VPX,VPY,VPZ,
     O     KIN,MM,NN,NLRLN,ICODE)
C
C***** START TRACKINGS
C
      NA1=N1
      NA2=N2
      NA3=N3
      NA4=N4
      DO I=1,NODE
        IF(NA1.EQ.IE(M,I))NNA(1)=I
        IF(NA2.EQ.IE(M,I))NNA(2)=I
        IF(NA3.EQ.IE(M,I))NNA(3)=I
      ENDDO
C
C ----- POINT P COINCIDES WITH ANY NODES OF WORKING ELEMENT MM
      IF(ICODE.EQ.0)GOTO 290
C
C----- CHECK IF POINT P IS ON ANY SIDE PLANE OF WORKING ELEMENT MM
C
      KN=0
      DO I=1,LS
        IF(KIN(I).EQ.1)THEN
          KN=KN+1
          KNC(KN)=I
        ENDIF
      ENDDO
      IF(KN.EQ.1 .OR. (KN.GT.1 .AND. ICALL.EQ.1))THEN
        I=KNC(1)
        GOTO 140
      ENDIF
      IF(KN.GT.1)THEN
        DO I=1,3
          II=NNA(I)
          XXX(I)=XX(II)
          YYY(I)=YY(II)
          ZZZ(I)=ZZ(II)
        ENDDO
        DO 100 J=1,KN
          KS=KNC(J)
          DO I=1,3
            II=KIT(I,KS)
            II=IEW(MM,II)
            CANG=FCOS(XXX,YYY,ZZZ,XW(II,1),XW(II,2),XW(II,3),0)
            IF(DABS(CANG).GT.EPSX)GOTO 100
          ENDDO
          I=KS
          GOTO 140
  100   CONTINUE
        I=KNC(1)
        GOTO 140
c       WRITE(16,*)'ERROR OCCURRED AT ELTRK8 --- AFTER MMLOC'
c       STOP
      ENDIF
C
C ***** POINT P IS NOT ON ANY SIDE PLANE OF ELEMENT MM, IT IS INSIDE
C ***** ELEMENT MM
C
      ID=0
      CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O            XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  131 CONTINUE
C
      CALL TRAK2H
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
      IF(SDT1.EQ.0.)THEN
        SDTT1=0.0D0
        RETURN
      ENDIF
      IF(SDT1.NE.SDT)GOTO 250
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 131
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 1 IN ELTRK8: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
  140 CONTINUE
      N1=KIT(1,I)
      N2=KIT(2,I)
      N3=KIT(3,I)
      N4=KIT(4,I)
      N1=IEW(MM,N1)
      N2=IEW(MM,N2)
      N3=IEW(MM,N3)
      N4=IEW(MM,N4)
C
C *** YES, POINT P IS ON A SIDE PLANE N1,N2,N3,N4
C *** CHECK IF THIS POINT IS ON THE SIDE OF TETRAGON N1,N2,N3,N4
C
      CALL ONPLAN 
     I    (XW(1,1),XW(1,2),XW(1,3),N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C
      CALL CKSIDE
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NODE,N1,N2,N3,N4,MXNPW,EPSX,
     O     KON,J1,J2)
      IF(KON.EQ.1)GOTO 150
C
C *** NO, POINT P IS NOT ON ANY SIDE
C
      ID=1
      CALL CKCNEL
     I            (LRLW,NLRLW,MM,N1,N2,N3,MXNPW,24,
     O             ME1,KOUNT)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  141 CONTINUE
C
      DO 145 I=1,KOUNT
        IF(I.EQ.2)MM=ME1
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2H
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M          N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  145 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 141
      ENDIF
      IF(KOUNT.EQ.1)THEN
        N1=NA1
        N2=NA2
        N3=NA3
        N4=NA4
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 2 AT ELTRK8: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
  150 CONTINUE
C
C *** YES, POINT P IS ON THE SIDE WITH NODES J1 AND J2
C
      CALL ONLINE
     I      (XW(1,1),XW(1,2),XW(1,3),J1,J2,MXNPW,
     M       XP,YP,ZP)
C
      J12=J1
  151 CONTINUE
      IF(DABS(XP-XW(J12,1)).LE.EPSX .AND. DABS(YP-XW(J12,2)).LE.EPSX
     >    .AND. DABS(ZP-XW(J12,3)).LE.EPSX)THEN
        CALL REPLAS(XW(J12,1),XW(J12,2),XW(J12,3),VXW(J12,1),VXW(J12,2),
     >              VXW(J12,3),XP,YP,ZP,VPX,VPY,VPZ)
        NLRLN=NLRLW(J12)
        NN=J12
        GOTO 290
      ELSEIF(J12.EQ.J1)THEN
        J12=J2
        GOTO 151
      ENDIF
C
      ID=2
      NLRL1=NLRLW(J1)
      NLRL2=NLRLW(J2)
      KOUNT=0
      DO 153 I1=1,NLRL1
        M1=LRLW(I1,J1)
        DO 152 I2=1,NLRL2
          M2=LRLW(I2,J2)
          IF(M1.EQ.M2)THEN
            KOUNT=KOUNT+1
            MK(KOUNT)=M1
          ENDIF
  152   CONTINUE
  153 CONTINUE
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  154 CONTINUE
C
      DO 155 I=1,KOUNT
        MM=MK(I)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2H
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M          N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  155 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 154
      ENDIF
C
      N1=NA1
      N2=NA2
      N3=NA3
      N4=NA4
      CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
      SDTT1=SDT
c     CALL IBWCK8
c    I    (IBW(J1),IBW(J2),XP,YP,ZP,VPX,VPY,VPZ,NXW,NYW,NZW,M,IE,SDT,
c    O     MAXEL,  N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
      RETURN
C
C $$$$$ START CONSEQUENT TRACKINGS
C
  250 CONTINUE
c
      if(idps.eq.1)then
        npfgs=npfgs+1
        xsfg(npfgs,1)=0.5d0*(xp+xq)
        xsfg(npfgs,2)=0.5d0*(yp+yq)
        xsfg(npfgs,3)=0.5d0*(zp+zq)
        mplocs(npfgs)=m
        dtreal=delt-0.5d0*(sdt+sdt1)
        xsi=(delt-dtreal)/delt
        cc=cpsp*xsi+cps*(1.0d0-xsi)
        csfg(npfgs)=cc*dexp(-ramada*dtreal)
      endif
      SDT=SDT1
      CALL REPLAS(XQ,YQ,ZQ,VQX,VQY,VQZ,XP,YP,ZP,VPX,VPY,VPZ)
C
C ***** CHECK IF POINT P COINSIDES WITH ANY GLOBAL NODES
C
      CALL CKCOIN
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NLRLW,NODE,
     O     N1,N2,N3,N4,MXNPW,EPSX,  NN,NLRLN)
      IF(NN.EQ.0)GOTO 400
C
C ***** YES, P COINSIDES WITH NODE NN, THEN PREPARE DATA FOR TRACK1
C
      CALL REPLAS(XW(NN,1),XW(NN,2),XW(NN,3),VXW(NN,1),VXW(NN,2),
     >            VXW(NN,3),  XP,YP,ZP,VPX,VPY,VPZ)
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
  290 CONTINUE
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  300 CONTINUE
C
      DO 350 I=1,NLRLN
        MM=LRLW(I,NN)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
        DO 305 J=1,NODE
          IF(NN.EQ.IEW(MM,J))GOTO 310
  305   CONTINUE
C
  310   CONTINUE
C
        CALL TRAK1H
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,J,MXELW, EPSX, idi,
     O          XQ,YQ,ZQ,VQX,VQY,VQZ,N1,N2,N3,N4,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  350 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 300
      ENDIF
C
      IF(IBW(NN).NE.0)THEN
        IF(IBW(NN).EQ.1 .OR. IBW(NN).EQ.11 .OR. IBW(NN).EQ.14 .OR.
     >     IBW(NN).EQ.15 .OR. IBW(NN).EQ.18 .OR. IBW(NN).EQ.21 .OR.
     >     IBW(NN).EQ.24.OR. IBW(NN).EQ.25 .OR. IBW(NN).EQ.29)THEN
          I=1
        ELSEIF(IBW(NN).EQ.2 .OR. IBW(NN).EQ.12 .OR. IBW(NN).EQ.16 .OR.
     >     IBW(NN).EQ.22 .OR. IBW(NN).EQ.26 .OR. IBW(NN).EQ.30)THEN
          I=2
        ELSEIF(IBW(NN).EQ.3 .OR. IBW(NN).EQ.13 .OR. IBW(NN).EQ.17 .OR.
     >     IBW(NN).EQ.23 .OR. IBW(NN).EQ.27 .OR. IBW(NN).EQ.31)THEN
          I=3
        ELSEIF(IBW(NN).EQ.4 .OR. IBW(NN).EQ.28 .OR. IBW(NN).EQ.32)THEN
          I=4
        ELSEIF(IBW(NN).EQ.5)THEN
          I=5
        ELSE
          I=6
        ENDIF
        N1=KIT(1,I)
        N2=KIT(2,I)
        N3=KIT(3,I)
        N4=KIT(4,I)
        N1=IE(M,N1)
        N2=IE(M,N2)
        N3=IE(M,N3)
        N4=IE(M,N4)
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
      WRITE(16,*)'ERROR MESSAGE 2 AT ELTRK8: POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C ***** NO, P DOESN'T COINSIDE WITH ANY NODES, THEN PREPARE DATA FOR
C       TRACK2
C
  400 CONTINUE
C
C ***** CHECK IF P IS ON ANY SIDES OF TETRAGON N1,N2,N3,N4
C
      CALL ONPLAN
     I    (XW(1,1),XW(1,2),XW(1,3),N1,N2,N3,MXNPW,
     M     XP,YP,ZP)
C
      CALL CKSIDE
     I    (XW(1,1),XW(1,2),XW(1,3),XP,YP,ZP,NODE,
     O     N1,N2,N3,N4,MXNPW,EPSX,  KON,J1,J2)
      IF(KON.EQ.1)GOTO 491
C
C *** NO, IT IS NOT ON ANY SIDE.
C
C === CHECK IF N1,N2,N3,N4 ARE ON A BOUNDARY PLANE
C
      CALL CKCNEL
     I            (LRLW,NLRLW,MM,N1,N2,N3,MXNPW,24,
     O             ME1,KOUNT)
      IF(KOUNT.EQ.1)GOTO 490
C
C --- N1,N2,N3,N4 ARE NOT ON A BOUNDARY PLANE. THEN CHECK WHICH PLANE
C --- POINT P IS ONTO
C
      ME=MM
      MM=ME1
C
  431 CONTINUE
C
      CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O            XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
  490 CONTINUE
C
        ID=1
C
C === DETERMINE HOW TO TRACK THIS POINT
C
  484 CONTINUE
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  485 CONTINUE
C
      CALL TRAK2H
     I            (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I             MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M             N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
      IF(SDT1.EQ.0.)THEN
        SDTT1=0.0D0
        RETURN
      ENDIF
      IF(SDT1.NE.SDT)GOTO 250
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 485
      ENDIF
      IF(KOUNT.EQ.1)THEN
        CALL BDYPLN
     I       (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,1,
     I        XW(1,1),XW(1,2),XW(1,3),
     O        N1,N2,N3,N4)
        CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
        SDTT1=SDT
        RETURN
      ENDIF
C
      IF(MM.EQ.ME1)THEN
        MM=ME
        GOTO 431
      ENDIF
C
      WRITE(16,*)'ERROR MESSAGE 3 AT ELTRK8 --- POINT',N,' CAN NOT BE',
     >          ' TRACKED'
      STOP
C
C ***** YES
C
  491 CONTINUE
C
      CALL ONLINE
     I      (XW(1,1),XW(1,2),XW(1,3),J1,J2,MXNPW,
     M       XP,YP,ZP)
C
      J12=J1
  495 CONTINUE
      IF(DABS(XP-XW(J12,1)).LE.EPSX .AND. DABS(YP-XW(J12,2)).LE.EPSX
     >    .AND. DABS(ZP-XW(J12,3)).LE.EPSX)THEN
        CALL REPLAS(XW(J12,1),XW(J12,2),XW(J12,3),VXW(J12,1),VXW(J12,2),
     >              VXW(J12,3),XP,YP,ZP,VPX,VPY,VPZ)
        NLRLN=NLRLW(J12)
        NN=J12
        GOTO 290
      ELSEIF(J12.EQ.J1)THEN
        J12=J2
        GOTO 495
      ENDIF
C
      ID=2
      NLRL1=NLRLW(J1)
      NLRL2=NLRLW(J2)
      KOUNT=0
      DO 493 I1=1,NLRL1
        M1=LRLW(I1,J1)
        DO 492 I2=1,NLRL2
          M2=LRLW(I2,J2)
          IF(M1.EQ.M2)THEN
            KOUNT=KOUNT+1
            MK(KOUNT)=M1
          ENDIF
  492   CONTINUE
  493 CONTINUE
C
C *** DETERMINE HOW TO TRACK THIS POINT
C
      IF(IDETQ.EQ.2)THEN
        IJUDGE=2
      ELSE
        IJUDGE=1
      ENDIF
C
  494 CONTINUE
C
      DO 499 I=1,KOUNT
        MM=MK(I)
        CALL WRKARY(IEW,XW,VXW,MXELW,8,MXNPW,MM,NODE,8,
     O              XXX,YYY,ZZZ,VXXX,VYYY,VZZZ)
C
        CALL TRAK2H
     I         (XP,YP,ZP,VPX,VPY,VPZ,XXX,YYY,ZZZ,VXXX,VYYY,VZZZ,SDT,
     I          MM,IBF,IEW,IJUDGE,ID,MXELW, EPSX,J1,J2, idi,
     M          N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDT1)
C
        IF(SDT1.EQ.0.)THEN
          SDTT1=0.0D0
          RETURN
        ENDIF
        IF(SDT1.NE.SDT)GOTO 250
  499 CONTINUE
C
      IF(IJUDGE.EQ.1)THEN
        IJUDGE=2
        GOTO 494
      ENDIF
C
      N1=J1
      N2=J2
      CALL BDYPLN
     I         (IBW,IE,KIT,XX,YY,ZZ,XP,YP,ZP,M,MXNPW,MAXEL,LS,EPSX,2,
     I          XW(1,1),XW(1,2),XW(1,3),
     O          N1,N2,N3,N4)
      CALL REPLAS(XP,YP,ZP,VPX,VPY,VPZ,XQ,YQ,ZQ,VQX,VQY,VQZ)
      SDTT1=SDT
c     CALL IBWCK8
c    I    (IBW(J1),IBW(J2),XP,YP,ZP,VPX,VPY,VPZ,NXW,NYW,NZW,M,IE,SDT,
c    O     MAXEL,  N1,N2,N3,N4,XQ,YQ,ZQ,VQX,VQY,VQZ,SDTT1)
C
      RETURN
      END
c
c
c
      SUBROUTINE MMLOC
     I    (XW,IEW,VXW,NLRLW,KIT,NELW,MXNPW,MXELW,LS,EPSX,
     I     NODE,  XP,YP,ZP,VPX,VPY,VPZ,
     O     KIN,MM,NN,NLRLN,ICODE)
C  11/17/93
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION XW(MXNPW,3),IEW(MXELW,8),NLRLW(MXNPW)
      DIMENSION VXW(MXNPW,3),KIN(LS),KIT(4,LS)
      DIMENSION XXX(8),YYY(8),ZZZ(8),XWW(8),YWW(8),ZWW(8)
C
      ICODE=-1
      XPP=XP
      YPP=YP
      ZPP=ZP
      DO 105 MP=1,NELW
        XP=XPP
        YP=YPP
        ZP=ZPP
        DO II=1,NODE
          IEM=IEW(MP,II)
          XXX(II)=XW(IEM,1)
          YYY(II)=XW(IEM,2)
          ZZZ(II)=XW(IEM,3)
          IF(DABS(XP-XXX(II)).LE.EPSX .AND.
     >       DABS(YP-YYY(II)).LE.EPSX .AND.
     >       DABS(ZP-ZZZ(II)).LE.EPSX)THEN
            CALL REPLAS(XXX(II),YYY(II),ZZZ(II),VXW(IEM,1),VXW(IEM,2),
     >                  VXW(IEM,3),XP,YP,ZP,VPX,VPY,VPZ)
            MM=MP
            NLRLN=NLRLW(IEM)
            NN=IEM
            ICODE=0
            RETURN
          ENDIF
        ENDDO
C
        DO 100 III=1,LS
          DO II1=1,3
            JJ=KIT(II1,III)
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
            KIN(III)=1
            CALL ONPLAN
     I         (XWW,YWW,ZWW,1,2,3,8,
     M          XP,YP,ZP)
          ENDIF
          IF(CANG.GT.0.0D0 .AND. DABS(CANG).GT.EPSX)KIN(III)=2
  100   CONTINUE
        MM=MP
        RETURN
  105 CONTINUE
C
      RETURN
      END
