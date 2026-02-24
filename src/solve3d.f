C $$$$$$$$$$$$ SOLVE3D.F
C       9/ 7/94     2:00 pm
C
      SUBROUTINE BLKITR(R,RW, CMTRXL,RLDL, CMTRXG,RLDG, GNLR,LNOJCN,
     1 NNPLR,LMAXDF, OMI,TOLB,NITER,IBUG,KPR,IDCS)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO SOLVE THE MATRIX EQUATION WITH BLOCK ITERATION.  FIRST,
C ------- THE BLOCK MATRIX EQUATION IS ASSEMBLED OUT OF THE GLOBAL
C ------- MATRIX EQUATION.  THEN THE INTRA-BOUNDARY CONDITIONS ARE
C ------- IMPLEMENTED.  FINALLY, THE BLOCK MATRIX EQUATION IS SOLVED
C ------- WITH DIRECT BAND MATRIX SOLVER.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: RW(NNP), CMTRXG(NNP,MXJBD), RLDG(NNP),
C -------        GNLR(NTNNP,NREGN), LNOJCN(MXJBD,NNPLR(K),NREGN),
C -------        NNPLR(NREGN), LMAXDF(NREGN), TOLB, NITER, IBUG, KPR,
C -------        OMI.
C
C ------- OUTPUT: R(NNP), RW(NNP).
C
C-------- INPUTING WORKING ARRAYS: CMTRXL(LMXNP,LMXBW), RLDL(LMXNP).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 GNLR
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /LGEOM/ LTMXNP,LMXNP,LMXBW,MXREGN,NREGN
C
      DIMENSION R(MXADNP),RW(MXADNP)
      DIMENSION CMTRXL(LMXNP,LMXBW),RLDL(LMXNP)
      DIMENSION CMTRXG(MXADNP,MXJBD),RLDG(MXADNP)
      DIMENSION GNLR(LTMXNP,MXREGN),LNOJCN(MXJBD,LMXNP,MXREGN)
      DIMENSION NNPLR(MXREGN),LMAXDF(MXREGN)
C
C ------ PUT ZERO-TH ITERATE INTO ARRAY R
C
      DO 110 NP=1,NNP
  110 R(NP)=RW(NP)
C
C ****** START INTERATION LOOP
C
      IF(IBUG.NE.0 .AND. KPR.NE.0) WRITE(16,1000)
C
      DO 690 IT=1,NITER
C
C ------ FOR EACH ITERATION, SOLVE FOR NREGN REGIONS
C
      DO 590 K=1,NREGN
C
      LHALFB=LMAXDF(K)
      LIHBP=LHALFB+1
      IBAND=2*LHALFB+1
      LNNP=NNPLR(K)
C
C ------- PUT GLOBAL LOAD VECTOR INTO CORRESPONDING LOCAL LOAD VECTOR
C ------- AND INITIATE THE LOCAL MATRIX CMTRXL(LNNP,IBAND)
C
      DO 210 LI=1,LNNP
      NP=GNLR(LI,K)
      RLDL(LI)=RLDG(NP)
      DO 210 J=1,IBAND
      CMTRXL(LI,J)=0.0
  210 CONTINUE
C
C ------ ASSEMBLE LOCAL COEFFICIENT MATRIX CMTRXL(LMXNP,LMXBW) FROM THE
C ------ GLOBAL COEFFICIENT MATRIX CMTRXG(MAXNP,MXJBD) AND INCORPORATE
C ------ INTERFACIAL DIRICHLET BOUNDARY CONDITIONS INTO RLDL(LMXNP).
C
      DO 490 LI=1,LNNP
      NI=GNLR(LI,K)
C
      DO 390 J=1,MXJBD
      LJ=LNOJCN(J,LI,K)
C
      IF(LJ.LE.0) GO TO 390
C
      NJ=GNLR(LJ,K)
      IF(LJ.GT.LNNP) GO TO 380
      LJB=LJ-LI+LIHBP
      CMTRXL(LI,LJB)=CMTRXG(NI,J)
      GO TO 390
C
  380 RLDL(LI)=RLDL(LI)-CMTRXG(NI,J)*R(NJ)
C
  390 CONTINUE
C
  490 CONTINUE
C
C ------ SOLVE THE BLOCK EQUATIONS
C
      CALL SOLVE(1,CMTRXL,RLDL,LNNP,LHALFB,LMXNP,LMXBW)
      CALL SOLVE(2,CMTRXL,RLDL,LNNP,LHALFB,LMXNP,LMXBW)
C
C ------ PUT THE NEWLY OBTAINED BLOCK SOLTUION INTO THE GLOBAL SOLUTION.
C
      DO 560 LI=1,LNNP
      NP=GNLR(LI,K)
      R(NP)=RLDL(LI)*OMI + (1.0D0-OMI)*R(NP)
  560 CONTINUE
C
  590 CONTINUE
C
C ------ CHECK IF THE CONVERGENT SOLUTION IS OBTAINED?
C
      IF(IDCS.EQ.1)THEN
        DIFMAX=0.0
        NOCCUR=1
        DO 660 NP=1,NNP
          DIF=R(NP)-RW(NP)
          DIF=DABS(DIF)
          IF(DIF.LE.DIFMAX) GO TO 660
          DIFMAX=DIF
          NOCCUR=NP
  660   CONTINUE
      ELSEIF(IDCS.EQ.2)THEN
        DIFMAX=-1.0D38
        NOCCUR=1
        DO 670 NP=1,NNP
          IF(RW(NP).EQ.0.0) GO TO 670
          DIF=(R(NP)-RW(NP))/RW(NP)
          DIF=DABS(DIF)
          IF(DIF.LE.DIFMAX) GO TO 670
          DIFMAX=DIF
          NOCCUR=NP
  670   CONTINUE
C
c       ERROR=0.0D0
c       DO 670 NP=1,NNP
c         ERROR=ERROR+(R(NP)-RW(NP))**2
c 670   CONTINUE
c       DIFMAX=DSQRT(ERROR/DBLE(NNP))
      ENDIF
C
C ------- UP DATA THE ITERATE
C
      DO 680 NP=1,NNP
      RW(NP)=R(NP)
  680 CONTINUE
C
C------- PRINT ITERATION INFORMATION
C
      IF(IBUG.NE.0 .AND. KPR.NE.0) WRITE(16,1100) IT,DIFMAX,TOLB,NOCCUR
C
      IF(IT.EQ.1) GO TO 690
C
      IF(DIFMAX.LT.TOLB) GO TO 990
C
  690 CONTINUE
C
      WRITE(16,2000) IT,NITER,DIFMAX,TOLB,NOCCUR
C
C
  990 CONTINUE
C
 1000 FORMAT('0'//15X,'   IT   DIFMAX       TOLB        NOCCUR'/15X,
     1 '   --   ------       ----        ------')
 1100 FORMAT(' ',15X,I5,2D12.4,I10)
 2000 FORMAT('0'/10X,' *** WARNING: NO CONVERGENCE IN BLKI AFTER  ',I4,
     1 '   ITERATIONS'/10X,'  NITER =',I4,'  DIFMAX =',D11.4,
     2 '   TOLB  =',D11.4,'  NOCCUR =',I4)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE SOLVE(KKK,C,R,NNP,IHALFB,MAXNP,MAXBW)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO SOVE A MATRIX EQUATION WITH BAND MATRIX SOLVER.
C
C********1*********2*********3*********4*********5*********5*********7**
C
C ------- INPUT: C(NEQ,NBAND), R(NEQ), NNP, IHALFB, KKK,
C -------         WHERE NNP=NEQ AND IHALFB=(NBAND-1)/2
C
C ------- OUTPUT: C(NEQ).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION C(MAXNP,MAXBW),R(MAXNP)
C
      IHBP=IHALFB+1
C
C  IF KKK = 1, THEN TRIANGULARIZE THE BAND MATRIX C(NP,IB), BUT
C  IF KKK = 2, THEN SIMPLY SOLVE WITH THE RIGHT-HAND SIDE R(NP)
C
      IF (KKK.EQ.2) GO TO 50
C
C  TRIANGULARIZE MATRIX C(NP,IB)
C
      NU=NNP-IHALFB
      DO 20 NI=1,NU
         PIVOTI=1.0D0/C(NI,IHBP)
         NJ=NI+1
         IB=IHBP
         NK=NI+IHALFB
         DO 10 NL=NJ,NK
            IB=IB-1
            A=-C(NL,IB)*PIVOTI
            C(NL,IB)=A
            JB=IB+1
            KB=IB+IHALFB
            LB=IHBP-IB
            DO 10 MB=JB,KB
               NB=LB+MB
   10          C(NL,MB)=C(NL,MB)+A*C(NI,NB)
   20    CONTINUE
      NR=NU+1
      NU=NNP-1
      NK=NNP
      IF(NR.GT.NU) RETURN
      DO 40 NI=NR,NU
         PIVOTI=1.0D0/C(NI,IHBP)
         NJ=NI+1
         IB=IHBP
         DO 30 NL=NJ,NK
            IB=IB-1
            A=-C(NL,IB)*PIVOTI
            C(NL,IB)=A
            JB=IB+1
            KB=IB+IHALFB
            LB=IHBP-IB
            DO 30 MB=JB,KB
               NB=LB+MB
   30          C(NL,MB)=C(NL,MB)+A*C(NI,NB)
   40    CONTINUE
      RETURN
C
C  MODIFY LOAD VECTOR R(NP)
C
   50 NU=NNP+1
      IBAND=2*IHALFB+1
      DO 70 NI=2,IHBP
         IB=IHBP-NI+1
         NJ=1
         SUM=0.0
         DO 60 JB=IB,IHALFB
            SUM=SUM+C(NI,JB)*R(NJ)
   60       NJ=NJ+1
   70    R(NI)=R(NI)+SUM
      IB=1
      NL=IHBP+1
      DO 90 NI=NL,NNP
         NJ=NI-IHBP+1
         SUM=0.0
         DO 80 JB=IB,IHALFB
            SUM=SUM+C(NI,JB)*R(NJ)
   80       NJ=NJ+1
   90    R(NI)=R(NI)+SUM
C
C  BACK SOLVE
C
      R(NNP)=R(NNP)/C(NNP,IHBP)
      DO 110 IB=2,IHBP
         NI=NU-IB
         NJ=NI
         MB=IHALFB+IB
         SUM=0.0
         DO 100 JB=NL,MB
            NJ=NJ+1
  100       SUM=SUM+C(NI,JB)*R(NJ)
  110    R(NI)=(R(NI)-SUM)/C(NI,IHBP)
      MB=IBAND
      DO 130 IB=NL,NNP
         NI=NU-IB
         NJ=NI
         SUM=0.0
         DO 120 JB=NL,MB
            NJ=NJ+1
  120       SUM=SUM+C(NI,JB)*R(NJ)
  130    R(NI)=(R(NI)-SUM)/C(NI,IHBP)
      RETURN
      END
C
C
C
C
      SUBROUTINE PISS
     I        (MXTNOD,MXJBD,MXADNP, NNP,NLRN,NITER, OMI,EPS,KTIM,
     I         IBUG, C,RI,RL,LRN,IDCS,
     O         R)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO SOLVE A MATRIX EQUATION WITH POINTWISE ITERATION SOLUTION
C ------- STRATEGIES.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: RI(NNP), RL(NNP), C(NNP,MXJBD), LRN(MXJBD,NNP),
C -------        OMI, EPS, NITER, IBUG, AND KTIM.
C
C ------- OUTPUT: R(NNP) - FINAL SOLUTION; AND RI(NNP)-ITERATE.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION R(MXADNP),RI(MXADNP),NLRN(MXTNOD)
      DIMENSION C(MXADNP,MXJBD),RL(MXADNP),LRN(MXJBD,MXADNP)
C
C ------- PRINT ITERATION INFORMATION IF DESIRED.
C
      IF(IBUG.NE.0 .AND. KTIM.NE.0) WRITE(16,1000)
C
      DO 290 IT=1,NITER
C
C ------- FOR EACH ITERATION, PUT THE LOAD VECTOR INTO R(MAXNP).
C
      DO 210 NP=1,NNP
  210 R(NP)=RL(NP)
C
C ------- THE MATRIX C = L + I + U, WHERE L IS THE LOWER TRIANGULAR
C ------- MATRIX AND U IS THE UPPER TRIANGULAR MATRIX AND I IS THE
C ------- DIAGONAL MATRIX.
C ------- NOW ADD U.RI TO THE RIGHT HAND SIDE, WHERE RI IS THE PREVIOUS
C ------- ITERATE.
C
      DO 230 NP=1,NNP
      DO 220 I=1,NLRN(NP)
      LNODE=LRN(I,NP)
      IF(LNODE.EQ.0) GO TO 220
      IF(LNODE.LE.NP) GO TO 220
      R(NP)=R(NP)-C(NP,I)*RI(LNODE)
  220 CONTINUE
  230 CONTINUE
C
C ------- START TO COMPUTE NEW ITERATE WITH POINTWISE ITERATIONS.
C
      DO 250 NP=1,NNP
      DO 240 I=1,NLRN(NP)
      LNODE=LRN(I,NP)
      IF(LNODE.EQ.0) GO TO 240
      IF(LNODE.GE.NP) GO TO 240
      R(NP)=R(NP)-C(NP,I)*R(LNODE)
  240 CONTINUE
      DO 245 I=1,NLRN(NP)
      LNODE=LRN(I,NP)
      IF(LNODE.NE.NP) GO TO 245
      R(NP)=OMI*(R(NP)/C(NP,I)) + (1.0D0-OMI)*RI(NP)
  245 CONTINUE
  250 CONTINUE
C
C ------- CHECK IF A CONVERGENT SOLUTION IS OBTAINED?
C
        NNCVN=0
        RELERR=-1.0D0
        ABSERR=-1.0D0
        DO 260 NP=1,NNP
          ABSDIF=DABS(R(NP)-RI(NP))
          ABSERR=DMAX1(ABSERR,ABSDIF)
          IF(RI(NP).EQ.0.0D0 .AND. IDCS.EQ.2)GOTO 260
          IF(RI(NP).NE.0.0D0) RELERR=DMAX1(RELERR,DABS(ABSDIF/RI(NP)))
          IF(IDCS.EQ.1 .AND. ABSERR.LE.EPS)GOTO 260
          IF(IDCS.EQ.2 .AND. RELERR.LE.EPS)GOTO 260
          NNCVN=NNCVN+1
  260   CONTINUE
C
c       ERROR=0.0D0
c       DO 270 NP=1,NNP
c         ERROR=ERROR+(R(NP)-RI(NP))**2
c 270   CONTINUE
c       ABSERR=DSQRT(ERROR/DBLE(NNP))
C
C ------- UPDATE THE ITERATE.
C
      DO 280 NP=1,NNP
  280 RI(NP)=R(NP)
C
C -------  PRINT ITERATION INFORMATION IF DESIRED.
C
      IF(IBUG.NE.0 .AND. KTIM.NE.0) WRITE(16,1100) IT,ABSERR,RELERR,
     1 NNCVN
C
      IF(IT.EQ.1) GO TO 290
      IF(ABSERR .LE. EPS .AND. IDCS.EQ.1) GO TO 990
      IF(RELERR .LE. EPS .AND. IDCS.EQ.2) GO TO 990
C
  290 CONTINUE
C
      WRITE(16,2000) IT,NITER,ABSERR,RELERR,NNCVN
C
  990 CONTINUE
C
 1000 FORMAT('0'/36X,'   IT   ABS ERR     REL ERR   NO. OF NODES'/36X,
     1 '   --   --- ---     --- ---   --- -- -----')
 1100 FORMAT(36X,I5,1PD12.4,1PD12.4,5X,I5)
 2000 FORMAT('0'/18X,' ::: WARNING: NO CONVERGENCE IN PISS AFTER ',
     1 I4,' ITERATIONS'/18X,' NPITER =',I4,'   ABSERR =',1PD12.4/18X,
     2 '   RELERR =',1PD12.4,'    NNCVN =',I4)
C
      RETURN
      END
C
C
C
C
      SUBROUTINE CONECT(nc,nd,nt)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C---- This subroutine generates nodal connections from element connectiv
C
C-------------------------------------------------
C Input Variables Needed for this Subroutine are:
C-------------------------------------------------
C   NN - Total number of nodes
C   NEL - Total number of elements
C   NODE(I,J) - Node number of the Jth node of element I
C
C----------------------------------------
C The Subroutine Generates the Following:
C----------------------------------------
C   NT(I) - Total number of nodes connected to node I
C   NC(J,I) - Node number of the Jth node connected to node I, in ascend
C             order of the node number
C   ND(I) - Position of the "diagonal" element in row I
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
C
      DIMENSION nc(mxjbd,maxnp),nd(maxnp),nt(maxnp)
C
C ******** initialization.   JPG ************
C
      nn=nnp
      mb=mxjbd
C
      DO 10 I=1,NN
        NT(I)=0
        ND(I)=0
  10  CONTINUE
C
C ******** calculate the number of elements each row
C
      DO np=1,nn
        DO icn=1,mb
          IF (nc(icn,np).EQ.0) GO TO 20
        ENDDO
   20   nt(np)=icn-1
      ENDDO
C
C---- Arrange the connected nodes in ascending order
C
      DO 50 I=1,NN
      DO 60 J=1,NT(I)-1
      DO 60 K=J+1,NT(I)
      IF(NC(k,i).LT.NC(j,i)) THEN
        MIN=NC(k,i)
        NC(k,i)=NC(j,i)
        NC(j,i)=MIN
      ENDIF
 60   CONTINUE
 50   CONTINUE
C
C---- Locate the diagonal element for each row
C
      DO 70 I=1,NN
      DO 80 J=1,NT(I)
      IF(NC(j,i).EQ.I) THEN
        ND(I)=J
        GO TO 70
      ENDIF
 80   CONTINUE
 70   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE PPCG
     I  (C,XX,LRN,NLRN, IEIGEN,GG,EPS,SQEPS,IBUG,KTIM, MAXNP,MXTNOD,
     M   MXJBD,NNP,ID,  SK,RK,ZK,PK,
     O   XK)
C                                                                       ppcg 015
C ------- Preconditioned conjugated gradient matrix solver              ppcg 020
C ------- Algorithm according to Preconditioned Conjugate-Gradient 2    ppcg 025
C ------- (PCG2), A Computer Program For Solving Groundwater Flow       ppcg 030
C ------- Equations, by Mary C. Hill, USGS Water-Resources              ppcg 035
C ------- Investigation Report 90-4048.                                 ppcg 040
C                                                                       ppcg 045
C ------- Parameters:                                                   ppcg 050
C ------- INPUT --                                                      ppcg 055
C         sk(neqn) = working array of dimension neqn                    ppcg 060
C         rk(neqn) = working array of dimension neqn                    ppcg 065
C         zk(neqn) = working array of dimension neqn                    ppcg 070
C         lrn(jband,nnp) = surounding nodes of node nnp                 ppcg 075
C         c(neqn,nbwd) = coefficient matrix of equation cx=b, note that ppcg 080
C                        c must be in compacted form with diagnal elemenppcg 085
C                        blocks as the first column                     ppcg 090
C         pk(neqn) = working array of dimension neqn                    ppcg 095
C         xx(neqn) = variable vector of matrix equation cx=b; RHS on inpppcg 100
C         nnp = # of nodes                                              ppcg 105
C         jband = band width of matrix lrn                              ppcg 110
C         nbwd = band width of coefficient matrix                       ppcg 115
C                                                                       ppcg 120
C ------- OUTPUT --                                                     ppcg 125
C         xk(neqn) = solution for the matrix equation cx=b              ppcg 130
C                                                                       ppcg 135
      IMPLICIT REAL*8(a-h,o-z)                                          ppcg 140
C                                                                       ppcg 145
      DIMENSION sk(maxnp),rk(maxnp),zk(maxnp),nlrn(mxtnod)              ppcg 165
      DIMENSION lrn(mxjbd,maxnp),c(maxnp,mxjbd),xx(maxnp)               ppcg 170
      DIMENSION pk(maxnp),xk(maxnp)                                     ppcg 175
C                                                                       ppcg 180
c     jband=mxjbd                                                       ppcg 185
c     nbwd=mxjbd                                                        ppcg 190
      npiter=nnp
c                                                                       ppcg 195
      IF (ibug.NE.0 .and. ktim.NE.0) WRITE (16,1000)                    ppcg 200
C                                                                       ppcg 210
C ------- put initial guess in xk(nnp)                                  ppcg 215
C                                                                       ppcg 220
      DO 90 np=1,nnp                                                    ppcg 225
        xk(np)=zk(np)                                                   ppcg 230
   90 CONTINUE                                                          ppcg 235
C                                                                       ppcg 240
C ------- initialization of computations                                ppcg 245
C ------- 1. calculate r0=b-cx                                          ppcg 250
C                                                                       ppcg 255
      DO 110 np=1,nnp                                                   ppcg 260
          rowsum = 0.0d0                                                ppcg 265
          DO 100 ibwd=1,nlrn(np)                                        ppcg 270
            IF (lrn(ibwd,np).EQ.0) GO TO 100                            ppcg 275
            rowsum = rowsum + c(np,ibwd)*xk(lrn(ibwd,np))               ppcg 280
  100     CONTINUE                                                      ppcg 285
          rk(np) = xx(np) - rowsum                                      ppcg 290
  110 CONTINUE                                                          ppcg 295
C                                                                       ppcg 300
C ------- calculate polynomial preconditioner                           ppcg 305
C                                                                       ppcg 310
       CALL POLYP
     I   (C,RK,LRN,NLRN, IEIGEN,GG, MXTNOD,MAXNP,MXJBD,NNP,
     M    ZK,
     O    SK)
C                                                                       ppcg 320
C ------- product of c.pk, for k=0, pk=sk                               ppcg 325
C                                                                       ppcg 330
      DO 280 np=1,nnp                                                   ppcg 335
        rowsum = 0.0d0                                                  ppcg 340
        DO 270 ibwd=1,nlrn(np)                                          ppcg 345
          IF (lrn(ibwd,np).EQ.0) GO TO 270                              ppcg 350
          rowsum = rowsum + c(np,ibwd)*sk(lrn(ibwd,np))                 ppcg 355
  270   CONTINUE                                                        ppcg 360
        zk(np) = rowsum                                                 ppcg 365
  280 CONTINUE                                                          ppcg 370
C                                                                       ppcg 375
C ------- finally, alphk                                                ppcg 380
C                                                                       ppcg 385
      banum = dotprd(sk,rk,nnp)
      alphk = banum / dotprd(zk,sk,nnp)
C                                                                       ppcg 400
C ------- calculate x1
C
      abserr = -1.0d0
      relerr= -1.0d0
      DO 310 np=1,nnp
        delxnp = alphk*sk(np)
        IF(ID.eq.1) THEN
          IF (DABS(delxnp).GT.abserr) THEN
            abserr = DABS(delxnp)
            ncveq = np
            if(xk(np).ne.0.0) relerr=dabs(delxnp/xk(np))
          ENDIF
        ELSEIF(ID.EQ.2) THEN
          if(xk(np).ne.0.0d0) then
            delrel=delxnp/xk(np)
            IF (DABS(delrel).GT.relerr) THEN
              relerr = DABS(delrel)
              ncveq = np
              abserr=dabs(delxnp)
            ENDIF
          endif
        ENDIF
C
        xk(np) = xk(np) + delxnp
        pk(np) = sk(np)
  310 CONTINUE
C
C ------- update rk
C
      resdmx=-1.0d38
      DO 320 np=1,nnp
        rknp = rk(np) - alphk*zk(np)
        resdmx=dmax1(resdmx,dabs(rknp))
        rk(np)=rknp
  320 CONTINUE
C                                                                       ppcg 500
      IF(ID.EQ.1) THEN
        IF(abserr.LT.eps .and. resdmx.LT.sqeps) RETURN
      ELSEIF(ID.EQ.2) THEN
        IF(relerr.LT.eps .and. resdmx.LT.sqeps) RETURN
      ENDIF
C                                                                       ppcg 470
C ------- print iteration residuals                                     ppcg 475
C                                                                       ppcg 480
      IF (ibug.NE.0 .and. ktim.NE.0)                                    ppcg 485
     .  WRITE (16,1001) 0,ncveq,abserr,relerr,xk(ncveq),rk(ncveq),eps
C                                                                       ppcg 540
        betad = banum
C                                                                       ppcg 550
C ------- start iteration loop                                          ppcg 555
C                                                                       ppcg 560
      DO 900 it=1,npiter                                                ppcg 565
C                                                                       ppcg 570
C ------- calculate polynomial preconditioner                           ppcg 575
C                                                                       ppcg 580
       CALL POLYP
     I   (C,RK,LRN,NLRN, IEIGEN,GG, MXTNOD,MAXNP,MXJBD,NNP,
     M    ZK,
     O    SK)
C                                                                       ppcg 590
C ------- calculate betak                                               ppcg 595
C                                                                       ppcg 600
        banum = dotprd(sk,rk,nnp)
        betak = banum / betad
C                                                                       ppcg 615
C ------- calculate pk                                                  ppcg 620
C                                                                       ppcg 625
        DO 400 np=1,nnp                                                 ppcg 630
          pk(np) = sk(np) + betak*pk(np)                                ppcg 635
  400   CONTINUE                                                        ppcg 640
C                                                                       ppcg 645
C ------- product of c.pk                                               ppcg 650
C                                                                       ppcg 655
      DO 560 np=1,nnp                                                   ppcg 660
        rowsum = 0.0d0                                                  ppcg 665
        DO 550 ibwd=1,nlrn(np)                                          ppcg 670
          IF (lrn(ibwd,np).EQ.0) GO TO 550                              ppcg 675
          rowsum = rowsum + c(np,ibwd)*pk(lrn(ibwd,np))                 ppcg 680
  550   CONTINUE                                                        ppcg 685
        zk(np) = rowsum                                                 ppcg 690
  560 CONTINUE                                                          ppcg 695
C                                                                       ppcg 700
C ------- alphk                                                         ppcg 705
C                                                                       ppcg 710
        alphk = banum / dotprd(zk,pk,nnp)
C                                                                       ppcg 720
C ------- calculate xk                                                  ppcg 725
C                                                                       ppcg 730
        abserr = -1.0d0
        relerr=-1.0d0
        DO 610 np=1,nnp
          delxnp = alphk*pk(np)
          IF(ID.eq.1) THEN
            IF (DABS(delxnp).GT.abserr) THEN
              abserr = DABS(delxnp)
              ncveq = np
              if(xk(np).ne.0.0) relerr=dabs(delxnp/xk(np))
            ENDIF
          ELSEIF(ID.EQ.2) THEN
            if(xk(np).ne.0.0d0) then
              delrel=delxnp/xk(np)
              IF (DABS(delrel).GT.relerr) THEN
                relerr = DABS(delrel)
                ncveq = np
                abserr=dabs(delxnp)
              ENDIF
            endif
          ENDIF
C
          xk(np) = xk(np) + delxnp
  610   CONTINUE
C
C ------- update rk, the residual array
C
        resdmx=-1.0d38
        DO 620 np=1,nnp
          rknp = rk(np) - alphk*zk(np)
          resdmx=dmax1(resdmx,dabs(rknp))
          rk(np)=rknp
  620   CONTINUE
c
      IF(ID.EQ.1) THEN
        IF(abserr.LT.eps .and. resdmx.LT.sqeps) RETURN
      ELSEIF(ID.EQ.2) THEN
        IF(relerr.LT.eps .and. resdmx.LT.sqeps) RETURN
      ENDIF
C                                                                       ppcg 785
C ------- print iteration residuals                                     ppcg 790
C                                                                       ppcg 795
        IF (ibug.NE.0 .and. ktim.NE.0)                                  ppcg 800
     .  WRITE (16,1001) it,ncveq,abserr,relerr,xk(ncveq),rk(ncveq),eps
C                                                                       ppcg 815
        betad = banum
  900 CONTINUE                                                          ppcg 865
      WRITE (16,1100)                                                   ppcg 870
C     PRINT 1100                                                        ppcg 875
C                                                                       ppcg 880
 1000 FORMAT(1X,'   it',' ncveq','      abserr','      relerr',
     .          '   xk(ncveq)','   rk(ncveq)','         eps'/
     .       1X,'   --',' -----',' -----------',' -----------',
     .          ' -----------',' -----------',' -----------')
 1001 FORMAT(1X,I5,1x,i5,5(1pd12.4))
 1100 FORMAT(1X,'===> WARNING: nonconvergence in ppcg occurs')          ppcg 910
C                                                                       ppcg 915
      END                                                               ppcg 920
C                                                                       ppcg 925
C                                                                       ppcg 930
C                                                                       ppcg 935
c                                                                       ppcg 940
      FUNCTION dotprd(vectr1,vectr2,neqn)
C                                                                       dotp 010
C ------- Calcuate dot product of two vectors                           dotp 015
C                                                                       dotp 020
      IMPLICIT REAL*8(a-h,o-z)                                          dotp 025
      DIMENSION vectr1(neqn),vectr2(neqn)
C                                                                       dotp 035
      sum = 0.0d0                                                       dotp 040
      DO ieq=1,neqn                                                     dotp 045
        sum = sum + vectr1(ieq)*vectr2(ieq)
      ENDDO                                                             dotp 055
C                                                                       dotp 060
      dotprd = sum
C                                                                       dotp 070
      END                                                               dotp 075
C                                                                       dotp 080
C                                                                       dotp 085
C                                                                       dotp 090
C                                                                       dotp 095
       SUBROUTINE POLYP
     I   (C,RK,LRN,NLRN, IEIGEN,GG, MXTNOD,MAXNP,MXJBD,NNP,
     M    ZK,
     O    SK)
C                                                                       poly 010
C ------- Polynomial preconditioner for conjugated gradient method;     poly 015
C ------- the input matrix c must in compacted form with nodal connectivpoly 020
C ------- given by array lrn(nnp,jband);  This routine is written       poly 025
C ------- especially for multi-region groundwater flow and transport    poly 030
C ------- modelings; The input and output are                           poly 035
C ------- INPUT:                                                        poly 040
C -------   lrn(jp,np) = nodes connected to node np, note that the      poly 045
C -------                diagonal nodes of lrn must be moved to the     poly 050
C ------                 first column of matrix lrn                     poly 055
C -------   c(neqn,jj) = the coefficient matrix of matrix equation cx=b poly 060
C -------   rk(neqn)   = residual of one iterate, namely, rk=b-cx       poly 065
C -------   zk(neqn)   = working array of dimension neqn                poly 070
C -------   nnp        = # of nodes                                     poly 075
C -------   jband      = the band width of matrix lrn                   poly 080
C -------   nbwd       = the band width of matrix c                     poly 085
C -------   gg         = the upper bound on the maximum eigenvalue of c poly 090
C -------   ieigen     = signal of parameter estimation for gg          poly 095
C -------                zero: not requested                            poly 100
C -------                non-zero: requested                            poly 105
C ------- OUTPUT:                                                       poly 110
C -------   sk(neqn)   = deviate of two successive iterates, namely     poly 115
C -------                sk=x(k+1)-x(k)                                 poly 120
C                                                                       poly 125
      IMPLICIT REAL*8 (a-h,o-z)                                         poly 130
C                                                                       poly 135
      DIMENSION sk(maxnp),zk(maxnp),nlrn(mxtnod)                        poly 155
      DIMENSION lrn(mxjbd,maxnp),c(maxnp,mxjbd),rk(maxnp)               poly 160
C                                                                       poly 165
      DATA c00,c10,c20/-0.46875d0,1.6875d0,-2.25d0/                     poly 170
C                                                                       poly 175
c     jband=mxjbd                                                       poly 180
      nbwd=mxjbd                                                        poly 185
C                                                                       poly 190
C ------- estimate of optimization coefficients requested               poly 195
C                                                                       poly 200
      IF (ieigen.NE.0) THEN                                             poly 205
        gg = -1.0d0                                                     poly 210
        DO 110 irow=1,nnp                                               poly 215
          rowsum = 0.0d0                                                poly 220
          DO 100 icol=1,nbwd                                            poly 225
            rowsum = rowsum + DABS(c(irow,icol))                        poly 230
  100     CONTINUE                                                      poly 235
          IF (rowsum.GT.gg) gg=rowsum                                   poly 240
  110   CONTINUE                                                        poly 245
      ENDIF                                                             poly 250
      c0 = c00*gg*gg*gg                                                 poly 255
      c1 = c10*gg*gg                                                    poly 260
      c2 = c20*gg                                                       poly 265
C                                                                       poly 270
C ------- calculate sk, iterate 1                                       poly 275
C                                                                       poly 280
      DO 210 np=1,nnp                                                   poly 285
          crprod = 0.0d0
          DO 200 jbwd=1,nlrn(np)                                        poly 295
            IF (lrn(jbwd,np).EQ.0) GO TO 200                            poly 300
            crprod = crprod + c(np,jbwd)*rk(lrn(jbwd,np))
  200     CONTINUE                                                      poly 310
          zk(np) = c2*rk(np) + crprod
  210 CONTINUE                                                          poly 320
C                                                                       poly 325
C ------- calculate sk, iterate 2                                       poly 330
C                                                                       poly 335
      DO 500 ical=2,3                                                   poly 340
        IF (ical.EQ.2) THEN                                             poly 345
          ck = c1                                                       poly 350
        ELSE                                                            poly 355
          DO 300 np=1,nnp                                               poly 360
            zk(np) = sk(np)                                             poly 365
  300     CONTINUE                                                      poly 370
          ck = c0                                                       poly 375
        ENDIF                                                           poly 380
C                                                                       poly 385
        DO 410 np=1,nnp                                                 poly 390
            crprod = 0.0d0
            DO 400 jbwd=1,nlrn(np)                                      poly 400
                IF (lrn(jbwd,np).EQ.0) GO TO 400                        poly 405
                crprod = crprod + c(np,jbwd)*zk(lrn(jbwd,np))
  400       CONTINUE                                                    poly 415
            sk(np) = ck*rk(np) + crprod
  410   CONTINUE                                                        poly 425
  500 CONTINUE                                                          poly 430
C                                                                       poly 435
      END                                                               poly 440
C                                                                       poly 445
C                                                                       poly 450
C                                                                       poly 455
C                                                                       poly 460
      SUBROUTINE ILUCG
     I  (A,B,NC,ND,NT, EPS,SQEPS,IBUG,KTIM, MAXNP,MXTNOD,MXJBD,
     M  NNP,ID,  S,R,Q,P,AA,IL,
     O  X)
C                                                                       iluc 015
      IMPLICIT REAL*8(A-H,O-Z)                                          iluc 020
C                                                                       iluc 025
C---- To solves sparse asymmetric system of equations.                  iluc 030
C  A - The coefficient matrix                                           iluc 035
C  B - The right hand side vector coming in; the solution vector going oiluc 040
C  AA,IL,P,Q,R,S,X - Dummy variables                                    iluc 045
C ---  EPS - The reduction in error from the initial guess error. A smaliluc 050
C            value will result in more accurate solution but more iteratiluc 055
C                                                                       iluc 060
      DIMENSION a(maxnp,mxjbd),aa(maxnp,mxjbd),b(maxnp),il(maxnp)       iluc 080
      DIMENSION nc(mxjbd,maxnp),nd(maxnp),nt(MXTNOD),p(maxnp)           iluc 085
      DIMENSION q(maxnp),s(maxnp),r(maxnp),x(maxnp)                     iluc 090
C                                                                       iluc 095
      mb=mxjbd                                                          iluc 100
      nn=nnp                                                            iluc 105
C                                                                       iluc 110
      KOUNT=0                                                           iluc 115
      ITMAX=nnp
C                                                                       iluc 135
C---- store the coefficient matrix in another location for LU decompositiluc 140
C                                                                       iluc 145
      DO 10 I=1,NN                                                      iluc 150
      DO 10 J=1,MB                                                      iluc 155
 10   AA(I,J)=A(I,J)                                                    iluc 160
C                                                                       iluc 165
C---- approximate LU decomposition of the sparse matrix AA              iluc 170
C                                                                       iluc 175
      DO 20 I=1,NN                                                      iluc 180
 20   IL(I)=0                                                           iluc 185
      DO 30 I=1,NN                                                      iluc 190
      IF(ND(I).EQ.1) GO TO 1                                            iluc 195
      DO 40 J=ND(I),NT(I)                                               iluc 200
      NJ=NC(j,i)                                                        iluc 205
      DO 50 K=1,ND(I)-1                                                 iluc 210
      NK=NC(k,i)                                                        iluc 215
      DO 60 JK=ND(NK)+1,NT(NK)                                          iluc 220
      IF(NC(jk,nk).EQ.NJ) THEN                                          iluc 225
        AA(I,J)=AA(I,J)-AA(I,K)*AA(NK,JK)                               iluc 230
        GO TO 50                                                        iluc 235
      ENDIF                                                             iluc 240
 60   CONTINUE                                                          iluc 245
 50   CONTINUE                                                          iluc 250
 40   CONTINUE                                                          iluc 255
  1   IF(ND(I).EQ.NT(I)) GO TO 30                                       iluc 260
      DO 70 J=ND(I)+1,NT(I)                                             iluc 265
      NJ=NC(j,i)                                                        iluc 270
      IL(NJ)=IL(NJ)+1                                                   iluc 275
      NL=NC(il(nj),nj)                                                  iluc 280
      IF(IL(NJ).EQ.1) GO TO 70                                          iluc 285
      DO 80 K=1,IL(NJ)-1                                                iluc 290
      NK=NC(k,nj)                                                       iluc 295
      DO 90 LK=ND(NK)+1,NT(NK)                                          iluc 300
      IF(NC(lk,nk).EQ.NL) THEN                                          iluc 305
        AA(NJ,IL(NJ))=AA(NJ,IL(NJ))-AA(NJ,K)*AA(NK,LK)                  iluc 310
        GO TO 80                                                        iluc 315
      ENDIF                                                             iluc 320
 90   CONTINUE                                                          iluc 325
 80   CONTINUE                                                          iluc 330
 70   AA(NJ,IL(NJ))=AA(NJ,IL(NJ))/AA(I,ND(I))                           iluc 335
 30   CONTINUE                                                          iluc 340
C                                                                       iluc 345
C---- solution of the equations (LU)X=B, to get an initial guess        iluc 350
C                                                                       iluc 355
      DO 100 I=1,NN                                                     iluc 360
      X(I)=B(I)                                                         iluc 365
      IF(ND(I).EQ.1) GO TO 100                                          iluc 370
      DO 110 J=1,ND(I)-1                                                iluc 375
      NJ=NC(j,i)                                                        iluc 380
110   X(I)=X(I)-AA(I,J)*X(NJ)                                           iluc 385
100   CONTINUE                                                          iluc 390
      DO 120 I=NN,1,-1                                                  iluc 395
      IF(ND(I).EQ.NT(I)) GO TO 120                                      iluc 400
      DO 130 J=ND(I)+1,NT(I)                                            iluc 405
      NJ=NC(j,i)                                                        iluc 410
130   X(I)=X(I)-AA(I,J)*X(NJ)                                           iluc 415
120   X(I)=X(I)/AA(I,ND(I))                                             iluc 420
C                                                                       iluc 425
C---- error corresponding to initial guess                              iluc 430
C                                                                       iluc 435
      ERR0=0.                                                           iluc 440
      SUM0=0.                                                           iluc 445
      DO 140 I=1,NN                                                     iluc 450
      SUM=0.                                                            iluc 455
      DO 150 J=1,NT(I)                                                  iluc 460
      NJ=NC(j,i)                                                        iluc 465
150   SUM=SUM+A(I,J)*X(NJ)                                              iluc 470
      R(I)=B(I)-SUM                                                     iluc 475
      SUM0=SUM0+DABS(X(I))                                              iluc 480
140   ERR0=ERR0+DABS(R(I))                                              iluc 485
      IF(SUM0.EQ.0.) SUM0=1.                                            iluc 490
C                                                                       iluc 495
c ****** The following block is commented out by us to force iteration.
c     if(ibug.ne.0 .and. ktim.ne.0) then                                iluc 500
c       WRITE(16,1000)                                                  iluc 505
c       WRITE(16,1010) 0,eps,err0/sum0                                  iluc 510
c     endif                                                             iluc 515
C                                                                       iluc 520
c     IF(ERR0/SUM0.LT.eps) GO TO 3                                      iluc 525
c ***************************************************
C                                                                       iluc 530
c ******* indices not passed. JPG *******                               iluc 535
C                                                                       iluc 540
       CALL LLTINV
     I   (AA,R,NC,ND,NT, MXTNOD,MAXNP,MXJBD,NNP,
     M    IL,
     O    S)
      DO 160 I=1,NN                                                     iluc 550
160   IL(I)=0                                                           iluc 555
      DO 170 I=1,NN                                                     iluc 560
      Q(I)=0.                                                           iluc 565
      DO 170 J=1,NT(I)                                                  iluc 570
      NJ=NC(j,i)                                                        iluc 575
      IL(NJ)=IL(NJ)+1                                                   iluc 580
170   Q(I)=Q(I)+A(NJ,IL(NJ))*S(NJ)                                      iluc 585
C                                                                       iluc 590
C---- iteration loop starts here                                        iluc 595
C                                                                       iluc 600
  2   KOUNT=KOUNT+1                                                     iluc 605
C                                                                       iluc 610
C---- solution of the equations (UtU)X=B                                iluc 615
C                                                                       iluc 620
      DO 180 I=1,NN                                                     iluc 625
      IL(I)=ND(I)                                                       iluc 630
      P(I)=Q(I)                                                         iluc 635
      IF(ND(I).EQ.1) GO TO 180                                          iluc 640
      DO 190 J=1,ND(I)-1                                                iluc 645
      NJ=NC(j,i)                                                        iluc 650
      IL(NJ)=IL(NJ)+1                                                   iluc 655
190   P(I)=P(I)-AA(NJ,IL(NJ))*P(NJ)                                     iluc 660
180   P(I)=P(I)/AA(I,ND(I))                                             iluc 665
      DO 200 I=NN,1,-1                                                  iluc 670
      IF(ND(I).EQ.NT(I)) GO TO 200                                      iluc 675
      DO 210 J=ND(I)+1,NT(I)                                            iluc 680
      NJ=NC(j,i)                                                        iluc 685
210   P(I)=P(I)-AA(I,J)*P(NJ)                                           iluc 690
200   P(I)=P(I)/AA(I,ND(I))                                             iluc 695
C                                                                       iluc 700
C---- find the value of Alpha for optimum step length                   iluc 705
C                                                                       iluc 710
      ALNUM=0.                                                          iluc 715
      ALDEN=0.                                                          iluc 720
      DO 220 I=1,NN                                                     iluc 725
      ALNUM=ALNUM+R(I)*S(I)                                             iluc 730
220   ALDEN=ALDEN+P(I)*Q(I)                                             iluc 735
      IF(ALDEN.EQ.0.) GO TO 3                                           iluc 740
      ALPHA=ALNUM/ALDEN                                                 iluc 745
c ******* The following block up to the statement
c ******* if(err/err0.lt.eps) go to 3 is the original.  This
c ******* block will be modified so that correct convergent
c ******* criterion is used.
c      DO 230 I=1,NN                                                    iluc 750
c230   X(I)=X(I)+ALPHA*P(I)                                             iluc 755
c      ERR=0.                                                           iluc 760
c     DO 240 I=1,NN                                                     iluc 765
c     SUM=0.                                                            iluc 770
c     DO 250 J=1,NT(I)                                                  iluc 775
c     NJ=NC(j,i)                                                        iluc 780
c  250   SUM=SUM+A(I,J)*P(NJ)                                           iluc 785
c     R(I)=R(I)-ALPHA*SUM                                               iluc 790
c  240   ERR=ERR+DABS(R(I))                                             iluc 795
C                                                                       iluc 800
c      if(ibug.ne.0 .and. ktim.ne.0) then                               iluc 805
c        WRITE(16,1010) kount,eps,err/err0                              iluc 810
c     endif                                                             iluc 815
C                                                                       iluc 820
c     IF(ERR/ERR0.LT.EPS) GO TO 3                                       iluc 825
c *****************************************************
c
      abserr=-1.0d0
      relerr=-1.0d0
      DO 230 I=1,NN
        delxnp=alpha*p(i)
        if(id.eq.1) then
          if(dabs(delxnp).gt.abserr) then
            abserr=dabs(delxnp)
            ncveq=i
            if(x(i).ne.0.0d0) relerr=dabs(delxnp/x(i))
          endif
        elseif(id.eq.2) then
          if(x(i).ne.0.0d0) then
            delrel=delxnp/x(i)
            if(dabs(delrel) .gt. relerr) then
              relerr=dabs(delrel)
              ncveq=i
              abserr=dabs(delxnp)
            endif
          endif
        endif
        X(I)=X(I)+ALPHA*P(I)
  230 continue
c
      resdmx=-1.0d38
      DO 240 I=1,NN
        SUM=0.
        DO J=1,NT(I)
          NJ=NC(j,i)
          SUM=SUM+A(I,J)*P(NJ)
        ENDDO
        rknp=r(i)-alpha*sum
        resdmx=dmax1(resdmx,rknp)
        r(i)=rknp
  240 continue
c
      if(id.eq.1) then
        if(abserr.LT.eps .and. resdmx.LT.sqeps) go to 3
      elseif(id.eq.2) then
        if(relerr.LT.eps .and. resdmx.LT.sqeps) go to 3
      endif
c
C ----- PRINT ITERATION RESIDUALS
C
      IF(IBUG.NE.0 .AND. KTIM.NE.0 .AND. KOUNT.EQ.1)WRITE(16,1000)
      if(ibug.ne.0 .and. ktim.ne.0)
     . write(16,1001) kount,ncveq,abserr,relerr,x(ncveq),r(ncveq),eps
C                                                                       iluc 830
       CALL LLTINV
     I   (AA,R,NC,ND,NT, MXTNOD,MAXNP,MXJBD,NNP,
     M    IL,
     O    S)
C                                                                       iluc 850
C---- find the value of Beta                                            iluc 855
C                                                                       iluc 860
      BETNUM=0.                                                         iluc 865
      DO 260 I=1,NN                                                     iluc 870
260   BETNUM=BETNUM+R(I)*S(I)                                           iluc 875
      BETA=BETNUM/ALNUM                                                 iluc 880
C     IF(ABS(BETA).LT.EPS) GO TO 3                                      iluc 885
      DO 270 I=1,NN                                                     iluc 890
270   IL(I)=0                                                           iluc 895
      DO 280 I=1,NN                                                     iluc 900
      SUM=0.                                                            iluc 905
      DO 290 J=1,NT(I)                                                  iluc 910
      NJ=NC(j,i)                                                        iluc 915
      IL(NJ)=IL(NJ)+1                                                   iluc 920
290   SUM=SUM+A(NJ,IL(NJ))*S(NJ)                                        iluc 925
280   Q(I)=BETA*Q(I)+SUM                                                iluc 930
C                                                                       iluc 935
      IF(KOUNT.GT.ITMAX) THEN                                           iluc 940
c       WRITE(16,2015) err                                              iluc 945
        write(16,1100)
        GO TO 3                                                         iluc 950
      ENDIF                                                             iluc 955
      GO TO 2                                                           iluc 960
C                                                                       iluc 965
  3   DO 300 I=1,NN                                                     iluc 970
300   B(I)=X(I)                                                         iluc 975
C                                                                       iluc 980
      RETURN                                                            iluc 985
C                                                                       iluc 990
c1000 FORMAT('  it','         eps','       error'/                      iluc 995
c    .       '  --','         ---','       -----')                      iluc1000
c1010 FORMAT(I4,2D12.4)                                                 iluc1005
 1000 FORMAT(1X,'   it',' ncveq','      abserr','      relerr',
     .          '   xk(ncveq)','   rk(ncveq)','         eps'/
     .       1X,'   --',' -----',' -----------',' -----------',
     .          ' -----------',' -----------',' -----------')
 1001 FORMAT(1X,I5,1x,i5,5(1pd12.4))
c2015 format(' '/1x,'*** Max. No. of Iterations Exceeded in ILUCG'/1x,  iluc1010
c    1 '*** error = ',d12.4)                                            iluc1015
 1100 format(1x,'WARNING: Nonconvergency in ILUCG occurs')
C                                                                       iluc1020
      END                                                               iluc1025
C                                                                       iluc1030
C                                                                       iluc1035
C                                                                       iluc1040
C                                                                       iluc1045
      SUBROUTINE LLTINV
     I   (A,B,NC,ND,NT, MXTNOD,MAXNP,MXJBD,NNP,
     M    IL,
     O    X)
C                                                                       llti 010
      IMPLICIT REAL*8(A-H,O-Z)                                          llti 015
C                                                                       llti 020
C---- For solving the set of equations (LLt)X=B                         llti 025
C                                                                       llti 030
      DIMENSION a(maxnp,mxjbd),b(maxnp),il(maxnp),nc(mxjbd,maxnp)       llti 050
      DIMENSION nd(maxnp),nt(mxtnod),x(maxnp)                           llti 055
c     DIMENSION A(NN,MB),B(NN),IL(NN),NC(mb,nn),ND(NN),NT(NN),X(NN)     llti 060
c                                                                       llti 065
      nn=nnp                                                            llti 070
c     mb=mxjbd                                                          llti 075
C                                                                       llti 080
      DO 10 I=1,NN                                                      llti 085
      X(I)=B(I)                                                         llti 090
      IF(ND(I).EQ.1) GO TO 10                                           llti 095
      DO 20 J=1,ND(I)-1                                                 llti 100
      NJ=NC(j,i)                                                        llti 105
 20   X(I)=X(I)-A(I,J)*X(NJ)                                            llti 110
 10   CONTINUE                                                          llti 115
      DO 30 I=NN,1,-1                                                   llti 120
      IL(I)=ND(I)                                                       llti 125
      IF(ND(I).EQ.NT(I)) GO TO 30                                       llti 130
      DO 40 J=ND(I)+1,NT(I)                                             llti 135
      NJ=NC(j,i)                                                        llti 140
      IL(NJ)=IL(NJ)-1                                                   llti 145
 40   X(I)=X(I)-A(NJ,IL(NJ))*X(NJ)                                      llti 150
 30   CONTINUE                                                          llti 155
      RETURN                                                            llti 160
      END                                                               llti 165
C
C
      SUBROUTINE ELENOD
     I    (IEM1,IEM2,
     O     NODE,NSIDE,IK)
C
      IF(IEM1.EQ.0)THEN
        NODE=4
        NSIDE=4
        IK=3
      ELSEIF(IEM2.EQ.0)THEN
        NODE=6
        NSIDE=5
        IK=2
      ELSE
        NODE=8
        NSIDE=6
        IK=1
      ENDIF
      RETURN
      END
C
      SUBROUTINE WARMSG
     I      (N,MAXN,SUBNAM,VARNAM,NO)
C
      CHARACTER*6 SUBNAM,VARNAM
C
      IF(N.GT.MAXN)THEN
        WRITE(16,10)N,VARNAM,SUBNAM,NO
  10    FORMAT(I8,' .GT.',A6,' IN SUB. ',A6,' ----->',I2)
        STOP
      ENDIF
      RETURN
      END
