C $$$$$$$$$$$$ CHEMI3D.F       4/22/95
C  5/9/95 ---  Q34ADB (VARIABLE BC)
c  2/1/95 ---  nlrnd=nlrn when idzoom=0
C      8/17/94         4:30 pm
C
      SUBROUTINE HMCTRN(C,CP,CS,DTI,F, RI,RL, X,IE,IB,LRL,NLRL,CMX,
c     . sk,rk,pk,aa,il,nd,NLRN, IRHO,AKHC,RHOMU,cnstkr,DINTS,
     . sk,rk,pk,aa,il,nd,NLRN, IRHO,AKHC,cnstkr,DINTS,
     1 LRN, CMATRX,RLD, NNPLR,LMAXDF,GNLR,LNOJCN,CMTRXL,RLDL, DCOSB,ISB,
     2 NPBB,BFLX, WETAB,V,VP,VEAVG,H,TH,THP, WWRK, THN,AKDC,
     3 LES,SOS,ISTYP, NPW,WSS,IWTYP, ILSV,IMSV, QCB,ICTYP,ISC,NPCB,
     > QNB,INTYP,ISN,NPNB, CVB,IVTYP,ISV,NPVB, CDB,IDTYP,NPDB, RHOTYP,
     > PROP,SPP,RKD,TRANC,NBDYB,IBDY,XW,VXW,CW, MWLOC,LRLW,NLRLW,IEW,
     > IBW,DL468, ieww,ibww,lrlww,nlrlww,xww,vxww,dl468w,
     > IBE,IBCHK,XPFG,CPFG,MPLOC,XSFG,CSFG,MPLOCS,XWFG,CWFG,MPLOCW,
     > XEFG,CEFG,MPLOCE,NFGM,NFGMB,NFGMBB,MAXFGW,MINFGW,CMAXFG,CMINFG,
     > ISE,NEPWN,NEPW,ISED,NDBD,MPLOCD,NLRND,NCFG,DTIFG,NPFG,NEFGS,
     7 KPR,KDSK,KDIG,KOUT,JTM,IRXN, IBUG,TITLE,NPROB,EPSX,SQEPS)
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 GNLR
      CHARACTER TITLE*70
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /LGEOM/ LTMXNP,LMXNP,LMXBW,MXREGN,NREGN
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
      COMMON /PCG/ GG,IEIGEN
      COMMON /NINTR/ KPR0,KDSK0,NSTRf,NSTRt,KSSf,KSS,IGEOM
      COMMON /CTIM/ DELT,CHNG,DELMAX,TMAX,DELT0,TIME
C
      COMMON /FINTE/ NCYLf,NITERf,NPITERf,KSP,KGRAV,IPNTSf
      COMMON /TINTE/ NCMt,KVIt,NITERt,NPITERt,IPNTSt,MICONF,IFLUX
      COMMON /TREAL/ W,WV,OME,OMI,TOLA,TOLB
      COMMON /TADP/ ADPEPS,ADPARM,IZOOM,IDZOOM,IEPC,NXG,NYG,NZG,
     >              NXW,NYW,NZW,NXD,NYD,NZD,IDETQ
C
      COMMON /CELS/ MXSEL,MXSPR,MXSDP,NSEL,NSPR,NSDP,KSAI
      COMMON /CNPS/ MXWNP,MXWPR,MXWDP,NWNP,NWPR,NWDP,KWAI
C
      COMMON /TCBC/ MXCNP,MXCES,MXCPR,MXCDP,NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /TNBC/ MXNNP,MXNES,MXNPR,MXNDP,NNNP,NNES,NNPR,NNDP,KNAI
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /TDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
      COMMON /TFLUX/ MXLSV,MXMSV
C
      COMMON /WETX/ APHA1,APHA2,APAH3,APHA4
      COMMON /WETY/ BETA1,BETA2,BETA3,BETA4
      COMMON /WETZ/ GAMA1,GAMA2,GAMA3,GAMA4
C
      COMMON /TFLOW/ FRATE(14),FLOW(14),TFLOW(14,7)
C
      COMMON /SAZFM/ MXNPFG,MXKGL,MXKGLD,MXNEP,MXEPW,MXNPW,MXELW,MXNDB,
     >               mxnpww,mxelww
      COMMON /CHEM/ MXNCC,NCC
C
      DIMENSION C(MAXNP,MXNCC),CP(MAXNP,MXNCC),CMX(MXNCC)
      DIMENSION CS(MAXNP,MXNCC),DTI(MAXNP,MXNCC)
      DIMENSION F(MAXNP,3,MXNCC)
      DIMENSION RI(MXADNP),RL(MXADNP)
      DIMENSION sk(mxadnp),rk(mxadnp),pk(mxadnp),aa(mxadnp,mxjbd)
      DIMENSION il(mxadnp),nd(mxadnp),NLRN(maxnp)
C
      DIMENSION X(MAXNP,3),IE(MAXEL,11),LRL(MXKBD,MAXNP),NLRL(MAXNP)
      DIMENSION IB(MAXNP),LRN(MXJBD,MXADNP)
      DIMENSION CMATRX(MXADNP,MXJBD),RLD(MXADNP)
C
      DIMENSION NNPLR(MXREGN),LMAXDF(MXREGN)
      DIMENSION GNLR(LTMXNP,MXREGN),LNOJCN(MXJBD,LMXNP,MXREGN)
      DIMENSION CMTRXL(LMXNP,LMXBW),RLDL(LMXNP)
C
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES),NPBB(MAXBNP)
      DIMENSION BFLX(MAXBNP,2,MXNCC)
C
      DIMENSION WETAB(12,MAXEL), V(MAXNP,3),VP(MAXNP,3)
      DIMENSION VEAVG(MAXNP,3,MXNCC),WWRK(MAXNP),AKHC(8,MAXEL,7)
      DIMENSION TH(MAXEL,8),THP(MAXEL,8),THN(MAXNP,2,MXNCC),
     >          AKDC(8,MXKGLD,8),H(MAXNP)
      DIMENSION RHOTYP(MAXMAT),SPP(MXSPPM,MAXMAT,4)
c      DIMENSION RHOMU(MXNCC),DINTS(MXNCC)
      DIMENSION DINTS(MXNCC)
C
C     DIMENSION SOSF(MXSDP,MXSPR,2),TSOSF(MXSDP,MXSPR)
      DIMENSION ISTYP(MXSEL,MXNCC),LES(MXSEL),SOS(MXSPR,2)
C     DIMENSION WSSF(MXWDP,MXWPR,2),TWSSF(MXWDP,MXWPR)
      DIMENSION IWTYP(MXWNP,MXNCC),NPW(MXWNP),WSS(MXWPR,2)
C
C     DIMENSION QCBF(MXCDP,MXCPR),TQCBF(MXCDP,MXCPR)
      DIMENSION ICTYP(MXCES,MXNCC),ISC(5,MXCES),NPCB(MXCNP),QCB(MXCPR)
C
C     DIMENSION CVBF(MXVDP,MXVPR),TCVBF(MXVDP,MXVPR)
      DIMENSION IVTYP(MXVES,MXNCC),ISV(5,MXVES),NPVB(MXVNP),CVB(MXVPR)
C
C     DIMENSION CDBF(MXDDP,MXDPR),TCDBF(MXDDP,MXDPR)
      DIMENSION IDTYP(MXDNP,MXNCC),NPDB(MXDNP),CDB(MXDPR)
C
C     DIMENSION QNBF(MXNDP,MXNPR),TQNBF(MXNDP,MXNPR)
      DIMENSION INTYP(MXNES,MXNCC),ISN(5,MXNES),NPNB(MXNNP),QNB(MXNPR)
C
      DIMENSION RKD(MAXMAT,MXNCC),TRANC(MAXMAT,MXNCC)
      DIMENSION PROP(MAXMAT,MXMPPM),NBDYB(MAXNP),IBDY(MXTUBS)
      DIMENSION XW(MXNPW,3),VXW(MXNPW,3),LRLW(24,MXNPW,3),MWLOC(MXNPW)
      DIMENSION NLRLW(MXNPW,3),IEW(MXELW,8,3),IBW(MXNPW,3),CW(MXNPW)
      DIMENSION DL468(8,MXNPW,3)
c
      DIMENSION XWw(MXNPWw,3),VXWw(MXNPWw,3),LRLWw(24,MXNPWw,3)
      DIMENSION NLRLWw(MXNPWw,3),IEWw(MXELWw,8,3),IBWw(MXNPWw,3)
      DIMENSION DL468w(8,MXNPWw,3)
C
      DIMENSION KPR(MXNTI),KDSK(MXNTI)
C
      DIMENSION IBE(MAXEL),IBCHK(MAXEL)
      DIMENSION XPFG(MXNPFG,3),CPFG(MXNPFG,MXNCC),MPLOC(MXNPFG)
      DIMENSION XSFG(MXNPFG,3),CSFG(MXNPFG,MXNCC),MPLOCS(MXNPFG)
      DIMENSION XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC),MPLOCW(MXNPFG)
      DIMENSION XEFG(MXNEP,3),CEFG(MXNEP,MXNCC),MPLOCE(MXNPFG)
      DIMENSION NFGM(MAXEL),NFGMB(MAXEL),MAXFGW(MXELW,MXNCC)
      DIMENSION MINFGW(MXELW,MXNCC),CMAXFG(MXELW,MXNCC),NEPWN(MXELW)
      DIMENSION CMINFG(MXELW,MXNCC),ISE(MXKGL,8),NEPW(MXELW,MXEPW)
      DIMENSION ISED(MXKGLD,9),NDBD(MXNDB),MPLOCD(MXNDB),NLRND(MXADNP)
      DIMENSION ILSV(MXLSV,5),IMSV(MXMSV,4),NFGMBB(MAXEL)
      DIMENSION NCFG(2*MXNCC),DTIFG(MXNPFG,MXNCC)
      DIMENSION CC(8,7),CQW(7),CQP(7),XX(8),YY(8),ZZ(8),DL(8),DNX(8)
C
      DATA RAMADA/0.0D0/
C
      IF(KSS.NE.0) GO TO 600
C
C $$$$$$$ PERFORM STEADY STATE COMPUTATION
C
      DO 160 K=1,NCC
      DO 160 NP=1,NNP
      C(NP,K)=OME*C(NP,K)+(1.0D0-OME)*CP(NP,K)
  160 CONTINUE
C
      DO 170 M=1,NEL
      DO 170 IQ=1,8
      THP(M,IQ)=TH(M,IQ)
  170 CONTINUE
C
      KDIG=KDIG+1
      IF(IBUG.NE.0) WRITE(16,4800) KDIG,TIME,DELT
      IF(IBUG.NE.0) WRITE(16,4850)
C
      CALL THNODE(THN(1,1,1),TH,THP,PROP,RKD(1,1),WWRK,IE,X)
      CALL DISPC
     I    (X,IE,V,VP,TH,THP,H,PROP,W,KSS,IQUAR,XWFG,CWFG,NEL,cnstkr,
c     I     ISED,1,MAXEL,7,KSP,KVIT,SPP,RHOMU,DINTS,IRHO,RHOTYP,
     I     ISED,1,MAXEL,7,KSP,KVIT,SPP,DINTS,IRHO,RHOTYP,
     O     AKHC)
      CALL AFABTA
     I       (X,IE,  V,VP,  W,KSS,IOPTIM,  THN(1,2,1),PROP,
     O        WETAB)
C
C ******* BEGIN THE NONLINEAR ITERATION LOOP
C
      EPS=0.5D0*TOLA
      ID=3
      IFLUX=0
      DO 540 ITER=1,NITERt
C
C ------- ASSEMBEL THE COEFFICIENT MATRIX AND CONSTRUCT THE LOAD VECTOR
C
        DO 300 K=1,NCC
          DO NP=1,NNP
            RI(NP)=C(NP,K)
          ENDDO
          CALL TASEMB
     >      (CMATRX,RLD,C,CP,CS(1,K),X,IE,LRN,NLRL,LRL,WETAB,V,VP,TH,
     >       THP,AKHC,AKDC, LES,ISTYP,SOS,NPW,IWTYP,WSS,IRHO,RHOTYP,
     >       PROP,DINTS,RKD,TRANC,DTI(1,K),W,WV,MICONF,ID,
     >       KSS,ISED,XPFG,MAXNP,MAXNP,NEL,K)
C
C ------- APPLY BOUNDARY CONDITIONS
C
          IF(MICONF.EQ.0 .OR. K.GT.3)THEN
            CALL TBC
     >         (CMATRX,RLD,CS(1,K),X,IE,LRN,NLRN, DCOSB,ISB,
     1          V,VP, QCB,ISC,ICTYP(1,K), QNB,ISN,INTYP(1,K),
     2          CVB,ISV,IVTYP(1,K), CDB,IDTYP(1,K),NPDB, W,KSS,ID,MAXNP)
          ENDIF
C
C ------- SOLVE THE MATRIX EQAUTION WITH BLOCK ITERATIONS
C
c         WRITE(16,*)'THE FOLLOWING HISTORY IS FOR TRANSPORT PART'
          IF(IPNTSt.EQ.0) THEN
            CALL BLKITR(RL,RI,CMTRXL,RLDL, CMATRX,RLD, GNLR,
     >           LNOJCN,NNPLR,LMAXDF,OMI,EPS,NPITERt,IBUG,KPR0,2)
          ELSE IF (IPNTSt.EQ.1) THEN
            CALL PISS
     I        (MAXNP,MXJBD,MXADNP, NNP,NLRN,NPITERt, OMI,EPS,
     I         KPR0,IBUG, CMATRX,RI,RLD,LRN,2,
     O         RL)
          ELSE IF (IPNTSt.EQ.2) THEN
            CALL PPCG
     I        (CMATRX,RLD,LRN,NLRN, IEIGEN,GG,EPS,SQEPS,IBUG,KPR0,
     M         MXADNP,MAXNP,MXJBD,NNP,2,  SK,RK,RI,PK,
     O         RL)
          ELSE IF (IPNTSt.EQ.3) THEN
            CALL ILUCG
     I       (CMATRX,RLD,LRN,ND,NLRN, EPS,SQEPS,IBUG,KPR0,
     M        MXADNP,MAXNP,MXJBD,NNP,2, SK,RK,RI,PK,AA,IL,
     O        RL)
          ENDIF
C
          DIFMAX=0.0
          NOCCUR=1
          DO 220 NP=1,NNP
            IF(C(NP,K).EQ.0.0) GO TO 220
            DIF=(C(NP,K)-RL(NP))/C(NP,K)
            DIF=DABS(DIF)
            IF(DIF.LE.DIFMAX) GO TO 220
            DIFMAX=DIF
            KMAX=K
            NOCCUR=NP
  220     CONTINUE
C
C ------- UPDATE NONLINEAR ITERATE
C
          DO 230 NP=1,NNP
            C(NP,K)=OME*RL(NP)+(1.0D0-OME)*C(NP,K)
  230     CONTINUE
  300   CONTINUE
C
        IF(NITERt.EQ.1) GO TO 550
        IF(IBUG.NE.0) WRITE(16,5400) ITER,KMAX,DIFMAX,TOLA,NOCCUR
        IF(ITER.EQ.1) GO TO 540
        IF(DIFMAX.LE.TOLA) GO TO 550
C
  540 CONTINUE
C
      WRITE(16,5500) ITER,NITERt,KMAX,DIFMAX,TOLA,NOCCUR
C
  550 CONTINUE
      DO 555 K=1,NCC
      DO 555 NP=1,NNP
  555 CP(NP,K)=C(NP,K)
C
      RETURN
C
C $$$$$$$ PERFORM TRANSIENT-STATE OR TRANSIENT COMPUTATION
C
C ------- COMPUTE SOURCE AND B. C. VALUE AT THE PRESENT TIME
C
  600 CONTINUE
      EPS = 0.5D0*TOLB
C
C ------- UPDATE CP AND ESTIMATE NONLINEAR ITERATE CW FOR COMPUTING
C ------- COEFFICIENT MATRIX AND LOAD VECTOR.
C
      DO 625 K=1,NCC
      DO 625 NP=1,NNP
        CP(NP,K)=C(NP,K)
  625 CONTINUE
C
        CALL DISPC
     I    (X,IE,V,VP,TH,THP,H,PROP,W,KSS,IQUAR,XWFG,CWFG,NEL,cnstkr,
c     I     ISED,1,MAXEL,7,KSP,KVIT,SPP,RHOMU,DINTS,IRHO,RHOTYP,
     I     ISED,1,MAXEL,7,KSP,KVIT,SPP,DINTS,IRHO,RHOTYP,
     O     AKHC)
C
C ******* COMPUTE ADVECTION CONCENTRATION
C
      IF(LGRN.NE.0) THEN
        NPFGW=0
        NPFGEP=0
        NPFGS=0
        DO M=1,NEL
          IBCHK(M)=0
        ENDDO
        DO I=1,MXNPFG
          MPLOCE(I)=0
        ENDDO
C
        DO 650 K=1,NCC
C
          CALL THNODE(THN(1,1,K),TH,THP,PROP,RKD(1,K),WWRK,IE,X)
C
          IF(MICONF.EQ.1 .AND. K.LE.3)THEN
            DO NP=1,NNP
              CS(NP,K)=CP(NP,K)
              DTI(NP,K)=1.0D0/DELT
              VEAVG(NP,1,K)=0.0D0
              VEAVG(NP,2,K)=0.0D0
              VEAVG(NP,3,K)=0.0D0
            ENDDO
            GOTO 650
          ENDIF
C
          DO I=1,NNP
            VEAVG(I,1,K)=0.5D0*(V(I,1)+VP(I,1))/THN(I,1,K)
            VEAVG(I,2,K)=0.5D0*(V(I,2)+VP(I,2))/THN(I,1,K)
            VEAVG(I,3,K)=0.5D0*(V(I,3)+VP(I,3))/THN(I,1,K)
          ENDDO
C
C +++++++ FOR THE CASE OF USING LAGRANGIAN APPROACH
C
C ------- DETERMINE NTAU AND DTAU FOR LAGRANGIAN STEP INTEGRATION
C
          CALL ADVBC
     I      (IE,X,V,VP,THN(1,1,K),CP(1,K),
     I       CVB,IVTYP(1,K),NPVB,ISV, QCB,ICTYP(1,K),NPCB,ISC,
     I       CDB,IDTYP(1,K),NPDB, DCOSB,ISB,NPBB,
     O       CS(1,K),
     M       RI,RL,RLD)
C
C ------- COMPUTE LAGRANGIAN CONCENTRATON CS
C
C %%%%% IDTI<0 --> Calculate concentrations CS and DTI
C %%%%% IDTI<0 & IBF=2 ----> Calculate cons. CSFG only by using CPFG
C %%%%% IDTI=0 --> Calculate concentrations CSFG only
C %%%%% IDTI=1 --> Calculate DTIFG only
C %%%%% IDTI=2 --> Calculate both DTIFG and CSFG arrays
C
          IBF=1
          IDTI=-1
ccc          print *,'k=',k,' call gntrak --->1'
          CALL GNTRAK
     I     (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,MXTUBS,
     I      MXNEP,NNP,NEL,NEFGS,NTUBS,NXW,NYW,NZW,IDTI,IZOOM,EPSX,
     I      IBF,IDETQ, DELT,RAMADA,CP(1,K),X,VEAVG(1,1,K),CPFG(1,K),
     I      XPFG,ISE,NFGMB, IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,ISB,DCOSB,
     I      NBDYB,IBDY,DL468, ieww,ibww,lrlww,nlrlww,dl468w,
     i      mxnpww,mxelww,nxg,nyg,nzg,nwnp,npw,iwtyp(1,k),wss(1,2),
     i      mxwnp,mxwpr,
     O      CS(1,K),DTI(1,K),DTIFG(1,K),NPFGS,XSFG,CSFG(1,K),MPLOCS,
     O      IBCHK,XEFG,CEFG(1,K),MPLOCE,
     M      XW,VXW,XWw,VXWw)
C
          IF(IZOOM.EQ.1)THEN
            NCFG(K)=NPFGW
            NCFG(K+NCC)=NPFGEP
C
C ------  forward tracking of activated fine grids
C
            IBF=2
            NINIT=1
ccc            print *,'k=',k,' call hptrak --->1'
            CALL HPTRAK
     I      (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,
     I      MXTUBS,MXNEP,MXNEP,MXNPFG,
     I      NNP,NEL,NEFGS,NTUBS,NXW,NYW,NZW, NPFG,NINIT,
     I      IBF,IDETQ, DELT,TMAX,RAMADA, IDTI,IZOOM,EPSX,
     I      CP,X,VEAVG(1,1,K),CPFG(1,K),XPFG,MPLOC,ISE,NFGMB,
     I      IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,ISB,DCOSB,NBDYB,IBDY,
     I      DL468,nwnp,npw,iwtyp(1,k),wss(1,2),mxwnp,mxwpr,
     O      CS(1,K),DTI(1,K),DTIFG(1,K),NPFGW,XWFG,CWFG(1,K),MPLOCW,
     O      NPFGEP,XEFG,CEFG(1,K),MPLOCE,
     M      XW,VXW)
C
C ------ forward tracking of global nodes
C
            IDTI=0
ccc            print *,'before call gntrak ---> 2'
            CALL GNTRAK
     I      (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,MXTUBS,
     I      MXNEP,NNP,NEL,NEFGS,NTUBS,NXW,NYW,NZW,IDTI,IZOOM,EPSX,
     I      IBF,IDETQ, DELT,RAMADA, CP(1,K),X,VEAVG(1,1,K),CPFG(1,K),
     I      XPFG,ISE,NFGMB,IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,ISB,DCOSB,
     I      NBDYB,IBDY,DL468, ieww,ibww,lrlww,nlrlww,dl468w,
     i      mxnpww,mxelww,nxg,nyg,nzg,nwnp,npw,iwtyp(1,k),wss(1,2),
     i      mxwnp,mxwpr,
     O      CS(1,K),DTI(1,K),DTIFG(1,K),NPFGW,XWFG,CWFG(1,K),MPLOCW,
     O      IBCHK,XEFG,CEFG(1,K),MPLOCE,
     M      XW,VXW,XWw,VXWw)
C
          ENDIF
C
C ------ INCORPORATE BOUNDARY CONDITION FOR THE LAGRANGIAN STEP
C
          CALL ADVBC
     I      (IE,X,V,VP,THN(1,1,K),CP(1,K),
     I       CVB,IVTYP(1,K),NPVB,ISV, QCB,ICTYP(1,K),NPCB,ISC,
     I       CDB,IDTYP(1,K),NPDB, DCOSB,ISB,NPBB,
     O       CS(1,K),
     M       RI,RL,RLD)
C
  650   CONTINUE
c
ccc        print *,'after backward and forward tracking, npfgw=',npfgw
ccc        print *,'after backward and forward tracking, npfgep=',npfgep
C
C ----- make up the concentrations of the other components at the point
C       generated by the K-th component
C ----- (1) the EPCOF points; (2) forward tracked nodes
C
        IF(IEPC.EQ.1)THEN
ccc          print *,'hptrak for epcof points'
          IBF=1
          IDTI=0
          DO K=1,NCC
            KWB=NCFG(K+NCC)+1
            IF(K.EQ.NCC)THEN
              KWE=NPFGEP
            ELSE
              KWE=NCFG(K+1+NCC)
            ENDIF
            DO IK=1,K-1
              CALL HPTRAK
     I          (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,
     I           MXTUBS,MXNEP,MXNPFG,MXNEP,
     I           NNP,NEL,NEFGS,NTUBS,NXW,NYW,NZW, NPFG,KWB,IBF,IDETQ,
     I           DELT,TMAX,RAMADA, IDTI,IZOOM,EPSX, CP(1,IK),X,
     I           VEAVG(1,1,IK),CPFG(1,IK),XPFG,MPLOC,ISE,NFGMB,IE,IB,
     I           LRL,NLRL,IEW,IBW,LRLW,NLRLW,ISB,DCOSB,NBDYB,IBDY,DL468,
     I           nwnp,npw,iwtyp(1,ik),wss(1,2),mxwnp,mxwpr,
     O           CS(1,IK),DTI(1,IK),DTIFG(1,IK),KWE,XEFG,CEFG(1,IK),
     O           MPLOCE,NPFGS,XSFG,CSFG(1,IK),MPLOCS,
     M           XW,VXW)
            ENDDO
            DO IK=K+1,NCC
              CALL HPTRAK
     I          (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,
     I           MXTUBS,MXNEP,MXNPFG,MXNEP,NNP,NEL,NEFGS,NTUBS,
     I           NXW,NYW,NZW, NPFG,KWB, IBF,IDETQ, DELT,TMAX,RAMADA,
     I           IDTI,IZOOM,EPSX,CP(1,IK),X,VEAVG(1,1,IK),CPFG(1,IK),
     I           XPFG,MPLOC,ISE,NFGMB,IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,
     I           ISB,DCOSB,NBDYB,IBDY,DL468,nwnp,npw,iwtyp(1,ik),
     I           wss(1,2),mxwnp,mxwpr,
     O           CS(1,IK),DTI(1,IK),DTIFG(1,IK),KWE,XEFG,CEFG(1,IK),
     O           MPLOCE,NPFGS,XSFG,CSFG(1,IK),MPLOCS,
     M           XW,VXW)
            ENDDO
          ENDDO
        ENDIF
C
        IF(IZOOM.EQ.1)THEN
ccc          print *,'hptrak for making up forward points'
          IBF=1
          IDTI=0
          DO K=1,NCC
            KWB=NCFG(K)+1
            IF(K.EQ.NCC)THEN
              KWE=NPFGW
            ELSE
              KWE=NCFG(K+1)
            ENDIF
            DO IK=1,K-1
              CALL HPTRAK
     I          (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,
     I           MXTUBS,MXNEP,MXNPFG,MXNPFG,NNP,NEL,NEFGS,NTUBS,
     I           NXW,NYW,NZW, NPFG,KWB, IBF,IDETQ, DELT,TMAX,RAMADA,
     I           IDTI,IZOOM,EPSX, CP(1,IK),X,VEAVG(1,1,IK),CPFG(1,IK),
     I           XPFG,MPLOC,ISE,NFGMB,IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,
     I           ISB,DCOSB,NBDYB,IBDY,DL468,nwnp,npw,iwtyp(1,ik),
     i           wss(1,2),mxwnp,mxwpr,
     O           CS(1,IK),DTI(1,IK),DTIFG(1,IK),KWE,XWFG,CWFG(1,IK),
     O           MPLOCW,NPFGEP,XSFG,CSFG(1,IK),MPLOCS,
     M           XW,VXW)
            ENDDO
            DO IK=K+1,NCC
              CALL HPTRAK
     I          (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,
     I           MXTUBS,MXNEP,MXNPFG,MXNPFG,NNP,NEL,NEFGS,NTUBS,
     I           NXW,NYW,NZW, NPFG,KWB, IBF,IDETQ, DELT,TMAX,RAMADA,
     I           IDTI,IZOOM,EPSX, CP(1,IK),X,VEAVG(1,1,IK),CPFG(1,IK),
     I           XPFG,MPLOC,ISE,NFGMB,IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,
     I           ISB,DCOSB,NBDYB,IBDY,DL468,nwnp,npw,iwtyp(1,ik),
     i           wss(1,2),mxwnp,mxwpr,
     O           CS(1,IK),DTI(1,IK),DTIFG(1,IK),KWE,XWFG,CWFG(1,IK),
     O           MPLOCW,NPFGEP,XSFG,CSFG(1,IK),MPLOCS,
     M           XW,VXW)
            ENDDO
          ENDDO
        ENDIF
C
C ----- solve Lagrangian concentrations due to reaction
C
C ----- the following block is for global nodes
C
ccc        print *,'advrx for cs'
        CALL ADVRX
     I    (CS, IRHO,CS,
c     I     MAXNP,MAXNP,NNP,IRXN,OME,EPS,EPSX,IBUG,X,IE,H,X(1,1),X(1,2),
     I     MAXNP,MAXNP,NNP,IRXN,OME,EPS,EPSX,X,IE,H,X(1,1),X(1,2),
c     I     X(1,3),NLRN,1,RKD,PROP,SPP,KSP,TH,THN,DTI,DELT,RHOMU,DINTS)
     I     X(1,3),NLRN,1,RKD,PROP,SPP,KSP,TH,THN,DTI,DELT,DINTS)
C
C ----- the following block is for forward tracked points
C
ccc        print *,'advrx for cwfg'
        CALL ADVRX
     I      (CWFG, IRHO,CS, MXNPFG,MAXNP,
c     I       NPFGW,IRXN,OME,EPS,EPSX,IBUG,X,IE,H,XWFG(1,1),XWFG(1,2),
     I       NPFGW,IRXN,OME,EPS,EPSX,X,IE,H,XWFG(1,1),XWFG(1,2),
     I       XWFG(1,3),MPLOCW,0,RKD,PROP,SPP,KSP,TH,THN,DTI,DELT,
c     I       RHOMU,DINTS)
     I       DINTS)
C
C ----- the following block is for EPCOF points
C
        CALL ADVRX
     I      (CEFG, IRHO,CS, MXNEP,MAXNP,
c     I       NPFGEP,IRXN,OME,EPS,EPSX,IBUG,X,IE,H,XEFG(1,1),XEFG(1,2),
     I       NPFGEP,IRXN,OME,EPS,EPSX,X,IE,H,XEFG(1,1),XEFG(1,2),
     I       XEFG(1,3),MPLOCE,0,RKD,PROP,SPP,KSP,TH,THN,DTI,DELT,
c     I       RHOMU,DINTS)
     I       DINTS)
C
C ------ determine sharp front elements
C
        IF(IZOOM.EQ.1)THEN
ccc          print *,'determining rough element'
          DO M=1,NEL
            IE(M,11)=0
          ENDDO
          DO K=1,NCC
            CALL SFDET
     I          (MAXNP,MAXEL,MXNPFG,MXNPFG,NPFGW,ADPARM,ADPEPS,CMX(K),
     I           X,CS(1,K),XWFG,CWFG(1,K),MPLOCW,
     M           IE)
          ENDDO
C
          DO K=1,NCC
            CALL SFDET
     I          (MAXNP,MAXEL,MXNEP,MXNPFG,NPFGEP,ADPARM,ADPEPS,CMX(K),
     I           X,CS(1,K),XEFG,CEFG(1,K),MPLOCE,
     M           IE)
          ENDDO
C
C ----- NOW, WE HAVE ALREADY DETERMINED ALL THE SHARP-FRONT ELEMENTS AT
C       THE CURRENT TIME STEP.  IN THE FOLLOWING, WE APPLY REGULAR
C       REFINEMENT FOR EVERY SHARP-FRONT ELEMENT AND FIND
C       OUT THE ADVECTION CONCENTRATION AT EACH HIDDEN POINT BY
C       BACKWARD TRACKING
C
ccc          print *,'before call fgdet'
          CALL FGDET
     I      (IE,IBCHK,X,IB,CP,CS,VEAVG,LRL,NLRL,ISE,NFGMB,XPFG,CPFG,
     I       MPLOC,IEW,IBW,LRLW,NLRLW,XW,CW,MWLOC,VXW,NEL,NEFGS,NPFG,
     I       NNP,NTUBS,NCC,MICONF,IEPC,NXG,NYG,NZG,NPFGEP,NXW,NYW,NZW,
     I       IDETQ,TMAX,DELT,RAMADA,ADPARM,ADPEPS,EPSX,IZOOM,CMX,DTI,
     I       DTIFG,XEFG,CEFG,MPLOCE,IBE,ISB,DCOSB,NBDYB,IBDY,MAXBES,
     i       MXTUBS,MAXNP,MAXEL,MXNPFG,MXKGL,MXNPW,MXELW,MXKBD,MXNEP,
     I       MXNCC,DL468,nwnp,npw,iwtyp,wss(1,2),mxwnp,mxwpr,
     O       NPFGS,XSFG,CSFG,MPLOCS)
c
ccc          print *,'after fgdet, npfgs=',npfgs
C
C ----- the following block is for imbedded regular fine grids
C
          CALL ADVRX
     I      (CSFG,  IRHO,CS, MXNPFG,MXNPFG,
c     I       NPFGS,IRXN,OME,EPS,EPSX,IBUG,X,IE,H,XSFG(1,1),XSFG(1,2),
     I       NPFGS,IRXN,OME,EPS,EPSX,X,IE,H,XSFG(1,1),XSFG(1,2),
     I       XSFG(1,3),MPLOCS,2,RKD,PROP,SPP,KSP,TH,THN,DTIFG,DELT,
c     I       RHOMU,DINTS)
     I       DINTS)
C
C ----- PREPARE INFORMATION FOR CHECKING THE HIGH AND LOW IN SUBELEMENTS
C       NOTE: ARRAYS MOLOC, NFGMB, AND NFGM ARE WORKING ARRAYS IN THE
C             FOLLOWING.
C
C ----- RECHECK EPCOF POINTS
C
          NPE=0
          DO 653 N=1,NPFGEP
            M=MPLOCE(N)
            IF(IE(M,11).EQ.0)GOTO 653
            NPE=NPE+1
            XEFG(NPE,1)=XEFG(N,1)
            XEFG(NPE,2)=XEFG(N,2)
            XEFG(NPE,3)=XEFG(N,3)
            DO K=1,NCC
              CEFG(NPE,K)=CEFG(N,K)
            ENDDO
            MPLOCE(NPE)=M
  653     CONTINUE
          NPFGEP=NPE
C
C ----- ADD NPFGEP TO NPFGW
C
          NPFGEP=NPE
          NPWW=NPFGW
          DO N=1,NPFGEP
            NPWW=NPWW+1
            CALL WARMSG(NPWW,MXNPFG,'HMCTRN','MXNPFG',1)
            XWFG(NPWW,1)=XEFG(N,1)
            XWFG(NPWW,2)=XEFG(N,2)
            XWFG(NPWW,3)=XEFG(N,3)
            DO K=1,NCC
              CWFG(NPWW,K)=CEFG(N,K)
            ENDDO
            MPLOCW(NPWW)=MPLOCE(N)
          ENDDO
c
ccc          print *,'the number for determining MXNPFG is',npww
C
C ---- STORE NFGMB INFORMATION AT THE PREVIOUS TIME STEP
C
          DO M=1,NEL
            NFGMBB(M)=NFGMB(M)
          ENDDO
C
          DO M=1,NEL
            NFGMB(M)=0
            NFGM(M)=0
          ENDDO
          DO 655 N=1,NPWW
            M=MPLOCW(N)
            IF(IE(M,11).EQ.0)GOTO 655
            NFGM(M)=NFGM(M)+1
  655     CONTINUE
          DO 660 M=1,NEL
            IF(M.EQ.1)THEN
              NFGMB(M)=0
              GOTO 660
            ENDIF
            NFGMB(M)=NFGMB(M-1)+NFGM(M-1)
  660     CONTINUE
c         NW=NFGMB(NEL)+NFGM(NEL)
          DO M=1,NEL
            NFGM(M)=0
          ENDDO
          DO 675 N=1,NPWW
            M=MPLOCW(N)
            IF(IE(M,11).EQ.0)GOTO 675
            NFGM(M)=NFGM(M)+1
            NM=NFGMB(M)+NFGM(M)
            MPLOC(NM)=N
  675     CONTINUE
C
C ----- DETERMINE ISE ARRAY AND HIGH-LOW, MAXFGW,MINFGW,CMAXFG,CMINFG
C
          DO I=1,MXNPFG
            MPLOCE(I)=0
          ENDDO
C
          NEFGS=0
ccc          print *,'before call isehil'
          CALL ISEHIL
     I    (MAXEL,MXKGL,MXNPFG,MXELW,MXNPW,MXEPW,MXNCC,NEL,NCC,
     I     NPFGW,EPSX,NFGM,NFGMB,MPLOCS,IE,XWFG,CWFG,
     I     MPLOC,XSFG,CSFG,NXG,NYG,NZG,
     M     NEPWN,NEPW,
     O     NPFGS,NEFGS,NLGELM,ISE,MAXFGW,MINFGW,CMAXFG,CMINFG,MPLOCE)
C
ccc          print *,'the total number of rough elements is ',nlgelm
ccc          print *,'after isehil, npfgs,nefgs=',npfgs,nefgs
C
C ----- END OF ADVECTION COMPUTATION
C
ccc          PRINT *,'THIS IS THE END OF ADVECTION COMPUTATION'
        ENDIF
C
C +++++++ FOR THE CASE OF NOT USING LAGRANGIAN APPOACH
C
      ELSEIF(LGRN.EQ.0)THEN
        NEFGS=0
        DO K=1,NCC
          DO NP=1,NNP
            CS(NP,K)=CP(NP,K)
          ENDDO
        ENDDO
        CALL THNODE(THN(1,1,1),TH,THP,PROP,RKD(1,1),WWRK,IE,X)
        CALL AFABTA
     I       (X,IE,  V,VP,  W,KSS,IOPTIM,  THN(1,2,1),PROP,
     O        WETAB)
      ENDIF
C
C +++++++ MERGE TWO CASES
C
      KDIG=KDIG+1
      KFLOW = 1
      IF(IBUG.NE.0 .AND. KPR(JTM).NE.0) WRITE(16,4800) KDIG,TIME,DELT
      IF(IBUG.NE.0 .AND. KPR(JTM).NE.0) WRITE(16,4850)
C
C %%%%%%% BEGIN THE NONLINEAR ITERATION LOOP
C
C ------- MAKE INITIAL GUESS FOR BLOCK ITERATION
C
C
C ----- PREPARE DIFFUSION INFORMATION FOR COMPOSING MATRIX EQUATIONS
C
      IF(IDZOOM.EQ.0)THEN
        NDFG=NNP
        DO NP=1,NNP
          NLRND(NP)=NLRN(NP)
        ENDDO
        DO K=1,NCC                                                       3/15/95
          DO NP=1,NNP                                                    3/15/95
            CWFG(NP,K)=CS(NP,K)                                          3/15/95
          ENDDO                                                          3/15/95
        ENDDO                                                            3/15/95
      ELSE
ccc        print *,'before call dfprep'
        CALL DFPREP
     I    (MAXNP,MAXEL,MXNPFG,MXKGL,MXKGLD,MXNPW,MXADNP,MXJBD,MXKBD,
     I     MXNCC,MXNDB,IDZOOM,NCC,NEL,NNP,NXD,NYD,NZD,EPSX,X,CS,IB,IE,
     I     NLRN,LRL,NLRL,ISE,NFGMB,XSFG,CSFG,NEFGS, MXMSV,MXLSV,IPNTSt,
     M     NFGM,XW,MWLOC,MPLOCW,XWFG,CWFG,
     O     NLS, ILSV,IMSV,
     O     MPLOC,ISED,NDBD,MPLOCD,NDFG,NEFGD,NDB,NLRND,LRN,ND)
C
ccc      print *,'after DFPREP, NDFG,NEFGD, and NDB are',NDFG,NEFGD,NDB
      ENDIF
C
C ---- PREPARE THE INITIAL GUESS FOR DIFFUSION ZOOMING ITERATION LOOP
C
      IF(NDFG.GT.NNP)THEN
        NNP1=NNP+1
        DO K=1,NCC
          DO NP=1,NNP
            DTIFG(NP,K)=DTI(NP,K)
          ENDDO
          IF(MICONF.NE.0 .AND. K.LE.3)THEN
            DO NP=NNP1,NDFG
              DTIFG(NP,K)=1.0D0/DELT
            ENDDO
          ELSE
            IDTI=1
            IBF=1
            CALL HPTRAK
     I          (MAXNP,MAXEL,MXNPFG,MXKGL,MXKBD,MXNPW,MXELW,MAXBES,
     I           MXTUBS,MXNEP,MXNEP,MXNPFG, NNP,NEL,NEFGS,NTUBS,
     I           NXW,NYW,NZW, NPFG,NNP1,IBF,IDETQ, DELT,TMAX,RAMADA,
     I           IDTI,IZOOM,EPSX,CP(1,K),X,VEAVG(1,1,K),CPFG(1,K),XPFG,
     I           MPLOC,ISE,NFGMBB,IE,IB,LRL,NLRL,IEW,IBW,LRLW,NLRLW,
     I           ISB,DCOSB,NBDYB,IBDY,DL468,nwnp,npw,iwtyp(1,k),
     i           wss(1,2),mxwnp,mxwpr,
     O           CS(1,K),DTI(1,K),DTIFG(1,K),NDFG,XWFG,CWFG(1,K),
     O           MPLOC,NPFGEP,XEFG,CEFG(1,K),MPLOCE,
     M           XW,VXW)
          ENDIF
        ENDDO
c
ccc        print *,' before call DISPC'
c
        CALL DISPC
     I    (X,IE,V,VP,TH,THP,H,PROP,W,KSS,IQUAR,XWFG,CWFG,NEFGD,cnstkr,
c     I     ISED,2,MXKGLD,8,KSP,KVIT,SPP,RHOMU,DINTS,IRHO,RHOTYP,
     I     ISED,2,MXKGLD,8,KSP,KVIT,SPP,DINTS,IRHO,RHOTYP,
     O     AKDC)
c
ccc        print *,' after call DISPC'
c
C
C ***** overwrite the Dirichlet boundary value
C
        DO I=1,NDNP
          NP=NPDB(I)
          K1=1
          IF(MICONF.NE.0)K1=4
          DO K=K1,NCC
            ITYP=IDTYP(I,K)
            CWFG(NP,K)=CDB(ITYP)
          ENDDO
        ENDDO
        DO K=1,NCC
          DO N=1,NNP
            C(N,K)=CWFG(N,K)
          ENDDO
          DO N=1,NDFG
            CPFG(N,K)=CWFG(N,K)
          ENDDO
        ENDDO
C
      ELSE
        DO I=1,NDNP
          NP=NPDB(I)
          K1=1
          IF(MICONF.NE.0)K1=4
          DO K=K1,NCC
            ITYP=IDTYP(I,K)
            CS(NP,K)=CDB(ITYP)
          ENDDO
        ENDDO
        DO K=1,NCC
          DO N=1,NNP
            C(N,K)=CS(N,K)
          ENDDO
        ENDDO
      ENDIF
C
      IF(LGRN.NE.0)THEN
        NNITER=1
      ELSE
        NNITER=NITERt
      ENDIF
      IFLUX=1
C
      DO 740 ITER=1,NNITER
C
        DO 730 K=1,NCC
C
C ------- ASSEMBLE COEFFICIENT MATRIX AND CONSTRUCT LOAD VECTOR
C
ccc          print *,'ID=',ID,' BEFORE CALL TASEMB'
          IF(IDZOOM.EQ.0 .OR. NDFG.LE.NNP)THEN
            ID=3
            CALL TASEMB
     >      (CMATRX,RLD,C,CP,CS(1,K),X,IE,LRN,NLRL,LRL, WETAB,V,VP,TH,
     >       THP,AKHC,AKDC,LES,ISTYP,SOS,NPW,IWTYP,WSS,IRHO,RHOTYP,
     >       PROP,DINTS,RKD,TRANC,DTI(1,K),W,WV,MICONF,ID,
     >       KSS,ISED,XPFG,MAXNP,MAXNP,NEL,K)
          ELSE
            ID=1
            CALL TASEMB
     >      (CMATRX,RLD,C,CP,CWFG(1,K),X,IE,LRN,NLRL,LRL, WETAB,V,VP,TH,
     >       THP,AKHC,AKDC,LES,ISTYP,SOS,NPW,IWTYP,WSS,IRHO,RHOTYP,
     >       PROP,DINTS,RKD,TRANC,DTI(1,K),W,WV,MICONF,ID,
     >       KSS,ISED,XPFG,MXNPFG,MAXNP,NEL,K)
          ENDIF
          IF(ID.NE.3)THEN
            ID=2
ccc            print *,'ID=',ID,' BEFORE CALL TASEMB'
            CALL TASEMB
     >      (CMATRX,RLD,C,CP,CWFG(1,K),X,IE,LRN,NLRL,LRL, WETAB,V,VP,TH,
     >       THP,AKHC,AKDC,LES,ISTYP,SOS,NPW,IWTYP,WSS,IRHO,RHOTYP,
     >       PROP,DINTS,RKD,TRANC,DTIFG(1,K),W,WV,MICONF,ID,
     >       KSS,ISED,XWFG,MXNPFG,MXNPFG,NEFGD,K)
           ENDIF
C
C ------- APPLY BOUNDARY CONDITIONS
C
          IF(MICONF.EQ.0 .OR. K.GT.3)THEN
            ID=1
ccc            print *,'ID=',ID,' BEFORE CALL TBC'
            CALL TBC
     >        (CMATRX,RLD,CS(1,K),X,IE,LRN,NLRND, DCOSB,ISB,
     1         V,VP, QCB,ISC,ICTYP(1,K), QNB,ISN,INTYP(1,K),
     2         CVB,ISV,IVTYP(1,K), CDB,IDTYP(1,K),NPDB, W,KSS,ID,MXADNP)
            IF(IDZOOM.NE.0 .AND. NDFG.GT.NNP)THEN
              ID=2
ccc              print *,'ID=',ID,' BEFORE CALL TBC1'
              CALL TBC1
     I      (NDB, NLRND,LRN,CPFG(1,K),IE,NDBD,MPLOCW,
     I       MPLOCD,IB,NLS,IFLUX, IMSV, ILSV, MXMSV,MXLSV,
     I       X,XWFG,EPSX, XPFG(1,1),NDFG,
     O       RLD,CMATRX)
            ENDIF
          ENDIF
C
C ------- SOLVE THE MATRIX EQUATION BY BLOCK ITERATIONS
C
c         WRITE(16,*)'THE FOLLOWING HISTORY IS FOR TRANSPORT PART'
          DO NP=1,NNP
            RI(NP)=C(NP,K)
          ENDDO
          IF(NDFG.GT.NNP)THEN
            DO NP=1,NDFG
              RI(NP)=CPFG(NP,K)
            ENDDO
          ENDIF
          IF(IPNTSt.EQ.0 .AND. IDZOOM.EQ.0) THEN
            CALL BLKITR(RL,RI,CMTRXL,RLDL, CMATRX,RLD, GNLR,
     >           LNOJCN,NNPLR, LMAXDF,OMI,EPS,NPITERT,IBUG,KPR(JTM),2)
          ELSE IF (IPNTSt.EQ.1 .OR. (IPNTSt.EQ.0 .AND. IDZOOM.NE.0))
     >           THEN
            CALL PISS
     I        (MXADNP,MXJBD,MXADNP, NDFG,NLRND,NPITERt,OMI,EPS,
     I         KPR(JTM),IBUG, CMATRX,RI,RLD,LRN,2,
     O         RL)
          ELSE IF (IPNTSt.EQ.2) THEN
            CALL PPCG
     I        (CMATRX,RLD,LRN,NLRND,IEIGEN,GG,EPS,SQEPS,IBUG,KPR(JTM),
     M         MXADNP,MXADNP,MXJBD,NDFG,2,  SK,RK,RI,PK,
     O         RL)
          ELSE IF (IPNTSt.EQ.3) THEN
            CALL ILUCG
     I       (CMATRX,RLD,LRN,ND,NLRND, EPS,SQEPS,IBUG,KPR(JTM),
     M        MXADNP,MXADNP,MXJBD,NDFG,2, SK,RK,RI,PK,AA,IL,
     O        RL)
          END IF
C
C ------- CHECK THE CONVERGENCY OF NONLINEAR LOOP
C
        IF(LGRN.EQ.0)THEN
          DIFMAX=0.0
          NOCR=1
          DO 720 NP=1,NNP
            IF(C(NP,K).EQ.0.0) GO TO 720
            DIF=(C(NP,K)-RL(NP))/C(NP,K)
            DIF=DABS(DIF)
            IF(DIF.LE.DIFMAX) GO TO 720
            DIFMAX=DIF
            KMAX=K
            NOCR=NP
  720     CONTINUE
        ENDIF
C
C ------- UPDATE NONLINEAR ITERATE CW FOR COMPUTING COEFFICIENT
C ------- MATRIX AND LOAD VECTOR.
C
        DO 725 NP=1,NNP
          C(NP,K)=OME*RL(NP)+(1.0D0-OME)*C(NP,K)
  725   CONTINUE
C
        IF(NDFG.GT.NNP)THEN
          DO NP=1,NDFG
            CPFG(NP,K)=RL(NP)
          ENDDO
        ENDIF
C
  730   CONTINUE
C
        IF(NNITER.EQ.1) GO TO 750
        IF(IBUG.NE.0.AND.KPR(JTM).NE.0) WRITE(16,5400)
     .    ITER,KMAX,DIFMAX,TOLB,NOCR
        IF(ITER.EQ.1) GO TO 740
        IF(DIFMAX.LT.TOLB) GO TO 750
C
  740 CONTINUE
C
C %%%%%%% END OF NONLINEAR LOOP
C
      WRITE(16,7500) JTM,ITER,NNITER,KMAX,DIFMAX,TOLB,NOCR
C
  750 CONTINUE
C
      IF(IMID.EQ.0) GO TO 830
C
      DO 810 K=1,NCC
      DO 810 NP=1,NNP
        C(NP,K)=2.0D0*C(NP,K)-CP(NP,K)
  810 CONTINUE
      DO 820 NPP=1,NDNP
        NP=NPDB(NPP)
        DO 820 K=1,NCC
          ITYP=IDTYP(NPP,K)
          C(NP,K)=CDB(ITYP)
  820 CONTINUE
C
  830 CONTINUE
C
      DO 850 K=1,NCC
C
C ------- CALCULATE MATERIAL FLUX FX AND FZ
C
        CALL FLUX(F(1,1,K),CMATRX,C(1,K),X,IE,V,AKHC,IQUAR)
C
C ------- DETERMINE FLUX THROUGH ALL BOUNDARIES
C
        CALL TSFLOW
     >     (BFLX(1,1,K),X,IE,C(1,K),F(1,1,K),TH,RKD(1,K),TRANC(1,K),
     >      DCOSB,ISB,NPBB, SOS,ISTYP(1,K),LES, WSS,IWTYP(1,K),
     >      NPW,NPVB,NPDB,NPCB,NPNB,PROP,DELT,KFLOW,K)
C
C ------- PRINT VARIABLES AT CURRENT TIME STEP
C
        KDIAG=0
        CALL TPRINT(C(1,K),F(1,1,K), TIME,DELT, KPR(JTM),KOUT,KDIAG,
     >              JTM,K)
C
  850 CONTINUE
C
      NPFG=NPFGS
C
C ----- COMPUTE CONCENTRATION OF FINE GRID POINTS BY INTERPOLATION
C
      DO M=1,NEL
        IE(M,10)=IE(M,11)
      ENDDO
      IF(IDZOOM.NE.0)THEN
        DO 900 NP=1,NPFGS
          MP=MPLOCS(NP)
          XP=XSFG(NP,1)
          YP=XSFG(NP,2)
          ZP=XSFG(NP,3)
          CALL ELENOD
     I        (IE(MP,5),IE(MP,7),
     O         NODE,K,K)
C
          CALL INTERP
     I        (MAXNP,MAXEL,MXNPFG,MXKGLD,9,MXNCC,NODE,MP,XP,YP,ZP,X,C,
     I         IE,XWFG,CPFG,NEL,NEFGD,IDZOOM,10,ISED,NFGM,1,NCC,EPSX,
     O         CQP)
C
          CALL INTERP
     I        (MAXNP,MAXEL,MXNPFG,MXKGLD,9,MXNCC,NODE,MP,XP,YP,ZP,X,CS,
     I         IE,XWFG,CWFG,NEL,NEFGD,IDZOOM,10,ISED,NFGM,1,NCC,EPSX,
     O         CQW)
C
          DO K=1,NCC
            CSFG(NP,K)=CSFG(NP,K)+CQP(K)-CQW(K)
          ENDDO
  900   CONTINUE
      ELSEIF(IZOOM.EQ.0)THEN
        DO 910 NP=1,NPFGS
          MP=MPLOCS(NP)
          XP=XSFG(NP,1)
          YP=XSFG(NP,2)
          ZP=XSFG(NP,3)
          CALL ELENOD
     I        (IE(MP,5),IE(MP,7),
     O         NODE,I,I)
          DO I=1,NODE
            IEMI=IE(MP,I)
            XX(I)=X(IEMI,1)
            YY(I)=X(IEMI,2)
            ZZ(I)=X(IEMI,3)
            DO K=1,NCC
              CC(I,K)=C(IEMI,K)-CS(IEMI,K)
            ENDDO
          ENDDO
          CALL BASE
c     I        (XX,YY,ZZ,XP,YP,ZP,M,NODE,1,
     I        (XX,YY,ZZ,XP,YP,ZP,NODE,1,
     O         DL,DNX,DNX,DNX)
          DO K=1,NCC
            DO IQ=1,NODE
              CSFG(NP,K)=CSFG(NP,K)+CC(IQ,K)*DL(IQ)
            ENDDO
          ENDDO
  910   CONTINUE
      ENDIF
C
      IF(KDSK(JTM).EQ.1) CALL TSTORE(X,IE,C,F,TITLE,NPROB,JTM,TIME)
C
      RETURN
C
 4000 FORMAT('0'/1X,'*** ITM = ',I4,' NTAU = ',I3,'  DTAU = ',1PD12.4)
 4800 FORMAT('1'//1X,'***********************************************',
     1 '******************************'///' DIAGNOSTIC TABLE',I4,
     2 '.. AT TIME =',1PD12.4,', (DELT =', 1PD12.4,')')
 4850 FORMAT(///' TABLE OF ITERATIVE PARAMETERS'// 6X,
     1 'ITERATION',7X,' MAX DIF',6X,'TOLERANCE',6X,
     > ' MAX OCCURANCE NODE')
 5400 FORMAT(5X,I10,3X,I2,2X,E12.4,3X,E12.4,15X,I10)
 5500 FORMAT('0'/5X,'** WARNING: NO CONVERGENCE AFTER',I4,' ITERATIONS',
     1       /8X,'NITER =',I4,' KMAX=',I2,' DIFMAX =',D12.4,
     >       ' TOLA =',D12.4,' NOCCUR =',I4)
 7500 FORMAT('0'/5X,'** WARNING: NO CONVERGENCE AT',I4,'-TH TIME STEP AF
     >      TER',I4,' ITERATIONS'/8X,'NITER =',I4,' KMAX=',I2,
     >      ' DIFMAX =',D12.4,' TOLB =',D12.4,' NOCCUR =',I4)
C
      END
C
C
C
      SUBROUTINE DISPC
     I   (X,IE,V,VP,TH,THP,H,PROP,W,KSS,IQUAR,XWFG,CWFG,NEL,cnstkr,ISED,
c     I    ID,MAXNEL,NO3,KSP,KVIT,SPP,RHOMU,DINTS,IRHO,RHOTYP,
     I    ID,MAXNEL,NO3,KSP,KVIT,SPP,DINTS,IRHO,RHOTYP,
     O    AKDC)
C
C $$$$$ TO COMPUTE THE DISPERSION COEFFICENTS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 N(8),L1(3),L2(3),L3(3),LL1(4),LL2(4),LL3(4),LL4(4)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /SAZFM/ MXNPFG,MXKGL,MXKGLD,MXNEP,MXEPW,MXNPW,MXELW,MXNDB,
     >               mxnpww,mxelww
      COMMON /CHEM/ MXNCC,NCC
C
      DIMENSION AKDC(8,MAXNEL,NO3),X(MAXNP,3)
      DIMENSION IE(MAXEL,11),XWFG(MXNPFG,3),CWFG(MXNPFG,MXNCC)
      DIMENSION V(MAXNP,3),VP(MAXNP,3),ISED(MXKGLD,9)
      DIMENSION TH(MAXEL,8),THP(MAXEL,8),H(MAXNP)
c      DIMENSION PROP(MAXMAT,MXMPPM),RHOMU(MXNCC),DINTS(MXNCC)
      DIMENSION PROP(MAXMAT,MXMPPM),DINTS(MXNCC)
      DIMENSION SPP(MXSPPM,MAXMAT,4),RHOTYP(MAXMAT)
      DIMENSION XQ(8),YQ(8),ZQ(8),DNX(8)
      DIMENSION VXQ(8),VYQ(8),VZQ(8),VXG(8),VYG(8),VZG(8),THG(8)
      DIMENSION S(8),T(8),U(8),HQ(8),THH(8)
      DIMENSION XX(8,3),VXX(8,3),VXXP(8,3),CKQ(8,7),CKG(7)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0, -1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA T/-1.0D0,-1.0D0,1.0D0,1.0D0, -1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA U/-1.0D0,-1.0D0,-1.0D0,-1.0D0, 1.0D0,1.0D0,1.0D0,1.0D0/
C
      IF(IQUAR.EQ.1 .OR. IQUAR.EQ.3)THEN
        P=0.577350269189626D0
        L1(1)=0.666666666666667D0
        L1(2)=0.166666666666667D0
        L1(3)=0.166666666666667D0
        L2(1)=0.166666666666667D0
        L2(2)=0.666666666666667D0
        L2(3)=0.166666666666667D0
        L3(1)=0.166666666666667D0
        L3(2)=0.166666666666667D0
        L3(3)=0.666666666666667D0
        LL1(1)=0.58541020D0
        LL1(2)=0.13819660D0
        LL1(3)=0.13819660D0
        LL1(4)=0.13819660D0
        LL2(1)=0.13819660D0
        LL2(2)=0.58541020D0
        LL2(3)=0.13819660D0
        LL2(4)=0.13819660D0
        LL3(1)=0.13819660D0
        LL3(2)=0.13819660D0
        LL3(3)=0.58541020D0
        LL3(4)=0.13819660D0
        LL4(1)=0.13819660D0
        LL4(2)=0.13819660D0
        LL4(3)=0.13819660D0
        LL4(4)=0.58541020D0
      ELSE
        P=1.0D0
        L1(1)=1.0D0
        L1(2)=0.0D0
        L1(3)=0.0D0
        L2(1)=0.0D0
        L2(2)=1.0D0
        L2(3)=0.0D0
        L3(1)=0.0D0
        L3(2)=0.0D0
        L3(3)=1.0D0
        LL1(1)=1.0D0
        LL1(2)=0.0D0
        LL1(3)=0.0D0
        LL1(4)=0.0D0
        LL2(1)=0.0D0
        LL2(2)=1.0D0
        LL2(3)=0.0D0
        LL2(4)=0.0D0
        LL3(1)=0.0D0
        LL3(2)=0.0D0
        LL3(3)=1.0D0
        LL3(4)=0.0D0
        LL4(1)=0.0D0
        LL4(2)=0.0D0
        LL4(3)=0.0D0
        LL4(4)=1.0D0
      ENDIF
C
      W1=W
      W2=1.0D0-W
      IF(KSS.NE.0) GO TO 100
      W1=1.0D0
      W2=0.0
  100 CONTINUE
C
C ******* COMPUTE DISPERSION COEFFICIENT
C
      DO 690 M=1,NEL
        IF(M.LE.0)GOTO 690
        IF(ID.EQ.1)THEN
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NODE,IQ,IQ)
          MTYP=IE(M,9)
C
          DO 510 IQ=1,NODE
            NP=IE(M,IQ)
            XQ(IQ)=X(NP,1)
            YQ(IQ)=X(NP,2)
            ZQ(IQ)=X(NP,3)
            VXQ(IQ)=V(NP,1)*W1+VP(NP,1)*W2
            VYQ(IQ)=V(NP,2)*W1+VP(NP,2)*W2
            VZQ(IQ)=V(NP,3)*W1+VP(NP,3)*W2
  510     CONTINUE
        ELSEIF(ID.EQ.2)THEN
          MM=ISED(M,9)
          CALL ELENOD
     I        (IE(MM,5),IE(MM,7),
     O         NODE,IQ,IQ)
          MTYP=IE(MM,9)
          IF(IRHO.EQ.1)RHO=RHOTYP(MTYP)
          DO IQ=1,NODE
            IEM=IE(MM,IQ)
            DO J=1,3
              XX(IQ,J)=X(IEM,J)
              VXX(IQ,J)=V(IEM,J)
              VXXP(IQ,J)=VP(IEM,J)
            ENDDO
            IF(IRHO.EQ.1)HQ(IQ)=H(IEM)
c           HPQ(IQ)=HP(IEM)
          ENDDO
          CALL ELENOD
     I        (ISED(M,5),ISED(M,7),
     O         NQ,IQ,IQ)
          DO 520 IQ=1,NQ
            NI=ISED(M,IQ)
            XQ(IQ)=XWFG(NI,1)
            YQ(IQ)=XWFG(NI,2)
            ZQ(IQ)=XWFG(NI,3)
            DO K=1,NCC
              CKQ(IQ,K)=CWFG(NI,K)
            ENDDO
            CALL BASE
c     I         (XX(1,1),XX(1,2),XX(1,3),XQ(IQ),YQ(IQ),ZQ(IQ),MM,NODE,1,
     I         (XX(1,1),XX(1,2),XX(1,3),XQ(IQ),YQ(IQ),ZQ(IQ),NODE,1,
     O          N,DNX,DNX,DNX)
            IF(KVIT.GT.0)THEN
              THH(IQ)=0.0D0
c             THHP(IQ)=0.0D0
              DO I=1,NODE
                THH(IQ)=THH(IQ)+N(I)*HQ(I)
c               THHP(IQ)=THHP(IQ)+N(IQ)*HPQ(I)
              ENDDO
            ELSE
              THH(IQ)=0.0D0
c             THHP(IQ)=0.0D0
              DO I=1,NODE
                THH(IQ)=THH(IQ)+N(I)*TH(MM,I)
c               THHP(IQ)=THHP(IQ)+N(I)*THP(MM,I)
              ENDDO
            ENDIF
            VXQ(IQ)=0.0D0
            VYQ(IQ)=0.0D0
            VZQ(IQ)=0.0D0
            DO I=1,NODE
              VXQ(IQ)=VXQ(IQ)+N(I)*(W1*VXX(I,1)+W2*VXXP(I,1))
              VYQ(IQ)=VYQ(IQ)+N(I)*(W1*VXX(I,2)+W2*VXXP(I,2))
              VZQ(IQ)=VZQ(IQ)+N(I)*(W1*VXX(I,3)+W2*VXXP(I,3))
            ENDDO
  520     CONTINUE
        ENDIF
C
        AL=PROP(MTYP,2)
        AT=PROP(MTYP,3)
        DD=PROP(MTYP,4)*PROP(MTYP,5)
C
C ------- EVALUATE VELOCITY COMPONENT AT GAUSSIAN POINTS
C
        DO 650 KG=1,NODE
          IF(NODE.EQ.8)THEN
            SS=P*S(KG)
            TT=P*T(KG)
            UU=P*U(KG)
            SM=1.0D0-SS
            SP=1.0D0+SS
            TM=1.0D0-TT
            TP=1.0D0+TT
            UM=1.0D0-UU
            UP=1.0D0+UU
            N(1)=0.125D0*SM*TM*UM
            N(2)=0.125D0*SP*TM*UM
            N(3)=0.125D0*SP*TP*UM
            N(4)=0.125D0*SM*TP*UM
            N(5)=0.125D0*SM*TM*UP
            N(6)=0.125D0*SP*TM*UP
            N(7)=0.125D0*SP*TP*UP
            N(8)=0.125D0*SM*TP*UP
          ELSEIF(NODE.EQ.6)THEN
            XSI=-P
            ND=KG
            IF(ND.GT.3)THEN
              XSI=P
              ND=KG-3
            ENDIF
            DL1=L1(ND)
            DL2=L2(ND)
            DL3=L3(ND)
            SM=1.0D0-XSI
            SP=1.0D0+XSI
            N(1)=0.5D0*SM*DL1
            N(2)=0.5D0*SM*DL2
            N(3)=0.5D0*SM*DL3
            N(4)=0.5D0*SP*DL1
            N(5)=0.5D0*SP*DL2
            N(6)=0.5D0*SP*DL3
          ELSE
            N(1)=LL1(KG)
            N(2)=LL2(KG)
            N(3)=LL3(KG)
            N(4)=LL4(KG)
          ENDIF
C
          VXG(KG)=0.0
          VYG(KG)=0.0
          VZG(KG)=0.0
          DO 620 IQ=1,NODE
            VXG(KG)=VXG(KG)+VXQ(IQ)*N(IQ)
            VYG(KG)=VYG(KG)+VYQ(IQ)*N(IQ)
            VZG(KG)=VZG(KG)+VZQ(IQ)*N(IQ)
  620     CONTINUE
C
          IF(ID.EQ.1)THEN
            THG(KG)=W1*TH(M,KG)+W2*THP(M,KG)
          ELSEIF(ID.EQ.2)THEN
            IF(KVIT.GT.0)THEN
              HNP=0.0D0
              DO K=1,NCC
                CKG(K)=0.0D0
              ENDDO
              DO I=1,NQ
                HNP=HNP+N(I)*THH(I)
                DO K=1,NCC
                  CKG(K)=CKG(K)+N(I)*CKQ(I,K)
                ENDDO
              ENDDO
              CALL SPFUNC
     I        (HNP,SPP,MTYP,KSP,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
c     I         CKG,1,IRHO,RHOMU,DINTS,RHO,cnstkr,
     I         CKG,1,IRHO,DINTS,RHO,cnstkr,
     O         THG(KG),0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     O         AKDC(KG,M,7))
              AKDC(KG,M,8)=THG(KG)
C             HNP=0.0D0
C             DO I=1,NQ
C               HNP=HNP+N(I)*THHP(I)
C             ENDDO
C    I        (HNP,SPP,MTYP,KSP,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
C    I         CKG(1,KG),1,RHOMU,DINTS,RHO,
C    O         AKHC(KG,M,9),0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
C    O         AKDC(KG,M,7))
            ELSE
              THG(KG)=0.0D0
              AKDC(KG,M,7)=1.0D0
              DO I=1,NQ
                THG(KG)=THG(KG)+N(I)*THH(I)
c               AKDC(KG,M,9)=AKDC(KG,M,9)+N(I)*THHP(I)
              ENDDO
              AKDC(KG,M,8)=THG(KG)
            ENDIF
          ENDIF
C
          VG=DSQRT(VXG(KG)*VXG(KG)+VYG(KG)*VYG(KG)+VZG(KG)*VZG(KG))
          IF(VG.NE.0.0)THEN
            VGI=1.0D0/VG
            AKDC(KG,M,1)=(AL*VXG(KG)*VXG(KG)+AT*(VYG(KG)*VYG(KG)+
     1                    VZG(KG)*VZG(KG)))*VGI+DD*THG(KG)
            AKDC(KG,M,2)=(AL*VYG(KG)*VYG(KG)+AT*(VZG(KG)*VZG(KG)+
     1                    VXG(KG)*VXG(KG)))*VGI+DD*THG(KG)
            AKDC(KG,M,3)=(AL*VZG(KG)*VZG(KG)+AT*(VXG(KG)*VXG(KG)+
     1                    VYG(KG)*VYG(KG)))*VGI+DD*THG(KG)
            AKDC(KG,M,4)=(AL-AT)*VXG(KG)*VYG(KG)*VGI
            AKDC(KG,M,5)=(AL-AT)*VXG(KG)*VZG(KG)*VGI
            AKDC(KG,M,6)=(AL-AT)*VYG(KG)*VZG(KG)*VGI
          ELSE
            AKDC(KG,M,1)=DD*THG(KG)
            AKDC(KG,M,2)=DD*THG(KG)
            AKDC(KG,M,3)=DD*THG(KG)
            AKDC(KG,M,4)=0.0
            AKDC(KG,M,5)=0.0
            AKDC(KG,M,6)=0.0
          ENDIF
C
  650   CONTINUE
  690 CONTINUE
C
      RETURN
      END
C
C
c
      SUBROUTINE AFABTA
     I    (X,IE,  V,VP,  W,KSS,IOPTIM,  THN,PROP,
     O     WETAB)
C
C ------- TO CALCULATE THE WEIGHTING FACTORS ON ALL THE SIDES OF EACH
C ------- ELEMENT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
C
      DIMENSION WETAB(12,MAXEL)
      DIMENSION X(MAXNP,3),IE(MAXEL,11)
      DIMENSION V(MAXNP,3),VP(MAXNP,3)
      DIMENSION THN(MAXNP)
      DIMENSION PROP(MAXMAT,MXMPPM)
C
      DIMENSION NODE(8)
C
      W1=W
      W2=1.0D0-W
      IF(KSS.NE.0)GOTO 100
      W1=1.0D0
      W2=0.0
  100 CONTINUE
C
      DO 490 M=1,NEL
        MTYP=IE(M,9)
        AL=PROP(MTYP,2)
        AT=PROP(MTYP,3)
        DD=PROP(MTYP,4)*PROP(MTYP,5)
        IF(IE(M,5).EQ.0)THEN
          NNQ=4
          NQ=6
        ELSEIF(IE(M,7).EQ.0)THEN
          NNQ=6
          NQ=9
        ELSE
          NNQ=8
          NQ=12
        ENDIF
C
        DO 210 IQ=1,NNQ
          NODE(IQ)=IE(M,IQ)
  210   CONTINUE
C
        DO 390 I=1,NQ
C
          IF(NQ.EQ.6)THEN
            GOTO(301,302,303,304,305,306),I
          ELSEIF(NQ.EQ.9)THEN
            GOTO(311,312,313,314,315,316,317,318,319),I
          ELSE
            GOTO(321,322,323,324,325,326,327,328,329,330,331,332),I
          ENDIF
C
  301     N1=NODE(1)
          N2=NODE(2)
          GOTO 350
C
  302     N1=NODE(1)
          N2=NODE(3)
          GOTO 350
C
  303     N1=NODE(1)
          N2=NODE(4)
          GOTO 350
C
  304     N1=NODE(2)
          N2=NODE(3)
          GOTO 350
C
  305     N1=NODE(2)
          N2=NODE(4)
          GOTO 350
C
  306     N1=NODE(3)
          N2=NODE(4)
          GOTO 350
C
  311     N1=NODE(1)
          N2=NODE(2)
          GOTO 350
C
  312     N1=NODE(2)
          N2=NODE(3)
          GOTO 350
C
  313     N1=NODE(1)
          N2=NODE(3)
          GOTO 350
C
  314     N1=NODE(4)
          N2=NODE(5)
          GOTO 350
C
  315     N1=NODE(5)
          N2=NODE(6)
          GOTO 350
C
  316     N1=NODE(4)
          N2=NODE(6)
          GOTO 350
C
  317     N1=NODE(1)
          N2=NODE(4)
          GOTO 350
C
  318     N1=NODE(2)
          N2=NODE(5)
          GOTO 350
C
  319     N1=NODE(3)
          N2=NODE(6)
          GOTO 350
C
  321     N1=NODE(1)
          N2=NODE(2)
          GO TO 350
C
  322     N1=NODE(4)
          N2=NODE(3)
          GO TO 350
C
  323     N1=NODE(5)
          N2=NODE(6)
          GO TO 350
C
  324     N1=NODE(8)
          N2=NODE(7)
          GO TO 350
C
  325     N1=NODE(1)
          N2=NODE(4)
          GO TO 350
C
  326     N1=NODE(2)
          N2=NODE(3)
          GO TO 350
C
  327     N1=NODE(5)
          N2=NODE(8)
          GO TO 350
C
  328     N1=NODE(6)
          N2=NODE(7)
          GO TO 350
C
  329     N1=NODE(1)
          N2=NODE(5)
          GO TO 350
C
  330     N1=NODE(2)
          N2=NODE(6)
          GO TO 350
C
  331     N1=NODE(3)
          N2=NODE(7)
          GO TO 350
C
  332     N1=NODE(4)
          N2=NODE(8)
C
  350     DISTX=X(N2,1)-X(N1,1)
          DISTY=X(N2,2)-X(N1,2)
          DISTZ=X(N2,3)-X(N1,3)
          DIST=DSQRT(DISTX*DISTX+DISTY*DISTY+DISTZ*DISTZ)
          DCSX=DISTX/DIST
          DCSY=DISTY/DIST
          DCSZ=DISTZ/DIST
          VXX=0.5D0*((V(N1,1)+V(N2,1))*W1+(VP(N1,1)+VP(N2,1))*W2)
          VYY=0.5D0*((V(N1,2)+V(N2,2))*W1+(VP(N1,2)+VP(N2,2))*W2)
          VZZ=0.5D0*((V(N1,3)+V(N2,3))*W1+(VP(N1,3)+VP(N2,3))*W2)
          THE=0.5D0*(THN(N1)+THN(N2))
          VV=DSQRT(VXX*VXX+VYY*VYY+VZZ*VZZ)
          IF(VV.EQ.0.0)GOTO 370
          VAL=VXX*DCSX+VYY*DCSY+VZZ*DCSZ
          VEL=DIST*VAL
          VVI=1.0D0/VV
C
          IF(VVI.NE.0.0) THEN
            DXX=(AL*VXX*VXX+AT*(VYY*VYY+VZZ*VZZ))*VVI + DD*THE
            DYY=(AL*VYY*VYY+AT*(VZZ*VZZ+VXX*VXX))*VVI + DD*THE
            DZZ=(AL*VZZ*VZZ+AT*(VXX*VXX+VYY*VYY))*VVI + DD*THE
            DXY=(AL-AT)*VXX*VYY*VVI
            DXZ=(AL-AT)*VXX*VZZ*VVI
            DYZ=(AL-AT)*VYY*VZZ*VVI
          ELSE
            DXX=DD*THE
            DYY=DD*THE
            DZZ=DD*THE
            DXY=0.0
            DXZ=0.0
            DYZ=0.0
          ENDIF
C
          DAL=DABS(DCSX*(DCSX*DXX+DCSY*DXY+DCSZ*DXZ)+DCSY*(DCSX*DXY+
     1        DCSY*DYY+DCSZ*DYZ)+DCSZ*(DCSX*DXZ+DCSY*DYZ+DCSZ*DZZ))
          DISP=2.0D0*DAL
          IF(IOPTIM.EQ.0)GOTO 380
          IF(DISP.EQ.0.0 .OR. VEL.EQ.0.0)GOTO 360
          WETAB(I,M)=1.0D0/DTANH(VEL/DISP)- DISP/VEL
          GOTO 390
C
  360     IF(VEL.EQ.0.0) WETAB(I,M)=0.0
          IF(VEL.GT.0.0) WETAB(I,M)=1.0D0
          IF(VEL.LT.0.0) WETAB(I,M)=-1.0D0
          GOTO 390
C
  370     WETAB(I,M)=0.0
          GOTO 390
C
  380     WETAB(I,M)=1.0D0
          IF(VEL.LT.0.0) WETAB(I,M)=-1.0D0
          IF(VEL.EQ.0.0) WETAB(I,M)=0.0
          IF(DISP.EQ.0.0)GOTO 390
          IF(DABS(VEL/DISP).LT.1.0D-10) WETAB(I,M)=0.0
  390   CONTINUE
  490 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE THNODE(THN,TH,THP,PROP,RKD,WWRK,IE,X)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO MOISTURE CONTENT AT NODAL POINTS.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: TH(MAXEL,8),THP(MAXEL,8),PROP(MAXMAT,MXPPM)
C -------        IE(MAXNP,11),X(MAXNP,3)
C
C ------- OUTPUT: THN(MAXNP,2)
C
C ------- WORKING ARRAY: WWRK(MAXNP)
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /TREAL/ W,WV,OME,OMI,TOLA,TOLB
C
      DIMENSION THN(MAXNP,2),TH(MAXEL,8),THP(MAXEL,8),
     >          PROP(MAXMAT,MXMPPM)
      DIMENSION WWRK(MAXNP),IE(MAXEL,11),X(MAXNP,3),RKD(MAXMAT)
      DIMENSION JV(4),XX(4),YY(4),ZZ(4),A(4)
C
      W1=W
      W2=1.0D0-W
C
      DO 110 NP=1,NNP
        WWRK(NP)=0.0
        THN(NP,1)=0.0
        THN(NP,2)=0.0
  110 CONTINUE
C
      DO 290 M=1,NEL
        MTYP=IE(M,9)
        RHOB=PROP(MTYP,1)
        RRKD=RKD(MTYP)
C
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NQ,NP1,NP1)
C
        NP1=IE(M,1)
        NP2=IE(M,2)
        NP3=IE(M,3)
        NP4=IE(M,4)
        IF(NQ.NE.4)THEN
          NP5=IE(M,5)
          NP6=IE(M,6)
          IF(NQ.EQ.8)THEN
            NP7=IE(M,7)
            NP8=IE(M,8)
          ENDIF
        ENDIF
C
C ----- COMPUTE THE VOLUME OF THE ELEMENT
C
        IF(NQ.EQ.8)THEN
          NV=6
        ELSEIF(NQ.EQ.6)THEN
          NV=3
        ELSE
          NV=1
        ENDIF
        VOL=0.0D0
        DO 120 IV=1,NV
          IF(NQ.EQ.4)THEN
            JV(1)=NP1
            JV(2)=NP2
            JV(3)=NP3
            JV(4)=NP4
          ELSEIF(NQ.EQ.6)THEN
            IF(IV.EQ.1)THEN
              JV(1)=NP4
              JV(2)=NP1
              JV(3)=NP2
              JV(4)=NP3
            ELSEIF(IV.EQ.2)THEN
              JV(1)=NP5
              JV(2)=NP4
              JV(3)=NP2
              JV(4)=NP3
            ELSE
              JV(1)=NP6
              JV(2)=NP4
              JV(3)=NP5
              JV(4)=NP3
            ENDIF
          ELSE
            IF(IV.EQ.1)THEN
              JV(1)=NP5
              JV(2)=NP1
              JV(3)=NP2
              JV(4)=NP4
            ELSEIF(IV.EQ.2)THEN
              JV(1)=NP6
              JV(2)=NP5
              JV(3)=NP2
              JV(4)=NP4
            ELSEIF(IV.EQ.3)THEN
              JV(1)=NP8
              JV(2)=NP5
              JV(3)=NP6
              JV(4)=NP4
            ELSEIF(IV.EQ.4)THEN
              JV(1)=NP7
              JV(2)=NP2
              JV(3)=NP3
              JV(4)=NP4
            ELSEIF(IV.EQ.5)THEN
              JV(1)=NP8
              JV(2)=NP6
              JV(3)=NP7
              JV(4)=NP4
            ELSE
              JV(1)=NP6
              JV(2)=NP7
              JV(3)=NP4
              JV(4)=NP2
            ENDIF
          ENDIF
          DO KK=1,4
            JVKK=JV(KK)
            XX(KK)=X(JVKK,1)
            YY(KK)=X(JVKK,2)
            ZZ(KK)=X(JVKK,3)
          ENDDO
C
          DO KK=1,4
            IF(KK.EQ.1)THEN
              K1=2
              K2=3
              K3=4
            ELSEIF(KK.EQ.2)THEN
              K1=1
              K2=3
              K3=4
            ELSEIF(KK.EQ.3)THEN
              K1=1
              K2=2
              K3=4
            ELSE
              K1=1
              K2=2
              K3=3
            ENDIF
            A(KK)=(-1.0D0)**(KK+1)*(XX(K1)*YY(K2)*ZZ(K3)+
     1            YY(K1)*ZZ(K2)*XX(K3)+ZZ(K1)*XX(K2)*YY(K3)-
     2            XX(K3)*YY(K2)*ZZ(K1)-YY(K3)*ZZ(K2)*XX(K1)-
     3            ZZ(K3)*XX(K2)*YY(K1))
          ENDDO
          VOL6=0.0D0
          DO KK=1,4
            VOL6=VOL6+A(KK)
          ENDDO
          VOL=VOL+DABS(VOL6)/6.0D0
  120   CONTINUE
C
        DO 250 IQ=1,NQ
          NP=IE(M,IQ)
          THN(NP,1)=THN(NP,1)+(W1*TH(M,IQ)+W2*THP(M,IQ)+RHOB*RRKD)*VOL
          THN(NP,2)=THN(NP,2)+(W1*TH(M,IQ)+W2*THP(M,IQ))*VOL
          WWRK(NP)=WWRK(NP)+VOL
  250   CONTINUE
  290 CONTINUE
C
      DO 390 NP=1,NNP
      THN(NP,1)=THN(NP,1)/WWRK(NP)
      THN(NP,2)=THN(NP,2)/WWRK(NP)
  390 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE VOLUME
     I    (MAXNP,MAXEL,X,Y,Z,IE,M,
     O     VOL)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION X(MAXNP),Y(MAXNP),Z(MAXNP),IE(MAXEL,11)
      DIMENSION KVB(4,6,3),XX(4),YY(4),ZZ(4),A(4)
C
      DATA KVB /5,1,2,4, 6,5,2,4, 8,5,6,4, 7,2,3,4, 8,6,7,4, 6,7,4,2,
     >          4,1,2,3, 5,4,2,3, 6,4,5,3, 0,0,0,0, 0,0,0,0, 0,0,0,0,
     >          4,1,2,3, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0/
C
      IF(IE(M,5).EQ.0)THEN
        ID=3
        NV=1
      ELSEIF(IE(M,7).EQ.0)THEN
        ID=2
        NV=3
      ELSE
        ID=1
        NV=6
      ENDIF
C
      VOL=0.0D0
      DO 120 IV=1,NV
        DO I=1,4
          II=KVB(I,IV,ID)
          JVKK=IE(M,II)
          XX(I)=X(JVKK)
          YY(I)=Y(JVKK)
          ZZ(I)=Z(JVKK)
        ENDDO
C
        DO KK=1,4
          IF(KK.EQ.1)THEN
            K1=2
            K2=3
            K3=4
          ELSEIF(KK.EQ.2)THEN
            K1=1
            K2=3
            K3=4
          ELSEIF(KK.EQ.3)THEN
            K1=1
            K2=2
            K3=4
          ELSE
            K1=1
            K2=2
            K3=3
          ENDIF
          A(KK)=(-1.0D0)**(KK+1)*(XX(K1)*YY(K2)*ZZ(K3)+
     1          YY(K1)*ZZ(K2)*XX(K3)+ZZ(K1)*XX(K2)*YY(K3)-
     2          XX(K3)*YY(K2)*ZZ(K1)-YY(K3)*ZZ(K2)*XX(K1)-
     3          ZZ(K3)*XX(K2)*YY(K1))
        ENDDO
        VOL6=0.0D0
        DO KK=1,4
          VOL6=VOL6+A(KK)
        ENDDO
        VOL=VOL+DABS(VOL6)/6.0D0
  120 CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE FLUX(F,CMATRX,C,X,IE,V,AKHC,IQUAR)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE THE MATERIAL FLUXES.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: C(NNP), X(NNP,3), IE(NEL,9), V(NNP,3), AKHC(8,MAXRL,7).
C
C ------- OUTPUT: F(NNP,3)
C
C ------- WORKING ARRAYS: CMATRX(NNP,1)
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /WETX/ APHA1,APHA2,APHA3,APHA4
      COMMON /WETY/ BETA1,BETA2,BETA3,BETA4
      COMMON /WETZ/ GAMA1,GAMA2,GAMA3,GAMA4
C
      DIMENSION F(MAXNP,3)
      DIMENSION CMATRX(MXADNP,MXJBD),C(MAXNP)
      DIMENSION X(MAXNP,3),IE(MAXEL,11)
      DIMENSION V(MAXNP,3), AKHC(8,MAXEL,7)
      DIMENSION QB(8,8),QRX(8),QRY(8),QRZ(8)
      DIMENSION XQ(8),YQ(8),ZQ(8),CQ(8)
      DIMENSION AKXYZQ(8,6)
C
C ------- INITIALIZE THE FLUX FX(NP), FY(NP), AND FZ(NP).
C
      DO 100 NP=1,NNP
      F(NP,1)=0.0
      F(NP,2)=0.0
      F(NP,3)=0.0
  100 CONTINUE
C
C ------- COMPUTE THE FLUX COMPONENTS BY APPLYING THE FINITE ELEMENT
C ------- METHOD TO THE DISPERSION TERMS
C
      DO 160 NP=1,NNP
  160 CMATRX(NP,1)=0.0
C
C ------- COMPUTE THE ELEMENT MATRIX QB AND QRX, QRY, & QRZ
C
      DO 290 M=1,NEL
C
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,IQ,IQ)
C
        APHA1=0.0
        APHA2=0.0
        APHA3=0.0
        APHA4=0.0
        BETA1=0.0
        BETA2=0.0
        BETA3=0.0
        BETA4=0.0
        GAMA1=0.0
        GAMA2=0.0
        GAMA3=0.0
        GAMA4=0.0
C
        DO 210 IQ=1,NODE
          NP=IE(M,IQ)
          XQ(IQ)=X(NP,1)
          YQ(IQ)=X(NP,2)
          ZQ(IQ)=X(NP,3)
          CQ(IQ)=C(NP)
          DO 205 I=1,6
            AKXYZQ(IQ,I)=AKHC(IQ,M,I)
  205     CONTINUE
  210   CONTINUE
C
        CALL TQ468DV(QB,QRX,QRY,QRZ,CQ,XQ,YQ,ZQ,AKXYZQ,NODE,IQUAR)
C
C ------- ASSEMBLE QB(IQ,JQ) INTO THE GLOBAL MATRIX C(NP,1) AND
C ------- FORM THE LOAD VECTOR FX(NP), FY(NY), AND FZ(NP)
C
        DO 280 IQ=1,NODE
          NI=IE(M,IQ)
          DO 240 JQ=1,NODE
            CMATRX(NI,1)=CMATRX(NI,1)+QB(IQ,JQ)
  240     CONTINUE
          F(NI,1)=F(NI,1)+QRX(IQ)
          F(NI,2)=F(NI,2)+QRY(IQ)
          F(NI,3)=F(NI,3)+QRZ(IQ)
  280   CONTINUE
  290 CONTINUE
C
C ------- SOLVE THE MATRIX EQUATION CX=B
C
      DO 370 NP=1,NNP
        F(NP,1)=F(NP,1)/CMATRX(NP,1)
        F(NP,2)=F(NP,2)/CMATRX(NP,1)
        F(NP,3)=F(NP,3)/CMATRX(NP,1)
  370 CONTINUE
C
C ------ ADD THE ADVECTION FLUX TO DISPERSION FLUX
C
      DO 380 NP=1,NNP
        F(NP,1)=F(NP,1)+V(NP,1)*C(NP)
        F(NP,2)=F(NP,2)+V(NP,2)*C(NP)
        F(NP,3)=F(NP,3)+V(NP,3)*C(NP)
  380 CONTINUE
C
      RETURN
      END
C
C
c
      SUBROUTINE FLUX1(F,CMATRX,CPFG,XWFG,X,ISED,IE,V,AKDC,NEFGD,
     >                 IQUAR,NDFG)
C 11/30/93
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE THE MATERIAL FLUXES AT DIFFUSION FINE GRIDS.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /WETX/ APHA1,APHA2,APHA3,APHA4
      COMMON /WETY/ BETA1,BETA2,BETA3,BETA4
      COMMON /WETZ/ GAMA1,GAMA2,GAMA3,GAMA4
      COMMON /SAZFM/ MXNPFG,MXKGL,MXKGLD,MXNEP,MXEPW,MXNPW,MXELW,MXNDB,
     >               mxnpww,mxelww
C
      DIMENSION F(MXNPFG,3),XWFG(MXNPFG,3)
      DIMENSION CMATRX(MXADNP,MXJBD),CPFG(MXNPFG)
      DIMENSION X(MAXNP,3),IE(MAXEL,11),ISED(MXKGLD,9)
      DIMENSION V(MAXNP,3),AKDC(8,MXKGLD,8)
      DIMENSION QB(8,8),QRX(8),QRY(8),QRZ(8)
      DIMENSION XQ(8),YQ(8),ZQ(8),CQ(8),VXQ(8),VYQ(8),VZQ(8)
      DIMENSION AKXYZQ(8,6),DL(8),DNX(8)
C
C ------- INITIALIZE THE FLUX FX(NP), FY(NP), AND FZ(NP).
C
      DO 100 NP=1,NDFG
      F(NP,1)=0.0D0
      F(NP,2)=0.0D0
      F(NP,3)=0.0D0
  100 CONTINUE
C
C ------- COMPUTE THE FLUX COMPONENTS BY APPLYING THE FINITE ELEMENT
C ------- METHOD TO THE DISPERSION TERMS
C
      DO 160 NP=1,NDFG
        CMATRX(NP,1)=0.0D0
        CMATRX(NP,2)=0.0D0
        CMATRX(NP,3)=0.0D0
        CMATRX(NP,4)=0.0D0
        CMATRX(NP,5)=0.0D0
  160 CONTINUE
C
C ------- COMPUTE THE ELEMENT MATRIX QB AND QRX, QRY, & QRZ
C
      DO 290 M=1,NEFGD
C
        MP=ISED(M,9)
        CALL ELENOD
     I      (ISED(M,5),ISED(M,7),
     O       NODE,IQ,IQ)
C
        DO 210 IQ=1,NODE
          NP=ISED(M,IQ)
          XQ(IQ)=XWFG(NP,1)
          YQ(IQ)=XWFG(NP,2)
          ZQ(IQ)=XWFG(NP,3)
          CQ(IQ)=CPFG(NP)
          DO 205 I=1,6
            AKXYZQ(IQ,I)=AKDC(IQ,M,I)
  205     CONTINUE
  210   CONTINUE
C
        CALL TQ468DV(QB,QRX,QRY,QRZ,CQ,XQ,YQ,ZQ,AKXYZQ,NODE,IQUAR)
C
C ------- ASSEMBLE QB(IQ,JQ) INTO THE GLOBAL MATRIX C(NP,1) AND
C ------- FORM THE LOAD VECTOR FX(NP), FY(NY), AND FZ(NP)
C
        DO 280 IQ=1,NODE
          NI=ISED(M,IQ)
          DO 240 JQ=1,NODE
            CMATRX(NI,1)=CMATRX(NI,1)+QB(IQ,JQ)
  240     CONTINUE
          F(NI,1)=F(NI,1)+QRX(IQ)
          F(NI,2)=F(NI,2)+QRY(IQ)
          F(NI,3)=F(NI,3)+QRZ(IQ)
  280   CONTINUE
C
C ----- COMPUTE ADVECTION FLUX
C
        DO IQ=1,NODE
          NN=IE(MP,IQ)
          XQ(IQ)=X(NN,1)
          YQ(IQ)=X(NN,2)
          ZQ(IQ)=X(NN,3)
          VXQ(IQ)=V(NN,1)
          VYQ(IQ)=V(NN,2)
          VZQ(IQ)=V(NN,3)
        ENDDO
        DO 285 IQ=1,NODE
          NP=ISED(M,IQ)
          IF(CMATRX(NP,2).EQ.1.0D0)GOTO 285
          CMATRX(NP,2)=1.0D0
          CALL BASE
c     I     (XQ,YQ,ZQ,XWFG(NP,1),XWFG(NP,2),XWFG(NP,3),MP,NODE,1,
     I     (XQ,YQ,ZQ,XWFG(NP,1),XWFG(NP,2),XWFG(NP,3),NODE,1,
     O      DL,DNX,DNX,DNX)
          VX=0.0D0
          VY=0.0D0
          VZ=0.0D0
          DO IIQ=1,NODE
            VX=VX+VXQ(IIQ)*DL(IIQ)
            VY=VY+VYQ(IIQ)*DL(IIQ)
            VZ=VZ+VZQ(IIQ)*DL(IIQ)
          ENDDO
          CMATRX(NP,3)=CMATRX(NP,3)+VX*CPFG(NP)
          CMATRX(NP,4)=CMATRX(NP,4)+VY*CPFG(NP)
          CMATRX(NP,5)=CMATRX(NP,5)+VZ*CPFG(NP)
  285   CONTINUE
  290 CONTINUE
C
C ------- SOLVE THE MATRIX EQUATION CX=B
C
      DO 370 NP=1,NDFG
        IF(CMATRX(NP,1).EQ.0.0D0)GOTO 370
        F(NP,1)=F(NP,1)/CMATRX(NP,1)+CMATRX(NP,3)
        F(NP,2)=F(NP,2)/CMATRX(NP,1)+CMATRX(NP,4)
        F(NP,3)=F(NP,3)/CMATRX(NP,1)+CMATRX(NP,5)
  370 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE TQ468DV(QB,QRX,QRY,QRZ, CQ,XQ,YQ,ZQ,AKXYZQ,NODE,IQUAR)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE THE INTEGRATION OF N(I)*N(J) AND -N(I)*D>GRAD(C)
C ------- OVER AN ELEMENT.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: CQ(8), XQ(8), YQ(8), ZQ(8), AKXYZQ(8,6)
C
C ------- OUTPUT: QB(8,8), QRX(8), QRZ(8), AND QRZ(8).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(8),L1(3),L2(3),L3(3),LL1(4),LL2(4),LL3(4),LL4(4)
C
      DIMENSION QB(8,8),QRX(8),QRY(8),QRZ(8),CQ(8),XQ(8),YQ(8),ZQ(8)
      DIMENSION AKXYZQ(8,6), AKXYZK(6)
      DIMENSION S(8),T(8),U(8), DNX(8),DNY(8),DNZ(8),W(8)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0, -1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA T/-1.0D0,-1.0D0,1.0D0,1.0D0, -1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA U/-1.0D0,-1.0D0,-1.0D0,-1.0D0, 1.0D0,1.0D0,1.0D0,1.0D0/
C
      IF(IQUAR.EQ.1 .OR. IQUAR.EQ.3)THEN
        P=0.577350269189626D0
        L1(1)=0.666666666666667D0
        L1(2)=0.166666666666667D0
        L1(3)=0.166666666666667D0
        L2(1)=0.166666666666667D0
        L2(2)=0.666666666666667D0
        L2(3)=0.166666666666667D0
        L3(1)=0.166666666666667D0
        L3(2)=0.166666666666667D0
        L3(3)=0.666666666666667D0
        LL1(1)=0.58541020D0
        LL1(2)=0.13819660D0
        LL1(3)=0.13819660D0
        LL1(4)=0.13819660D0
        LL2(1)=0.13819660D0
        LL2(2)=0.58541020D0
        LL2(3)=0.13819660D0
        LL2(4)=0.13819660D0
        LL3(1)=0.13819660D0
        LL3(2)=0.13819660D0
        LL3(3)=0.58541020D0
        LL3(4)=0.13819660D0
        LL4(1)=0.13819660D0
        LL4(2)=0.13819660D0
        LL4(3)=0.13819660D0
        LL4(4)=0.58541020D0
      ELSE
        P=1.0D0
        L1(1)=1.0D0
        L1(2)=0.0D0
        L1(3)=0.0D0
        L2(1)=0.0D0
        L2(2)=1.0D0
        L2(3)=0.0D0
        L3(1)=0.0D0
        L3(2)=0.0D0
        L3(3)=1.0D0
        LL1(1)=1.0D0
        LL1(2)=0.0D0
        LL1(3)=0.0D0
        LL1(4)=0.0D0
        LL2(1)=0.0D0
        LL2(2)=1.0D0
        LL2(3)=0.0D0
        LL2(4)=0.0D0
        LL3(1)=0.0D0
        LL3(2)=0.0D0
        LL3(3)=1.0D0
        LL3(4)=0.0D0
        LL4(1)=0.0D0
        LL4(2)=0.0D0
        LL4(3)=0.0D0
        LL4(4)=1.0D0
      ENDIF
C
C ------- INITIATE MATRICES QB(IQ,JQ), QRX(IQ), QRY(IQ) & QRZ(IQ)
C
      DO 100 IQ=1,8
      QRX(IQ)=0.0
      QRY(IQ)=0.0
      QRZ(IQ)=0.0
      DO 100 JQ=1,8
  100 QB(IQ,JQ)=0.0
C
C ------- SUMMATION OF THE INTEGRAND OVER THE GAUSSIAN POINTS
C
      DO 490 KG=1,NODE
C
C ------- DETERMINE LOACAL COORDINATE OF GAUSSIAN POINT KG
C
C ------- CALCULATE VALUES OF BASIS FUNCTIONS N(IQ) AND THEIR
C ------- DERIVATIVES DNX(IQ), DNY(IQ), AND DNZ(IQ), W.R.T. TO
C ------- X, Y, AND Z, RESPECTIVELY, AT THE GAUSSIAN POINT KG.
C
        IF(NODE.EQ.8)THEN
          SS=P*S(KG)
          TT=P*T(KG)
          UU=P*U(KG)
          CALL SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,0,
     O       N,DNX,DNY,DNZ,W,DJAC)
        ELSEIF(NODE.EQ.6)THEN
          IF(KG.LE.3)THEN
            XSI=-P
            KKG=KG
          ELSE
            XSI=P
            KKG=KG-3
          ENDIF
          DL1=L1(KKG)
          DL2=L2(KKG)
          DL3=L3(KKG)
          CALL SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,0,
     O       N,DNX,DNY,DNZ,W,DJAC)
          DJAC=DJAC/3.0D0
        ELSEIF(NODE.EQ.4)THEN
          D1=LL1(KG)
          D2=LL2(KG)
          D3=LL3(KG)
          D4=LL4(KG)
          CALL SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,0,
     O       N,DNX,DNY,DNZ,W,DJAC)
          DJAC=0.25D0*DJAC
        ENDIF
C
        DO 270 I=1,6
          AKXYZK(I)=AKXYZQ(KG,I)
  270   CONTINUE
C
        AKXK=AKXYZK(1)*DJAC
        AKYK=AKXYZK(2)*DJAC
        AKZK=AKXYZK(3)*DJAC
        AKXYK=AKXYZK(4)*DJAC
        AKXZK=AKXYZK(5)*DJAC
        AKYZK=AKXYZK(6)*DJAC
C
C ------- ACCUMULATE THE SUMS TO OBTAIN THE MATRIX INTEGRALS QB(IQ,JQ)
C ------- AND QRX(IQ), QRY(IQ), AND QRZ(IQ)
C
        DO 390 IQ=1,NODE
          DO 350 JQ=1,NODE
            QB(IQ,JQ)=QB(IQ,JQ)+ N(IQ)*N(JQ)*DJAC
            QRX(IQ)=QRX(IQ)-N(IQ)*CQ(JQ)*(AKXK*DNX(JQ)+AKXYK*DNY(JQ)+
     1              AKXZK*DNZ(JQ))
            QRY(IQ)=QRY(IQ)-N(IQ)*CQ(JQ)*(AKXYK*DNX(JQ)+AKYK*DNY(JQ)+
     1              AKYZK*DNZ(JQ))
            QRZ(IQ)=QRZ(IQ)-N(IQ)*CQ(JQ)*(AKXZK*DNX(JQ)+AKYZK*DNY(JQ)+
     1              AKZK*DNZ(JQ))
  350     CONTINUE
  390   CONTINUE
C
  490 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE ADVBC
     I    (IE,X,V,VP,THN,CP,
     I     CVB,IVTYP,NPVB,ISV,
     I     QCB,ICTYP,NPCB,ISC, CDB,IDTYP,NPDB, DCOSB,ISB,NPBB,
     O     CS,
     M     RI,RL,RLD)
C
C $$$$$ TO APPLY BOUNDARY CONDITIONS TO LAGRANGIAN STEP.
C
C NOTE: ONLY ONE COMPONENT IS CONSIDERED IN THIS ROUTINE.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
C
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /TCBC/ MXCNP,MXCES,MXCPR,MXCDP,NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /TDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
C
      DIMENSION IE(MAXEL,11),X(MAXNP,3),NPBB(MAXBNP)
      DIMENSION CS(MAXNP),THN(MAXNP,2),RI(MXADNP),RL(MXADNP),RLD(MXADNP)
      DIMENSION V(MAXNP,3),VP(MAXNP,3),CP(MAXNP)
C
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES)
      DIMENSION CVB(MXVPR),IVTYP(MXVES),NPVB(MXVNP),ISV(5,MXVES)
      DIMENSION QCB(MXCPR),ICTYP(MXCES),NPCB(MXCNP),ISC(5,MXCES)
      DIMENSION CDB(MXDPR),IDTYP(MXDNP),NPDB(MXDNP)
C
      DIMENSION RQI(4),RQL(4)
      DIMENSION XQ(4),YQ(4),ZQ(4),VXQ(4),VYQ(4),VZQ(4),THQ(4),COFTH(3)
C
      DIMENSION KGB(4,6,3)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
C ******* APPLY CAUCHY CONDITION: QC=V.N.C - N.(THETA)D.GRAD(C)
C
  100 IF(NCES.EQ.0)GOTO 500
C
      DO 110 NI=1,NNP
        RI(NI)=0.0
        RL(NI)=0.0
  110 CONTINUE
C
      DO 150 MP=1,NCES
        ITYP=ICTYP(MP)
        QCBMP=QCB(ITYP)
        MPB=ISC(5,MP)
        LS=ISB(5,MPB)
        M=ISB(6,MPB)
        IF(IE(M,5).EQ.0)THEN
          ID=3
        ELSEIF(IE(M,7).EQ.0)THEN
          ID=2
        ELSE
          ID=1
        ENDIF
C
        NODE=4
        DO 130 IQ=1,4
          I=KGB(IQ,LS,ID)
          IF(I.EQ.0)THEN
            NODE=3
            GOTO 130
          ENDIF
          NI=IE(M,I)
          XQ(IQ)=X(NI,1)
          YQ(IQ)=X(NI,2)
          ZQ(IQ)=X(NI,3)
          VXQ(IQ)=(V(NI,1)+VP(NI,1))*0.5D0
          VYQ(IQ)=(V(NI,2)+VP(NI,2))*0.5D0
          VZQ(IQ)=(V(NI,3)+VP(NI,3))*0.5D0
  130   CONTINUE
C
        CALL Q34ADB(RQI,RQL,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB(1,MPB),QCBMP,THQ,
     >              1,NODE,1,IQUAR)
C
        DO 140 IQ=1,NODE
          I=KGB(IQ,LS,ID)
          NI=IE(M,I)
          RL(NI)=RL(NI)-RQL(IQ)
          RI(NI)=RI(NI)-RQI(IQ)
  140   CONTINUE
C
  150 CONTINUE
C
      DO 190 NPP=1,NCNP
        NP=NPCB(NPP)
        NP=NPBB(NP)
        IF(RL(NP).NE.0.0) THEN
          CS(NP)=RI(NP)/RL(NP)
          RLD(NP)=RL(NP)
        ELSE
          RLD(NP)=0.0D0
        ENDIF
  190 CONTINUE
C
c     IF(KSORP.NE.1)GOTO 500
C
c     DO 210 NI=1,NNP
c       RI(NI)=0.0
c       RL(NI)=0.0
c 210 CONTINUE
C
c     DO 250 MP=1,NCES
c       MPB=ISC(5,MP)
c       LS=ISB(5,MPB)
c       M=ISB(6,MPB)
c       IF(IE(M,5).EQ.0)THEN
c         ID=3
c       ELSEIF(IE(M,7).EQ.0)THEN
c         ID=2
c       ELSE
c         ID=1
c       ENDIF
C
c       NODE=4
c       DO 230 IQ=1,4
c         I=KGB(IQ,LS,ID)
c         IF(I.EQ.0)THEN
c           NODE=3
c           GOTO 230
c         ENDIF
c         NI=IE(M,I)
c         XQ(IQ)=X(NI,1)
c         YQ(IQ)=X(NI,2)
c         ZQ(IQ)=X(NI,3)
c         VXQ(IQ)=THN(NI,1)
c         THQ(IQ)=THN(NI,2)
c         VYQ(IQ)=0.0D0
c         VZQ(IQ)=0.0D0
c         IF(IQ.LE.3)COFTH(IQ)=1.0D0
c 230   CONTINUE
C
c       CALL Q34ADB(RQI,RQL,XQ,YQ,ZQ,VXQ,VYQ,VZQ,COFTH,QCBMP,THQ,1,
c    >            NODE,2,IQUAR)
C
c       DO 240 IQ=1,NODE
c         I=KGB(IQ,LS,ID)
c         NI=IE(M,I)
c         IF(RLD(NI).NE.0.0D0)THEN
c           RL(NI)=RL(NI)+RQL(IQ)
c           RI(NI)=RI(NI)+RQI(IQ)*CS(NI)+(RQL(IQ)-RQI(IQ))*CP(NI)
c         ENDIF
c 240   CONTINUE
C
c 250 CONTINUE
C
c     DO 290 NPP=1,NCNP
c       NP=NPCB(NPP)
c       NP=NPBB(NP)
c       IF(RL(NP).NE.0.0)CS(NP)=RI(NP)/RL(NP)
c 290 CONTINUE
C
C ******* APPLY VARIABLE BOUNDARY CONDITIONS
C
  500 IF(NVES.EQ.0)GOTO 700
C
      DO 510 NP=1,NNP
        RI(NP)=0.0
        RL(NP)=0.0
  510 CONTINUE
C
      DO 550 MP=1,NVES
C
        ITYP=IVTYP(MP)
        CINMP=CVB(ITYP)
        MPB=ISV(5,MP)
        LS=ISB(5,MPB)
        M=ISB(6,MPB)
        IF(IE(M,5).EQ.0)THEN
          ID=3
        ELSEIF(IE(M,7).EQ.0)THEN
          ID=2
        ELSE
          ID=1
        ENDIF
C
        NODE=4
        DO 530 IQ=1,4
          I=KGB(IQ,LS,ID)
          IF(I.EQ.0)THEN
            NODE=3
            GOTO 530
          ENDIF
          NI=IE(M,I)
          XQ(IQ)=X(NI,1)
          YQ(IQ)=X(NI,2)
          ZQ(IQ)=X(NI,3)
          VXQ(IQ)=(V(NI,1)+VP(NI,1))*0.5D0
          VYQ(IQ)=(V(NI,2)+VP(NI,2))*0.5D0
          VZQ(IQ)=(V(NI,3)+VP(NI,3))*0.5D0
  530   CONTINUE
C
        CALL Q34ADB(RQI,RQL,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB(1,MPB),CINMP,
     >              THQ,3,NODE,1,IQUAR)
C
        DO 540 IQ=1,NODE
          I=KGB(IQ,LS,ID)
          NI=IE(M,I)
          RL(NI)=RL(NI)-RQL(IQ)
          RI(NI)=RI(NI)-RQI(IQ)
  540   CONTINUE
C
  550 CONTINUE
C
      DO 590 NPP=1,NVNP
        NP=NPVB(NPP)
        NP=NPBB(NP)
        IF(RL(NP).NE.0.0) THEN
          CS(NP)=RI(NP)/RL(NP)
          RLD(NP)=RL(NP)
        ELSE
          RLD(NP)=0.0D0
        ENDIF
  590 CONTINUE
C
c     IF(KSORP.NE.1)GOTO 700
C
c     DO 610 NP=1,NNP
c       RI(NP)=0.0
c       RL(NP)=0.0
c 610 CONTINUE
C
c     DO 650 MP=1,NVES
C
c       MPB=ISV(5,MP)
c       LS=ISB(5,MPB)
c       M=ISB(6,MPB)
c       IF(IE(M,5).EQ.0)THEN
c         ID=3
c       ELSEIF(IE(M,7).EQ.0)THEN
c         ID=2
c       ELSE
c         ID=1
c       ENDIF
C
c       NODE=4
c       DO 630 IQ=1,4
c         I=KGB(IQ,LS,ID)
c         IF(I.EQ.0)THEN
c           NODE=3
c           GOTO 630
c         ENDIF
c         NI=IE(M,I)
c         XQ(IQ)=X(NI,1)
c         YQ(IQ)=X(NI,2)
c         ZQ(IQ)=X(NI,3)
c         VXQ(IQ)=THN(NI,1)
c         THQ(IQ)=THN(NI,2)
c         VYQ(IQ)=0.0D0
c         VZQ(IQ)=0.0D0
c         IF(IQ.LE.3)COFTH(IQ)=1.0D0
c 630   CONTINUE
C
c       CALL Q34ADB(RQI,RQL,XQ,YQ,ZQ,VXQ,VYQ,VZQ,COFTH,CINMP,
c    >              THQ,3,NODE,1,IQUAR)
C
c       DO 640 IQ=1,NODE
c         I=KGB(IQ,LS,ID)
c         NI=IE(M,I)
c         IF(RLD(NI).NE.0.0D0)THEN
c           RL(NI)=RL(NI)+RQL(IQ)
c           RI(NI)=RI(NI)+RQI(IQ)*CS(NI)+(RQL(IQ)-RQI(IQ))*CP(NI)
c         ENDIF
c 640   CONTINUE
C
c 650 CONTINUE
C
c     DO 690 NPP=1,NVNP
c       NP=NPVB(NPP)
c       NP=NPBB(NP)
c       IF(RL(NP).NE.0.0)CS(NP)=RI(NP)/RL(NP)
c 690 CONTINUE
C
C ******* APPLY DIRICHLET BOUNDARY CONDITION
C
  700 CONTINUE
      IF(NDNP.NE.0) THEN
        DO 790 NPP=1,NDNP
          NP=NPDB(NPP)
          ITYP=IDTYP(NPP)
          BB=CDB(ITYP)
          CS(NP)=BB
  790   CONTINUE
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE Q34ADB(RQI,RQL,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB,CINMP,THN,
     >                  IBC,NODE,IADV,IQUAR)
C
C $$$$$ TO COMPUTE BOUNDARY-SURFACE VOLUME FLUXES AND MATERIAL FLUXES
C       OVER A BOUNDARY SURFACE IN LAGRANGIAN STEP.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(4)
C
      DIMENSION RQL(4),RQI(4),XQ(4),YQ(4),ZQ(4),THN(4)
      DIMENSION VXQ(4),VYQ(4),VZQ(4),DCOSB(3)
      DIMENSION S(4),T(4),DNSS(4),DNTT(4)
      DIMENSION DL1(3),DL2(3),DL3(3),WG(3)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0/, T/-1.0D0,-1.0D0,1.0D0,1.0D0/
C
      DATA WG/0.333333333333333D0, 0.333333333333333D0,
     >        0.333333333333333D0/
C
      IF(IQUAR.EQ.3 .OR. IQUAR.EQ.4)THEN
        P=0.577350269189626D0
        DL1(1)=0.666666666666667D0
        DL1(2)=0.166666666666667D0
        DL1(3)=0.166666666666667D0
        DL2(1)=0.166666666666667D0
        DL2(2)=0.666666666666667D0
        DL2(3)=0.166666666666667D0
        DL3(1)=0.166666666666667D0
        DL3(2)=0.166666666666667D0
        DL3(3)=0.666666666666667D0
      ELSE
        P=1.0D0
        DL1(1)=1.0D0
        DL1(2)=0.0D0
        DL1(3)=0.0D0
        DL2(1)=0.0D0
        DL2(2)=1.0D0
        DL2(3)=0.0D0
        DL3(1)=0.0D0
        DL3(2)=0.0D0
        DL3(3)=1.0D0
      ENDIF
C
C ------- INITIATE VECTOR RQL(IQ) AND RQI(IQ)
C
      DO 100 IQ=1,4
        RQL(IQ)=0.0
        RQI(IQ)=0.0
  100 CONTINUE
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS IF NODE.EQ.3
C
      IF(NODE.EQ.3)THEN
        DXDDL2=XQ(2)-XQ(1)
        DYDDL2=YQ(2)-YQ(1)
        DZDDL2=ZQ(2)-ZQ(1)
        DXDDL3=XQ(3)-XQ(1)
        DYDDL3=YQ(3)-YQ(1)
        DZDDL3=ZQ(3)-ZQ(1)
        DETX=DYDDL2*DZDDL3-DYDDL3*DZDDL2
        DETY=DXDDL2*DZDDL3-DXDDL3*DZDDL2
        DETZ=DXDDL2*DYDDL3-DXDDL3*DYDDL2
        DET1=DSQRT(DETX*DETX+DETY*DETY+DETZ*DETZ)*0.5D0
      ENDIF
C
C ------- SUMMATION OF THE INTEGRAND OVER THE GAUSSIAN POINTS
C
      DO 690 KG=1,NODE
C
C ------- DETERMINE LOACAL COORDINATE OF GAUSSIAN POINT KG
C
        IF(NODE.EQ.4)THEN
C
          SS=P*S(KG)
          TT=P*T(KG)
          SM=1.0D0-SS
          SP=1.0D0+SS
          TM=1.0D0-TT
          TP=1.0D0+TT
          N(1)=0.25D0*SM*TM
          N(2)=0.25D0*SP*TM
          N(3)=0.25D0*SP*TP
          N(4)=0.25D0*SM*TP
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS
C
          DNSS(1)=-0.25D0*TM
          DNSS(2)= 0.25D0*TM
          DNSS(3)= 0.25D0*TP
          DNSS(4)=-0.25D0*TP
          DNTT(1)=-0.25D0*SM
          DNTT(2)=-0.25D0*SP
          DNTT(3)= 0.25D0*SP
          DNTT(4)= 0.25D0*SM
          DXDSS=0.0D0
          DYDSS=0.0D0
          DZDSS=0.0D0
          DXDTT=0.0D0
          DYDTT=0.0D0
          DZDTT=0.0D0
          DO 290 IQ=1,4
            DXDSS=DXDSS+XQ(IQ)*DNSS(IQ)
            DYDSS=DYDSS+YQ(IQ)*DNSS(IQ)
            DZDSS=DZDSS+ZQ(IQ)*DNSS(IQ)
            DXDTT=DXDTT+XQ(IQ)*DNTT(IQ)
            DYDTT=DYDTT+YQ(IQ)*DNTT(IQ)
            DZDTT=DZDTT+ZQ(IQ)*DNTT(IQ)
  290     CONTINUE
          DETZ=DXDSS*DYDTT-DYDSS*DXDTT
          DETY=-DXDSS*DZDTT+DZDSS*DXDTT
          DETX=DYDSS*DZDTT-DZDSS*DYDTT
          DET=DSQRT(DETX*DETX+DETY*DETY+DETZ*DETZ)
C
        ELSE
C
          N(1)=DL1(KG)
          N(2)=DL2(KG)
          N(3)=DL3(KG)
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS
C
          DET=DET1*WG(KG)
C
        ENDIF
C
C ------- ACCUMULATE THE SUMS TO OBTAIN THE FLUX INTEGRALS RQL AND RQI
C
        GOTO(310,410,510)IBC
C
C ******* CAUCHY CONDITIONS
C
  310   VXK=0.0
        VYK=0.0
        VZK=0.0
        DO 320 IQ=1,NODE
          VXK=VXK+VXQ(IQ)*N(IQ)
          VYK=VYK+VYQ(IQ)*N(IQ)
          VZK=VZK+VZQ(IQ)*N(IQ)
  320   CONTINUE
        VNK=VXK*DCOSB(1)+VYK*DCOSB(2)+VZK*DCOSB(3)
        QBMP=CINMP                                                       5/10/95
        IF(IADV.EQ.2)THEN
          THK=0.0D0
          DO 330 IQ=1,NODE
            THK=THK+THN(IQ)*N(IQ)
  330     CONTINUE
          QBMP=THK
        ENDIF
        DO 390 IQ=1,NODE
          RQI(IQ)=RQI(IQ)+N(IQ)*QBMP*DET
          RQL(IQ)=RQL(IQ)+N(IQ)*VNK*DET
  390   CONTINUE
        GOTO 690
C
C ****** NEUMANN CONDITIONS
C
  410   CONTINUE
        GOTO 690
C
C ******* VARIABLE CONDITIONS
C
  510   VXK=0.0
        VYK=0.0
        VZK=0.0
        DO 520 IQ=1,NODE
          VXK=VXK+VXQ(IQ)*N(IQ)
          VYK=VYK+VYQ(IQ)*N(IQ)
          VZK=VZK+VZQ(IQ)*N(IQ)
  520   CONTINUE
        VNK=VXK*DCOSB(1)+VYK*DCOSB(2)+VZK*DCOSB(3)
        IF(VNK.GE.0.0 .AND. IADV.EQ.1)GOTO 690
        QBMP=CINMP*VNK                                                    5/9/95
        IF(IADV.EQ.2)THEN
          THK=0.0D0
          DO 530 IQ=1,NODE
            THK=THK+THN(IQ)*N(IQ)
          CONTINUE
  530     QBMP=THK
        ENDIF
        DO 590 IQ=1,NODE
          RQI(IQ)=RQI(IQ)+N(IQ)*QBMP*DET
          RQL(IQ)=RQL(IQ)+N(IQ)*VNK*DET
  590   CONTINUE
C
  690 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE TASEMB
     > (CMATRX,RLD,CW,CP,CSTAR,X,IE,LRN,NLRL,LRL, WETAB,
     >  V,VP,TH,THP,AKHC,AKDC, LES,ISTYP,SOS,NPW,IWTYP,WSS,IRHO,
     >  RHOTYP,PROP,DINTS,RKD,TRANC, DTI,W,WV,MICONF,ID,KSS,
     >  ISED,XPFG,
     >  MXNOD1,MXNOD2,NNEL, KKK)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO ASSEMBLE THE GLOBAL COEFFICIENT MATRIX AND GLOBAL LOAD
C ------- VECTOR.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: CW(MAXNP),CP(MAXNP),CSTAR(MAXNP), X(NNP,3),
C -------        IE(NEL,9),LRN(MXJBD,NNP),NLRL(MAXNP),LRL(MXKBD,NNP),
C -------        WETAB(12,NEL),V(MAXNP,3),
C -------        VP(NNP,3),TH(NEL,8),THP(NEL,8),
C -------        AKDC(8,NEL,6), LES(NSEL),ISTYP(NSEL),
C -------        SOS(NSPR,2), NPW(NWNP),IWTYP(NWNP),WSS(NWPR,2),
C -------        PROP(NMAT,NMPPM), DELT, AND KSS.
C
C ------- OUTPUT: CMATRX(NNP,MXJBD) AND RLD(NNP).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KD,LAMBDA
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /SAZFM/ MXNPFG,MXKGL,MXKGLD,MXNEP,MXEPW,MXNPW,MXELW,MXNDB,
     >               mxnpww,mxelww
      COMMON /CELS/ MXSEL,MXSPR,MXSDP,NSEL,NSPR,NSDP,KSAI
      COMMON /CNPS/ MXWNP,MXWPR,MXWDP,NWNP,NWPR,NWDP,KWAI
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /TADP/ ADPEPS,ADPARM,IZOOM,IDZOOM,IEPC,NXG,NYG,NZG,
     >              NXW,NYW,NZW,NXD,NYD,NZD,IDETQ
C
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
      COMMON /WETX/ APHA1,APHA2,APHA3,APHA4
      COMMON /WETY/ BETA1,BETA2,BETA3,BETA4
      COMMON /WETZ/ GAMA1,GAMA2,GAMA3,GAMA4
C
      COMMON /CHEM/ MXNCC,NCC
      COMMON /MICROB/ GRATE,YCOEFF,RTARDS,RTARDO,RTARDN,SCOEFF,
     >                ECOEFF,DCOEFF,SATURC,PCOEFF,COFK
C
      DIMENSION CMATRX(MXADNP,MXJBD),RLD(MXADNP)
      DIMENSION CP(MAXNP,MXNCC),CW(MAXNP,MXNCC),CSTAR(MXNOD1)
      DIMENSION X(MAXNP,3),IE(MAXEL,11),LRN(MXJBD,MXADNP),NLRL(MAXNP)
      DIMENSION RKD(MAXMAT,MXNCC),TRANC(MAXMAT,MXNCC),LRL(MXKBD,MAXNP)
      DIMENSION WETAB(12,MAXEL),V(MAXNP,3),VP(MAXNP,3),DTI(MXNOD2)
      DIMENSION TH(MAXEL,8),THP(MAXEL,8),AKDC(8,MXKGLD,8),DINTS(MXNCC)
      DIMENSION RHOTYP(MAXMAT),PROP(MAXMAT,MXMPPM),AKHC(8,MAXEL,7)
      DIMENSION SOS(MXSPR,2),LES(MXSEL),ISTYP(MXSEL,MXNCC)
      DIMENSION WSS(MXWPR,2),NPW(MXWNP),IWTYP(MXWNP,MXNCC)
      DIMENSION ISED(MXKGLD,9),XPFG(MXNPFG,3)
C
      DIMENSION QA(8,8),QAA(8,8),QB(8,8),QC(8,8),QD(8,8),QV(8,8),QR(8)
      DIMENSION QE(8,8),XG(8),YG(8),ZG(8),VXG(8),VYG(8),VZG(8)
      DIMENSION XQ(8),YQ(8),ZQ(8),VXQ(8),VYQ(8),VZQ(8),THG(8),DL(8)
      DIMENSION AKXYZQ(8,6),RHOQ(8),CCG(8,7),DNX(8)
      DIMENSION DSDCQ(8,4),CWQ(8,3),RSQ(4,8),KD(4),CNP(7),CCQ(8,7)
      DIMENSION GRATE(4),YCOEFF(4),RTARDS(4),RTARDO(4),RTARDN(4),
     >          SCOEFF(4),ECOEFF(4),DCOEFF(4),SATURC(4),PCOEFF(4)
C
      W1=W
      W2=1.0D0-W
      W1V=WV
      W2V=1.0D0-WV
      IF(KSS.NE.0) GO TO 100
      DO 90 NP=1,NNP
        DTI(NP)=0.0D0
 90   CONTINUE
      W1=1.0D0
      W2=0.0
      W1V=1.0D0
      W2V=0.0
  100 IF(IMID.NE.0 .AND. ID.EQ.1) THEN
        W1=1.0D0
        W2=0.0
        W1V=1.0D0
        W2V=0.0
      ENDIF
      VTERM=0.0
      IF(KSS.EQ.0) VTERM=1.0D0
      IF(LGRN.EQ.0) VTERM=1.0D0
C
C ------- INITIALIZE MATRICES CMATRX(NP,IB) AND RLD(NP)
C
      IF(ID.NE.2)THEN
        DO 150 NP=1,MXADNP
          RLD(NP)=0.0
          DO 150 I=1,MXJBD
            CMATRX(NP,I)=0.0
  150   CONTINUE
      ENDIF
C
C ******* LOOP OVER ALL ELEMENTS TO FORM THE GLOBAL MATRIX EQUATION
C
      DO 690 MM=1,NNEL
C
        IF(ID.EQ.2)THEN
          M=ISED(MM,9)
          CALL ELENOD
     I        (ISED(MM,5),ISED(MM,7),
     O         NQ,IQ,IQ)
        ELSEIF(ID.EQ.1)THEN
          M=MM
          IF(IE(M,11).NE.0)GOTO 690
        ELSE
          M=MM
        ENDIF
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,IQ,IQ)
        IF(ID.NE.2) NQ=NODE
        DO IQ=1,NODE
          NP=IE(M,IQ)
          XQ(IQ)=X(NP,1)
          YQ(IQ)=X(NP,2)
          ZQ(IQ)=X(NP,3)
          DO K=1,NCC
            CCG(IQ,K)=CP(NP,K)
          ENDDO
          IF(IMID.EQ.0)THEN
            VXG(IQ)=W1V*V(NP,1)+W2V*VP(NP,1)
            VYG(IQ)=W1V*V(NP,2)+W2V*VP(NP,2)
            VZG(IQ)=W1V*V(NP,3)+W2V*VP(NP,3)
          ELSE
            VXG(IQ)=0.5D0*(V(NP,1)+VP(NP,1))
            VYG(IQ)=0.5D0*(V(NP,2)+VP(NP,2))
            VZG(IQ)=0.5D0*(V(NP,3)+VP(NP,3))
          ENDIF
        ENDDO
C
        MTYP=IE(M,9)
        IF(IRHO.EQ.1)RHOW=1.0D0/RHOTYP(MTYP)
        RHOSTR=1.0D0
        SOSQ=0.0
        SOSC=0.0
        SOSQC=0.0D0
        IF(NSEL.NE.0) THEN
          IDO=1
          I=0
  210     CONTINUE
            I=I+1
            MP=LES(I)
            IF(MP.NE.M) IDO=0
          IF(IDO.EQ.1 .AND. I.LT.NSEL)GOTO 210
          IF(IDO.EQ.0)THEN
            ITYP=ISTYP(I,KKK)
            SOSQ=SOS(ITYP,1)
            SOSC=SOS(ITYP,2)
            SOSQC=0.0D0
            IF(SOSQ.LE.0.0D0)GOTO 220
            SOSQC=SOSQ
C
C ***** IF DENSITY EFFECT IS TAKEN INTO ACCOUNT
C
            IF(IRHO.EQ.1)THEN
              DO K=1,NCC
                ITYP=ISTYP(I,K)
                CNP(K)=SOS(ITYP,2)
              ENDDO
c              CALL RHOFUN(RHOSTR,AMU,RHOW,DINTS,CNP,RHOMU)
              CALL RHOFUN(RHOSTR,AMU,RHOW,DINTS,CNP)
              SOSQC=SOSQ*RHOSTR
            ENDIF
          ENDIF
C
        ENDIF
  220   CONTINUE
C
        IF(IWET.NE.0 .AND. LGRN.EQ.0)THEN
          IF(NQ.EQ.4)THEN
            APHA1=WETAB(1,M)
            APHA2=WETAB(2,M)
            APHA3=WETAB(3,M)
            BETA1=WETAB(4,M)
            BETA2=WETAB(5,M)
            BETA3=WETAB(6,M)
          ELSEIF(NQ.EQ.6)THEN
            APHA1=WETAB(1,M)
            APHA2=WETAB(2,M)
            APHA3=WETAB(3,M)
            BETA1=WETAB(4,M)
            BETA2=WETAB(5,M)
            BETA3=WETAB(6,M)
            GAMA1=WETAB(7,M)
            GAMA2=WETAB(8,M)
            GAMA3=WETAB(9,M)
          ELSE
            APHA1=WETAB(1,M)
            APHA2=WETAB(2,M)
            APHA3=WETAB(3,M)
            APHA4=WETAB(4,M)
            BETA1=WETAB(5,M)
            BETA2=WETAB(6,M)
            BETA3=WETAB(7,M)
            BETA4=WETAB(8,M)
            GAMA1=WETAB(9,M)
            GAMA2=WETAB(10,M)
            GAMA3=WETAB(11,M)
            GAMA4=WETAB(12,M)
          ENDIF
        ENDIF
C
        KD(4)=RKD(MTYP,KKK)
        DO I=1,MIN0(3,NCC)
          KD(I)=RKD(MTYP,I)
        ENDDO
        RHOB=PROP(MTYP,1)
        LAMBDA=TRANC(MTYP,KKK)
C
        DO 250 IQ=1,NQ
          NP=IE(M,IQ)
          IF(ID.EQ.2)THEN
            NP=ISED(MM,IQ)
            XG(IQ)=XPFG(NP,1)
            YG(IQ)=XPFG(NP,2)
            ZG(IQ)=XPFG(NP,3)
          ENDIF
C
          VXQ(IQ)=0.0D0
          VYQ(IQ)=0.0D0
          VZQ(IQ)=0.0D0
          IF(ID.NE.2)THEN
            DO K=1,NCC
              CCQ(IQ,K)=CCG(IQ,K)
            ENDDO
          ELSE
            DO K=1,NCC
              CCQ(IQ,K)=0.0D0
            ENDDO
          ENDIF
          DO I=1,6
            AKXYZQ(IQ,I)=0.0D0
          ENDDO
          IF(MICONF.NE.1 .OR. KKK.GT.3)THEN
            IF(ID.NE.2)THEN
              VXQ(IQ)=VXG(IQ)
              VYQ(IQ)=VYG(IQ)
              VZQ(IQ)=VZG(IQ)
              DO I=1,6
                AKXYZQ(IQ,I)=AKHC(IQ,M,I)
              ENDDO
            ELSE
              DO I=1,6
                AKXYZQ(IQ,I)=AKDC(IQ,MM,I)
              ENDDO
              CALL BASE
c     I            (XQ,YQ,ZQ,XG(IQ),YG(IQ),ZG(IQ),M,NODE,1,
     I            (XQ,YQ,ZQ,XG(IQ),YG(IQ),ZG(IQ),NODE,1,
     O             DL,DNX,DNX,DNX)
              DO I=1,NODE
                VXQ(IQ)=VXQ(IQ)+DL(I)*VXG(I)
                VYQ(IQ)=VYQ(IQ)+DL(I)*VYG(I)
                VZQ(IQ)=VZQ(IQ)+DL(I)*VZG(I)
                DO K=1,NCC
                  CCQ(IQ,K)=CCQ(IQ,K)+DL(I)*CCG(IQ,K)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
C
          DO I=1,3
            RSQ(I,IQ)=0.0D0
            CWQ(IQ,I)=0.0D0
            DSDCQ(IQ,I)=0.0D0
          ENDDO
          DSDCQ(IQ,4)=0.0D0
C
          DO I=1,MIN0(4,NCC)
            DSDCQ(IQ,I)=KD(I)
          ENDDO
C
          IF(KSS.EQ.0)THEN
            IF(IMID.EQ.0)THEN
              DO I=1,NCC
                CNP(I)=W1*CW(NP,I)+W2*CCQ(IQ,I)
              ENDDO
            ELSE
              DO I=1,NCC
                CNP(I)=0.5D0*(CW(NP,I)+CCQ(IQ,I))
              ENDDO
            ENDIF
C
C ***** Assign the current concentration of microbes to CSQ array for
C       calculating removal rate
C
            DO I=1,MIN0(3,NCC)
              CWQ(IQ,I)=CNP(I)
            ENDDO
C
            CALL RXRATE
     I          (CNP,KKK,
     O           RSQ(1,IQ))
          ELSE
            DO I=1,MIN0(3,NCC)
              CWQ(IQ,I)=0.0D0
            ENDDO
          ENDIF
C
  250   CONTINUE
C
        DO 260 KG=1,NQ
          IF(ID.NE.2)THEN
            IF(IMID.EQ.0)THEN
              THG(KG)=W1*TH(M,KG)+W2*THP(M,KG)
            ELSE
              THG(KG)=0.5D0*(TH(M,KG)+THP(M,KG))
            ENDIF
            RHOQ(KG)=AKHC(KG,M,7)
          ELSE
            THG(KG)=AKDC(KG,MM,8)
            RHOQ(KG)=AKDC(KG,MM,7)
            XQ(KG)=XG(KG)
            YQ(KG)=YG(KG)
            ZQ(KG)=ZG(KG)
          ENDIF
  260   CONTINUE
C
C ------- COMPUTE MATRICES QA(IQ,JQ), QAA(IQ,JQ), QB(IQ,JQ), QV(IQ,JQ),
C ------- AND QC(IQ,JQ) AND THE LOAD VECTOR QR(IQ) FOR EACH ELEMENT M.
C
        CALL TQ468
     >      (QA,QAA,QB,QC,QD,QE,QV,QR,XQ,YQ,ZQ,VXQ,VYQ,VZQ,THG,CCQ,RHOQ,
     >       AKXYZQ,RHOB,LAMBDA,RHOW,DINTS, SOSQ,SOSC,SOSQC,DSDCQ,CWQ,
     >       RSQ,KKK,IRHO, NQ)
C
C ------- ASSEMBLE QA(IQ,JQ), QAA(IQ,JQ), QB(IQ,JQ)/QV(IQ,JQ), AND
C ------- QC(IQ,JQ) INTO THE GLOBAL MATRIX CMATRX(NP,I).
C ------- CMATRX(NP,I) =QB+QC+(QA+QAA)/DELT.
C ------- FORM THE LOAD VECTOR, RLD(NP) = QR + QA/DELT*CSTAR +
C ------- QAA/DELT*CP.
C
        DO 390 IQ=1,NQ
          IF(ID.NE.2)THEN
            NI=IE(M,IQ)
          ELSE
            NI=ISED(MM,IQ)
          ENDIF
          RLD(NI)=RLD(NI)+QR(IQ)
          DO 340 JQ=1,NQ
            IF(ID.NE.2)THEN
              NJ=IE(M,JQ)
            ELSE
              NJ=ISED(MM,JQ)
            ENDIF
            CPNJ=CSTAR(NJ)
C
            IF(IMID.NE.0) GO TO 305
C
C ------- FOR THE CASE OF NON MID-DIFFERENCE
C
            QA(IQ,JQ)=QA(IQ,JQ)*DTI(NI)
            QAA(IQ,JQ)=QAA(IQ,JQ)*DTI(NI)
            RLD(NI)=RLD(NI)+QA(IQ,JQ)*CSTAR(NJ)+QAA(IQ,JQ)*CPNJ
            GO TO 310
C
C ------- FOR THE CASE OF MID-DIFFERENCE
C
  305       QA(IQ,JQ)=QA(IQ,JQ)*DTI(NI)*2.0D0
            QAA(IQ,JQ)=QAA(IQ,JQ)*DTI(NI)*2.0D0
            RLD(NI)=RLD(NI) + QA(IQ,JQ)*(CSTAR(NJ)+CP(NJ,KKK))*0.5D0 +
     1              QAA(IQ,JQ)*(CPNJ+CP(NJ,KKK))*0.5D0
C
C ------- MERGE NON MID-DIFFERENCE AND MID-DIFFERENCE CASES
C
  310       CONTINUE
            IF(VTERM.EQ.0.0D0)THEN
              ADVECT=0.0D0
            ELSE
              ADVECT=VTERM*W2V*QV(IQ,JQ)*CP(NJ,KKK)
            ENDIF
            RLD(NI)=RLD(NI)-W2*(QB(IQ,JQ)+QC(IQ,JQ)-QD(IQ,JQ)-QE(IQ,JQ))
     1             *CSTAR(NJ)-ADVECT
C
            DO 325 I=1,MXJBD
              LNODE=LRN(I,NI)
              IF(NJ.EQ.LNODE) GO TO 330
  325       CONTINUE
C
            WRITE(16,1000) NI,M,JQ
            STOP
C
  330       CMATRX(NI,I)=CMATRX(NI,I)+QA(IQ,JQ)+QAA(IQ,JQ)+
     >                 W1*(QB(IQ,JQ)+QC(IQ,JQ)-QD(IQ,JQ)-QE(IQ,JQ))+
     >                 VTERM*W1V*QV(IQ,JQ)
C
  340     CONTINUE
  390   CONTINUE
  690 CONTINUE
C
C ------- INCORPORATE WELL SOURCE/SINK CONDITIONS
C
      IF(NWNP.EQ.0 .OR. ID.EQ.2) GO TO 990
      DO 790 I=1,NWNP
        NI=NPW(I)
        ITYP=IWTYP(I,KKK)
        WSSQ=WSS(ITYP,1)
        WSSC=WSS(ITYP,2)
        IF(WSSQ.LE.0.0) GO TO 790
C
        IF(IRHO.EQ.1)THEN
          RHO=0.0D0
          VOLT=0.0D0
          DO J=1,NLRL(NI)
            MP=LRL(J,NI)
            MTYP=IE(MP,9)
            CALL VOLUME
     I          (MAXNP,MAXEL,X(1,1),X(1,2),X(1,3),IE,MP,
     O           VOL)
            RHO=RHO+PROP(MTYP,7)*VOL
            VOLT=VOLT+VOL
          ENDDO
          RHOW=RHO/VOLT
C
          DO K=1,NCC
            ITYP=IWTYP(I,K)
            CNP(K)=WSS(ITYP,2)
            QR(K)=CW(NI,K)   
          ENDDO
c          CALL RHOFUN(RHOSTR,AMU,RHOW,DINTS,CNP,RHOMU)
          CALL RHOFUN(RHOSTR,AMU,RHOW,DINTS,CNP)
c          CALL RHOFUN(RHO,AMU,RHOW,DINTS,QR,RHOMU)
          CALL RHOFUN(RHO,AMU,RHOW,DINTS,QR)
          WSSQC=WSSQ*RHOSTR/RHO
        ENDIF
C
        RLD(NI)=RLD(NI)+WSSQ*WSSC
        RLD(NI)=RLD(NI)-W2*WSSQC*CSTAR(NI)
  750   DO 760 J=1,MXJBD
          LNODE=LRN(J,NI)
          IF(LNODE.EQ.NI) CMATRX(NI,J)=CMATRX(NI,J)+W1*WSSQC
  760   CONTINUE
  790 CONTINUE
C
  990 CONTINUE
C
 1000 FORMAT('1'/5X,'*** WARNING: NONE OF THE LOWER-LEFT NODE IN EQUATIO
     1N',I3,/5X,'***  IS CORRESPONDING TO ',I5,'-TH ELEMENT-S',I2,
     2'-TH NODE; STOP  ****')
C
      RETURN
      END
C
C
c
      SUBROUTINE ADVRX
c     I     (CS,IRHO,CSTAR,MAXNOD,MXDTIN,NODE,IRXN,OMET,EPS,EPSX,IBUG,X,
     I     (CS,IRHO,CSTAR,MAXNOD,MXDTIN,NODE,IRXN,OMET,EPS,EPSX,X,
     I      IE,H,XW1,XW2,XW3,MPLOCW,IGLOBL,RKD,PROP,SPP,KSP,THM,THN,
c     I      DTI,DELT,RHOMU,DINTS)
     I      DTI,DELT,DINTS)
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 KD(4)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /TINTE/ NCMt,KVIt,NITERt,NPITERt,IPNTSt,MICONF,IFLUX
      COMMON /CHEM/  MXNCC,NCC
C
      DIMENSION CS(MAXNOD,MXNCC),IE(MAXEL,11),CSTAR(MAXNP,MXNCC)
      DIMENSION RKD(MAXMAT,MXNCC),H(MAXNP),PROP(MAXMAT,MXMPPM)
c      DIMENSION RHOMU(MXNCC),DINTS(MXNCC)
      DIMENSION DINTS(MXNCC)
      DIMENSION THN(MAXNP,2,MXNCC),CNP(7),DTI(MXDTIN,MXNCC),THM(MAXEL,8)
      DIMENSION SPP(MXSPPM,MAXMAT,4),CW(7),XW3(MAXNOD)
      DIMENSION X(MAXNP,3),XW1(MAXNOD),XW2(MAXNOD),MPLOCW(MAXNOD)
      DIMENSION RSQ(4),XX(8),YY(8),ZZ(8),HQ(8),DL(8),DNX(8)
C     DIMENSION CSQ(4)
C
C ----- if iglobl=1,2 ==> use variable time interval,i.e. dti
C ----- if iglobl=0   ==> use fixed time interval, i.e. delt, only for
C ----- forward tracked fine grids and EPCOF point
C
C
C     IF(IBUG.NE.0)THEN
C       WRITE(16,7000)ITM,DELT
C     ENDIF
C
      DO 500 NP=1,NODE
C
        IF(IRXN.EQ.-1)THEN
C -------- for special cases
          IF(CS(NP,2)-2.0D0*CS(NP,1) .GE. 0.0D0)THEN
            CS(NP,2)=CS(NP,2)-2.0D0*CS(NP,1)
            CS(NP,1)=0.0D0
          ELSE
            CS(NP,1)=CS(NP,1)-0.5D0*CS(NP,2)
            CS(NP,2)=0.0D0
          ENDIF
          GOTO 500
        ENDIF
C
        DO I=1,NCC
          CW(I)=CS(NP,I)
        ENDDO
C
        XP=XW1(NP)
        YP=XW2(NP)
        ZP=XW3(NP)
        IF(IGLOBL.NE.1)THEN
          DT=DELT
          M=MPLOCW(NP)
          MTYP=IE(M,9)
          DO I=1,MIN0(3,NCC)
            KD(I)=RKD(MTYP,I)
          ENDDO
          RHOB=PROP(MTYP,1)
C
          CALL ELENOD
     I        (IE(M,5),IE(M,7),
     O         NQ,I,I)
          DO I=1,NQ
            IEM=IE(M,I)
            XX(I)=X(IEM,1)
            YY(I)=X(IEM,2)
            ZZ(I)=X(IEM,3)
C
C ***** IF THE FINE GRIDS COINSIDES WITH A GLOBAL NODE, THE CONC.
C       IS OBTAINED FROM THE CSTAR ARRAY DIRECTLY.
C
            IF(DABS(XP-XX(I)).LE.EPSX.AND. DABS(YP-YY(I)).LE.EPSX
     >         .AND. DABS(ZP-ZZ(I)).LE.EPSX)THEN
              DO K=1,NCC
                CS(NP,K)=CSTAR(IEM,K)
              ENDDO
              GOTO 500
            ENDIF
C
            IF(IRHO.EQ.1)HQ(I)=H(IEM)
          ENDDO
          CALL BASE
c     I              (XX,YY,ZZ,XP,YP,ZP,M,NQ,1,
     I              (XX,YY,ZZ,XP,YP,ZP,NQ,1,
     O               DL,DNX,DNX,DNX)
          IF(IRHO.EQ.1)THEN
            HNP=0.0D0
            DO I=1,NQ
              HNP=HNP+DL(I)*HQ(I)
            ENDDO
            CALL SPFUNC
     I        (HNP,SPP,MTYP,KSP,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
c     I         CW,1,0,RHOMU,DINTS,RHO,cnstkr,
     I         CW,1,0,DINTS,RHO,cnstkr,
     O         TH,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,RHOM)
          ELSE
            TH=0.0D0
            DO I=1,NQ
              TH=TH+DL(I)*THM(M,I)
            ENDDO
          ENDIF
        ENDIF
C
        DO 400 K=1,NITERt
c          KMAX=1
c          NPMAX=1
          DIFMAX=-1.0D38
          DO 300 KKK=1,NCC
            IF(IGLOBL.NE.1)KD(4)=RKD(MTYP,KKK)
            IF(IGLOBL.EQ.1 .OR. IGLOBL.EQ.2)DT=1.0D0/DTI(NP,KKK)
C
            DO I=1,NCC
              CNP(I)=CW(I)
            ENDDO
            DO I=1,3
              RSQ(I)=0.0D0
c             CSQ(I)=0.0D0
            ENDDO
            RSQ(4)=0.0D0
C
C ****** Assign the current concentration of microbes to CSQ array for
C        calculating removal rate
C
c           DO I=1,MIN0(3,NCC)
c             CSQ(I)=CNP(I)
c           ENDDO
C
            CALL RXRATE
     I             (CNP,KKK,
     O              RSQ)
C
c           DO I=1,MIN0(3,NCC)
c             RSQ(I)=RSQ(I)*CNP(I)
c           ENDDO
C
            IF(KKK.LE.3)THEN
c             CW(KKK)=(RSQ(KKK)*DT+CS(NP,KKK))/(1.0D0+RSQ(4)*DT)
              CW(KKK)=CS(NP,KKK)/(1.0D0-(RSQ(KKK)-RSQ(4))*DT)        
            ELSEIF(IGLOBL.EQ.1)THEN
              CW(KKK)=0.0D0
              DO I=1,3
                CW(KKK)=CW(KKK)-THN(NP,1,I)*RSQ(I)*CNP(I)
              ENDDO
              CW(KKK)=CW(KKK)*DT/THN(NP,1,KKK)
              CW(KKK)=CS(NP,KKK)/(1.0D0-CW(KKK))
            ELSE
              CW(KKK)=0.0D0
              DO I=1,3
                CW(KKK)=CW(KKK)-(TH+RHOB*KD(I))*RSQ(I)*CNP(I)
              ENDDO
              CW(KKK)=CW(KKK)*DT/(TH+RHOB*KD(4))
              CW(KKK)=CS(NP,KKK)/(1.0D0-CW(KKK))
            ENDIF
C ----- find the max. relative erros
C
            IF(CNP(KKK).NE.0.0D0)THEN
              DIF=DABS((CNP(KKK)-CW(KKK))/CNP(KKK))
              IF(DIF.GT.DIFMAX)THEN
                DIFMAX=DIF
c                KMAX=KKK
c                NPMAX=NP
              ENDIF
            ENDIF
C
C ***** UPDATE THE CONCENTRATION AFTER THE ABOVE CALCULATION
C
            CW(KKK)=OMET*CW(KKK)+(1.0D0-OMET)*CNP(KKK)
  300     CONTINUE
C
C ***** check nonlinear loop convergence
C
C         IF(IBUG.NE.0)THEN
C           WRITE(16,8010)NP,KMAX,NPMAX,DIFMAX,EPS
C         ENDIF
          IF(NITERt.EQ.1)GOTO 450
          IF(DIFMAX.LE.EPS)GOTO 450
C
  400   CONTINUE
        WRITE(16,7500)NP,K,NITERt,DIFMAX,EPS
C
  450   CONTINUE
        DO I=1,NCC
          CS(NP,I)=CW(I)
        ENDDO
C
  500 CONTINUE
 7000 FORMAT('1',' ITERATION INFORMATION IN TIME STEP =',I5,
     1 ', (DELT =',1PD12.4,')'///' TABLE OF ITERATION PARAMETERS'//6X,
     2 'ITERATION',7X,' MAX DIF',6X,'TOLERANCE')
 7500 FORMAT('0'/5X,'** WARNING: NO CONVERGENCE AT',I4,'-TH POINT AFT',
     1'ER',I4,' ITERATIONS'/8X,'NITER =',I4,' DIFMAX =',D12.4,' TOLB =',
     2 D12.4)
 8010 FORMAT(1X,3I6,1PD12.4,1PD12.4)
      RETURN
      END
C
c
c
      SUBROUTINE RXRATE
     I                 (CNP,KKK,
     O                  RSQ)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /MICROB/ GRATE,YCOEFF,RTARDS,RTARDO,RTARDN,SCOEFF,
     >                ECOEFF,DCOEFF,SATURC,PCOEFF,COFK
      COMMON /CHEM/MXNCC,NCC
C
      DIMENSION CNP(7),RSQ(4)
      DIMENSION GRATE(4),YCOEFF(4),RTARDS(4),RTARDO(4),RTARDN(4),
     >          SCOEFF(4),ECOEFF(4),DCOEFF(4),SATURC(4),PCOEFF(4)
C
C ***** calculate inhibition coefficient
C
        IF(KKK.EQ.3 .OR. KKK.EQ.4 .OR. KKK.EQ.7)THEN
          IF(COFK.EQ.0.0 .AND. CNP(5).EQ.0.0)THEN
            HICOQ=1.0D0
          ELSE
            HICOQ=COFK/(COFK+CNP(5))
          ENDIF
        ENDIF
C
C ***** Calculate removal rate of substrate under aerobic & anaerobic
C       conditions
C
C
C ------  FOR ONE COMPONENT TEST
          IF(KKK.LE.3 .AND. NCC.EQ.1)THEN
            RSQ(KKK)=0.0D0
            RETURN
          ENDIF
C
          IF(KKK.LE.3)THEN
C           for microbe #1, #2, and #3
            IF(RTARDS(KKK).EQ.0.0D0 .AND. CNP(4).EQ.0.0D0)THEN
              T1=1.0D0
            ELSE
              T1=CNP(4)/(RTARDS(KKK)+CNP(4))
            ENDIF
            IF(KKK.EQ.2)THEN
              C2=CNP(6)
            ELSE
              C2=CNP(5)
            ENDIF
            IF(RTARDO(KKK).EQ.0.0 .AND. C2.EQ.0.0)THEN
              T2=1.0D0
            ELSE
              T2=C2/(RTARDO(KKK)+C2)
            ENDIF
            IF(RTARDN(KKK).EQ.0.0 .AND. CNP(7).EQ.0.0)THEN
c           IF(RTARDN(KKK).EQ.0.0 .or. CNP(7).EQ.0.0)THEN
              T3=1.0D0
            ELSE
              T3=CNP(7)/(RTARDN(KKK)+CNP(7))
            ENDIF
            RSQ(KKK)=GRATE(KKK)*T1*T2*T3
            RSQ(4)=DCOEFF(KKK)
            IF(KKK.EQ.3)THEN
              IF(RTARDS(KKK+1).EQ.0.0 .AND. CNP(4).EQ.0.0)THEN
                T1=1.0D0
              ELSE
                T1=CNP(4)/(RTARDS(KKK+1)+CNP(4))
              ENDIF
              IF(RTARDO(KKK+1).EQ.0.0 .AND. CNP(6).EQ.0.0)THEN
                T2=1.0D0
              ELSE
                T2=CNP(6)/(RTARDO(KKK+1)+CNP(6))
              ENDIF
              IF(RTARDN(KKK+1).EQ.0.0 .AND. CNP(7).EQ.0.0)THEN
                T3=1.0D0
              ELSE
                T3=CNP(7)/(RTARDN(KKK+1)+CNP(7))
              ENDIF
              RSQ(KKK)=HICOQ*GRATE(KKK+1)*T1*T2*T3+RSQ(KKK)
              RSQ(4)=RSQ(4)+HICOQ*DCOEFF(KKK+1)
            ENDIF
          ENDIF
C
          IF(KKK.EQ.4)THEN
C           for substrate
            DO I=1,4
              IF(I.EQ.2 .OR. I.EQ.4)THEN
                C2=CNP(6)
              ELSE
                C2=CNP(5)
              ENDIF
              IF(RTARDS(I).EQ.0.0 .AND. CNP(4).EQ.0.0)THEN
c               WRITE(16,*)'ERROR WITH RETARDED SUBSTRATE SATURATION ',
c    >                     'CONSTANT = 0.0'
c               STOP
                T1=0.0D0
              ELSE
                T1=1.0D0/(RTARDS(I)+CNP(4))
              ENDIF
              IF(RTARDO(I).EQ.0.0 .AND. C2.EQ.0.0)THEN
                T2=1.0D0
              ELSE
                T2=C2/(RTARDO(I)+C2)
              ENDIF
              IF(RTARDN(I).EQ.0.0 .AND. CNP(7).EQ.0.0)THEN
c             IF(RTARDN(I).EQ.0.0 .or. CNP(7).EQ.0.0)THEN
                T3=1.0D0
              ELSE
                T3=CNP(7)/(RTARDN(I)+CNP(7))
              ENDIF
              IF(YCOEFF(I).EQ.0.0)THEN
                WRITE(16,*)'ERROR WITH YIELD COEFF = 0.0'
                STOP
              ENDIF
              RSQ(I)=GRATE(I)/YCOEFF(I)*T1*T2*T3
            ENDDO
            RSQ(3)=RSQ(3)+RSQ(4)*HICOQ
            RSQ(4)=0.0D0
          ENDIF
C
          IF(KKK.EQ.5 .OR. KKK.EQ.6)THEN
            DO I=KKK-4,3,7-KKK
              IF(I.EQ.3 .AND. KKK.EQ.6)THEN
                II=4
              ELSE
                II=I
              ENDIF
              IF(RTARDS(II).EQ.0.0 .AND. CNP(4).EQ.0.0)THEN
                T1=1.0D0
              ELSE
                T1=CNP(4)/(RTARDS(II)+CNP(4))
              ENDIF
              IF(RTARDO(II).EQ.0.0 .AND. CNP(KKK).EQ.0.0)THEN
c               WRITE(16,*)'ERROR WITH RETARDED OXYGEN OR NITRATE ',
c    >                     'SATURATION COEFF = 0.0'
c               STOP
                T2=0.0D0
              ELSE
                T2=1.0D0/(RTARDO(II)+CNP(KKK))
c               T2=1.0D0
              ENDIF
c             IF(RTARDN(II).EQ.0.0 .or. CNP(7).EQ.0.0)THEN
              IF(RTARDN(II).EQ.0.0 .AND. CNP(7).EQ.0.0)THEN
                T3=1.0D0
              ELSE
                T3=CNP(7)/(RTARDN(II)+CNP(7))
              ENDIF
              IF(SATURC(II).EQ.0.0 .AND. CNP(KKK).EQ.0.0)THEN
                IF(ECOEFF(II)*DCOEFF(II).EQ.0.0)THEN
                  T4=1.0D0
                ELSE
                  WRITE(16,*)'ERROR WITH OXYGEN OR NITRATE SATURATION ',
     >                       'CONST FOR DECAY = 0.0'
                  STOP
                ENDIF
              ELSE
                T4=ECOEFF(II)*DCOEFF(II)/(SATURC(II)+CNP(KKK))
              ENDIF
              RSQ(I)=SCOEFF(II)*GRATE(II)*T1*T2*T3+T4
            ENDDO
          ENDIF
          IF(KKK.EQ.7)THEN
C           for neutrient
            DO I=1,4
              IF(RTARDS(I).EQ.0.0 .AND. CNP(4).EQ.0.0)THEN
                T1=1.0D0
              ELSE
                T1=CNP(4)/(RTARDS(I)+CNP(4))
              ENDIF
              IF(I.EQ.2 .OR. I.EQ.4)THEN
                C2=CNP(6)
              ELSE
                C2=CNP(5)
              ENDIF
              IF(RTARDO(I).EQ.0.0  .AND. C2.EQ.0.0)THEN
                T2=1.0D0
              ELSE
                T2=C2/(RTARDO(I)+C2)
              ENDIF
              IF(RTARDN(I).EQ.0.0 .AND. CNP(7).EQ.0.0)THEN
c               WRITE(16,*)'ERROR WITH RETARDED NUTRIEN SATURATION ',
c    >                     'CONSTANT = 0.0'
c               STOP
                T3=0.0D0
              ELSE
                T3=1.0D0/(RTARDN(I)+CNP(7))
c               T3=1.0D0
              ENDIF
              IF(YCOEFF(I).EQ.0.0)THEN
                WRITE(16,*)'ERROR WITH YIELD COEFF = 0.0'
                STOP
              ELSE
                RSQ(I)=PCOEFF(I)*GRATE(I)/YCOEFF(I)*T1*T2*T3
              ENDIF
            ENDDO
            RSQ(3)=RSQ(3)+RSQ(4)*HICOQ
            RSQ(4)=0.0D0
          ENDIF
          RETURN
          END
C
c
c
      SUBROUTINE TQ468
     >     (QA,QAA,QB,QC,QD,QE,QV,QR, XQ,YQ,ZQ,VXQ,VYQ,VZQ,THG,CCQ,RHOQ,
     >      AKXYZQ,RHOB,LAMBDA,RHOW,DINTS,SOSQ,SOSC,SOSQC,DSDCQ,CWQ,RSQ,
     >      KKK,IRHO, NODE)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE ELEMENT MATRICES AND ELEMENT LOAD VECTORS.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: XQ(8), YQ(8), ZQ(8), VXQ(8), VYQ(8), VZQ(8), THG(8),
C -------        AKXYZQ(8,6), RHOB, LAMBDA, SOSQ,SOSC,
C -------        DSDCQ(8,4), RSQ(8,4), AND CWQ(8,3).
C
C ------- OUTPUT: QA(8,8), QAA(8,8), QB(8,8), QC(8,8), QV(8,8),
C -------         AND QR(8).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(8),LAMBDA,L1(3),L2(3),L3(3),LL1(4),LL2(4),LL3(4),LL4(4)
C
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
      COMMON /CHEM/ MXNCC,NCC
C
      DIMENSION QA(8,8),QAA(8,8),QB(8,8),QC(8,8),QD(8,8),QE(8,8),
     >          QV(8,8),QR(8)
      DIMENSION XQ(8),YQ(8),ZQ(8),VXQ(8),VYQ(8),VZQ(8),THG(8)
      DIMENSION AKXYZQ(8,6), AKXYZK(6),RHOQ(8),DINTS(MXNCC)
      DIMENSION DSDCQ(8,4),RSQ(4,8),CWQ(8,3),DSDCK(4),CCQ(8,7)
C
      DIMENSION DNX(8),DNY(8),DNZ(8), W(8),WG(8,3)
      DIMENSION S(8),T(8),U(8),RSS(4),CC(3)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0, -1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA T/-1.0D0,-1.0D0,1.0D0,1.0D0, -1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA U/-1.0D0,-1.0D0,-1.0D0,-1.0D0, 1.0D0,1.0D0,1.0D0,1.0D0/
C
      DATA WG/1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,
     >        0.333333333333333D0,0.333333333333333D0,
     >        0.333333333333333D0,0.333333333333333D0,
     >        0.333333333333333D0,0.333333333333333D0,0.0D0,0.0D0,
     >        0.25D0,0.25D0,0.25D0,0.25D0,0.0D0,0.0D0,0.0D0,0.0D0/
C
      IF(IQUAR.EQ.1 .OR. IQUAR.EQ.3)THEN
        P=0.577350269189626D0
        L1(1)=0.666666666666667D0
        L1(2)=0.166666666666667D0
        L1(3)=0.166666666666667D0
        L2(1)=0.166666666666667D0
        L2(2)=0.666666666666667D0
        L2(3)=0.166666666666667D0
        L3(1)=0.166666666666667D0
        L3(2)=0.166666666666667D0
        L3(3)=0.666666666666667D0
        LL1(1)=0.58541020D0
        LL1(2)=0.13819660D0
        LL1(3)=0.13819660D0
        LL1(4)=0.13819660D0
        LL2(1)=0.13819660D0
        LL2(2)=0.58541020D0
        LL2(3)=0.13819660D0
        LL2(4)=0.13819660D0
        LL3(1)=0.13819660D0
        LL3(2)=0.13819660D0
        LL3(3)=0.58541020D0
        LL3(4)=0.13819660D0
        LL4(1)=0.13819660D0
        LL4(2)=0.13819660D0
        LL4(3)=0.13819660D0
        LL4(4)=0.58541020D0
      ELSE
        P=1.0D0
        L1(1)=1.0D0
        L1(2)=0.0D0
        L1(3)=0.0D0
        L2(1)=0.0D0
        L2(2)=1.0D0
        L2(3)=0.0D0
        L3(1)=0.0D0
        L3(2)=0.0D0
        L3(3)=1.0D0
        LL1(1)=1.0D0
        LL1(2)=0.0D0
        LL1(3)=0.0D0
        LL1(4)=0.0D0
        LL2(1)=0.0D0
        LL2(2)=1.0D0
        LL2(3)=0.0D0
        LL2(4)=0.0D0
        LL3(1)=0.0D0
        LL3(2)=0.0D0
        LL3(3)=1.0D0
        LL3(4)=0.0D0
        LL4(1)=0.0D0
        LL4(2)=0.0D0
        LL4(3)=0.0D0
        LL4(4)=1.0D0
      ENDIF
C
C
C ------- INITIATE MATRICES QA, QAA, QB, QV, QC, AND QR
C
      DO 110 IQ=1,8
      QR(IQ)=0.0
      DO 110 JQ=1,8
      QA(IQ,JQ)=0.0
      QAA(IQ,JQ)=0.0
      QB(IQ,JQ)=0.0
      QC(IQ,JQ)=0.0
      QD(IQ,JQ)=0.0
      QE(IQ,JQ)=0.0D0
      QV(IQ,JQ)=0.0
  110 CONTINUE
C
      DO 490 KG=1,NODE
C
C ------- DETERMINE LOACAL COORDINATE  OF GAUSSIAN POINT KG
C
        IF(NODE.EQ.8)THEN
          ID=1
          SS=P*S(KG)
          TT=P*T(KG)
          UU=P*U(KG)
        ELSEIF(NODE.EQ.6)THEN
          ID=2
          XSI=-P
          ND=KG
          IF(ND.GT.3)THEN
            ND=KG-3
            XSI=P
          ENDIF
          DL1=L1(ND)
          DL2=L2(ND)
          DL3=L3(ND)
        ELSE
          ID=3
          D1=LL1(KG)
          D2=LL2(KG)
          D3=LL3(KG)
          D4=LL4(KG)
        ENDIF
C
C ------- CALCULATE VALUES OF N(IQ), DNX(IQ), DNY(IQ), DNZ(IQ), W(IQ),
C ------- DWX(IQ), DWY(IQ), DWZ(IQ), AND DJAC.
C
        CALL SHAPE
     >      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,1,
     >       N,DNX,DNY,DNZ,W,DJAC)
C
        DO 210 I=1,6
          AKXYZK(I)=AKXYZQ(KG,I)
  210   CONTINUE
C
        VXK=0.0
        VYK=0.0
        VZK=0.0
        DO I=1,3
          CC(I)=0.0D0
        ENDDO
        DO I=1,4
          RSS(I)=0.0D0
          DSDCK(I)=0.0
        ENDDO
C
        DO 215 IQ=1,NODE
          VXK=VXK+VXQ(IQ)*N(IQ)
          VYK=VYK+VYQ(IQ)*N(IQ)
          VZK=VZK+VZQ(IQ)*N(IQ)
          DO I=1,4
            DSDCK(I)=DSDCK(I)+DSDCQ(IQ,I)*N(IQ)
          ENDDO
          DO I=1,3
            RSS(I)=RSS(I)+N(IQ)*RSQ(I,IQ)
            CC(I)=CC(I)+N(IQ)*CWQ(IQ,I)
          ENDDO
  215   CONTINUE
C
        THK=THG(KG)
        RHOK=RHOQ(KG)
C
        IF(IRHO.EQ.1)THEN
C ***** COMPUTE V.GRAD(RHO/RHOW)
C
c         COFRHO=1.0D0
C ----- the following block is for testing salt water intrusion problem
c         COFRHO=0.0245D0
C ------------------------
C
C ------ GRAD(RHO/RHOW) IS BASED ON GRAD(CP)
C
          EEX=0.0D0
          EEY=0.0D0
          EEZ=0.0D0
          DO 440 K=1,NCC
            COFRHO=RHOW-1.0D0/DINTS(K)
            DO IQ=1,NODE
              EEX=EEX+DNX(IQ)*CCQ(IQ,K)
              EEY=EEY+DNY(IQ)*CCQ(IQ,K)
              EEZ=EEZ+DNZ(IQ)*CCQ(IQ,K)
            ENDDO
            EEX=EEX*COFRHO
            EEY=EEY*COFRHO
            EEZ=EEZ*COFRHO
  440     CONTINUE
          E=(EEX*VXK+EEY*VYK+EEZ*VZK)/RHOK
        ELSE
          E=0.0D0
        ENDIF
C
        DJAC=DJAC*WG(KG,ID)
C
        DXX=DJAC*AKXYZK(1)
        DYY=DJAC*AKXYZK(2)
        DZZ=DJAC*AKXYZK(3)
        DXY=DJAC*AKXYZK(4)
        DXZ=DJAC*AKXYZK(5)
        DYZ=DJAC*AKXYZK(6)
C
        VXK=VXK*DJAC
        VYK=VYK*DJAC
        VZK=VZK*DJAC
C
        SOSQK=SOSQC*DJAC
        IF(SOSQ.LT.0.0) THEN
          SOSQK=0.0
          SOSQ=0.0D0
        ENDIF
C
        A=DJAC*THK
        AA=DJAC*RHOB*DSDCK(4)
        B=SOSQ*SOSC*DJAC
        C=DJAC*(LAMBDA*(THK+RHOB*DSDCK(4)))+SOSQK/RHOK
        E=E*DJAC
        DO I=1,3
          RSS(I)=CC(I)*RSS(I)
        ENDDO
C
        DO 390 IQ=1,NODE
          QR(IQ)=QR(IQ)+B*N(IQ)
          IF(KKK.LE.3)THEN
            QR(IQ)=QR(IQ)+DJAC*(THK+RHOB*DSDCK(KKK))*RSS(KKK)*N(IQ)
          ENDIF
          DO 350 JQ=1,NODE
            WN=N(IQ)*N(JQ)
            DWXDNX=DNX(IQ)*DNX(JQ)
            DWXDNY=DNX(IQ)*DNY(JQ)
            DWXDNZ=DNX(IQ)*DNZ(JQ)
            DWYDNX=DNY(IQ)*DNX(JQ)
            DWYDNY=DNY(IQ)*DNY(JQ)
            DWYDNZ=DNY(IQ)*DNZ(JQ)
            DWZDNX=DNZ(IQ)*DNX(JQ)
            DWZDNY=DNZ(IQ)*DNY(JQ)
            DWZDNZ=DNZ(IQ)*DNZ(JQ)
            WDNX=W(IQ)*DNX(JQ)
            WDNY=W(IQ)*DNY(JQ)
            WDNZ=W(IQ)*DNZ(JQ)
            QA(IQ,JQ)=QA(IQ,JQ) + A*WN
            QAA(IQ,JQ)=QAA(IQ,JQ) + AA*WN
            QB(IQ,JQ)=QB(IQ,JQ)+DWXDNX*DXX+(DWXDNY+DWYDNX)*DXY+
     >                DWYDNY*DYY+(DWYDNZ+DWZDNY)*DYZ+DWZDNZ*DZZ+
     >                (DWXDNZ+DWZDNX)*DXZ
            QC(IQ,JQ)=QC(IQ,JQ) + C*WN
            QD(IQ,JQ)=QD(IQ,JQ)-DJAC*(THK+RHOB*DSDCK(4))*RSS(4)*WN
            QE(IQ,JQ)=QE(IQ,JQ) + E*WN
            IF(KKK.GT.3)THEN
              DO I=1,3
                QD(IQ,JQ)=QD(IQ,JQ)-DJAC*(THK+RHOB*DSDCK(I))*
     >                  RSS(I)*WN
              ENDDO
            ENDIF
            QV(IQ,JQ)=QV(IQ,JQ) + (VXK*WDNX+VYK*WDNY+VZK*WDNZ)
  350     CONTINUE
  390   CONTINUE
C
  490   CONTINUE
C
      IF(ILUMP.EQ.0) RETURN
C
      DO 940 I=1,NODE
        SUM=0.0
        SUMAA=0.0
        SUMC=0.0
        SUMD=0.0
        SUME=0.0D0
        DO 920 J=1,NODE
          SUM=SUM+QA(I,J)
          SUMAA=SUMAA+QAA(I,J)
          SUMC=SUMC+QC(I,J)
          SUMD=SUMD+QD(I,J)
          SUME=SUME+QE(I,J)
          QA(I,J)=0.0
          QAA(I,J)=0.0
          QC(I,J)=0.0
          QD(I,J)=0.0D0
          QE(I,J)=0.0D0
  920   CONTINUE
        QA(I,I)=SUM
        QAA(I,I)=SUMAA
        QC(I,I)=SUMC
        QD(I,I)=SUMD
        QE(I,I)=SUME
  940 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,ID,
     O       N,DNX,DNY,DNZ,W,DJAC)
C
C ----- 1/26/93(TESTED)
C
C $$$$$ TO COMPUTE THE BASIS AND WEIGHTING FUNCTIONS, THEIR
C       DERIVATIVES WITH RESPECT TO X, Y, Z, AND THE JACOBIAN
C       AT A GAUSSIAN POINT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(8)
      REAL*8 J1,J2,J3,J4,J5,J6,J7,J8,J9
      REAL*8 J11,J12,J13,J21,J22,J23,J31,J32,J33
C
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
      COMMON /WETX/ APHA1,APHA2,APHA3,APHA4
      COMMON /WETY/ BETA1,BETA2,BETA3,BETA4
      COMMON /WETZ/ GAMA1,GAMA2,GAMA3,GAMA4
C
      DIMENSION DNX(8),DNY(8),DNZ(8),W(8), XQ(8),YQ(8),ZQ(8)
      DIMENSION DNSS(8),DNTT(8),DNUU(8),A(4),B(4),C(4),D(4)
C
C ------- DEFINE UPSTREAM WEIGHTING FUNCTIONS AND THEIR DERIVATIVES
C
      WF1(X,AL)=0.25D0*((2.0D0-3.0D0*AL)-2.0D0*X+3.0D0*AL*X*X)
      WF2(X,AL)=0.25D0*((2.0D0+3.0D0*AL)+2.0D0*X-3.0D0*AL*X*X)
C
      WF3(D1,D2,D3,D4,AL2,AL3,AL4)=D1+3.0D0*D1*(D2*AL2+D3*AL3+D4*AL4)
C
      WF4(D1,D2,D3,AL2,AL3)=D1+3.0D0*D1*(D2*AL2+D3*AL3)
C
C ******* COMPUTE BASIS FUNCTIONS AND THEIR DERIVATIVES
C
      IF(NODE.EQ.8)THEN
C
C ----- FOR HEXAHEDRAL ELEMENTS
C
        SM=1.0D0-SS
        SP=1.0D0+SS
        TM=1.0D0-TT
        TP=1.0D0+TT
        UM=1.0D0-UU
        UP=1.0D0+UU
C
        N(1)=0.125D0*SM*TM*UM
        N(2)=0.125D0*SP*TM*UM
        N(3)=0.125D0*SP*TP*UM
        N(4)=0.125D0*SM*TP*UM
        N(5)=0.125D0*SM*TM*UP
        N(6)=0.125D0*SP*TM*UP
        N(7)=0.125D0*SP*TP*UP
        N(8)=0.125D0*SM*TP*UP
C
        DNSS(1)=-.125D0*TM*UM
        DNSS(2)=.125D0*TM*UM
        DNSS(3)=.125D0*TP*UM
        DNSS(4)=-.125D0*TP*UM
        DNSS(5)=-.125D0*TM*UP
        DNSS(6)=.125D0*TM*UP
        DNSS(7)=.125D0*TP*UP
        DNSS(8)=-.125D0*TP*UP
C
        DNTT(1)=-.125D0*SM*UM
        DNTT(2)=-.125D0*SP*UM
        DNTT(3)=.125D0*SP*UM
        DNTT(4)=.125D0*SM*UM
        DNTT(5)=-.125D0*SM*UP
        DNTT(6)=-.125D0*SP*UP
        DNTT(7)=.125D0*SP*UP
        DNTT(8)=.125D0*SM*UP
C
        DNUU(1)=-.125D0*SM*TM
        DNUU(2)=-.125D0*SP*TM
        DNUU(3)=-.125D0*SP*TP
        DNUU(4)=-.125D0*SM*TP
        DNUU(5)= .125D0*SM*TM
        DNUU(6)= .125D0*SP*TM
        DNUU(7)= .125D0*SP*TP
        DNUU(8)= .125D0*SM*TP
C
        SUM1=0.0
        SUM2=0.0
        SUM3=0.0
        SUM4=0.0
        SUM5=0.0
        SUM6=0.0
        SUM7=0.0
        SUM8=0.0
        SUM9=0.0
C
        DO 290 I=1,8
          SUM1=SUM1+XQ(I)*DNSS(I)
          SUM2=SUM2+YQ(I)*DNSS(I)
          SUM3=SUM3+ZQ(I)*DNSS(I)
          SUM4=SUM4+XQ(I)*DNTT(I)
          SUM5=SUM5+YQ(I)*DNTT(I)
          SUM6=SUM6+ZQ(I)*DNTT(I)
          SUM7=SUM7+XQ(I)*DNUU(I)
          SUM8=SUM8+YQ(I)*DNUU(I)
          SUM9=SUM9+ZQ(I)*DNUU(I)
  290   CONTINUE
C
        DJAC=SUM1*(SUM5*SUM9-SUM6*SUM8)+SUM2*(SUM6*SUM7-SUM4*SUM9)+
     1       SUM3*(SUM4*SUM8-SUM5*SUM7)
C
        DJACI=1.0D0/DJAC
C
        SUMI1=DJACI*(SUM5*SUM9-SUM6*SUM8)
        SUMI2=DJACI*(SUM3*SUM8-SUM2*SUM9)
        SUMI3=DJACI*(SUM2*SUM6-SUM3*SUM5)
        SUMI4=DJACI*(SUM6*SUM7-SUM4*SUM9)
        SUMI5=DJACI*(SUM1*SUM9-SUM3*SUM7)
        SUMI6=DJACI*(SUM3*SUM4-SUM1*SUM6)
        SUMI7=DJACI*(SUM4*SUM8-SUM5*SUM7)
        SUMI8=DJACI*(SUM2*SUM7-SUM1*SUM8)
        SUMI9=DJACI*(SUM1*SUM5-SUM2*SUM4)
C
        DO 390 I=1,8
          DNX(I)=SUMI1*DNSS(I)+SUMI2*DNTT(I)+SUMI3*DNUU(I)
          DNY(I)=SUMI4*DNSS(I)+SUMI5*DNTT(I)+SUMI6*DNUU(I)
          DNZ(I)=SUMI7*DNSS(I)+SUMI8*DNTT(I)+SUMI9*DNUU(I)
  390   CONTINUE
C
      ELSEIF(NODE.EQ.6)THEN
C
C ----- FOR PENTAHEDRAL ELEMENTS
C
        SM=1.0D0-XSI
        SP=1.0D0+XSI
        N(1)=0.5D0*SM*DL1
        N(2)=0.5D0*SM*DL2
        N(3)=0.5D0*SM*DL3
        N(4)=0.5D0*SP*DL1
        N(5)=0.5D0*SP*DL2
        N(6)=0.5D0*SP*DL3
C
        DNSS(1)=-.5D0*DL1
        DNSS(2)=-.5D0*DL2
        DNSS(3)=-.5D0*DL3
        DNSS(4)=.5D0*DL1
        DNSS(5)=.5D0*DL2
        DNSS(6)=.5D0*DL3
C
        DNTT(1)=.5D0*SM
        DNTT(2)=.0D0
        DNTT(3)=-.5D0*SM
        DNTT(4)=.5D0*SP
        DNTT(5)=.0D0
        DNTT(6)=-.5D0*SP
C
        DNUU(1)=.0D0
        DNUU(2)=.5D0*SM
        DNUU(3)=-.5D0*SM
        DNUU(4)=.0D0
        DNUU(5)=.5D0*SP
        DNUU(6)=-.5D0*SP
C
        X3416=XQ(3)+XQ(4)-XQ(1)-XQ(6)
        X3526=XQ(3)+XQ(5)-XQ(2)-XQ(6)
        X1436=XQ(1)+XQ(4)-XQ(3)-XQ(6)
        X2536=XQ(2)+XQ(5)-XQ(3)-XQ(6)
        Y3416=YQ(3)+YQ(4)-YQ(1)-YQ(6)
        Y3526=YQ(3)+YQ(5)-YQ(2)-YQ(6)
        Y1436=YQ(1)+YQ(4)-YQ(3)-YQ(6)
        Y2536=YQ(2)+YQ(5)-YQ(3)-YQ(6)
        Z3416=ZQ(3)+ZQ(4)-ZQ(1)-ZQ(6)
        Z3526=ZQ(3)+ZQ(5)-ZQ(2)-ZQ(6)
        Z1436=ZQ(1)+ZQ(4)-ZQ(3)-ZQ(6)
        Z2536=ZQ(2)+ZQ(5)-ZQ(3)-ZQ(6)
C
        J11=0.5D0*(XQ(6)-XQ(3)+DL1*X3416+DL2*X3526)
        J12=0.5D0*(YQ(6)-YQ(3)+DL1*Y3416+DL2*Y3526)
        J13=0.5D0*(ZQ(6)-ZQ(3)+DL1*Z3416+DL2*Z3526)
        J21=0.5D0*(X1436+XSI*X3416)
        J22=0.5D0*(Y1436+XSI*Y3416)
        J23=0.5D0*(Z1436+XSI*Z3416)
        J31=0.5D0*(X2536+XSI*X3526)
        J32=0.5D0*(Y2536+XSI*Y3526)
        J33=0.5D0*(Z2536+XSI*Z3526)
C
        DJAC=J11*J22*J33+J12*J23*J31+J13*J21*J32-
     >       J11*J32*J23-J12*J21*J33-J13*J31*J22
C
        DJACI=1.0D0/DJAC
C
        J1=DJACI*(J22*J33-J23*J32)
        J2=-DJACI*(J12*J33-J13*J32)
        J3=DJACI*(J12*J23-J13*J22)
        J4=-DJACI*(J21*J33-J23*J31)
        J5=DJACI*(J11*J33-J13*J31)
        J6=-DJACI*(J11*J23-J13*J21)
        J7=DJACI*(J21*J32-J22*J31)
        J8=-DJACI*(J11*J32-J12*J31)
        J9=DJACI*(J11*J22-J12*J21)
C
        DO 490 I=1,6
          DNX(I)=J1*DNSS(I)+J2*DNTT(I)+J3*DNUU(I)
          DNY(I)=J4*DNSS(I)+J5*DNTT(I)+J6*DNUU(I)
          DNZ(I)=J7*DNSS(I)+J8*DNTT(I)+J9*DNUU(I)
  490   CONTINUE
        DJAC=DABS(DJAC)*0.5D0
C
      ELSE
C
C ----- FOR TETRAHEDRAL ELEMENTS
C
        DJAC=0.0
        DO 550 KK=1,4
          IF(KK.EQ.1)THEN
            K1=2
            K2=3
            K3=4
          ELSEIF(KK.EQ.2)THEN
            K1=1
            K2=3
            K3=4
          ELSEIF(KK.EQ.3)THEN
            K1=1
            K2=2
            K3=4
          ELSE
            K1=1
            K2=2
            K3=3
          ENDIF
          A(KK)=(-1.0D0)**(KK+1)*(XQ(K1)*YQ(K2)*ZQ(K3)+
     1          YQ(K1)*ZQ(K2)*XQ(K3)+ZQ(K1)*XQ(K2)*YQ(K3)-
     2          XQ(K3)*YQ(K2)*ZQ(K1)-YQ(K3)*ZQ(K2)*XQ(K1)-
     3          ZQ(K3)*XQ(K2)*YQ(K1))
          B(KK)=(-1.0D0)**KK*(YQ(K1)*ZQ(K2)+
     1          YQ(K2)*ZQ(K3)+YQ(K3)*ZQ(K1)-
     2          YQ(K3)*ZQ(K2)-YQ(K2)*ZQ(K1)-
     3          YQ(K1)*ZQ(K3))
          C(KK)=(-1.0D0)**(KK+1)*(XQ(K1)*ZQ(K2)+
     1          XQ(K2)*ZQ(K3)+XQ(K3)*ZQ(K1)-
     2          XQ(K3)*ZQ(K2)-XQ(K2)*ZQ(K1)-
     3          XQ(K1)*ZQ(K3))
          D(KK)=(-1.0D0)**KK*(XQ(K1)*YQ(K2)+
     1          XQ(K2)*YQ(K3)+XQ(K3)*YQ(K1)-
     2          XQ(K3)*YQ(K2)-XQ(K2)*YQ(K1)-
     3          XQ(K1)*YQ(K3))
          DJAC=DJAC+A(KK)
  550   CONTINUE
        N(1)=D1
        N(2)=D2
        N(3)=D3
        N(4)=D4
C
        DO 590 KK=1,4
          DNX(KK)=B(KK)/DJAC
          DNY(KK)=C(KK)/DJAC
          DNZ(KK)=D(KK)/DJAC
  590   CONTINUE
        DJAC=DABS(DJAC)/6.0D0
      ENDIF
C
C ----- CHECK IF WEIGHTING FUNCTIONS ARE TO BE COMPUTED
C
      IF(ID.EQ.0)GOTO 800
C
C ******* COMPUTE WEIGHTING FUNCTINS AND THERI DERIVATIVES
C
      IF(IWET.EQ.0) GO TO 700
C
C ------- FOR THE CASE OF UPSTREAM WEIGHTING
C
      IF(NODE.EQ.8)THEN
        W(1)=WF1(SS,APHA1)*WF1(TT,BETA1)*WF1(UU,GAMA1)
        W(2)=WF2(SS,APHA1)*WF1(TT,BETA2)*WF1(UU,GAMA2)
        W(3)=WF2(SS,APHA2)*WF2(TT,BETA2)*WF1(UU,GAMA3)
        W(4)=WF1(SS,APHA2)*WF2(TT,BETA1)*WF1(UU,GAMA4)
        W(5)=WF1(SS,APHA3)*WF1(TT,BETA3)*WF2(UU,GAMA1)
        W(6)=WF2(SS,APHA3)*WF1(TT,BETA4)*WF2(UU,GAMA2)
        W(7)=WF2(SS,APHA4)*WF2(TT,BETA4)*WF2(UU,GAMA3)
        W(8)=WF1(SS,APHA4)*WF2(TT,BETA3)*WF2(UU,GAMA4)
C       WRITE(16,999)(W(I),I=1,8)
  999   FORMAT(8E10.3)
      ELSEIF(NODE.EQ.6)THEN
        W(1)=WF4(DL1,DL2,DL3,-APHA1,-APHA3)*WF1(XSI,GAMA1)
        W(2)=WF4(DL2,DL1,DL3,APHA1,-APHA2)*WF1(XSI,GAMA2)
        W(3)=WF4(DL3,DL1,DL2,APHA3,APHA2)*WF1(XSI,GAMA3)
        W(4)=WF4(DL1,DL2,DL3,-BETA1,-BETA3)*WF2(XSI,GAMA1)
        W(5)=WF4(DL2,DL1,DL3,BETA1,-BETA2)*WF2(XSI,GAMA2)
        W(6)=WF4(DL3,DL1,DL2,BETA3,BETA2)*WF2(XSI,GAMA3)
      ELSE
        D1=N(1)
        D2=N(2)
        D3=N(3)
        D4=N(4)
        W(1)=WF3(D1,D2,D3,D4,-APHA1,-APHA2,-APHA3)
        W(2)=WF3(D2,D1,D3,D4,APHA1,-BETA1,-BETA2)
        W(3)=WF3(D3,D1,D2,D4,APHA2,BETA1,-BETA3)
        W(4)=WF3(D4,D1,D2,D3,APHA3,BETA2,BETA3)
      ENDIF
C
      GOTO 800
C
C ------- FOR THE CASE OF GALERKIN WEIGHTING.
C
  700 DO 790 I=1,NODE
        W(I)=N(I)
  790 CONTINUE
C
  800 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE TBC
     >     (CMATRX,RLD,CSTAR,X,IE,LRN,NLRN,DCOSB,ISB,
     1      V,VP, QCB,ISC,ICTYP, QNB,ISN,INTYP,
     2      CVB,ISV,IVTYP, CDB,IDTYP,NPDB, W,KSS,ID,MXTNOD)
C 11/24/93
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO APPLY CAUCHY, NEUMANN, VARIABLE, AND DIRICHLET BOUNDARY
C ------- CONDITIONS.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: CSTAR(NNP),X(NNP,3),IE(NEL,9),
C -------        LRN(MXJBD,NNP),DCOSB(3,NBES),ISB(6,NBES),
C -------        V(NNP,3),VP(NNP,3),
C -------        QCB(NCPR),ISC(5,NCES),ICTYP(NCES),
C -------        QNB(NNPR),ISN(5,NNES),INTYP(NNES),  CVB(NVPR),
C -------        ISV(5,NVES),IVTYP(NVES),  CDB(NDPR),IDTYP(NDNP),
C -------        NPDB(NDNP).
C
C ------- OUTPUT: CMATRX(NNP,MXJBD) AND RLD(NNP).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /LGEOM/ LTMXNP,LMXNP,LMXBW,MXREGN,NREGN
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
C
      COMMON /TCBC/ MXCNP,MXCES,MXCPR,MXCDP,NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /TNBC/ MXNNP,MXNES,MXNPR,MXNDP,NNNP,NNES,NNPR,NNDP,KNAI
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /TDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
C
      DIMENSION CMATRX(MXADNP,MXJBD),RLD(MXADNP), CSTAR(MAXNP)
      DIMENSION X(MAXNP,3),IE(MAXEL,11),LRN(MXJBD,MXADNP),NLRN(MXTNOD)
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES)
C
      DIMENSION V(MAXNP,3),VP(MAXNP,3)
C
      DIMENSION QCB(MXCPR),ICTYP(MXCES),ISC(5,MXCES)
      DIMENSION QNB(MXNPR),INTYP(MXNES),ISN(5,MXNES)
      DIMENSION CVB(MXVPR),IVTYP(MXVES),ISV(5,MXVES)
      DIMENSION CDB(MXDPR),IDTYP(MXDNP),NPDB(MXDNP)
C
      DIMENSION BQ(4,4),RQ(4)
      DIMENSION XQ(4),YQ(4),ZQ(4),VXQ(4),VYQ(4),VZQ(4)
      DIMENSION KGB(4,6,3)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
      W1=W
      W2=1.0D0-W
      IF(KSS.EQ.0 .OR. IMID.EQ.1) THEN
        W1=1.0D0
        W2=0.0
      ENDIF
C
C ******* APPLY CAUCHY CONDITION: QC=V.N.C - N.(THETA)D.GRAD(C)
C
      IF(NCES.EQ.0) GO TO 300
C
      DO 190 MP=1,NCES
C
        ITYP=ICTYP(MP)
        QCBMP=QCB(ITYP)
        MPB=ISC(5,MP)
        LS=ISB(5,MPB)
        M=ISB(6,MPB)
        CALL ELENOD(IE(M,5),IE(M,7),NODE,IQ,ID)
C
        NODE=4
C
        DO 130 IQ=1,4
          I=KGB(IQ,LS,ID)
          IF(I.EQ.0)THEN
            NODE=3
            GOTO 130
          ENDIF
          NI=IE(M,I)
          XQ(IQ)=X(NI,1)
          YQ(IQ)=X(NI,2)
          ZQ(IQ)=X(NI,3)
          IF(IMID.EQ.0)THEN
            VXQ(IQ)=W1*V(NI,1)+W2*VP(NI,1)
            VYQ(IQ)=W1*V(NI,2)+W2*VP(NI,2)
            VZQ(IQ)=W1*V(NI,3)+W2*VP(NI,3)
          ELSE
            VXQ(IQ)=0.5D0*(V(NI,1)+VP(NI,1))
            VYQ(IQ)=0.5D0*(V(NI,2)+VP(NI,2))
            VZQ(IQ)=0.5D0*(V(NI,3)+VP(NI,3))
          ENDIF
  130   CONTINUE
C
        CALL Q34CNV(BQ,RQ,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB(1,MPB),QCBMP,1,
     >              NODE,IQUAR)
C
        DO 180 IQ=1,NODE
          I=KGB(IQ,LS,ID)
          NI=IE(M,I)
          RLD(NI)=RLD(NI) + RQ(IQ)
          DO 160 JQ=1,NODE
            J=KGB(JQ,LS,ID)
            NJ=IE(M,J)
            RLD(NI)=RLD(NI)-W2*CSTAR(NJ)
            DO 140 JJ=1,MXJBD
              LNODE=LRN(JJ,NI)
              IF(LNODE.EQ.NJ) GO TO 150
  140       CONTINUE
C
            WRITE (6,1000) MP,IQ,NI,JQ,NJ
            STOP
C
  150       CMATRX(NI,JJ)=CMATRX(NI,JJ) + W1*BQ(IQ,JQ)
  160     CONTINUE
  180   CONTINUE
C
  190 CONTINUE
C
C ******* APPLY NEUMANN CONDITION: QN= - N.(THETA)D.GRAD(C)
C
  300 IF(NNES.EQ.0) GO TO 500
C
      DO 390 MP=1,NNES
C
        ITYP=INTYP(MP)
        QNBMP=QNB(ITYP)
        MPB=ISN(5,MP)
        LS=ISB(5,MPB)
        M=ISB(6,MPB)
        CALL ELENOD(IE(M,5),IE(M,7),NODE,IQ,ID)
C
        NODE=4
C
        DO 330 IQ=1,4
          I=KGB(IQ,LS,ID)
          IF(I.EQ.0)THEN
            NODE=3
            GOTO 330
          ENDIF
          NI=IE(M,I)
          XQ(IQ)=X(NI,1)
          YQ(IQ)=X(NI,2)
          ZQ(IQ)=X(NI,3)
          IF(IMID.EQ.0)THEN
            VXQ(IQ)=W1*V(NI,1)+W2*VP(NI,1)
            VYQ(IQ)=W1*V(NI,2)+W2*VP(NI,2)
            VZQ(IQ)=W1*V(NI,3)+W2*VP(NI,3)
          ELSE
            VXQ(IQ)=0.5D0*(V(NI,1)+VP(NI,1))
            VYQ(IQ)=0.5D0*(V(NI,2)+VP(NI,2))
            VZQ(IQ)=0.5D0*(V(NI,3)+VP(NI,3))
          ENDIF
  330   CONTINUE
C
        CALL Q34CNV
     >    (BQ,RQ,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB(1,MPB),QNBMP,2,NODE,IQUAR)
C
        DO 380 IQ=1,NODE
          I=KGB(IQ,LS,ID)
          NI=IE(M,I)
          RLD(NI)=RLD(NI)+RQ(IQ)
  380   CONTINUE
C
  390 CONTINUE
C
C ******* APPLY VARIABLE BOUNDARY CONDITIONS
C
  500 IF(NVES.EQ.0) GO TO 700
C
      DO 590 MP=1,NVES
C
        ITYP=IVTYP(MP)
        CINMP=CVB(ITYP)
        MPB=ISV(5,MP)
        LS=ISB(5,MPB)
        M=ISB(6,MPB)
        CALL ELENOD(IE(M,5),IE(M,7),NODE,IQ,ID)
C
        NODE=4
C
        DO 530 IQ=1,4
          I=KGB(IQ,LS,ID)
          IF(I.EQ.0)THEN
            NODE=3
            GOTO 530
          ENDIF
          NI=IE(M,I)
          XQ(IQ)=X(NI,1)
          YQ(IQ)=X(NI,2)
          ZQ(IQ)=X(NI,3)
          IF(IMID.EQ.0)THEN
            VXQ(IQ)=W1*V(NI,1)+W2*VP(NI,1)
            VYQ(IQ)=W1*V(NI,2)+W2*VP(NI,2)
            VZQ(IQ)=W1*V(NI,3)+W2*VP(NI,3)
          ELSE
            VXQ(IQ)=0.5D0*(V(NI,1)+VP(NI,1))
            VYQ(IQ)=0.5D0*(V(NI,2)+VP(NI,2))
            VZQ(IQ)=0.5D0*(V(NI,3)+VP(NI,3))
          ENDIF
  530   CONTINUE
C
        CALL Q34CNV
     >    (BQ,RQ,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB(1,MPB),CINMP,3,NODE,IQUAR)
C
        DO 580 IQ=1,NODE
          I=KGB(IQ,LS,ID)
          NI=IE(M,I)
          RLD(NI)=RLD(NI) + RQ(IQ)
          DO 560 JQ=1,NODE
            J=KGB(JQ,LS,ID)
            NJ=IE(M,J)
            RLD(NI)=RLD(NI)-W2*CSTAR(NJ)
            DO 540 JJ=1,MXJBD
              LNODE=LRN(JJ,NI)
              IF(LNODE.EQ.NJ) GO TO 550
  540       CONTINUE
C
            WRITE (6,5000) MP,IQ,NI,JQ,NJ
            STOP
C
  550       CMATRX(NI,JJ)=CMATRX(NI,JJ)+W1*BQ(IQ,JQ)
  560     CONTINUE
  580   CONTINUE
C
  590 CONTINUE
C
C ******* APPLY DIRICHLET BOUNDARY CONDITION
C
  700 IF(NDNP.EQ.0) GO TO 900
C
      DO 740 NPP=1,NDNP
        NI=NPDB(NPP)
        ITYP=IDTYP(NPP)
C ------- put the Dirichlet concentration on the right-hand side
        BB=CDB(ITYP)
        RLD(NI)=BB
C ------- modify the row corresponding to Dirichlet node.
        DO 710 I=1,MXJBD
          CMATRX(NI,I)=0.0
          IB=LRN(I,NI)
          IF(IB.EQ.NI) CMATRX(NI,I)=1.0D0
  710   CONTINUE
C ------- modify the column corresponding to the Dirichlet node.  The
C ------- reason of this is to make the coefficient matrix symmetric.
        DO 720 INP=1,NLRN(NI)
          NP=LRN(INP,NI)
          IF(NP.EQ.NI .OR. NP.EQ.0) GO TO 720
          DO 715 IP=1,MXJBD
            IB=LRN(IP,NP)
            IF(IB.EQ.0) GO TO 715
            IF(IB.EQ.NI) THEN
              RLD(NP)=RLD(NP)-CMATRX(NP,IP)*RLD(NI)
              CMATRX(NP,IP)=0.0D0
              GO TO 720
            ENDIF
  715     CONTINUE
  720   CONTINUE
C
  740 CONTINUE
C
  900 CONTINUE
C
      RETURN
C
 1000 FORMAT('0','FOR',I4,'-TH CAUCHY SIDE',I2,'-TH NODE EQUATION NO.',
     1 I4/1X,'  WE CANNOT FIND THE COEFFICIENT FOR THE',I2,'-TH NODE',
     2 ' UNKNOWN NO.',I4,'   STOP')
 5000 FORMAT('0','FOR',I4,'-TH VB SIDE',I2,'-TH NODE EQUATION NO.',
     1 I4/1X,'  WE CANNOT FIND THE COEFFICIENT FOR THE',I2,'-TH NODE',
     2 ' UNKNOWN NO.',I4,'   STOP')
      END
C
c
c
      SUBROUTINE TBC1
     I    (NDB, NLRND,LRN,CPFG,IE,NDBD,MPLOCW,MPLOCD,IB,
     I     NLS,IFLUX, IMSV, ILSV, MXMSV,MXLSV,
     I     X,XWFG,EPSX, XPFG,NDFG,
     O     RLD,CMATRX)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /SAZFM/ MXNPFG,MXKGL,MXKGLD,MXNEP,MXEPW,MXNPW,MXELW,MXNDB,
     >               mxnpww,mxelww
C
      DIMENSION CMATRX(MXADNP,MXJBD),RLD(MXADNP),LRN(MXJBD,MXADNP)
      DIMENSION CPFG(MXNPFG),NLRND(MXADNP),X(MAXNP,3)
      DIMENSION XWFG(MXNPFG,3),IE(MAXEL,11),IB(MAXNP),MPLOCW(MXNPFG)
c     DIMENSION LRL(MXKBD,MAXNP),NLRL(MAXNP)
      DIMENSION NDBD(MXNDB),MPLOCD(MXNDB),XPFG(MXNPFG)
      DIMENSION IMSV(MXMSV,4),ILSV(MXLSV,5),KGB(4,6,3)
      DIMENSION XX(8),YY(8),ZZ(8),CC(8),DL(4)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
C ----- INITIATE XPFG(NP) AS AN INDEX ARRAY
C
      DO NP=1,NDFG
        XPFG(NP)=0.0D0
      ENDDO
C
C ----- APPLY THE INTRA-BOUNDARY CONDITIONS
C
      IF(IFLUX.NE.0)THEN
        IF(NLS.EQ.0)GOTO 100
        DO 90 MP=1,NLS
          NODE=4
          MM=ILSV(MP,5)
          DO 10 IQ=1,4
            NI=IMSV(MM,IQ)
            IF(IQ.EQ.1)N1=NI
            IF(IQ.EQ.2)N2=NI
            IF(IQ.EQ.3)N3=NI
            IF(NI.EQ.0)THEN
              NODE=3
              GOTO 10
            ENDIF
            IF(IQ.EQ.4)N4=NI
            XX(IQ)=X(NI,1)
            YY(IQ)=X(NI,2)
            ZZ(IQ)=X(NI,3)
  10      CONTINUE
          DO 30 IQ=1,NODE
            NI=ILSV(MP,IQ)
            IF(XPFG(NI).NE.0.0D0)GOTO 30
            IF(NI.LE.NNP)GOTO 30
            XPFG(NI)=1.0D0
            RLD(NI)=0.0D0
            DO 26 K=1,NODE
              NK=IMSV(MM,K)
              DO 25 JJ=1,NLRND(NI)
                LNODE=LRN(JJ,NI)
                IF(LNODE.EQ.NK)GOTO 26
  25          CONTINUE
              NLRND(NI)=NLRND(NI)+1
              NLD=NLRND(NI)
              LRN(NLD,NI)=NK
  26        CONTINUE
            XQ=XWFG(NI,1)
            YQ=XWFG(NI,2)
            ZQ=XWFG(NI,3)
            CALL LOCPLN
     I        (XQ,YQ,ZQ,XX,YY,ZZ,1,2,3,4,NODE,
     O         XSI,ETA,DL1,DL2,DL3,DL)
            CALL SLAVPT
     I          (XSI,ETA,DL1,DL2,DL3,NODE,NI,NLRND,LRN,N1,N2,N3,N4,
     O           CMATRX)
   30     CONTINUE
C
   90   CONTINUE
      ENDIF
C
  100 CONTINUE
C
      DO 600 N=1,NDB
        NP=NDBD(N)
C
        IF(NP.LE.NNP)GOTO 600
C
        MP=MPLOCW(N)
        LP=MPLOCD(N)
        CALL ELENOD(IE(MP,5),IE(MP,7), NQ,I,IL)
C
        NODE=4
        DO I=1,4
          IQ=KGB(I,LP,IL)
          IF(IQ.EQ.0)THEN
            NODE=3
            N4=0
            GOTO 550
          ENDIF
          NPP=IE(MP,IQ)
          IF(I.EQ.1)N1=NPP
          IF(I.EQ.2)N2=NPP
          IF(I.EQ.3)N3=NPP
          IF(I.EQ.4)N4=NPP
          XX(I)=X(NPP,1)
          YY(I)=X(NPP,2)
          ZZ(I)=X(NPP,3)
          CC(I)=CPFG(NPP)
        ENDDO
C
C ***** CHECK THE ELEMENT SIDE
C
 550    CONTINUE
        XQ=XWFG(NP,1)
        YQ=XWFG(NP,2)
        ZQ=XWFG(NP,3)
C
C ----- CHECK IF THE PLANE IS AN INTERIOR PLANE OF THE DOMAIN
C
c       CALL CKCNEL
c    I    (LRL,NLRL,MP,N1,N2,N3,MAXNP,MXKBD,
c    O     ME1,KOUNT)
C
C ----- IF KOUNT=2 AND IFLUX=1, THEN THE DIRICHLET INTRA-BOUNDARY
C       CONDITIONS ARE NOT USED
C
c       IF(KOUNT.EQ.2 .AND. IFLUX.NE.0)GOTO 600
C
C ----- IF KOUNT=1, THEN CHECK IF THIS PLANE IS SPECIFIED FROM EULERIAN
C       POINT OF VIEW
C
c       IF(KOUNT.EQ.1)THEN
C
C ----- CHECK IF THIS POINT IS ON ANY SIDE LINE OF THE PLANE
C
          CALL CKSIDE
     I     (X(1,1),X(1,2),X(1,3),XQ,YQ,ZQ,NQ,N1,N2,N3,N4,MAXNP,EPSX,
     O      KON,J1,J2)
          IF(KON.NE.0)THEN
            IF(IB(J1).EQ.5 .OR. IB(J2).EQ.5)GOTO 600
            RLD(NP)=0.0D0
            IF(IB(J1).EQ.1 .AND. IB(J2).EQ.1)THEN
              XS1=DSQRT((XQ-X(J1,1))**2+(YQ-X(J1,2))**2+(ZQ-X(J1,3))**2)
              XS2=DSQRT((X(J2,1)-X(J1,1))**2+(X(J2,2)-X(J1,2))**2+
     >                  (X(J2,3)-X(J1,3))**2)
              XSI=XS1/XS2
              CCFG=CPFG(J1)*(1.0D0-XSI)+CPFG(J2)*XSI
              GOTO 560
            ELSE
              DO 126 K=1,2
                IF(K.EQ.1)NK=J1
                IF(K.EQ.2)NK=J2
                DO 125 JJ=1,NLRND(NP)
                  LNODE=LRN(JJ,NP)
                  IF(LNODE.EQ.NK)GOTO 126
  125           CONTINUE
                NLD=NLRND(NP)+1
                NLRND(NP)=NLD
                LRN(NLD,NP)=NK
  126         CONTINUE
              IF(X(J2,1).NE.X(J1,1))THEN
                XSI=(2.0D0*XQ-X(J1,1)-X(J2,1))/(X(J2,1)-X(J1,1))
              ELSEIF(X(J2,2).NE.X(J1,2))THEN
                XSI=(2.0D0*YQ-X(J1,2)-X(J2,2))/(X(J2,2)-X(J1,2))
              ELSE
                XSI=(2.0D0*ZQ-X(J1,3)-X(J2,3))/(X(J2,3)-X(J1,3))
              ENDIF
              DL1=0.5D0*(1.0D0-XSI)
              DL2=0.5D0*(1.0D0+XSI)
              DO JJ=1,NLRND(NP)
                LNODE=LRN(JJ,NP)
                IF(LNODE.EQ.J1)THEN
                  CMATRX(NP,JJ)=-DL1
                ELSEIF(LNODE.EQ.J2)THEN
                  CMATRX(NP,JJ)=-DL2
                ELSEIF(LNODE.EQ.NP)THEN
                  CMATRX(NP,JJ)=1.0D0
                ELSE
                  CMATRX(NP,JJ)=0.0D0
                ENDIF
              ENDDO
              GOTO 600
            ENDIF
          ENDIF
C
          IF(IB(N1).EQ.5)GOTO 600
          IF(IB(N2).EQ.5)GOTO 600
          IF(IB(N3).EQ.5)GOTO 600
          IF(NODE.EQ.4)THEN
            IF(IB(N4).EQ.5)GOTO 600
          ENDIF
c       ENDIF
C
C ---- STARTING DOING INTERPOLATION FOR THIS DIFFUSION BOUNDARY POINT
C
        CALL LOCPLN
     I    (XQ,YQ,ZQ,XX,YY,ZZ,1,2,3,4,NODE,
     O     XSI,ETA,DL1,DL2,DL3,DL)
C
        RLD(NP)=0.0D0
        IF(NODE.EQ.3)THEN
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1 .AND. IB(N3).EQ.1)THEN
            CCFG=CC(1)*DL1+CC(2)*DL2+CC(3)*DL3
          ELSE
            CALL SLAVPT
     I          (XSI,ETA,DL1,DL2,DL3,NODE,NP,NLRND,LRN,N1,N2,N3,N4,
     O           CMATRX)
            GOTO 600
          ENDIF
        ELSE
          IF(IB(N1).EQ.1 .AND. IB(N2).EQ.1 .AND. IB(N3).EQ.1 .AND.
     >       IB(N4).EQ.1)THEN
            CCFG=0.25D0*(CC(1)*(1.0D0-XSI)*(1.0D0-ETA)+
     >                   CC(2)*(1.0D0+XSI)*(1.0D0-ETA)+
     >                   CC(3)*(1.0D0+XSI)*(1.0D0+ETA)+
     >                   CC(4)*(1.0D0-XSI)*(1.0D0+ETA))
          ELSE
            CALL SLAVPT
     I          (XSI,ETA,DL1,DL2,DL3,NODE,NP,NLRND,LRN,N1,N2,N3,N4,
     O           CMATRX)
            GOTO 600
          ENDIF
        ENDIF
C
  560   CONTINUE
C
C ------ FOR DIRICHLET BOUNDARY
C
        RLD(NP)=CCFG
        DO J=1,NLRND(NP)
          NQ=LRN(J,NP)
          CMATRX(NP,J)=0.0D0
          IF(NQ.EQ.NP)CMATRX(NP,J)=1.0D0
        ENDDO
C
  600 CONTINUE
C
      RETURN
 1000 FORMAT('0','FOR',I4,'-TH CAUCHY SIDE',I2,'-TH NODE EQUATION NO.',
     1 I4/1X,'  WE CANNOT FIND THE COEFFICIENT FOR THE',I2,'-TH NODE',
     2 ' UNKNOWN NO.',I4,'   STOP')
      END
c
c
c
      SUBROUTINE SLAVPT
     I    (XSI,ETA,DL1,DL2,DL3,NODE,NI,NLRND,LRN,N1,N2,N3,N4,
     O     CMATRX)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
C
      DIMENSION CMATRX(MXADNP,MXJBD),LRN(MXJBD,MXADNP),NLRND(MXADNP)
C
      n1i=0
      n2i=0
      n3i=0
      n4i=0
            IF(NODE.EQ.3)THEN
              DO JJ=1,NLRND(NI)
                LNODE=LRN(JJ,NI)
                IF(LNODE.EQ.N1)THEN
                  CMATRX(NI,JJ)=-DL1
                  n1i=1
                ELSEIF(LNODE.EQ.N2)THEN
                  CMATRX(NI,JJ)=-DL2
                  n2i=1
                ELSEIF(LNODE.EQ.N3)THEN
                  CMATRX(NI,JJ)=-DL3
                  n3i=1
                ELSEIF(LNODE.EQ.NI)THEN
                  CMATRX(NI,JJ)=1.0D0
                ELSE
                  CMATRX(NI,JJ)=0.0D0
                ENDIF
              ENDDO
              if(n1i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n1
                cmatrx(ni,ndd)=-dl1
              endif
              if(n2i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n2
                cmatrx(ni,ndd)=-dl2
              endif
              if(n3i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n3
                cmatrx(ni,ndd)=-dl3
              endif
            ELSE
              DO JJ=1,NLRND(NI)
                LNODE=LRN(JJ,NI)
                IF(LNODE.EQ.N1)THEN
                  CMATRX(NI,JJ)=-0.25D0*(1.0D0-XSI)*(1.0D0-ETA)
                  n1i=1
                ELSEIF(LNODE.EQ.N2)THEN
                  CMATRX(NI,JJ)=-0.25D0*(1.0D0+XSI)*(1.0D0-ETA)
                  n2i=1
                ELSEIF(LNODE.EQ.N3)THEN
                  CMATRX(NI,JJ)=-0.25D0*(1.0D0+XSI)*(1.0D0+ETA)
                  n3i=1
                ELSEIF(LNODE.EQ.N4)THEN
                  CMATRX(NI,JJ)=-0.25D0*(1.0D0-XSI)*(1.0D0+ETA)
                  n4i=1
                ELSEIF(LNODE.EQ.NI)THEN
                  CMATRX(NI,JJ)=1.0D0
                ELSE
                  CMATRX(NI,JJ)=0.0D0
                ENDIF
              ENDDO
              if(n1i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n1
                CMATRX(NI,ndd)=-0.25D0*(1.0D0-XSI)*(1.0D0-ETA)
              endif
              if(n2i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n2
                CMATRX(NI,ndd)=-0.25D0*(1.0D0+XSI)*(1.0D0-ETA)
              endif
              if(n3i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n3
                CMATRX(NI,ndd)=-0.25D0*(1.0D0+XSI)*(1.0D0+ETA)
              endif
              if(n4i.eq.0)then
                nlrnd(ni)=nlrnd(ni)+1
                ndd=nlrnd(ni)
                lrn(ndd,ni)=n4
                CMATRX(NI,ndd)=-0.25D0*(1.0D0+XSI)*(1.0D0+ETA)
              endif
            ENDIF
      RETURN
      END
C
C
c
      SUBROUTINE Q34CNV(BQ,RQ,XQ,YQ,ZQ,VXQ,VYQ,VZQ,DCOSB,QBMP,
     >                  IBC,NODE,IQUAR)
C
C $$$$$ TO COMPUTE BOUNDARY-SURFACE MATRIX AND BOUNDARY-SURFACE LOAD
C       VECTOR OVER A BOUNDARY SURFACE.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(4)
C
      DIMENSION BQ(4,4),RQ(4),XQ(4),YQ(4),ZQ(4)
      DIMENSION VXQ(4),VYQ(4),VZQ(4),DCOSB(3)
      DIMENSION S(4),T(4),DNSS(4),DNTT(4)
      DIMENSION DL1(3),DL2(3),DL3(3),WG(3)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0/, T/-1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA WG/0.333333333333333D0, 0.333333333333333D0,
     >        0.333333333333333D0/
C
      IF(IQUAR.EQ.3 .OR. IQUAR.EQ.4)THEN
        P=0.577350269189626D0
        DL1(1)=0.666666666666667D0
        DL1(2)=0.166666666666667D0
        DL1(3)=0.166666666666667D0
        DL2(1)=0.166666666666667D0
        DL2(2)=0.666666666666667D0
        DL2(3)=0.166666666666667D0
        DL3(1)=0.166666666666667D0
        DL3(2)=0.166666666666667D0
        DL3(3)=0.666666666666667D0
      ELSE
        P=1.0D0
        DL1(1)=1.0D0
        DL1(2)=0.0D0
        DL1(3)=0.0D0
        DL2(1)=0.0D0
        DL2(2)=1.0D0
        DL2(3)=0.0D0
        DL3(1)=0.0D0
        DL3(2)=0.0D0
        DL3(3)=1.0D0
      ENDIF
C
C ------- INITIATE MATRICES BQ(IQ,JQ) AND RQ(IQ)
C
      DO 100 IQ=1,4
        RQ(IQ)=0.0
        DO 100 JQ=1,4
          BQ(IQ,JQ)=0.0
  100 CONTINUE
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS IF NODE.EQ.3
C
      IF(NODE.EQ.3)THEN
        DXDDL2=XQ(2)-XQ(1)
        DYDDL2=YQ(2)-YQ(1)
        DZDDL2=ZQ(2)-ZQ(1)
        DXDDL3=XQ(3)-XQ(1)
        DYDDL3=YQ(3)-YQ(1)
        DZDDL3=ZQ(3)-ZQ(1)
        DETX=DYDDL2*DZDDL3-DYDDL3*DZDDL2
        DETY=DXDDL2*DZDDL3-DXDDL3*DZDDL2
        DETZ=DXDDL2*DYDDL3-DXDDL3*DYDDL2
        DET1=DSQRT(DETX*DETX+DETY*DETY+DETZ*DETZ)*0.5D0
      ENDIF
C
C ------- SUMMATION OF THE INTEGRAND OVER THE GAUSSIAN POINTS
C
      DO 690 KG=1,NODE
C
C ------- DETERMINE LOACAL COORDINATE OF GAUSSIAN POINT KG
C
        IF(NODE.EQ.4)THEN
C
          SS=P*S(KG)
          TT=P*T(KG)
          SM=1.0D0-SS
          SP=1.0D0+SS
          TM=1.0D0-TT
          TP=1.0D0+TT
          N(1)=0.25D0*SM*TM
          N(2)=0.25D0*SP*TM
          N(3)=0.25D0*SP*TP
          N(4)=0.25D0*SM*TP
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS
C
          DNSS(1)=-0.25D0*TM
          DNSS(2)= 0.25D0*TM
          DNSS(3)= 0.25D0*TP
          DNSS(4)=-0.25D0*TP
          DNTT(1)=-0.25D0*SM
          DNTT(2)=-0.25D0*SP
          DNTT(3)= 0.25D0*SP
          DNTT(4)= 0.25D0*SM
          DXDSS=0.0D0
          DYDSS=0.0D0
          DZDSS=0.0D0
          DXDTT=0.0D0
          DYDTT=0.0D0
          DZDTT=0.0D0
          DO 290 IQ=1,4
            DXDSS=DXDSS+XQ(IQ)*DNSS(IQ)
            DYDSS=DYDSS+YQ(IQ)*DNSS(IQ)
            DZDSS=DZDSS+ZQ(IQ)*DNSS(IQ)
            DXDTT=DXDTT+XQ(IQ)*DNTT(IQ)
            DYDTT=DYDTT+YQ(IQ)*DNTT(IQ)
            DZDTT=DZDTT+ZQ(IQ)*DNTT(IQ)
  290     CONTINUE
          DETZ=DXDSS*DYDTT-DYDSS*DXDTT
          DETY=-DXDSS*DZDTT+DZDSS*DXDTT
          DETX=DYDSS*DZDTT-DZDSS*DYDTT
          DET=DSQRT(DETX*DETX+DETY*DETY+DETZ*DETZ)
C
        ELSE
C
          N(1)=DL1(KG)
          N(2)=DL2(KG)
          N(3)=DL3(KG)
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS
C
          DET=DET1*WG(KG)
C
        ENDIF
C
C ------- ACCUMULATE THE SUMS TO OBTAIN THE MATRIX INTEGRALS BQ AND RQ
C
        GOTO(310,410,510)IBC
C
C ******* CAUCHY CONDITIONS
C
  310   VXK=0.0
        VYK=0.0
        VZK=0.0
        DO 320 IQ=1,NODE
          VXK=VXK+VXQ(IQ)*N(IQ)
          VYK=VYK+VYQ(IQ)*N(IQ)
          VZK=VZK+VZQ(IQ)*N(IQ)
  320   CONTINUE
        VNK=VXK*DCOSB(1)+VYK*DCOSB(2)+VZK*DCOSB(3)
        DO 390 IQ=1,NODE
          RQ(IQ)=RQ(IQ)-N(IQ)*QBMP*DET
          DO 350 JQ=1,NODE
            BQ(IQ,JQ)=BQ(IQ,JQ)-N(IQ)*VNK*N(JQ)*DET
  350     CONTINUE
  390   CONTINUE
        GOTO 690
C
C ******* NEUMANN CONDITIONS
C
  410   DO 490 IQ=1,NODE
          RQ(IQ)=RQ(IQ)-N(IQ)*QBMP*DET
  490   CONTINUE
        GOTO 690
C
C ******* VARIABLE CONDITIONS
C
  510   VXK=0.0
        VYK=0.0
        VZK=0.0
        DO 520 IQ=1,NODE
          VXK=VXK+VXQ(IQ)*N(IQ)
          VYK=VYK+VYQ(IQ)*N(IQ)
          VZK=VZK+VZQ(IQ)*N(IQ)
  520   CONTINUE
        VNK=VXK*DCOSB(1)+VYK*DCOSB(2)+VZK*DCOSB(3)
        IF(VNK.GE.0.0)GOTO 690
        DO 590 IQ=1,NODE
          RQ(IQ)=RQ(IQ)-N(IQ)*VNK*QBMP*DET
          DO 550 JQ=1,NODE
            BQ(IQ,JQ)=BQ(IQ,JQ)-N(IQ)*VNK*N(JQ)*DET
  550     CONTINUE
  590   CONTINUE
C
  690 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE TSFLOW(BFLX,X,IE,C,F,TH,RKD,TRANC,
     1 DCOSB,ISB,NPBB, SOS,ISTYP,LES,WSS,IWTYP,NPW, NPVB,NPDB,NPCB,NPNB,
     2 PROP, DELT, KFLOW, KKK)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE MATERIAL FLUXES, INCREMENTAL MASS FLOW, AND
C ------- ACCUMULATED MASS FLOW THROUGH ALL TYPES OF BOUNDARIES: AND
C ------- CHANGE OF MATERIALS IN THE REGION OF INTEREST.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: X(NNP,3), IE(NEL,9),
C -------        C(NNP), F(NNP,3), TH(NEL,8),
C -------        DCOSB(3,NBES), ISB(6,NBES), NPBB(NBNP), SOS(NSPR,2),
C -------        ISTYP(NSEL), LES(NSEL), WSS(NWPR,2), IWTYP(NWNP),
C -------        NPW(NWNP), NPVB(NVNP), NPDB(NDNP), NPCB(NCNP),
C -------        NPNB(NNNP), PROP(NMNNP,NMAT), DELT, AND KFLOW.
C
C ------- OUTPUT: FRATE(14), FLOW(14), AND TFLOW(14).
C
C ------- WORKING ARRAYS: BFLX(NBNP,2)
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KD,LAMBDA
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
C
      COMMON /CELS/ MXSEL,MXSPR,MXSDP,NSEL,NSPR,NSDP,KSAI
      COMMON /CNPS/ MXWNP,MXWPR,MXWDP,NWNP,NWPR,NWDP,KWAI
C
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /TDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
      COMMON /TCBC/ MXCNP,MXCES,MXCPR,MXCDP,NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /TNBC/ MXNNP,MXNES,MXNPR,MXNDP,NNNP,NNES,NNPR,NNDP,KNAI
C
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
C
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
C
      COMMON /TFLOW/ FRATE(14),FLOW(14),TFLOW(14,7)
C
      DIMENSION BFLX(MAXBNP,2)
      DIMENSION X(MAXNP,3),IE(MAXEL,11)
      DIMENSION C(MAXNP),F(MAXNP,3),TH(MAXEL,8)
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES),NPBB(MAXBNP)
C
      DIMENSION SOS(MXSPR,2),ISTYP(MXSEL),LES(MXSEL)
      DIMENSION WSS(MXWPR,2),IWTYP(MXWNP),NPW(MXWNP)
C
      DIMENSION NPVB(MXVNP),NPDB(MXDNP),NPCB(MXCNP),NPNB(MXNNP)
C
      DIMENSION PROP(MAXMAT,MXMPPM),RKD(MAXMAT)
      DIMENSION TRANC(MAXMAT)
C
      DIMENSION CQ(8),CSQ(8),XQ(8),YQ(8),ZQ(8),THG(8)
      DIMENSION RRQ(4),FFQ(4),XXQ(4),YYQ(4),ZZQ(4)
      DIMENSION QRY(7),QDY(7),QLY(7),SRCY(7)
C
      DIMENSION KGB(4,6,3)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
      DATA QRY /0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
      DATA QDY /0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
      DATA QLY /0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
      DATA SRCY /0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
C
C ******* CALCULATE NODAL FLOW RATES THROUGH ALL BOUNDARY NODES
C
      DO 110 NP=1,NBNP
        BFLX(NP,1)=BFLX(NP,2)
        BFLX(NP,2)=0.0
  110 CONTINUE
C
      DO 170 MP=1,NBES
C
        LS=ISB(5,MP)
        M=ISB(6,MP)
        IF(IE(M,5).EQ.0)THEN
          ID=3
        ELSEIF(IE(M,7).EQ.0)THEN
          ID=2
        ELSE
          ID=1
        ENDIF
        NODE=4
C
        DO 120 IQ=1,4
          I=KGB(IQ,LS,ID)
          IF(I.EQ.0)THEN
            NODE=3
            GOTO 120
          ENDIF
          NI=IE(M,I)
          XXQ(IQ)=X(NI,1)
          YYQ(IQ)=X(NI,2)
          ZZQ(IQ)=X(NI,3)
          FFQ(IQ)=DCOSB(1,MP)*F(NI,1)+DCOSB(2,MP)*F(NI,2)+
     >            DCOSB(3,MP)*F(NI,3)
  120   CONTINUE
C
        CALL Q34BB(RRQ,FFQ,XXQ,YYQ,ZZQ,NODE,IQUAR)
C
        DO 140 IQ=1,NODE
          NII=ISB(IQ,MP)
          IF(NII.EQ.0)GOTO 140
          BFLX(NII,2)=BFLX(NII,2)+RRQ(IQ)
  140   CONTINUE
C
  170 CONTINUE
C
      IF(KFLOW.GT.0) GO TO 200
      DO 180 NP=1,NBNP
  180 BFLX(NP,1)=BFLX(NP,2)
      DO 190 I=1,14
  190 TFLOW(I,KKK)=0.0
C
C ******* DETERMINE FLOWS AND FLOW RATES THROUGH VARIOUS TYPES OF NODES,
C ******* STARTING WITH THE NET FLOWS THROUGH ALL BOUNDARY NODES
C
  200 S=0.
      SP=0.
      DO 210 NP=1,NBNP
        S=S+BFLX(NP,2)
        SP=SP+BFLX(NP,1)
  210 CONTINUE
C
      FRATE(7)=S
      FLOW(7)=.5D0*(S+SP)*DELT
C
C ******* THROUGH DIRICHLET BOUNDARY NODES
C
      FRATE(1)=0.
      FLOW(1)=0.
      IF(NDNP.LE.0) GO TO 340
      S=0.
      SP=0.
      DO 330 NPP=1,NDNP
        NP=NPDB(NPP)
        DO 310 I=1,NBNP
          IJ=NPBB(I)
          IF(IJ.NE.NP) GO TO 310
          NII=I
          GO TO 320
  310   CONTINUE
  320   CONTINUE
        S=S+BFLX(NII,2)
        SP=SP+BFLX(NII,1)
  330 CONTINUE
      FRATE(1)=S
      FLOW(1)=0.5D0*(S+SP)*DELT
C
C ******* THROUGH CAUCHY NODES
C
  340 FRATE(2)=0.
      FLOW(2)=0.
      IF(NCNP.LE.0) GO TO 380
      S=0.
      SP=0.
      DO 370 NPP=1,NCNP
        NII=NPCB(NPP)
        S=S+BFLX(NII,2)
        SP=SP+BFLX(NII,1)
  370 CONTINUE
      FRATE(2)=S
      FLOW(2)=0.5D0*(S+SP)*DELT
C
C ******* THROUGH NEUMANN NODES
C
  380 FRATE(3)=0.0
      FLOW(3)=0.0
      IF(NNNP.LE.0) GO TO 400
      S=0.0
      SP=0.0
      DO 395 NPP=1,NNNP
        NII=NPNB(NPP)
        S=S+BFLX(NII,2)
        SP=SP+BFLX(NII,1)
  395 CONTINUE
      FRATE(3)=S
      FLOW(3)=0.5D0*(S+SP)*DELT
C
C ******* THROUGH RUN-IN SEEP-OUT NODES
C
  400 FRATE(4)=0.
      FLOW(4)=0.
      FRATE(5)=0.0
      FLOW(5)=0.0
      IF(NVNP.LE.0) GO TO 500
      S=0.
      SP=0.
      SM=0.0
      SMP=0.0
      DO 490 NPP=1,NVNP
        NII=NPVB(NPP)
        IF(BFLX(NII,2).LT.0.0) GO TO 460
        S=S+BFLX(NII,2)
        GO TO 470
  460   SM=SM+BFLX(NII,2)
  470   IF(BFLX(NII,1).LT.0.0) GO TO 480
        SP=SP+BFLX(NII,1)
        GO TO 490
  480   SMP=SMP+BFLX(NII,1)
  490 CONTINUE
      FRATE(4)=S
      FLOW(4)=0.5D0*(S+SP)*DELT
      FRATE(5)=SM
      FLOW(5)=0.5D0*(SM+SMP)*DELT
C
C ******* NUMERICAL FLOW THROUGH UNSPECIFIED BOUNDARY NODES
C
  500 S=0.
      SP=0.
      DO 510 I=1,5
        S=S+FRATE(I)
  510 SP=SP+FLOW(I)
      FRATE(6)=FRATE(7)-S
      FLOW(6)=FLOW(7)-SP
C
C ******* CALCULATE INCREASES OF INTEGRATED MATERIAL CONTENTS IN THE
C ******* FLUID QR AND IN SOLID QD; DETERMINE LOSSES DUE TO RADIOACTIVE
C ******* DECAY QL; AND COMPUTE INTEGRATE THE SOURCE/SINK SOURCES
C
      QRP=QRY(KKK)
      QDP=QDY(KKK)
      QLP=QLY(KKK)
      SOURSP=SRCY(KKK)
      QR=0.
      QD=0.
      QL=0.
      SOURCE=0.0
C
      DO 690 M=1,NEL
C
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,MP,MP)
C
        SOSQP=0.0
        SOSCP=0.0
        IF(NSEL.EQ.0) GO TO 605
        DO 603 MP=1,NSEL
          MS=LES(MP)
          IF(MS.NE.M) GO TO 603
          ITYP=ISTYP(MP)
          SOSQP=SOS(ITYP,1)
          SOSCP=SOS(ITYP,2)
          GO TO 605
  603   CONTINUE
C
  605   MTYP=IE(M,9)
        KD=RKD(MTYP)
        RHOB=PROP(MTYP,1)
        LAMBDA=TRANC(MTYP)
C
        DO 650 IQ=1,NODE
          NP=IE(M,IQ)
          XQ(IQ)=X(NP,1)
          YQ(IQ)=X(NP,2)
          ZQ(IQ)=X(NP,3)
          CQ(IQ)=C(NP)
C
          CSQ(IQ)=KD*CQ(IQ)
  650   CONTINUE
C
        DO 670 KG=1,NODE
          THG(KG)=TH(M,KG)
  670   CONTINUE
C
        CALL Q468R(QRM,QDM,SOSM,CQ,CSQ,THG,XQ,YQ,ZQ,SOSQP,SOSCP,
     >             NODE,IQUAR)
C
        QR=QR+QRM
        QD=QD+QDM*RHOB
        QLM=QRM+QDM*RHOB
        QL=QL+LAMBDA*QLM
        SOURCE=SOURCE-SOSM
  690 CONTINUE
C
C ******** INCORPORATE WELL SOURCE/SINK
C
      IF(NWNP.LE.0) GO TO 705
      DO 700 I=1,NWNP
        ITYP=IWTYP(I)
        WSSQ=WSS(ITYP,1)
        WSSC=WSS(ITYP,2)
        NI=NPW(I)
        IF(WSSQ.LT.0.0) SOURCE=SOURCE-WSSQ*C(NI)
  700 IF(WSSQ.GE.0.0) SOURCE=SOURCE-WSSQ*WSSC
C
  705 IF(KFLOW.GT.0) GO TO 710
      QRP=0.0
      QDP=0.0
      QLP=QL
      SOURSP=SOURCE
      SUM=QR+QD
      S=FRATE(7)+QL+SOURCE
C
  710 FLOW(8)=QR-QRP
      FRATE(8)=FLOW(8)/DELT
      IF(KFLOW.LE.0 .AND. SUM.NE.0.0D0) FRATE(8)=-S*QR/SUM
      FLOW(9)=QD-QDP
      FRATE(9)=FLOW(9)/DELT
      IF(KFLOW.LE.0 .AND. SUM.NE.0.0D0) FRATE(9)=-S*QD/SUM
C
      FRATE(10)=QL
      FRATE(11)=0.0
      FRATE(12)=0.0
      FRATE(13)=0.0
      FRATE(14)=SOURCE
      FLOW(10)=0.5D0*(QL+QLP)*DELT
      FLOW(11)=0.0
      FLOW(12)=0.0
      FLOW(13)=0.0
      FLOW(14)=0.5D0*(SOURCE+SOURSP)*DELT
C
      DO 720 I=1,14
        TFLOW(I,KKK)=TFLOW(I,KKK)+FLOW(I)
  720 CONTINUE
C
      QRY(KKK)=QR
      QDY(KKK)=QD
      QLY(KKK)=QL
      SRCY(KKK)=SOURCE
      RETURN
      END
C
C
C
      SUBROUTINE Q34BB(RQ,FQ,XQ,YQ,ZQ,NODE,IQUAR)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE THE NORMAL FLOW RATES (M/T) BY INTEGRATING THE
C ------- NORMAL FLUXES (M/L**2/T) OVER A BOUNDARY SURFACE.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: FQ(4), XQ(4), YQ(4), AND ZQ(4).
C
C ------- OUTPUT: RQ(4).
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(4),L1(3),L2(3),L3(3)
C
      DIMENSION RQ(4),FQ(4),XQ(4),YQ(4),ZQ(4)
      DIMENSION S(4),T(4),DNSS(4),DNTT(4),WG(3)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0/, T/-1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA WG/0.333333333333333D0, 0.333333333333333D0,
     1        0.333333333333333D0/
C
      IF(IQUAR.EQ.3 .OR. IQUAR.EQ.4)THEN
        P=0.577350269189626D0
        L1(1)=0.666666666666667D0
        L1(2)=0.166666666666667D0
        L1(3)=0.166666666666667D0
        L2(1)=0.166666666666667D0
        L2(2)=0.666666666666667D0
        L2(3)=0.166666666666667D0
        L3(1)=0.166666666666667D0
        L3(2)=0.166666666666667D0
        L3(3)=0.666666666666667D0
      ELSE
        P=1.0D0
        L1(1)=1.0D0
        L1(2)=0.0D0
        L1(3)=0.0D0
        L2(1)=0.0D0
        L2(2)=1.0D0
        L2(3)=0.0D0
        L3(1)=0.0D0
        L3(2)=0.0D0
        L3(3)=1.0D0
      ENDIF
C
C ------- INITIATE MATRICES RQ(IQ)
C
      DO 100 IQ=1,4
        RQ(IQ)=0.0
  100 CONTINUE
C
      IF(NODE.EQ.3)THEN
        DXDDL2=XQ(2)-XQ(1)
        DYDDL2=YQ(2)-YQ(1)
        DZDDL2=ZQ(2)-ZQ(1)
        DXDDL3=XQ(3)-XQ(1)
        DYDDL3=YQ(3)-YQ(1)
        DZDDL3=ZQ(3)-ZQ(1)
        DETX=DYDDL2*DZDDL3-DYDDL3*DZDDL2
        DETY=DXDDL2*DZDDL3-DXDDL3*DZDDL2
        DETZ=DXDDL2*DYDDL3-DXDDL3*DYDDL2
        DET1=DSQRT(DETX*DETX+DETY*DETY+DETZ*DETZ)*0.5D0
      ENDIF
C
C ------- SUMMATION OF THE INTEGRAND OVER THE GAUSSIAN POINTS
C
      DO 490 KG=1,NODE
C
C ------- DETERMINE LOACAL COORDINATE OF GAUSSIAN POINT KG
C
        IF(NODE.EQ.4)THEN
          SS=P*S(KG)
          TT=P*T(KG)
          SM=1.0D0-SS
          SP=1.0D0+SS
          TM=1.0D0-TT
          TP=1.0D0+TT
C
          N(1)=0.25D0*SM*TM
          N(2)=0.25D0*SP*TM
          N(3)=0.25D0*SP*TP
          N(4)=0.25D0*SM*TP
          DNSS(1)=-0.25D0*TM
          DNSS(2)= 0.25D0*TM
          DNSS(3)= 0.25D0*TP
          DNSS(4)=-0.25D0*TP
          DNTT(1)=-0.25D0*SM
          DNTT(2)=-0.25D0*SP
          DNTT(3)= 0.25D0*SP
          DNTT(4)= 0.25D0*SM
C
          DXDSS=0.0D0
          DYDSS=0.0D0
          DZDSS=0.0D0
          DXDTT=0.0D0
          DYDTT=0.0D0
          DZDTT=0.0D0
          DO 290 IQ=1,4
            DXDSS=DXDSS+XQ(IQ)*DNSS(IQ)
            DYDSS=DYDSS+YQ(IQ)*DNSS(IQ)
            DZDSS=DZDSS+ZQ(IQ)*DNSS(IQ)
            DXDTT=DXDTT+XQ(IQ)*DNTT(IQ)
            DYDTT=DYDTT+YQ(IQ)*DNTT(IQ)
            DZDTT=DZDTT+ZQ(IQ)*DNTT(IQ)
  290     CONTINUE
C
          DETZ=DXDSS*DYDTT-DYDSS*DXDTT
          DETY=-DXDSS*DZDTT+DZDSS*DXDTT
          DETX=DYDSS*DZDTT-DZDSS*DYDTT
          DET=DSQRT(DETX*DETX+DETY*DETY+DETZ*DETZ)
        ELSE
          N(1)=L1(KG)
          N(2)=L2(KG)
          N(3)=L3(KG)
C
C ----- COMPUTE JACOBIAN AT GAUSSIAN POINTS
C
          DET=DET1*WG(KG)
C
        ENDIF
C
C ------- ACCUMULATE THE SUMS TO OBTAIN THE MATRIX INTEGRALS RQ(IQ)
C
        FK=0.0D0
        DO 350 IQ=1,NODE
          FK=FK+FQ(IQ)*N(IQ)
  350   CONTINUE
C
        DO 390 IQ=1,NODE
          RQ(IQ)=RQ(IQ)+N(IQ)*FK*DET
  390   CONTINUE
C
  490 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE Q468R
     >     (QRM,QDM,SOSM,CQ,CSQ,THG,XQ,YQ,ZQ,SOSQP,SOSCP,NODE,IQUAR)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO COMPUTE THE MATERIAL INTEGRATION AND ELEMENT SOURCE
C ------- INTEGRATION OVER AN ELEMENT.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: CQ(8), CSQ(8), THG(8), XQ(8), YQ(8), ZQ(8), SOSQP, AND
C -------        SOSCP.
C
C ------- OUTPUT: QRM, QDM, AND SOSM.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 N(8),L1(3),L2(3),L3(3),LL1(4),LL2(4),LL3(4),LL4(4)
C
      DIMENSION CQ(8),CSQ(8),THG(8), XQ(8),YQ(8),ZQ(8)
      DIMENSION S(8),T(8),U(8),DNX(8),DNY(8),DNZ(8),W(8)
C
      DATA S/-1.0D0,1.0D0,1.0D0,-1.0D0, -1.0D0,1.0D0,1.0D0,-1.0D0/
      DATA T/-1.0D0,-1.0D0,1.0D0,1.0D0, -1.0D0,-1.0D0,1.0D0,1.0D0/
      DATA U/-1.0D0,-1.0D0,-1.0D0,-1.0D0, 1.0D0,1.0D0,1.0D0,1.0D0/
C
      IF(IQUAR.EQ.1 .OR. IQUAR.EQ.3)THEN
        P=0.577350269189626D0
        L1(1)=0.666666666666667D0
        L1(2)=0.166666666666667D0
        L1(3)=0.166666666666667D0
        L2(1)=0.166666666666667D0
        L2(2)=0.666666666666667D0
        L2(3)=0.166666666666667D0
        L3(1)=0.166666666666667D0
        L3(2)=0.166666666666667D0
        L3(3)=0.666666666666667D0
        LL1(1)=0.58541020D0
        LL1(2)=0.13819660D0
        LL1(3)=0.13819660D0
        LL1(4)=0.13819660D0
        LL2(1)=0.13819660D0
        LL2(2)=0.58541020D0
        LL2(3)=0.13819660D0
        LL2(4)=0.13819660D0
        LL3(1)=0.13819660D0
        LL3(2)=0.13819660D0
        LL3(3)=0.58541020D0
        LL3(4)=0.13819660D0
        LL4(1)=0.13819660D0
        LL4(2)=0.13819660D0
        LL4(3)=0.13819660D0
        LL4(4)=0.58541020D0
      ELSE
        P=1.0D0
        L1(1)=1.0D0
        L1(2)=0.0D0
        L1(3)=0.0D0
        L2(1)=0.0D0
        L2(2)=1.0D0
        L2(3)=0.0D0
        L3(1)=0.0D0
        L3(2)=0.0D0
        L3(3)=1.0D0
        LL1(1)=1.0D0
        LL1(2)=0.0D0
        LL1(3)=0.0D0
        LL1(4)=0.0D0
        LL2(1)=0.0D0
        LL2(2)=1.0D0
        LL2(3)=0.0D0
        LL2(4)=0.0D0
        LL3(1)=0.0D0
        LL3(2)=0.0D0
        LL3(3)=1.0D0
        LL3(4)=0.0D0
        LL4(1)=0.0D0
        LL4(2)=0.0D0
        LL4(3)=0.0D0
        LL4(4)=1.0D0
      ENDIF
C
      QRM=0.0
      QDM=0.0
      SOSM=0.0
C
      DO 490 KG=1,NODE
C
C ------- DETERMINE LOACAL COORDINATE OF GAUSSIAN POINT KG.
C
        IF(NODE.EQ.8)THEN
          SS=P*S(KG)
          TT=P*T(KG)
          UU=P*U(KG)
          CALL SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,0,
     O       N,DNX,DNY,DNZ,W,DJAC)
        ELSEIF(NODE.EQ.6)THEN
          IF(KG.LE.3)THEN
            XSI=-P
            KKG=KG
          ELSE
            XSI=P
            KKG=KG-3
          ENDIF
          DL1=L1(KKG)
          DL2=L2(KKG)
          DL3=L3(KKG)
          CALL SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,0,
     O       N,DNX,DNY,DNZ,W,DJAC)
          DJAC=DJAC/3.0D0
        ELSEIF(NODE.EQ.4)THEN
          D1=LL1(KG)
          D2=LL2(KG)
          D3=LL3(KG)
          D4=LL4(KG)
          CALL SHAPE
     I      (XQ,YQ,ZQ, SS,TT,UU, XSI,DL1,DL2,DL3, D1,D2,D3,D4, NODE,0,
     O       N,DNX,DNY,DNZ,W,DJAC)
          DJAC=0.25D0*DJAC
        ENDIF
C
C ------- INTERPOLATE TO OBTAIN WATER CONTENT AT THE GAUSSIAN POINT KG
C
        CQP=0.0
        CSQP=0.0
        SOSMK=0.0
        DO 390 IQ=1,NODE
          CQP=CQP+CQ(IQ)*N(IQ)
          CSQP=CSQP+CSQ(IQ)*N(IQ)
          IF(SOSQP.LT.0.0) SOSMK=SOSMK+SOSQP*CQ(IQ)*N(IQ)
  390   CONTINUE
        IF(SOSQP.GE.0.0) SOSMK=SOSQP*SOSCP
C
        THQP=THG(KG)
C
C ------- ACCUMULATE THE SUM TO EVALUATE THE INTEGRAL
C
        QRM=QRM+THQP*CQP*DJAC
        QDM=QDM+CSQP*DJAC
        SOSM=SOSM+SOSMK*DJAC
  490 CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE TPRINT(C,F, TIME,DELT,KPR,KOUT,KDIAG,ITIM,KKK)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO OUTPUT MATERIAL FLOWS, CONCENTRATION, AND MATERIAL FLUXES
C ------- AS SPECIFIED BY THE PARAMETER KPR.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: FRATE(14), FLOW(14), TFLOW(14), C(NNP), F(NNP,3),
C -------        TIME, DELT, KPR, KOUT, KDIAG, MAXNP,
C -------        NNP, AND ITIM.
C
C ------- OUTPUT: TO PRINT ALL INPUTS EXCEPT FOR KPR.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /TFLOW/ FRATE(14),FLOW(14),TFLOW(14,7)
C
      DIMENSION C(MAXNP),F(MAXNP,3)
      CHARACTER*10 NAME(7)
C
      DATA NAME/'MICROBE 1','MICROBE 2','MICROBE 3',
     >          'SUBSTRATE','OXYGEN   ','NITRATE  ','NUTRIENT  '/
C
      IF(KPR.GT.0) THEN
C ----- print flow through all types of boundaries
        KDIAG=KDIAG+1
        WRITE(16,1000) KKK,NAME(KKK),KDIAG,TIME,DELT,ITIM
        WRITE(16,1100) (FRATE(I),FLOW(I),TFLOW(I,KKK),I=1,14)
      END IF
C
      IF(KPR.GT.1) THEN
C ----- print concentrations
        KOUT=KOUT+1
        LINE=0
        DO 200 NI=1,NNP,4
          NJMN=NI
          NJMX=MIN0(NI+3,NNP)
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2000) KOUT,TIME,DELT,ITIM
          WRITE(16,2100) (NJ,C(NJ),NJ=NJMN,NJMX)
  200   CONTINUE
      END IF
C
      IF(KPR.GT.2) THEN
C ----- print material fluxes
        KOUT=KOUT+1
        LINE=0
        DO 300 NI=1,NNP,2
          NJMN=NI
          NJMX=MIN0(NI+1,NNP)
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0)  WRITE(16,3000) KOUT,TIME,DELT,ITIM
          WRITE(16,3100) (NJ,(F(NJ,I),I=1,3),NJ=NJMN,NJMX)
  300   CONTINUE
      END IF
C
      RETURN
C
 1000 FORMAT('1','**********',I3,'-TH CHEMICAL (',A10,') *************',
     2 '**************************'///
     > ' SYSTEM-FLOW TABLE',I4,'.. AT TIME =',1PD10.3,
     3 ' ,(DELT =',1PD10.3,')',' ITIM=',I5)
 1100 FORMAT(//1X,'TYPE OF FLOW',25X,'RATE',4X,'INC. FLOW',2X,
     1 'TOTAL FLOW'/1X,
     1 ' 1. THROUGH DIRICHLET BOUNDARY NODES . . . ',3(1PD11.3)/1X,
     2 ' 2. THROUGH CAUCHY BOUNDARY NODES . .  . . ',3(1PD11.3)/1X,
     3 ' 3. THROUGH NEUMANN BOUNDARY NODES . . . . ',3(1PD11.3)/1X,
     4 ' 4. THROUGH SEEPAGE NODES . . . . . .  . . ',3(1PD11.3)/1X,
     5 ' 5. THROUGH INFILTRATION NODES .. . .  . . ',3(1PD11.3)/1X,
     6 ' 6. THROUGH UNSPECIFIED NODES(NUMERICAL) . ',3(1PD11.3)/1X,
     7 ' 7. NET FLOW THROUGH ENTIRE BOUNDARY NODES ',3(1PD11.3)/1X,
     8 ' 8. INCREASE IN MATERIAL CONTENT (LIQUID) .',3(1PD11.3)/1X,
     9 ' 9. INCREASE IN MATERIAL CONTENT (SOLID) . ',3(1PD11.3)/1X,
     A '10. RADIOACTIVE LOSSES (LIQUID AND SOLID) .',3(1PD11.3)/1X,
     B '11. LOSS TO COMP. OF SKELTON . . . . . . ..',3(1PD11.3)/1X,
     C '12. LOSS THROUGH DISSOLVED PHASE . . . . ..',3(1PD11.3)/1X,
     D '13. LOSS THROUGH ADSORBED PHASE . . . . . .',3(1PD11.3)/1X,
     E '14. ARTIFICIAL SOURCES/SINKS . . . . . . ..',3(1PD11.3)/1X,
     F ' *** NOTE: (+) = OUT FROM, (-) = INTO THE REGION. '/1X,
     G ' *** RATE (M/T/L), INC. FLOW (M/L), TOTAL FLOW (M/L)'/)
 2000 FORMAT('1OUTPUT TABLE',I4,'.. CONCENTRATIONS(M/L**3) AT TIME =',
     > 1PD12.4   /5X,' (DELT =',1PD12.4,')'    ,' *** ITIME =',I6/1X,
     > 4(' NODE   C(M/L**3)  ')/1X,4(' ---- ------------ ')/)
 2100 FORMAT(1H ,4(I5,1PD13.4,1X))
 3000 FORMAT('1OUTPUT TABLE',I4,'.. MATERIAL FLUX (M/L**2/T) AT TIME =',
     > 1PD12.4/5X,' (DELT =',1PD12.4,')'    ,' *** ITIME =',I6//1X,
     > 2(' NODE     FX         FY         FZ     ')/1X,
     > 2(' ------------------------------------- ')/)
 3100 FORMAT(' ',2(I5,1PD11.3,1PD11.3,1PD11.3,1X))
C
      END
C
C
C
      SUBROUTINE TSTORE(X,IE,C,F,TITLE,NPROB,JTM,TIME)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO STORE PERTINENT QUANTITIES ON AUXILIARY DEVICE FOR FUTURE
C ------- USES, E.G., FOR PLOTTING.  WHAT DEVICE IS TO BE USED MUST BE
C ------- SPECIFIED IN THE JCL.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: X(NNP,3), IE(NEL,9), C(NNP), F(NNP,3),
C -------        TITLE(9), NPROB, NNP, NEL, NTI, AND TIME.
C
C ------- OUTPUT: STORE ALL INPUTS IN LOGICAL UNIT 2.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*70 TITLE
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /CHEM/ MXNCC,NCC
C
      DIMENSION X(MAXNP,3),IE(MAXEL,11)
      DIMENSION C(MAXNP,MXNCC),F(MAXNP,3,MXNCC)
C
      DATA NPPROB/-1/
C
      IF(NPPROB.EQ.(-1)) REWIND 12
      IF(NPPROB.NE.NPROB) THEN
        WRITE(12) TITLE,NPROB,NNP,NEL,NTI
        WRITE(12) ((X(N,I),I=1,3),N=1,NNP),((IE(M,I),M=1,NEL),I=1,11)
        WRITE(12) JTM,TIME,((C(N,K),N=1,NNP),K=1,NCC),(((F(N,I,K),
     >            I=1,3),N=1,NNP),K=1,NCC),(IE(M,11),M=1,NEL)
        NPPROB=NPROB
C
      ELSE 
C
        WRITE(12)JTM,TIME,((C(N,K),N=1,NNP),K=1,NCC),(((F(N,I,K),I=1,3),
     >          N=1,NNP),K=1,NCC),(IE(M,11),M=1,NEL)
C
      ENDIF
      RETURN
      END
