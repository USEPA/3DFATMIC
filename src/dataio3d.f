C $$$$$$$$$$$ DATAIO3D.F
C      9/ 2/94      9:30 am
C
      SUBROUTINE RDATIO
     >     (PROPf,SPP,DINTS,RHOMU,PROPt,RKD,TRANC,X,NNPLR,GNLR,IE,DCOSB,
     >      ISB,NPBB,IMOD,IRXN)
C
C***  READ MATERIAL PROPERTIES, SOIL PRPOPERTIY PARAMETERS, COORDINATE,
C***   SUBREGION DATA, ELEMENT CONNECTIVITY, AND MATERIAL TYPES.
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 GNLR
      CHARACTER DATNAM*1,TITLE*70
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /LGEOM/ LTMXNP,LMXNP,LMXBW,MXREGN,NREGN
      COMMON /NINTR/ KPR0,KDSK0,NSTRf,NSTRt,KSSf,KSSt,IGEOM
      COMMON /FINTE/ NCYLf,NITERf,NPITERf,KSP,KGRAV,IPNTSf
      COMMON /FREAL/ TOLA,TOLB,W,OME,OMI,GRAV,cnstkr
      COMMON /TINTE/ NCMt,KVIt,NITERt,NPITERt,IPNTSt,MICONF,IFLUX
      COMMON /SCMTL/ MAXMAT,MXSPPM,MXMPPM,NMAT,NMPPM,NSPPM
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /CHEM/ MXNCC,NCC
      COMMON /MICROB/ GRATE,YCOEFF,RTARDS,RTARDO,RTARDN,SCOEFF,
     >                ECOEFF,DCOEFF,SATURC,PCOEFF,COFK
C
      DIMENSION PROPf(MAXMAT,MXMPPM),SPP(MXSPPM,MAXMAT,4)
      DIMENSION DINTS(MXNCC),RHOMU(mxncc),PROPt(MAXMAT,MXMPPM)
      DIMENSION RKD(MAXMAT,MXNCC),TRANC(MAXMAT,MXNCC)
      DIMENSION GRATE(4),YCOEFF(4),RTARDS(4),RTARDO(4),RTARDN(4),
     >          SCOEFF(4),ECOEFF(4),DCOEFF(4),SATURC(4),PCOEFF(4)
      DIMENSION X(MAXNP,3),NNPLR(MXREGN),GNLR(LTMXNP,MXREGN)
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES), NPBB(MAXBNP)
      DIMENSION IE(MAXEL,11)
      DIMENSION IEM(8)
C
C ******* DATA SET 5: MATERIAL PROPERTIES
C
      READ(15,1000) DATNAM
      READ(15,*) NMAT,NMPPM,NCC,IRXN
      WRITE(16,2000) NMAT,NMPPM,NCC,IRXN
      IF(IMOD.NE.1)THEN
        WRITE(16,2100)
        DO I=1,NMAT
          READ(15,*) (PROPf(I,J),J=1,NMPPM)
          WRITE(16,2110) I,(PROPf(I,J),J=1,NMPPM)
        ENDDO
        WRITE(16,2120)
        READ(15,*) (DINTS(I),I=1,NCC)
        WRITE(16,2130) (DINTS(I),I=1,NCC)
        WRITE(16,2140)
        READ(15,*) (RHOMU(I),I=1,NCC)
        WRITE(16,2130) (RHOMU(I),I=1,NCC)
      ENDIF
C
      IF(IMOD.NE.10)THEN
        WRITE(16,2200)
        DO I=1,NMAT
          READ(15,*) (PROPt(I,J),J=1,NMPPM)
          WRITE(16,2210) I,(PROPt(I,J),J=1,NMPPM)
          READ(15,*) (RKD(I,J),J=1,NCC)
          WRITE(16,2215)(RKD(I,J),J=1,NCC)
C
C ----- read transformation rate constant
C
          READ(15,*) (TRANC(I,J),J=1,NCC)
          WRITE(16,2220)(TRANC(I,J),J=1,NCC)
        ENDDO
C ******* MICROBE-CHEMICAL INTERACTION CONSTANTS
C
C ----- read maximum specific oxygen-based growth rate for microbe #1,
C -----      maximum specific nitrate-based growth rate for microbe #2,
C -----      maximum specific oxygen-based growth rate for microbe #3,
C ----- and  maximum specific nitrate-based growth rate for microbe #3
C
        READ(15,*) (GRATE(I),I=1,4)
C
C ----- read yield coefficient for microbe #1 utilizing oxygen,
C -----      yield coefficient for microbe #2 utilizing nitrate,
C -----      yield coefficient for microbe #3 utilizing oxygen,
C ----- and  yield coefficient for microbe #3 utilizing nitrate
C
        READ(15,*) (YCOEFF(I),I=1,4)
C
C ----- read retarded substrate saturation constant under aerobic
C            conditions w.r.t. microbes #1,
C            retarded substrate saturation constane under anaerobic
C            conditions w.r.t. microbes #2,
C            retarded substrate saturation constant under aerobic
C            conditions w.r.t. microbes #3,
C       and  retarded substrate saturation constant under anaerobic
C            conditions w.r.t. microbes #3
C
        READ(15,*) (RTARDS(I),I=1,4)
C
C ----- read retarded oxygen saturation constant under aerobic
C            conditions w.r.t. microbes #1,
C -----      retarded oxygen saturation constant under anaerobic
C            conditions w.r.t. microbes #2,
C -----      retarded oxygen saturation constant under aerobic
C            conditions w.r.t. microbes #3,
C ----- and  retarded oxygen saturation constant under anaerobic
C            conditions w.r.t. microbes #3
C
        READ(15,*) (RTARDO(I),I=1,4)
C
C ----- read retarded nutrient saturation constant under aerobic
C            conditions w.r.t. microbes #1,
C -----      retarded nutrient saturation constant under anaerobic
C            conditions w.r.t. microbes #2,
C -----      retarded nutrient saturation constant under aerobic
C            conditions w.r.t. microbes #3,
C ----- and  retarded nutrient saturation constant under anaerobic
C            conditions w.r.t. microbes #3
C
        READ(15,*) (RTARDN(I),I=1,4)
C
C ----- read oxygen-use coefficient for syntheses by microbes #1,
C -----      nitrate-use coefficient for syntheses by microbes #2,
C -----      oxygen-use coefficient for syntheses by microbes #3,
C ----- and  nitrate-use coefficient for syntheses by microbes #3
C
        READ(15,*) (SCOEFF(I),I=1,4)
C
C ----- read oxygen-use coefficient for energy by microbes #1,
C -----      nitrate-use coefficient for energy by microbes #2,
C -----      oxygen-use coefficient for energy by microbes #3,
C ----- and  nitrate-use coefficient for energy by microbes #3
C
        READ(15,*) (ECOEFF(I),I=1,4)
C
C ----- read microbial decay coefficient of aerobic respiration of
C            microbe #1,
C -----      microbial decay coefficient of anaerobic respiration of
C            microbe #2,
C -----      microbial decay coefficient of aerobic respiration of
C            microbe #3,
C ----- and  microbial decay coefficient of anaerobic respiration of
C            microbe #3
C
        READ(15,*) (DCOEFF(I),I=1,4)
C
C ----- read oxygen saturation constant for decay w.r.t. microbe #1,
C -----      nitrate saturation constant for decay w.r.t. microbe #2,
C -----      oxygen saturation constant for decay w.r.t. microbe #3,
C ----- and  nitrate saturation constant for decay w.r.t. microbe #3
C
        READ(15,*) (SATURC(I),I=1,4)
C
C ----- read nutrien-use coefficient for the production of microbe #1
C -----      with aerobic respiration,
C -----      nutrien-use coefficient for the production of microbe #2
C -----      with anaerobic respiration,
C -----      nutrien-use coefficient for the production of microbe #3
C -----      with aerobic respiration,
C ----- and  nutrien-use coefficient for the production of microbe #3
C -----      with anaerobic respiration
C
        READ(15,*) (PCOEFF(I),I=1,4)
C
C ----- read inhibition coefficint
C
        READ(15,*) COFK
C
C ----- write microbe-chemical interaction constant
C
        WRITE(16,2225)(GRATE(I),I=1,4),(YCOEFF(I),I=1,4),(RTARDS(I),
     >              I=1,4),(RTARDO(I),I=1,4),(RTARDN(I),I=1,4),
     >              (SCOEFF(I),I=1,4),(ECOEFF(I),I=1,4),(DCOEFF(I),
     >              I=1,4),(SATURC(I),I=1,4),(PCOEFF(I),I=1,4),COFK
C
      ENDIF
C
C ******* DATA SET 6: SOIL PROPERTIES
C
      READ(15,1000) DATNAM
      READ(15,*) KSP,NSPPM,KCP,RHO,GRAV,VISC
      WRITE(16,2230) KSP,NSPPM,KCP,RHO,GRAV,VISC
C
      IF(KSP.EQ.0) THEN
C
C $$$$$ SOIL PROPERTIES ARE READ BY ANALYTICAL FUNCTIONS
C
        IF(NSPPM.NE.0) THEN
C
C  ------ READ SOIL PROPERTY PARAMETERS
C
          WRITE(16,2240)
          DO 205 I=1,NMAT
            READ(15,*) (SPP(J,I,1),J=1,NSPPM)
            WRITE(16,2245) I,(SPP(J,I,1),J=1,NSPPM)
  205     CONTINUE
C
C ------- PRINT SOIL PROPERTY PARAMETERS
C
          WRITE(16,2250)
          DO 210 I=1,NMAT
            READ(15,*) (SPP(J,I,2),J=1,NSPPM)
            WRITE(16,2255) I,(SPP(J,I,2),J=1,NSPPM)
  210     CONTINUE
C
        ELSE
C
C ------- ERROR IN READING NUMBER OF SOIL PROPERTY PARAMETERS
C
          WRITE(16,2260) NSPPM
          STOP
C
        END IF
C
      ELSE
C
C $$$$$ SOIL PROPERTIES ARE READ IN TABULAR FORM
C
        IF(NSPPM.GT.0) THEN
C
C ------- READ SOIL PROPERTY TABLE
C
          DO 230 I=1,NMAT
             READ(15,*) (SPP(J,I,4),J=1,NSPPM)
  230     CONTINUE
          DO 240 I=1,NMAT
            READ(15,*) (SPP(J,I,1),J=1,NSPPM)
  240     CONTINUE
          DO 250 I=1,NMAT
            READ(15,*) (SPP(J,I,2),J=1,NSPPM)
  250     CONTINUE
          DO 260 I=1,NMAT
            READ(15,*) (SPP(J,I,3),J=1,NSPPM)
  260     CONTINUE
C
C ------- PRINT SOIL PROPERTY TABLE
C
          WRITE(16,2270)
          DO 270 I=1,NMAT
            WRITE(16,2275) I,(SPP(J,I,4),(SPP(J,I,K),K=1,3),J=1,NSPPM)
  270     CONTINUE
C
        ELSE
C
C ------- ERROR IN NUMBER OF SOIL PROPERTY TABLE DATA POINTS
C
          WRITE(16,2280)
          STOP
C
        END IF
C
      END IF
C
      IF(KCP.NE.0) THEN
C
C ####  CONVERT SATURATED PERMEABILITY TO SATURATED CONDUCTIVITY
C
        DO 290 I=1,NMAT
          PKCF=RHO*GRAV/VISC
          DO 285 J=1,6
            PROPf(I,J)=PROPf(I,J)*PKCF
  285     CONTINUE
  290   CONTINUE
      END IF
C
C ******* DATA SET 7: NODE COORDINATES
C
      IF(NSTRf.GT.0)THEN
        REWIND 11
        IF(IPNTSf.EQ.0) THEN
          READ(11) TITLE,NPROB,NNP,NEL,NBNP,NBES,NTIF,NREGN,NVNP
          READ(11) ((X(N,I),I=1,3),N=1,NNP),
     >     ((IE(M,I),M=1,NEL),I=1,11),((DCOSB(I,M),I=1,3),M=1,NBES),
     >     ((ISB(I,M),I=1,6),M=1,NBES),(NPBB(N),N=1,NBNP),(NNPLR(N),N=1,
     >     NREGN),((GNLR(N,I),N=1,LTMXNP),I=1,NREGN)
        ELSE
          READ(11) TITLE,NPROB,NNP,NEL,NBNP,NBES,NTIF,NVNP
          READ(11) ((X(N,I),I=1,3),N=1,NNP),
     1     ((IE(M,I),M=1,NEL),I=1,11),((DCOSB(I,M),I=1,3),M=1,NBES),
     2     ((ISB(I,M),I=1,6),M=1,NBES),(NPBB(N),N=1,NBNP)
        END IF
      ENDIF
      IF(NSTRt.GT.0)THEN
        REWIND 12
        READ(12) TITLE,NPROB,NNP,NEL,NTI
        READ(12) ((X(N,I),I=1,3),N=1,NNP),((IE(M,I),M=1,NEL),I=1,11)
      ENDIF
C
      IF(NSTRF.GT.0 .OR. NSTRT.GT.0)THEN
        PRINT *,'TITLE IN RDATIO =',TITLE
        PRINT *,'NPROB, NTIF=',NPROB,NTIF
      ENDIF
c
      IF(NSTRf.EQ.0 .AND. NSTRt.EQ.0)THEN
        READ(15,1000) DATNAM
        READ(15,*) NNP
        WRITE(16,2300) NNP
C
        NPI=0
  310   READ(15,*) NI,NSEQ,NIAD,XI,YI,ZI,XIAD,YIAD,ZIAD
C
        IF(NI.GT.0) THEN
          NJ=NI+NSEQ
          DO 320 NP=NI,NJ
            I=NI+NIAD*(NP-NI)
            X(I,1)=XI+XIAD*  dble(NP-NI)
            X(I,2)=YI+YIAD*  dble(NP-NI)
            X(I,3)=ZI+ZIAD*  dble(NP-NI)
            NPI=NPI+1
  320     CONTINUE
          GO TO 310
C
        ELSE
          IF(NPI.NE.NNP) THEN
            WRITE(16,2315) NPI,NNP
            STOP
          END IF
C
        END IF
      ENDIF
C
C ------- PRINT NODAL POINT COORDINATES
C
      IF(MOD(IGEOM,2).NE.0) THEN
        LINE=0
        DO 365 NP=1,NNP,2
          NJMN=NP
          NJMX=MIN0(NP+1,NNP)
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2320)
          WRITE(16,2325) (NJ,(X(NJ,I),I=1,3),NJ=NJMN,NJMX)
  365   CONTINUE
      END IF
C
C ******* DATA SET 8: SUBREGION DATA
C
      IF(IPNTSf.EQ.0 .OR. IPNTSt.EQ.0) THEN
        IF(NSTRf.EQ.0 .AND. NSTRt.EQ.0)THEN
          READ(15,1000) DATNAM
          READ(15,*) NREGN
          WRITE(16,2400) NREGN
          CALL READN(NNPLR,MXREGN,NREGN)
C
          IF(MOD(IGEOM,2).NE.0) WRITE(16,2410)
          DO 480 K=1,NREGN
            LNNP=NNPLR(K)
            CALL READN(GNLR(1,K),LTMXNP,LNNP)
            IF(MOD(IGEOM,2).NE.0) WRITE(16,2415) K,LNNP,
     >       (GNLR(I,K),I=1,LNNP)
  480     CONTINUE
        ENDIF
      END IF
C
C ******* DATA SET 9: ELEMENT DATA
C
      IF(NSTRf.EQ.0 .AND. NSTRt.EQ.0)THEN
        READ(15,1000) DATNAM
        READ(15,*) NEL
        WRITE(16,2500) NEL
C
        MMP=0
  500   READ(15,*) MI,NSEQ,MIAD,(IEM(I),I=1,8),IEMAD
C
        IF(MI.GT.0) THEN
          CALL ELENOD
     I        (IEM(5),IEM(7),
     O         NODE,MP,MP)
C
C ####  READ ELEMENT CONNECTIVITY
C
          MJ=MI+NSEQ
          DO 520 MP=MI,MJ
            M=MI+(MP-MI)*MIAD
            DO 510 IQ=1,NODE
              NI=IEM(IQ)+(MP-MI)*IEMAD
              IE(M,IQ)=NI
  510       CONTINUE
            DO IQ=NODE+1,8
              IE(M,IQ)=0
            ENDDO
            MMP=MMP+1
  520     CONTINUE
          GO TO 500
C
        ELSE
C
          IF(MMP.NE.NEL) THEN
            WRITE(16,2510) MMP,NEL
            STOP
          ELSE
            DO M=1,NEL
              IE(M,9)=1
              IE(M,10)=0
              IE(M,11)=0
            ENDDO
          END IF
        ENDIF
C
      END IF
C
      IF(IE(1,5).EQ.0)THEN
        ISHAPE=4
      ELSEIF(IE(1,7).EQ.0)THEN
        ISHAPE=6
      ELSE
        ISHAPE=8
      ENDIF
      DO 560 M=2,NEL
        IF(ISHAPE.EQ.0)GOTO 560
        IF(ISHAPE.EQ.8)THEN
          IF(IE(M,5).EQ.0 .OR. IE(M,7).EQ.0)ISHAPE=0
        ELSEIF(ISHAPE.EQ.6)THEN
          IF(IE(M,5).EQ.0 .OR. IE(M,7).NE.0)ISHAPE=0
        ELSE
          IF(IE(M,5).NE.0)ISHAPE=0
        ENDIF
  560 CONTINUE
C
C ******* DATA SET 10: MATERIAL CORRECTIONS
C
      IF(NSTRf.EQ.0 .AND. NSTRt.EQ.0)THEN
        READ(15,1000) DATNAM
        READ(15,*) NCM
        WRITE(16,2600) NCM
        IF(NCM.NE.0) THEN
C
C ##### READ MATERIAL CORRECTIONS
C
          CALL READN(IE(1,9),MAXEL,NCM)
        END IF
      ENDIF
C
C
      IF(MOD(IGEOM,2).NE.0) THEN
C
C ##### PRINT ELEMENT INCIDENCES AND MATERIAL TYPES
C
        LINE=0
        DO 710 NI=1,NEL
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2710)
          WRITE(16,2715) NI,(IE(NI,K),K=1,9)
  710   CONTINUE
      END IF
C
C ------- CHECK IF MATERIAL TYPE FOR EACH ELEMENT IS CORRECT
C
      DO 820 M=1,NEL
        MTYP=IE(M,9)
        IF(MTYP.GT.0 .AND. MTYP.LE.NMAT) GO TO 820
        WRITE(16,2810) M
        STOP
  820 CONTINUE
C
      RETURN
C
 1000 FORMAT(A1)
C
 2000 FORMAT(' *** MATERIAL PROPERTIES *** '/5X,
     > ' NUMBER OF DIFFERENT MATERIALS, NMAT. . . . . . . . .',I5/5X,
     > ' NUMBER OF MATERIAL PROPERTIES PER MATERIAL, NMPPM  .',I5/5X,
     > ' NUMBER OF COMPONENTS, NCC  . . . . . . . . . . . . .',I5/5X,
     > ' INDEX OF TYPE OF KINETIC REACTION, IRXN. . . . . . .',I5/)
 2100 FORMAT('MAT NO.  SAT KXX    SAT KYY    SAT KZZ    SAT KXY  ',
     > '  SAT KXZ    SAT KXZ  '/5X,
     > '-------  -------    -------    -------    -------  ',
     > '  -------    -------  ')
 2110 FORMAT(5X,I7,6(1PD11.4))
 2120 FORMAT(/' *** INTRINSIC DENSITY OF EACH CHEMICAL ***'/
     > '   RHO1   ','   RHO2   ','   RHO3   ','  RHO4  ',
     > '   RHO5   ','   RHO6   ','   RHO7   ')
 2130 FORMAT(7G10.4)
 2140 FORMAT(/' *** DENSITY AND VISCOSITY FUNCTION COEFF.S ***'/
     > '    B1    ','    B2    ','    B3    ','   B4   ',
     > '    B5    ','    B6    ','    B7    ')
 2200 FORMAT('MAT. NO.   RHOB       AL        AT        AM      TAU    '
     4 /1X,'-------- --------  --------  --------  --------  -------- '
     6 //)
 2210 FORMAT(1X,I8, 5(1PD10.3)/1X,8X, 3(1PD10.3))
 2215 FORMAT(1X,'Kd',2X,7(1PD9.3,2X))
 2220 FORMAT('LAMDA',7(1PD9.3,2X))
 2225 FORMAT(//'1  **** MICROBE-CHEMICAL INTERACTION CONSTANTS ****'/
     >     10X,'MICROBE#1   MICROBE#2   MICROBE#3   MICROBE#3'/
     >     10X,'---------   ---------   ---------   ---------'/
     >     ' GRATE',4X,4(1PD9.3,3X)/' YCOEFF',3X,4(1PD9.3,3X)/
     >     ' RTARDS',3X,4(1PD9.3,3X)/' RTARDO',3X,4(1PD9.3,3X)/
     >     ' RTARDN',3X,4(1PD9.3,3X)/' SCOEFF',3X,4(1PD9.3,3X)/
     >     ' ECOEFF',3X,4(1PD9.3,3X)/' DCOEFF',3X,4(1PD9.3,3X)/
     >     ' SATURC',3X,4(1PD9.3,3X)/' PCOEFF',3X,4(1PD9.3,3X)/
     >     ' COFK',5X,1PD9.3)
C
 2230 FORMAT('1','**** SOIL PROPERTY PARAMETERS ****'//5X,
     > 'SOIL PROPERTY TABULAR INPUT CONTROL, KSP . . . . . .',I5/5X,
     > 'NUMBER OF SOIL PROPERTY PARAMETERS, NSPPM  . . . . .',I5/5X,
     > 'PERMEABILITY INPUT CONTROL, KCP  . . . . . . . . . .',I5/5X,
     > ' DENSITY OF WATER, RHO . . . . . . . . . . . . . .',1PD15.6/5X,
     > ' ACCELERATION OF GRAVITY, GRAV . . . . . . . . . .',1PD15.6/5X,
     > ' VISCOSITY OF WATER, VISC. . . . . . . . . . . . .',1PD15.6/)
 2240 FORMAT('0',' *** MOISTURE-CONTENT PARAMETERS'/5X,
     > 'MAT NO.   PARM1      PARM2      PARM3      PARM4      PARM5'/5X,
     > '-------   -----      -----      -----      -----      -----'/)
 2245 FORMAT(5X,I7,6(1PD11.4))
 2250 FORMAT('0','*** CONDUCTIVITY PARAMETERS'//5X,
     > 'MAT NO.   PARM1      PARM2      PARM3      PARM4      PARM5'/5X,
     > '-------   -----      -----      -----      -----      -----'/)
 2255 FORMAT(5X,I7,6(1PD11.4))
 2260 FORMAT('1','*** ERROR IN READING SOIL PROPERTIY PARAMETERS:',
     > ' NSPPM =',I5,'  STOP ***')
 2270 FORMAT('0 *** SOIL PROPERTY INTERPOLATION VALUES'//1X,
     > 'MAT. NO.    PRESSURE     MOIST. CONT.   REL. CONDUCT.',
     > '   WATER CAP.  '/1X,
     > '-------     --------     ------------   -------------',
     > '   ----------  ')
 2275 FORMAT(1X,I7,4(1PD15.6))
 2280 FORMAT('1','*** ERROR IN READING SOIL PROPERTY DATA:',
     > ' NSPPM =',I5,'  STOP ***')
C
 2300 FORMAT('1  **** NODAL COORDINATE DATA **** '//1X,
     > 'NO. OF NODAL POINTS, NNP . . . . . . . . . . . . . .',I5/)
 2315 FORMAT('1','*** ERROR IN READING NODE COORDINATE: ',
     > ' NPI =',I5,'  .NE.  NNP =',I5,':  STOP')
 2320 FORMAT('0'/1X,
     > 2(' NODE     X          Y          Z      ')/1X,
     > 2(' ---- ---------  ---------  ---------  '))
 2325 FORMAT(1X,2(I5,1PD11.3,1PD11.3,1PD11.3,1X))
C
 2400 FORMAT('1'/5X,'**** SUBREGION DATA ****'/5X,
     > 'NO. OF SUBRETION, NREGN  . . . . . . . . . . . . . .',I5//)
 2410 FORMAT('1'/5X,' $$$ OUTPUT GNLR(I,K) $$$ '/)
 2415 FORMAT('0',5X,' ---- SUBREGION NUMBER K =',I4,'  LNNP =',I5/
     > (6X,10I5))
C
 2500 FORMAT('1',' ***** ELEMENT DATA *****'/5X,
     > 'NO. OF ELEMENTS, NEL . . . . . . . . . . . . . . . .',I5//)
 2510 FORMAT(////'ERROR IN READING IE, MMP =',I5,' NEL =',I5,' STOP')
 2600 FORMAT('0',' ***** MATERIAL CORRECTION DATA ***'/5X,
     > 'NO. OF ELEMENTS REQUIRED MATERIAL CORRECTION, NCM  .',I5//)
 2710 FORMAT('0  **** GLOBAL INDICES OF ELEMENT NODES ****'//5X,
     1 '  ELM NOD1 NOD2 NOD3 NOD4 NOD5 NOD6 NOD7 NOD8 MTYP'/5X,
     2 '  --- ---- ---- ---- ---- ---- ---- ---- ---- ----')
 2715 FORMAT(5X,10I5)
C
 2810 FORMAT(////'   ERROR IN MATERIAL TYPE CODE FOR ELEMENT',I5///)
C
      END
C
C
      SUBROUTINE ELMAKE(NEL,IE,MAXEL,NO2)
C
C CCCC THE FOLLOWING SECTION IS USED FOR CREATING TRIANGULAR PRISM AND
C      TETRAHEDRAL ELEMENTS
C
      DIMENSION IE(MAXEL,NO2)
C
      DO 530 M=1,NEL
C ***** for triangular prism
c       MM=NEL+M
c       IE(MM,1)=IE(M,2)
c       IE(MM,2)=IE(M,3)
c       IE(MM,3)=IE(M,4)
c       IE(MM,4)=IE(M,6)
c       IE(MM,5)=IE(M,7)
c       IE(MM,6)=IE(M,8)
c       IE(MM,7)=0
c       IE(MM,8)=0
c       IE(MM,9)=IE(M,9)
c       MM=M
c       IE(MM,1)=IE(M,1)
c       IE(MM,2)=IE(M,2)
c       IE(MM,3)=IE(M,4)
c       IE(MM,4)=IE(M,5)
c       IE(MM,5)=IE(M,6)
c       IE(MM,6)=IE(M,8)
c       IE(MM,7)=0
c       IE(MM,8)=0
c       IE(MM,9)=IE(M,9)
C
C ***** for tetrahedral
        MM=NEL+M
        IE(MM,1)=IE(M,1)
        IE(MM,2)=IE(M,2)
        IE(MM,3)=IE(M,4)
        IE(MM,4)=IE(M,5)
        IE(MM,5)=0
        IE(MM,6)=0
        IE(MM,7)=0
        IE(MM,8)=0
        IE(MM,9)=IE(M,9)
        MM=2*NEL+M
        IE(MM,1)=IE(M,2)
        IE(MM,2)=IE(M,6)
        IE(MM,3)=IE(M,4)
        IE(MM,4)=IE(M,5)
        IE(MM,5)=0
        IE(MM,6)=0
        IE(MM,7)=0
        IE(MM,8)=0
        IE(MM,9)=IE(M,9)
        MM=3*NEL+M
        IE(MM,1)=IE(M,6)
        IE(MM,2)=IE(M,8)
        IE(MM,3)=IE(M,4)
        IE(MM,4)=IE(M,5)
        IE(MM,5)=0
        IE(MM,6)=0
        IE(MM,7)=0
        IE(MM,8)=0
        IE(MM,9)=IE(M,9)
        MM=4*NEL+M
        IE(MM,1)=IE(M,2)
        IE(MM,2)=IE(M,4)
        IE(MM,3)=IE(M,6)
        IE(MM,4)=IE(M,3)
        IE(MM,5)=0
        IE(MM,6)=0
        IE(MM,7)=0
        IE(MM,8)=0
        IE(MM,9)=IE(M,9)
        MM=5*NEL+M
        IE(MM,1)=IE(M,8)
        IE(MM,2)=IE(M,7)
        IE(MM,3)=IE(M,6)
        IE(MM,4)=IE(M,3)
        IE(MM,5)=0
        IE(MM,6)=0
        IE(MM,7)=0
        IE(MM,8)=0
        IE(MM,9)=IE(M,9)
        MM=M
        IE(MM,1)=IE(M,3)
        IE(MM,2)=IE(M,4)
        IE(MM,3)=IE(M,6)
        IE(MM,4)=IE(M,8)
        IE(MM,5)=0
        IE(MM,6)=0
        IE(MM,7)=0
        IE(MM,8)=0
        IE(MM,9)=IE(M,9)
  530 CONTINUE
c     NEL=2*NEL
      NEL=6*NEL
C
      RETURN
      END
C
C
      SUBROUTINE TSSDAT(SOSF,TSOSF,ISTYP,LES, WSSF,TWSSF,IWTYP,NPW)
C
C********1*********2*********3*********4*********5*********6*********7**
C ------- TO READ AND PRINT SOURCE/SINK
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER DATNAM*1
C
      COMMON /CELS/ MXSEL,MXSPR,MXSDP,NSEL,NSPR,NSDP,KSAI
      COMMON /CNPS/ MXWNP,MXWPR,MXWDP,NWNP,NWPR,NWDP,KWAI
      COMMON /CHEM/ MXNCC,NCC
C
      DIMENSION ISTYP(MXSEL,MXNCC),LES(MXSEL)
      DIMENSION SOSF(MXSDP,MXSPR,2),TSOSF(MXSDP,MXSPR)
      DIMENSION IWTYP(MXWNP,MXNCC),NPW(MXWNP)
      DIMENSION WSSF(MXWDP,MXWPR,2),TWSSF(MXWDP,MXWPR)
C
C ******* DATA SET 11: ELEMNET SOURCE/SINK
C
      READ(15,1000) DATNAM
      PRINT *, 'DATNAM in TSSDAT =',DATNAM
      READ(15,*) NSEL,NSPR,NSDP,KSAI
      WRITE(16,2200) NSEL,NSPR,NSDP,KSAI
C
      IF(NSEL.NE.0) THEN
C
C ------- READ AND WRITE ELEMENT SOURCE/SINK PROFILES
C
        DO 510 I=1,NSPR
        READ(15,*) (TSOSF(J,I),SOSF(J,I,1),SOSF(J,I,2),J=1,NSDP)
        WRITE(16,2210) I
        WRITE(16,2220) (TSOSF(J,I),SOSF(J,I,1),SOSF(J,I,2),J=1,NSDP)
  510   CONTINUE
C
C ------- READ SOURCE TYPE ASSIGNED TO EACH ELEMENT
C
        READ(15,*) (LES(I),I=1,NSEL)
        DO K=1,NCC
          CALL READN(ISTYP(1,K),MXSEL,NSEL)
          WRITE(16,5400)K
          WRITE(16,2230)
          WRITE(16,2235) (I,LES(I),ISTYP(I,K),I=1,NSEL)
        ENDDO
C
      END IF
C
C ******* DATA SET 12: WELL SOURCE/SINK
C
      READ(15,1000) DATNAM
      READ(15,*) NWNP,NWPR,NWDP,KWAI
      WRITE(16,2300) NWNP,NWPR,NWDP,KWAI
C
      IF(NWNP.NE.0) THEN
C
C -------  READ AND WRITE WELL SOURCE/SINK PROFILES
C
      DO 580 I=1,NWPR
      READ(15,*) (TWSSF(J,I),WSSF(J,I,1),WSSF(J,I,2),J=1,NWDP)
      WRITE(16,2310) I
      WRITE(16,2320) (TWSSF(J,I),WSSF(J,I,1),WSSF(J,I,2),J=1,NWDP)
  580 CONTINUE
C
C ------- READ WELL SOURCE/SINK NODES AND TYPE OF PROFILES ASSIGNED TO
C -------- EACH OF NWNP NODES.
C
      READ(15,*) (NPW(I),I=1,NWNP)
      DO K=1,NCC
        CALL READN(IWTYP(1,K),MXWNP,NWNP)
C
C ------- PRINT GLOBAL WELL NODE NUMBERS AND PROFILE TYPE OF WELL NODE
C
        WRITE(16,5400)K
        WRITE(16,2330)
        WRITE(16,2335) (I,NPW(I),IWTYP(I,K),I=1,NWNP)
      ENDDO
      ENDIF
C
      RETURN
C
 1000 FORMAT(A1)
C
 2200 FORMAT('1 *** ELEMENT SOURCE/SINK DATA ***'//5X,
     1 'NO. OF ELEMENT-SOURCE/SINK ELEMENTS . . . . . . . . . .',I5/5X,
     2 'NO. OF ELEMENT-SOURCE/SINK PROFILES . . . . . . . . . .',I5/5X,
     3 'NO. OF DATA POINTS ON ELEMENT-SOURCE/SINK PROFILES  . .',I5/5X,
     4 'ANALYTICAL ELEMENT-SOURCE/SINK INPUT CONTROL  . . . . .',I5/)
 2210 FORMAT('0'/5X,' PROFILE NO.',I2,//2('    TIME        SOSQ  ',
     1 '     SOSC  ')/2('     -----      SOSQ       SOSQC  '))
 2220 FORMAT(' ',2(3D11.3))
 2230 FORMAT('0'//9X,'ELEMENT NUMBER AND PROFILE TYPES OF SOURCE'//5X,
     1 3('    I  LES STYP',5X)/5X,3('    -  --- ----',5X))
 2235 FORMAT(' ',4X,3(3I5,5X))
 2300 FORMAT('1 *** WELL SOURCE/SINK DATA ***'//5X,
     1 'NO. OF WELL-SOURCE/SINK NODAL POINTS  . . . . . . . . .',I5/5X,
     2 'NO. OF WELL-SOURCE/SINK PROFILES  . . . . . . . . . . .',I5/5X,
     3 'NO. OF DATA POINTS ON WELL-SOURCE/SINK PROFILES . . . .',I5/5X,
     4 'ANALYTICAL WELL-SOURCE/SINK INPUT CONTROL . . . . . . .',I5/)
 2310 FORMAT('0'/5X,' PROFILE NO.',I2//2('     TIME       WSSQ  ',
     1 '     WSSC  ')/2('     ----       ----       ----  '))
 2320 FORMAT(' ',2(3D11.3))
 2330 FORMAT('0'//10X,' NODEL NUMBER AND PROFILE TYPE OF SOURCE'//5X,
     1 3('    I  NPW  TYP',5X)/5X,3('    -  ---  ----',5X))
 2335 FORMAT(' ',4X,3(3I5,5X))
 5400 FORMAT('0'//' ',' *** FOR CHEMICAL NO.',I2,' ***')
C
      END
C
C
C
C
      SUBROUTINE FSSDAT(SOSF,TSOSF,ISTYP,LES, WSSF,TWSSF,IWTYP,NPW)
C
C********1*********2*********3*********4*********5*********6*********7**
C ------- TO READ AND PRINT SOURCE/SINK
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER DATNAM*1
C
      COMMON /FCS/ MXSEL,MXSPR,MXSDP,NSEL,NSPR,NSDP,KSAI
      COMMON /FCW/ MXWNP,MXWPR,MXWDP,NWNP,NWPR,NWDP,KWAI
C
      DIMENSION LES(MXSEL),ISTYP(MXSEL)
      DIMENSION SOSF(MXSDP,MXSPR),TSOSF(MXSDP,MXSPR)
      DIMENSION NPW(MXWNP),IWTYP(MXWNP)
      DIMENSION WSSF(MXWDP,MXWPR),TWSSF(MXWDP,MXWPR)
C
C ******* DATA SET 12:  ELEMENT (DISTRIBUTED) SOURCE/SINK
C
      READ(15,1000) DATNAM
      PRINT *, 'DATNAM in FSSDAT =',DATNAM
      READ(15,*) NSEL,NSPR,NSDP,KSAI
      WRITE(16,2200) NSEL,NSPR,NSDP,KSAI
C
      IF(NSEL.NE.0) THEN
C
C $$$$$ READ AND PRINT ELEMENT SOURCE/SINK
C
        DO 210 I=1,NSPR
          READ(15,*) (TSOSF(J,I),SOSF(J,I),J=1,NSDP)
          WRITE(16,2210) I
          WRITE(16,2220) (TSOSF(J,I),SOSF(J,I),J=1,NSDP)
  210   CONTINUE
C
C ------- READ SOURCE TYPE ASSIGNED TO EACH ELEMENT
C
        READ(15,*) (LES(M),M=1,NSEL)
        CALL READN(ISTYP,MXSEL,NSEL)
C
C ------- PRINT ELEMENT SOURCE/SINK PROFILES AND TYPE
C
        LINE=0
        DO 220 I=1,NSEL,4
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2230)
          NJMN=I
          NJMX=MIN0(I+3,NSEL)
          WRITE(16,2240) (J,LES(J),ISTYP(J),J=NJMN,NJMX)
  220   CONTINUE
C
      END IF
C
C ******* DATA SET 13:  POINT (WELL) SOURCE/SINK
C
      READ(15,1000) DATNAM
      READ(15,*) NWNP,NWPR,NWDP,KWAI
      WRITE(16,2300) NWNP,NWPR,NWDP,KWAI
C
      IF(NWNP.NE.0) THEN
C
C $$$$$ READN AND WRITE POINT (WELL) SOURCE/SINK
C
        DO 310 I=1,NWPR
          READ(15,*) (TWSSF(J,I),WSSF(J,I),J=1,NWDP)
          WRITE(16,2310) I
          WRITE(16,2320) (TWSSF(J,I),WSSF(J,I),J=1,NWDP)
  310   CONTINUE
C
C ------- READ WELL SOURCE/SINK NODES AND TYPE OF PROFILES ASSIGNED TO
C ------- EACH OF NWNP NODES.
C
        READ(15,*) (NPW(I),I=1,NWNP)
        CALL READN(IWTYP,MXWNP,NWNP)
C
C ------- PRINT GLOBAL WELL NODE NUMBERS AND PROFILE TYPE OF WELL NODE
C
        LINE=0
        DO 320 I=1,NWNP,4
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2330)
          NJMN=I
          NJMX=MIN0(I+3,NWNP)
          WRITE(16,2340) (J,NPW(J),IWTYP(J),J=NJMN,NJMX)
  320   CONTINUE
C
      END IF
C
      RETURN
C
 1000 FORMAT(A1)
C
 2200 FORMAT('1 *** ELEMENT (DISTRIBUTED) SOURCE/SINK ***'/5X,
     > 'NO. OF ELEMENT-SOURCE/SINK ELEMENTS, NSEL . . . . .',I5/5X,
     > 'NO. OF ELEMENT-SOURCE/SINK PROFILES, NSPR . . . . .',I5/5X,
     > 'NO. OF DATA POINTS ON ELEMENT-SOURCE/SINK PROFILE .',I5/5X,
     > 'ANALYTICAL ELEMENT-SOURCE/SINK INPUT CONTROL  . . .',I5/)
 2210 FORMAT('0'/5X,'ELEMENT SOURCE/SINK PROFILE NO.',I2/1X,
     > 3(4X,'TIME      SOURCE  ')/1X,3(4X,'----      ------  '))
 2220 FORMAT(' ',3(1PD11.3,1PD11.3))
 2230 FORMAT('0'//10X,'GLOBAL ELEMENT NUMBER AND PROFILE TYPES OF',
     > ' DISTRIBUTED SOURCE/SINK ELEMENT'/1X,
     > 4('    I  LES STYP    ')/1X,4('    -  --- ----    '))
 2240 FORMAT(' ',4(3I5,4X))
C
 2300 FORMAT('0'//5X,' *** WELL (POINT) SOURCE/SINK ***'/5X,
     > 'NO. OF WELL-SOURCE/SINK NODES, NWNP . . . . . . . .',I5/5X,
     > 'NO. OF WELL-SOURCE/SINK PROFILES, NWPR  . . . . . .',I5/5X,
     > 'NO. OF DATA POINTS ON WELL-SOURCE/SINK PROFILE  . .',I5/5X,
     > 'ANALYTICAL WELL-SOURCE/SINK INPUT CONTROL . . . . .',I5/)
 2310 FORMAT('0'/5X,'POINT SOURCE/SINK PROFILE NO.',I2/1X,
     > 3(4X,'TIME      SOURCE  ')/1X,3(4X,'----      -----  '))
 2320 FORMAT(' ',3(1PD11.3,1PD11.3))
 2330 FORMAT('0'//10X,'GLOBAL NODAL NUMBER AND PROFILE TYPE OF WELLS',
     > ' SOUCE/SINK NODES'/1X,
     > 4('    I  NPW WTYP    ')/1X,4('    -  --- ----    '))
 2340 FORMAT(' ',4(3I5,4X))
C
      END
C
C
C
C
      SUBROUTINE TBCDAT(ISB,NPBB,IB,
     > CVBF,TCVBF,IVTYP,ISV,NPVB,CDBF,TCDBF,IDTYP,NPDB,
     > QCBF,TQCBF,ICTYP,ISC,NPCB,QNBF,TQNBF,INTYP,ISN,NPNB)
C
C********1*********2*********3*********4*********5*********6*********7**
C ------- TO READ AND PRINT BOUNDARY CONDITIONS
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER DATNAM*1
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /NINTR/ KPR0,KDSK0,NSTRf,NSTRt,KSSf,KSSt,IGEOM
      COMMON /TVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /TDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
      COMMON /TCBC/ MXCNP,MXCES,MXCPR,MXCDP,NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /TNBC/ MXNNP,MXNES,MXNPR,MXNDP,NNNP,NNES,NNPR,NNDP,KNAI
      COMMON /CHEM/ MXNCC,NCC
C
      DIMENSION ISB(6,MAXBES),NPBB(MAXBNP),IB(MAXNP)
      DIMENSION CVBF(MXVDP,MXVPR),TCVBF(MXVDP,MXVPR),IVTYP(MXVES,MXNCC)
      DIMENSION ISV(5,MXVES),NPVB(MXVNP)
      DIMENSION CDBF(MXDDP,MXDPR),TCDBF(MXDDP,MXDPR),IDTYP(MXDNP,MXNCC)
      DIMENSION NPDB(MXDNP)
      DIMENSION QCBF(MXCDP,MXCPR),TQCBF(MXCDP,MXCPR),ICTYP(MXCES,MXNCC)
      DIMENSION ISC(5,MXCES),NPCB(MXCNP)
      DIMENSION QNBF(MXNDP,MXNPR),TQNBF(MXNDP,MXNPR),INTYP(MXNES,MXNCC)
      DIMENSION ISN(5,MXNES),NPNB(MXNNP)
      DIMENSION NIMI(4),NJMJ(4)
C
C ******* DATA SET 20: VARIABLE BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      PRINT *, 'DATNAM in TBCDAT =',DATNAM
      READ(15,*) NVES,NVNP,NVPR,NVDP,KVAI
      WRITE(16,2100) NVES,NVNP,NVPR,NVDP,KVAI
C
      IF(NVES.GT.0) THEN
C
C ------- READ AND WRITE RUNIN CONCENTRATION PROFILES
C
        WRITE(16,2110)
        DO 610 I=1,NVPR
        READ(15,*) (TCVBF(J,I),CVBF(J,I),J=1,NVDP)
        WRITE(16,2120) I
        WRITE(16,2125) (TCVBF(J,I),CVBF(J,I),J=1,NVDP)
  610   CONTINUE
C
C ------- READ INCOMING CONCENTRATION TYPE ASSIGNED TO EACH OF VB SIDES
C
        DO K=1,NCC
          CALL READN(IVTYP(1,K),MXVES,NVES)
        ENDDO
C
C ------- READ FOUR GLOBAL NODE NUMBER OF ALL VB SIDES
C
        MPI=0
  620   READ(15,*) MI,NSEQ,MIAD,I1,I2,I3,I4,I1AD,I2AD,I3AD,I4AD
        IF(MI.EQ.0) GO TO 630
        MJ=MI+NSEQ
        DO 625 MP=MI,MJ
          I=MI+(MP-MI)*MIAD
          ISV(1,I)=I1+I1AD*(MP-MI)
          ISV(2,I)=I2+I2AD*(MP-MI)
          ISV(3,I)=I3+I3AD*(MP-MI)
          ISV(4,I)=I4+I4AD*(MP-MI)
          MPI=MPI+1
  625   CONTINUE
        GO TO 620
  630   IF(MPI.EQ.NVES) GO TO 635
        WRITE(16,2130)
        STOP
C
C ------- PRINT INPUTTED GLOBAL NODAL NUMBER AND RAINFALL TYPES OF ALL
C ------- VARIABLE BOUNDARY ELEMENT SIDES
C
  635   CONTINUE
        DO 645 K=1,NCC
        WRITE(16,5400)K
        LINE=0
        DO 640 MP=1,NVES,2
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2140)
          NJMN=MP
          NJMX=MIN0(MP+1,NVES)
          WRITE(16,2145) (J,(ISV(I,J),I=1,4),IVTYP(J,K),J=NJMN,NJMX)
  640   CONTINUE
  645   CONTINUE
C
C ------- READ GLOABAL NODAL NUMBER FOR EACH OF ALL VB NODES
C
        CALL READN(NPVB,MXVNP,NVNP)
C ------ ASSIGN BOUNDARY NODE INFORMATION TO IB ARRAY
        DO I=1,NVNP
          NP=NPVB(I)
          IB(NP)=2
        ENDDO
C
C ------ PRINT GLOBAL NODE NUMBER FOR EACH OF ALL VB NODES
C
        LINE=0
        DO 648 I=1,NVNP,6
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2150)
          NJMN=I
          NJMX=MIN0(I+5,NVNP)
          WRITE(16,2155) (J,NPVB(J),J=NJMN,NJMX)
  648   CONTINUE
C
C ------- COMPUTE BOUNDARY SIDE NUMBER FOR ALL VARIABLE BOUNDARY SIDES
C
        DO 659 MI=1,NVES
          NODE1=4
          IF(ISV(4,MI).EQ.0)NODE1=3
          DO 651 IQ=1,NODE1
            NIMI(IQ)=ISV(IQ,MI)
  651     CONTINUE
C
          DO 657 MJ=1,NBES
            NODE2=4
            IF(ISB(4,MJ).EQ.0)NODE2=3
            DO 652 JQ=1,NODE2
              IJ=ISB(JQ,MJ)
              NJMJ(JQ)=NPBB(IJ)
  652       CONTINUE
            IEQ=0
            DO 656 IQ=1,NODE1
              NI=NIMI(IQ)
              DO 653 JQ=1,NODE2
                NJ=NJMJ(JQ)
                IF(NJ.EQ.NI) GO TO 655
  653         CONTINUE
              GO TO 657
  655         IEQ=IEQ+1
  656       CONTINUE
            IF(IEQ.EQ.NODE1 .AND. IEQ.EQ.NODE2) GO TO 658
  657     CONTINUE
C
          WRITE(16,2160) MI
          STOP
C
  658     ISV(5,MI)=MJ
  659   CONTINUE
C
C ------- CHANGE NPVB FROM CONTAINING GLOBAL NODE NUMBER TO
C ------- CONTAINING BOUNDARY NODE NUMBER
C
        DO 669 NP=1,NVNP
          NI=NPVB(NP)
C
          DO 665 I=1,NBNP
            NJ=NPBB(I)
            IF(NJ.NE.NI) GO TO 665
            NII=I
            GO TO 667
  665     CONTINUE
C
        WRITE(16,2170) NP
        STOP
C
  667   NPVB(NP)=NII
  669   CONTINUE
C
C ------- PRINT COMPUTED BOUNDARY NODAL NUMBER FOR ALL VB NODES
C
        LINE=0
        DO 670 I=1,NVNP,6
        LINE=LINE+1
        IF(MOD(LINE-1,50).EQ.0) WRITE(16,2180)
        NJMN=I
        NJMX=MIN0(I+5,NVNP)
        WRITE(16,2185) (J,NPVB(J),J=NJMN,NJMX)
  670   CONTINUE
C
      END IF
C
C ******* DATA SET 21: DIRICHLET BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      READ(15,*) NDNP,NDPR,NDDP,KDAI
      WRITE(16,2200) NDNP,NDPR,NDDP,KDAI
C
      IF(NDNP.GT.0) THEN
C
C ------- READ AND PRINT DIRICHLET CONCENTRATION PROFILE
C
        WRITE(16,2210)
        DO 710 I=1,NDPR
          READ(15,*) (TCDBF(J,I),CDBF(J,I),J=1,NDDP)
          WRITE(16,2220) I
          WRITE(16,2225) (TCDBF(J,I),CDBF(J,I),J=1,NDDP)
  710   CONTINUE
C
C ------- READ GLOBAL NODE NUMBER OF ALL DIRICHLET NODES AND THE TYPE
C
        CALL READN(NPDB,MXDNP,NDNP)
        DO K=1,NCC
          CALL READN(IDTYP(1,K),MXDNP,NDNP)
        ENDDO
C ------ ASSIGN BOUNDARY NODE INFORMATION TO IB ARRAY
        DO I=1,NDNP
          NP=NPDB(I)
          IB(NP)=1
        ENDDO
C
C ------- PRINT GLOBAL NODE NUMBER AND PROFILE TYPES OF DIRICHLET NODES
C
      DO 725 K=1,NCC
        WRITE(16,5400)K
        LINE=0
        DO 720 I=1,NDNP,3
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2230)
          NJMN=I
          NJMX=MIN0(I+2,NDNP)
          WRITE(16,2235) (J,NPDB(J),IDTYP(J,K),J=NJMN,NJMX)
  720   CONTINUE
  725 CONTINUE
C
      END IF
C
C ******* DATA SET 21: CAUCHY BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      READ(15,*) NCES,NCNP,NCPR,NCDP,KCAI
      WRITE(16,2300) NCES,NCNP,NCPR,NCDP,KCAI
C
      IF(NCES.GT.0) THEN
C
        WRITE(16,2310)
        DO 810 I=1,NCPR
          READ(15,*) (TQCBF(J,I),QCBF(J,I),J=1,NCDP)
          WRITE(16,2320) I
          WRITE(16,2325) (TQCBF(J,I),QCBF(J,I),J=1,NCDP)
  810   CONTINUE
C
C ------- READ TYPE OF CAUCHY FLUX ASSIGNED TO CAUCHY SIDES
C
        DO K=1,NCC
          CALL READN(ICTYP(1,K),MXCES,NCES)
        ENDDO
C
C ------- READ FOUR GLOBAL NODE NUMBER OF ALL CAUCHY SIDES
C
        MPI=0
  820   READ(15,*) MI,NSEQ,MIAD,I1,I2,I3,I4,I1AD,I2AD,I3AD,I4AD
        IF(MI.EQ.0) GO TO 830
        MJ=MI+NSEQ
        DO 825 MP=MI,MJ
          I=MI+(MP-MI)*MIAD
          ISC(1,I)=I1+I1AD*(MP-MI)
          ISC(2,I)=I2+I2AD*(MP-MI)
          ISC(3,I)=I3+I3AD*(MP-MI)
          ISC(4,I)=I4+I4AD*(MP-MI)
          MPI=MPI+1
  825   CONTINUE
        GO TO 820
  830   IF(MPI.EQ.NCES) GO TO 835
        WRITE(16,2330)
        STOP
C
C ------- PRINT INPUTTED GLOBAL NODAL NUMBER AND FLUX TYPES OF ALL
C ------- CAUCHY BOUNDARY ELEMENT SIDES
C
  835 CONTINUE
      DO 845 K=1,NCC
        WRITE(16,5400)K
        LINE=0
        DO 840 MP=1,NCES,2
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2340)
          NJMN=MP
          NJMX=MIN0(MP+1,NCES)
          WRITE(16,2345) (J,(ISC(I,J),I=1,4),ICTYP(J,K),J=NJMN,NJMX)
  840   CONTINUE
  845 CONTINUE
C
C ------- READ GLOABAL NODAL NUMBER FOR EACH OF ALL CAUCHY NODES
C
        CALL READN(NPCB,MXCNP,NCNP)
C ------ ASSIGN BOUNDARY NODE INFORMATION TO IB ARRAY
        DO I=1,NCNP
          NP=NPCB(I)
          IB(NP)=3
        ENDDO
C
C ------ PRINT GLOBAL NODE NUMBER FOR EACH OF ALL CAUCHY NODES
C
        LINE=0
        DO 848 I=1,NCNP,6
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2350)
          NJMN=I
          NJMX=MIN0(I+5,NCNP)
          WRITE(16,2355) (J,NPCB(J),J=NJMN,NJMX)
  848   CONTINUE
C
C ------- COMPUTE BOUNDARY SIDE NUMBER FOR ALL CAUCHY SIDES
C
        DO 859 MI=1,NCES
          NODE1=4
          IF(ISC(4,MI).EQ.0)NODE1=3
          DO 851 IQ=1,NODE1
  851     NIMI(IQ)=ISC(IQ,MI)
C
          DO 857 MJ=1,NBES
            NODE2=4
            IF(ISB(4,MJ).EQ.0)NODE2=3
            DO 852 JQ=1,NODE2
              IJ=ISB(JQ,MJ)
              NJMJ(JQ)=NPBB(IJ)
  852       CONTINUE
            IEQ=0
            DO 856 IQ=1,NODE1
              NI=NIMI(IQ)
              DO 853 JQ=1,NODE2
                NJ=NJMJ(JQ)
                IF(NJ.EQ.NI) GO TO 855
  853         CONTINUE
              GO TO 857
  855         IEQ=IEQ+1
  856       CONTINUE
            IF(IEQ.EQ.NODE1 .AND. IEQ.EQ.NODE2) GO TO 858
  857     CONTINUE
C
          WRITE(16,2360) MI
          STOP
C
  858     ISC(5,MI)=MJ
  859   CONTINUE
C
C ------- CHANGE NPCB FROM CONTAINING GLOBAL NODE NUMBER TO
C ------- CONTAINING BOUNDARY NODE NUMBER
C
        DO 869 NP=1,NCNP
          NI=NPCB(NP)
C
          DO 865 I=1,NBNP
            NJ=NPBB(I)
            IF(NJ.NE.NI) GO TO 865
            NII=I
            GO TO 867
  865     CONTINUE
C
          WRITE(16,2370) NP
          STOP
C
  867     NPCB(NP)=NII
  869   CONTINUE
C
C ------- PRINT COMPUTED BOUNDARY NODAL NUMBER FOR ALL CAUCHY NODES
C
        LINE=0
        DO 870 I=1,NCNP,6
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2380)
          NJMN=I
          NJMX=MIN0(I+5,NCNP)
          WRITE(16,2385) (J,NPCB(J),J=NJMN,NJMX)
  870   CONTINUE
C
      END IF
C
C ******* DATA SET 23:  NEUMANN BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      READ(15,*) NNES,NNNP,NNPR,NNDP,KNAI
      WRITE(16,2400) NNES,NNNP,NNPR,NNDP,KNAI
C
      IF(NNES.GT.0) THEN
C
        WRITE(16,2410)
        DO 910 I=1,NNPR
          READ(15,*) (TQNBF(J,I),QNBF(J,I),J=1,NNDP)
          WRITE(16,2420) I
          WRITE(16,2425) (TQNBF(J,I),QNBF(J,I),J=1,NNDP)
  910   CONTINUE
C
C ------- READ NEUMANN FLUX TYPE ASSIGNED TO EACH NEUMANN SIDE
C
        DO K=1,NCC
           CALL READN(INTYP(1,K),MXNES,NNES)
        ENDDO
C
C ------- READ FOUR GLOBAL NODE NUMBER OF ALL NEUMANN SIDES
C
        MPI=0
  920   READ(15,*) MI,NSEQ,MIAD,I1,I2,I3,I4,I1AD,I2AD,I3AD,I4AD
        IF(MI.EQ.0) GO TO 930
        MJ=MI+NSEQ
        DO 925 MP=MI,MJ
          I=MI+(MP-MI)*MIAD
          ISN(1,I)=I1+I1AD*(MP-MI)
          ISN(2,I)=I2+I2AD*(MP-MI)
          ISN(3,I)=I3+I3AD*(MP-MI)
          ISN(4,I)=I4+I4AD*(MP-MI)
          MPI=MPI+1
  925   CONTINUE
        GO TO 920
  930   IF(MPI.EQ.NNES) GO TO 935
        WRITE(16,2430)
        STOP
C
C ------- PRINT INPUTTED GLOBAL NODAL NUMBER AND FLUX TYPES OF ALL
C ------- NEUMANN BOUNDARY ELEMENT SIDES
C
  935 CONTINUE
      DO 945 K=1,NCC
        WRITE(16,5400)K
        LINE=0
        DO 940 MP=1,NNES,2
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2440)
          NJMN=MP
          NJMX=MIN0(MP+1,NNES)
          WRITE(16,2445) (J,(ISN(I,J),I=1,4),INTYP(J,K),J=NJMN,NJMX)
  940   CONTINUE
  945 CONTINUE
C
C ------- READ GLOABAL NODAL NUMBER FOR EACH OF ALL NEUMANN NODES
C
        CALL READN(NPNB,MXNNP,NNNP)
C
        DO NP=1,NNNP
          NPN=NPNB(NP)
          IB(NPN)=4
        ENDDO
C
C ------ PRINT GLOBAL NODE NUMBER FOR EACH OF ALL NEUMANN NODES
C
        LINE=0
        DO 950 I=1,NNNP,6
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2450)
          NJMN=I
          NJMX=MIN0(I+5,NNNP)
          WRITE(16,2455) (J,NPNB(J),J=NJMN,NJMX)
  950   CONTINUE
C
C ------- COMPUTE BOUNDARY SIDE NUMBER FOR ALL NEUMANN SIDES
C
        DO 959 MI=1,NNES
          NODE1=4
          IF(ISN(4,MI).EQ.0)NODE1=3
          DO 951 IQ=1,NODE1
  951     NIMI(IQ)=ISN(IQ,MI)
C
          DO 957 MJ=1,NBES
            NODE2=4
            IF(ISB(4,MJ).EQ.0)NODE2=3
            DO 952 JQ=1,NODE2
              IJ=ISB(JQ,MJ)
              NJMJ(JQ)=NPBB(IJ)
  952       CONTINUE
            IEQ=0
            DO 956 IQ=1,NODE1
              NI=NIMI(IQ)
              DO 953 JQ=1,NODE2
                NJ=NJMJ(JQ)
                IF(NJ.EQ.NI) GO TO 955
  953         CONTINUE
              GO TO 957
  955         IEQ=IEQ+1
  956       CONTINUE
            IF(IEQ.EQ.NODE1 .AND. IEQ.EQ.NODE2) GO TO 958
  957     CONTINUE
C
          WRITE(16,2460) MI
          STOP
C
  958     ISN(5,MI)=MJ
  959   CONTINUE
C
C ------- CHANGE NPNB FROM CONTAINING GLOBAL NODE NUMBER TO
C ------- CONTAINING BOUNDARY NODE NUMBER
C
        DO 969 NP=1,NNNP
          NI=NPNB(NP)
C
          DO 965 I=1,NBNP
            NJ=NPBB(I)
            IF(NJ.NE.NI) GO TO 965
            NII=I
            GO TO 967
  965     CONTINUE
C
        WRITE(16,2470) NP
        STOP
C
  967   NPNB(NP)=NII
  969   CONTINUE
C
C ------- PRINT COMPUTED BOUNDARY NODAL NUMBER FOR ALL CAUCHY NODES
C
        LINE=0
        DO 970 I=1,NNNP,6
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,2480)
          NJMN=I
          NJMX=MIN0(I+5,NNNP)
          WRITE(16,2485) (J,NPNB(J),J=NJMN,NJMX)
  970   CONTINUE
C
      END IF
C
      RETURN
C
 1000 FORMAT(A1)
C
 2100 FORMAT('1 *** VARIABLE BOUNDARY CONDITION DATA ***'/5X,
     1 'NO. OF VARIABLE BOUNDARY ELEMENT SIDES  . . . . . . . .',I5/5X,
     2 'NO. OF VARIABLE BOUNDARY NODAL POINTS . . . . . . . . .',I5/5X,
     3 'NO. OF VARIABLE FLUX PROFILES . . . . . . . . . . . . .',I5/5X,
     4 'NO. OF DATA POINTS ON VARIABLE FLUX PROFILES. . . . . .',I5/5X,
     5 'ANALYTICAL VARIABLE FLUX INPUT CONTROL  . . . . . . . .',I5/)
 2110 FORMAT('0'///10X,'--- VARIABLE-FLUID CONCENTRATION PROFILE ---')
 2120 FORMAT('0'/5X,' PROFILE NO.',I2//3(4X,'TIME',6X,' RCIN ',2X)/
     > 3(4X,' ---- ',6X,' ---- ',2X))
 2125 FORMAT(' ',3(2D11.3))
 2130 FORMAT('1'/5X,' ERROR IN READING VB-ELEMENT-SIDES')
 2140 FORMAT('0'///10X,' --- INPUTTED VB SIDE INFORMATION ---'//5X,
     > 2('   MP  GN1  GN2  GN3  GN4 CTYP',5X)/5X,
     > 2('   --  ---  ---  ---  --- ----',5X))
 2145 FORMAT(' ',4X,2(6I5,5X))
 2150 FORMAT('0'/10X,' --- INPUTTED VARIABLE NODE DATA ---'//5X,
     1 6('    I NPVB  ')/5X,6('    - ----  '))
 2155 FORMAT(' ',4X,6(2I5,2X))
 2160 FORMAT('1'/5X,' CANNOT FIND A BOUNDARY SIDE COINCIDING WITH',
     1 I3,'-TH VARIABLE BOUNDARY SIDE: STOP ***')
 2170 FORMAT('1'/5X,' *** CANNOT FIND A BOUNDARY NODAL NUMBER FOR',
     1 I3,'-TH VARIABLE BOUNDARY NODE: STOP')
 2180 FORMAT('0'//10X,' --- COMPUTTED BOUNDARY NODE NUMBER OF VB NODE'//
     1 5X,6('    I NPVB',2X)/5X,6('    - ----',2X))
 2185 FORMAT(' ',4X,6(2I5,2X))
C
 2200 FORMAT('1 *** DIRICHLET BOUNDARY CONDITIONS ***'/5X,
     1 'NO. OF DIRICHLET NODAL POINTS . . . . . . . . . . . . .',I5/5X,
     2 'NO. OF DIRICHELT CONCENTRATION PROFILES . . . . . . . .',I5/5X,
     3 'NO. OF DATA POINTS ON DIRICHELT CONCENTRATION PROFILES.',I5/5X,
     4 'ANALYTICAL DIRICHLET BV INPUT CONTROL . . . . . . . . .',I5//)
 2210 FORMAT('0'///10X,'--- DIRICHLET CONCENTRATION PROFILE ---')
 2220 FORMAT('0'/5X,' PROFILE NO.',I2,//3(4X,'TIME',6X,' CONC ',2X)/
     > 3(4X,' ---- ',6X,' ---- ',2X))
 2225 FORMAT(' ',3(2D11.3))
 2230 FORMAT('0'//10X,' GLOBAL NODAL NUMBER AND PROFILE TYPE OF ',
     1 'DIRICHLET BOUNDARY NODES'//5X,3('    I NPDB DTYP',5X)/5X,
     2 3('    - ---- ----',5X))
 2235 FORMAT(' ',4X,3(3I5,5X))
C
 2300 FORMAT('1 *** CAUCHY BOUNDARY CONDITIONS ***'/5X,
     1 'NO. OF CAUCHY BOUNDARY ELEMENT SIDES  . . . . . . . . .',I5/5X,
     2 'NO. OF CAUCHY BOUNDARY NODAL POINTS . . . . . . . . . .',I5/5X,
     3 'NO. OF CAUCHY FLUX PROFILES   . . . . . . . . . . . . .',I5/5X,
     4 'NO. OF DATA POINTS ON CAUCHY FLUX PROFILES  . . . . . .',I5/5X,
     5 'ANALYTICAL CAUCHY FLUX INPUT CONTROL  . . . . . . . . .',I5/)
 2310 FORMAT('0'///10X,'--- CAUCHY FLUX PROFILE ---')
 2320 FORMAT('0'/5X,' PROFILE NO.',I2//3(4X,'TIME',6X,' FLUX ',2X)/
     > 3(4X,'----',6X,' ---- ',2X))
 2325 FORMAT(' ',3(2D11.3))
 2330 FORMAT('0',10X,' *** ERROR IN READING CAUCHY BOUNDARY ELEMENT',
     1 ' SIDE: STOP ***')
 2340 FORMAT('0'/10X,' --- INPUTTED CAUCHY SIDE DATA ---'//5X,
     1 2('   MP  GN1  GN2  GN3  GN4 CTYP',5X)/5X,
     2 2('   --  ---  ---  ---  --- ----',5X))
 2345 FORMAT(' ',4X,2(6I5,5X))
 2350 FORMAT('0'/10X,' --- INPUTTED CAUCHY NODE DATA ---'//5X,
     1 6('    I NPCB',2X)/5X,6('    - ----',2X))
 2355 FORMAT(' ',4X,6(2I5,2X))
 2360 FORMAT('1'/5X,' CANNOT FIND A BOUNDARY SIDE COINCIDING WITH',
     1 I3,'-TH CAUCHY SIDE: STOP ***')
 2370 FORMAT('1'/5X,' *** CANNOT FIND A BOUNDARY NODAL NUMBER FOR',
     1 I3,'-TH CAUCHY BOUNDARY NODE:  STOP')
 2380 FORMAT('0'/10X,' --- COMPUTED CAUCHY NODE DATA ---'//5X,
     1 6('    I NPCB',2X)/5X,6('    - ----',2X))
 2385 FORMAT(' ',4X,6(2I5,2X))
C
 2400 FORMAT('1 *** NEUMANN BOUNDARY CONDITIONS ***'/5X,
     1 'NO. OF NEUMANN BOUNDARY ELEMENT SIDES . . . . . . . . .',I5/5X,
     2 'NO. OF NEUMANN BOUNDARY NODAL POINTS  . . . . . . . . .',I5/5X,
     3 'NO. OF NEUMANN FLUX PROFILES  . . . . . . . . . . . . .',I5/5X,
     4 'NO. OF DATA POINTS ON NEUMANN FLUX PROFILES . . . . . .',I5/5X,
     5 'ANALYTICAL NEUMANN FLUX INPUT CONTROL . . . . . . . . .',I5/)
 2410 FORMAT('0'///10X,'--- NEUMANN FLUX PROFILE ---')
 2420 FORMAT('0'/5X,' PROFILE NO.',I2//3(4X,'TIME',6X,' FLUX ',2X)/
     > 3(4X,'----',6X,' ---- ',2X))
 2425 FORMAT(' ',3(2D11.3))
 2430 FORMAT('0',10X,' *** ERROR IN READING NEUMANN BOUNDARY ELEMENT',
     1 ' SIDE: STOP ***')
 2440 FORMAT('0'/10X,' --- INPUTTED NEUMANN SIDE DATA ---'//5X,
     1 2('   MP  GN1  GN2  GN3  GN4 CTYP',5X)/5X,
     2 2('   --  ---  ---  ---  --- ----',5X))
 2445 FORMAT(' ',4X,2(6I5,5X))
 2450 FORMAT('0'/10X,' --- INPUTTED NEUMANN NODE DATA ---'//5X,
     1 6('    I NPCB',2X)/5X,6('    - ----',2X))
 2455 FORMAT(' ',4X,6(2I5,2X))
 2460 FORMAT('1'/5X,' CANNOT FIND A BOUNDARY SIDE COINCIDING WITH',
     1 I3,'-TH NEUMANN SIDE: STOP ***')
 2470 FORMAT('1'/5X,' *** CANNOT FIND A BOUNDARY NODAL NUMBER FOR',
     1 I3,'-TH NEUMANN BOUNDARY NODE: STOP')
 2480 FORMAT('0'/10X,' --- COMPUTED NEUMANN NODE DATA ---'//5X,
     1 6('    I NPCB',2X)/5X,6('    - ----',2X))
 2485 FORMAT(' ',4X,6(2I5,2X))
 5400 FORMAT('0'//' ',' *** FOR CHEMICAL NO.',I2,' ***')
C
      END
C
C
C
C
      SUBROUTINE FBCDAT(ISB,NPBB,
     3 QCBF,TQCBF,ICTYP,ISC,NPCB,  QNBF,TQNBF,INTYP,ISN,NPNB,
     4 QVBF,TQVBF,IVTYP,ISV,NPVB, RSVAB,  HDBF,THDBF,IDTYP,NPDB)
C
C********1*********2*********3*********4*********5*********6*********7**
C ------- TO READ BOUNDARY CONDITIONS
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER DATNAM*1
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /FCBC/ MXCNP,MXCES,MXCPR,MXCDP,NCNP,NCES,NCPR,NCDP,KCAI
      COMMON /FNBC/ MXNNP,MXNES,MXNPR,MXNDP,NNNP,NNES,NNPR,NNDP,KNAI
      COMMON /FVBC/ MXVES,MXVNP,MXVPR,MXVDP,NVES,NVNP,NVPR,NVDP,KVAI
      COMMON /FDBC/ MXDNP,MXDPR,MXDDP,NDNP,NDPR,NDDP,KDAI
C
      DIMENSION ISB(6,MAXBES), NPBB(MAXBNP)
      DIMENSION QCBF(MXCDP,MXCPR),TQCBF(MXCDP,MXCPR),ICTYP(MXCES)
      DIMENSION ISC(5,MXCES),NPCB(MXCNP)
      DIMENSION QNBF(MXNDP,MXNPR),TQNBF(MXNDP,MXNPR),INTYP(MXNES)
      DIMENSION ISN(5,MXNES),NPNB(MXNNP)
      DIMENSION QVBF(MXVDP,MXVPR),TQVBF(MXVDP,MXVPR), IVTYP(MXVES)
      DIMENSION ISV(5,MXVES),NPVB(MXVNP)
      DIMENSION RSVAB(MXVNP,4)
      DIMENSION HDBF(MXDDP,MXDPR),THDBF(MXDDP,MXDPR),IDTYP(MXDNP)
      DIMENSION NPDB(MXDNP)
      DIMENSION NIMI(4),NJMJ(4)
C
C ******* DATA SET 16: RAINFALL/EVAPORATION-SEEPAGE BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      PRINT *, 'DATNAM in FBCDAT =',DATNAM
      READ(15,*) NVES,NVNP,NVPR,NVDP,KVAI
      WRITE(16,6000) NVES,NVNP,NVPR,NVDP,KVAI
C
      IF(NVES.GT.0) THEN
      WRITE(16,6100)
      DO 610 I=1,NVPR
         READ(15,*) (TQVBF(J,I),QVBF(J,I),J=1,NVDP)
      WRITE(16,6150) I
      WRITE(16,6155) (TQVBF(J,I),QVBF(J,I),J=1,NVDP)
  610 CONTINUE
C
C ------- READ RAINFALL/EVAPORATION TYPE ASSIGNED TO EACH RS SIDE
C
      CALL READN(IVTYP,MXVES,NVES)
C
C ------- READ FOUR GLOBAL NODE NUMBER FOR EACH OF ALL VARIABLE
C ------- BOUNDARY ELEMENT SIDES.
C
      MPI=0
  620 READ(15,*) MI,NSEQ,MIAD,I1,I2,I3,I4,I1AD,I2AD,I3AD,I4AD
      IF(MI.EQ.0) GO TO 630
      MJ=MI+NSEQ
      DO 625 MP=MI,MJ
      I=MI+(MP-MI)*MIAD
      ISV(1,I)=I1+(MP-MI)*I1AD
      ISV(2,I)=I2+(MP-MI)*I2AD
      ISV(3,I)=I3+(MP-MI)*I3AD
      ISV(4,I)=I4+(MP-MI)*I4AD
      MPI=MPI+1
  625 CONTINUE
      GO TO 620
  630 IF(MPI.EQ.NVES) GO TO 635
      WRITE(16,6300)
      STOP
C
C ------- PRINT INPUTTED GLOBAL NODAL NUMBER AND RAINFALL TYPES OF ALL
C ------- VARIABLE BOUNDARY ELEMENT SIDES.
C
  635 LINE=0
      DO 640 MP=1,NVES,2
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,6400)
      NJMN=MP
      NJMX=MIN0(MP+1,NVES)
      WRITE(16,6450) (J,(ISV(I,J),I=1,4),IVTYP(J),J=NJMN,NJMX)
  640 CONTINUE
C
C ------- READ GLOBAL NODAL NUMBER FOR EACH OF ALL VARIABLE NODES.
C
      CALL READN(NPVB,MXVNP,NVNP)
C
C ------- READ PONDING DEPTH AND MINIMUM HEAD FOR EACH OF ALL RS NODES
C
      CALL READR(RSVAB(1,1),MXVNP,NVNP)
      CALL READR(RSVAB(1,2),MXVNP,NVNP)
C
C ------- PRINT GLOBAL NODAL NUMBER, PONDING DEPTH AND MINIMUM PRESURE
C ------- RESSURE HEAD FOR ALL VARIABLE BOUNDARY NODES
C
      LINE=0
      DO 645 I=1,NVNP,2
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,6500)
      NJMN=I
      NJMX=MIN0(I+1,NVNP)
      WRITE(16,6550) (J,NPVB(J),(RSVAB(J,K),K=1,2),J=NJMN,NJMX)
  645 CONTINUE
C
C ------- COMPUTE BOUNDARY SIDE NUMBER FOR EACH OF ALL VARIABLE
C ------- BOUNDARY SIDES.
C
      DO 659 MI=1,NVES
      NODE1=4
      IF(ISV(4,MI).EQ.0)NODE1=3
      DO 651 IQ=1,NODE1
  651 NIMI(IQ)=ISV(IQ,MI)
C
      DO 657 MJ=1,NBES
      NODE2=4
      IF(ISB(4,MJ).EQ.0)NODE2=3
      DO 652 JQ=1,NODE2
      IJ=ISB(JQ,MJ)
  652 NJMJ(JQ)=NPBB(IJ)
      IEQ=0
      DO 656 IQ=1,NODE1
      NI=NIMI(IQ)
      DO 653 JQ=1,NODE2
      NJ=NJMJ(JQ)
      IF(NJ.EQ.NI) GO TO 655
  653 CONTINUE
      GO TO 657
  655 IEQ=IEQ+1
  656 CONTINUE
      IF(IEQ.EQ.NODE1 .AND. IEQ.EQ.NODE2) GO TO 658
  657 CONTINUE
C
      WRITE(16,6570) MI
      STOP
  658 ISV(5,MI)=MJ
C
  659 CONTINUE
C
C ------- CHANGE NPVB FROM CONTAINING GLOBAL NODAL NUMBER TO
C ------- CONTAINING BOUNDARY NODAL NUMBER.
C
      DO 669 NP=1,NVNP
      NI=NPVB(NP)
C
      DO 665 I=1,NBNP
      NJ=NPBB(I)
      IF(NJ.NE.NI) GO TO 665
      NII=I
      GO TO 667
  665 CONTINUE
C
      WRITE(16,6670) NP
      STOP
  667 NPVB(NP)=NII
C
  669 CONTINUE
C
C --------- PRINT COMPUTED BOUNDARY NODAL NUMBER FOR ALL VB NODES
C
      LINE=0
      DO 670 I=1,NVNP,6
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,6700)
      NJMN=I
      NJMX=MIN0(I+5,NVNP)
      WRITE(16,6750) (J,NPVB(J),J=NJMN,NJMX)
  670 CONTINUE
C
C ------- CHANGE ISV(I,MP) I=1,4 FROM CONTAINING GLOBAL NODAL NUMBER
C ------- TO CONTAINING COMPRESSED VARIABLE BOUNDARY NODAL NUMBER.
C
      DO 690 MP=1,NVES
      MPB=ISV(5,MP)
      NODE=4
      IF(ISB(4,MPB).EQ.0)NODE=3
      DO 685 IQ=1,NODE
      NB=ISB(IQ,MPB)
      DO 675 I=1,NVNP
      NI=NPVB(I)
      IF(NI.NE.NB) GO TO 675
      NII=I
      GO TO 680
  675 CONTINUE
      WRITE(16,6751) IQ,MP
      STOP
  680 ISV(IQ,MP)=NII
  685 CONTINUE
  690 CONTINUE
C
C ------- PRINT COMPUTED BOUNDARY NODAL NUMBER & SIDE NUMBER AND
C ------- RAINFALL TYPES FOR ALL VB SIDES
C
      LINE=0
      DO 695 MP=1,NVES,2
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,6900)
      NJMN=MP
      NJMX=MIN0(MP+1,NVES)
      WRITE(16,6950) (J,(ISV(I,J),I=1,5),IVTYP(J),J=NJMN,NJMX)
  695 CONTINUE
C
      END IF
C
C ******* DATA SET 17: DIRICHLET BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      READ(15,*) NDNP,NDPR,NDDP,KDAI
      WRITE(16,7000) NDNP,NDPR,NDDP,KDAI
C
      IF(NDNP.GT.0) THEN
      DO 710 I=1,NDPR
      READ(15,*) (THDBF(J,I),HDBF(J,I),J=1,NDDP)
      WRITE(16,7100) I
      WRITE(16,7155) (THDBF(J,I),HDBF(J,I),J=1,NDDP)
  710 CONTINUE
C
C ------- READ GLOBAL NODAL NUMBER OF ALL DIRICHLET NODES AND
C ------- THE TYPE OF TOTAL HEAD ASSIGNED TO EACH OF THEM.
C
      CALL READN(NPDB,MXDNP,NDNP)
      CALL READN(IDTYP,MXDNP,NDNP)
C
C --------- PRINT GLOBAL NODAL NUMBER AND PROFILE OF DIRICHLET NODES
C
      LINE=0
      DO 720 I=1,NDNP,4
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,7200)
      NJMN=I
      NJMX=MIN0(I+3,NDNP)
      WRITE(16,7250) (J,NPDB(J),IDTYP(J),J=NJMN,NJMX)
  720 CONTINUE
C
      END IF
C
C ******* DATA SET 18: CAUCHY BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      READ(15,*) NCES,NCNP,NCPR,NCDP,KCAI
      WRITE(16,8000) NCES,NCNP,NCPR,NCDP,KCAI
C
      IF(NCES.GT.0) THEN
      DO 810 I=1,NCPR
      READ(15,*) (TQCBF(J,I),QCBF(J,I),J=1,NCDP)
      WRITE(16,8100) I
      WRITE(16,8155) (TQCBF(J,I),QCBF(J,I),J=1,NCDP)
  810 CONTINUE
C
C ------- READ CAUCHY FLUX TYPE ASSIGNED TO EACH CAUCHY SIDE
C
      CALL READN(ICTYP,MXCES,NCES)
C
C ------- READ FOUR GLOBAL NODE NUMBER FOR EACH OF ALL CAUCHY SIDES
C
      MPI=0
  820 READ(15,*) MI,NSEQ,MIAD,I1,I2,I3,I4,I1AD,I2AD,I3AD,I4AD
      IF(MI.EQ.0) GO TO 830
      MJ=MI+NSEQ
      DO 825 MP=MI,MJ
      I=MI+(MP-MI)*MIAD
      ISC(1,I)=I1+(MP-MI)*I1AD
      ISC(2,I)=I2+(MP-MI)*I2AD
      ISC(3,I)=I3+(MP-MI)*I3AD
      ISC(4,I)=I4+(MP-MI)*I4AD
      MPI=MPI+1
  825 CONTINUE
      GO TO 820
  830 IF(MPI.EQ.NCES) GO TO 835
      WRITE(16,8300)
      STOP
C
C ------- PRINT INPUTTED GLOBAL NODAL NUMBER AND CAUCHY FLUX TYPES
C ------- FOR ALL CAUCHY BOUNDARY ELEMENT SIDES.
C
  835 LINE=0
      DO 840 MP=1,NCES,2
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,8400)
      NJMN=MP
      NJMX=MIN0(MP+1,NCES)
      WRITE(16,8450) (J,(ISC(I,J),I=1,4),ICTYP(J),J=NJMN,NJMX)
  840 CONTINUE
C
C ------- READ GLOBAL NODAL NUMBER FOR EACH OF ALL CAUCHY NODES.
C
      CALL READN(NPCB,MXCNP,NCNP)
C
C ------- PRINT GLOBAL NODAL NUMBER FOR ALL CAUCHY NODES
C
      LINE=0
      DO 845 I=1,NCNP,6
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,8500)
      NJMN=I
      NJMX=MIN0(I+5,NCNP)
      WRITE(16,8550) (J,NPCB(J),J=NJMN,NJMX)
  845 CONTINUE
C
C ------- COMPUTE BOUNDARY SIDE NUMBER FOR ALL CAUSHY SIDES
C
      DO 859 MI=1,NCES
      NODE1=4
      IF(ISC(4,MI).EQ.0)NODE1=3
      DO 851 IQ=1,NODE1
  851 NIMI(IQ)=ISC(IQ,MI)
C
      DO 857 MJ=1,NBES
      NODE2=4
      IF(ISB(4,MJ).EQ.0)NODE2=3
      DO 852 JQ=1,NODE2
      IJ=ISB(JQ,MJ)
  852 NJMJ(JQ)=NPBB(IJ)
      IEQ=0
      DO 856 IQ=1,NODE1
      NI=NIMI(IQ)
      DO 853 JQ=1,NODE2
      NJ=NJMJ(JQ)
      IF(NJ.EQ.NI) GO TO 855
  853 CONTINUE
      GO TO 857
  855 IEQ=IEQ+1
  856 CONTINUE
      IF(IEQ.EQ.NODE1 .AND. IEQ.EQ.NODE2) GO TO 858
  857 CONTINUE
C
      WRITE(16,8570) MI
      STOP
  858 ISC(5,MI)=MJ
C
  859 CONTINUE
C
C ------- CHANGE NPCB FROM CONTAINING GLOBAL NODAL NUMBER TO
C ------- CONTAINING BOUNDARY NODAL NUMBER.
C
      DO 869 NP=1,NCNP
      NI=NPCB(NP)
C
      DO 865 I=1,NBNP
      NJ=NPBB(I)
      IF(NJ.NE.NI) GO TO 865
      NII=I
      GO TO 867
  865 CONTINUE
C
      WRITE(16,8670) NP
      STOP
  867 NPCB(NP)=NII
C
  869 CONTINUE
C
C --------- PRINT COMPUTED BOUNDARY NODAL NUMBER FOR ALL CAUCHY NODES
C
      LINE=0
      DO 870 I=1,NCNP,6
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,8700)
      NJMN=I
      NJMX=MIN0(I+5,NCNP)
      WRITE(16,8750) (J,NPCB(J),J=NJMN,NJMX)
  870 CONTINUE
C
      END IF
C
C ******* DATA SET 19:  NEUMANN BOUNDARY CONDITIONS
C
      READ(15,1000) DATNAM
      READ(15,*) NNES,NNNP,NNPR,NNDP,KNAI
      WRITE(16,9000) NNES,NNNP,NNPR,NNDP,KNAI
C
      IF(NNES.GT.0) THEN
      DO 910 I=1,NNPR
      READ(15,*) (TQNBF(J,I),QNBF(J,I),J=1,NNDP)
      WRITE(16,9100) I
      WRITE(16,9155) (TQNBF(J,I),QNBF(J,I),J=1,NNDP)
  910 CONTINUE
C
C ------- READ NEUMANN FLUX TYPE ASSIGNED TO EACH NEUMANN SIDE
C
      CALL READN(INTYP,MXNES,NNES)
C
C ------- READ FOUR GLOBAL NODE NUMBER FOR EACH OF ALL NEUMANN SIDES
C
      MPI=0
  920 READ(15,*) MI,NSEQ,MIAD,I1,I2,I3,I4,I1AD,I2AD,I3AD,I4AD
      IF(MI.EQ.0) GO TO 930
      MJ=MI+NSEQ
      DO 925 MP=MI,MJ
      I=MI+(MP-MI)*MIAD
      ISN(1,I)=I1+(MP-MI)*I1AD
      ISN(2,I)=I2+(MP-MI)*I2AD
      ISN(3,I)=I3+(MP-MI)*I3AD
      ISN(4,I)=I4+(MP-MI)*I4AD
      MPI=MPI+1
  925 CONTINUE
      GO TO 920
  930 IF(MPI.EQ.NNES) GO TO 935
      WRITE(16,9300)
      STOP
C
C ------- PRINT INPUTTED GLOBAL NODAL NUMBER AND NEUMANN FLUX TYPES
C ------- FOR ALL NEUMANN BOUNDARY ELEMENT SIDES.
C
  935 LINE=0
      DO 940 MP=1,NNES,2
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,9400)
      NJMN=MP
      NJMX=MIN0(MP+1,NNES)
      WRITE(16,9450) (J,(ISN(I,J),I=1,4),INTYP(J),J=NJMN,NJMX)
  940 CONTINUE
C
C ------- READ GLOBAL NODAL NUMBER FOR EACH OF ALL NEUMANN NODES.
C
      CALL READN(NPNB,MXNNP,NNNP)
C
C ------- PRINT GLOBAL NODAL NUMBER FOR ALL NEUMANN NODES
C
      LINE=0
      DO 945 I=1,NNNP,6
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,9500)
      NJMN=I
      NJMX=MIN0(I+5,NNNP)
      WRITE(16,9550) (J,NPNB(J),J=NJMN,NJMX)
  945 CONTINUE
C
C ------- COMPUTE BOUNDARY SIDE NUMBER FOR EACH OF NEUMANN
C ------- BOUNDARY SIDES.
C
      DO 959 MI=1,NNES
      NODE1=4
      IF(ISN(4,MI).EQ.0)NODE1=3
      DO 951 IQ=1,NODE1
  951 NIMI(IQ)=ISN(IQ,MI)
C
      DO 957 MJ=1,NBES
      NODE2=4
      IF(ISB(4,MJ).EQ.0)NODE2=3
      DO 952 JQ=1,NODE2
      IJ=ISB(JQ,MJ)
  952 NJMJ(JQ)=NPBB(IJ)
      IEQ=0
      DO 956 IQ=1,NODE1
      NI=NIMI(IQ)
      DO 953 JQ=1,NODE2
      NJ=NJMJ(JQ)
      IF(NJ.EQ.NI) GO TO 955
  953 CONTINUE
      GO TO 957
  955 IEQ=IEQ+1
  956 CONTINUE
      IF(IEQ.EQ.NODE1 .AND. IEQ.EQ.NODE2) GO TO 958
  957 CONTINUE
C
      WRITE(16,9570) MI
      STOP
  958 ISN(5,MI)=MJ
C
  959 CONTINUE
C
C ------- CHANGE NPNB FROM CONTAINING GLOBAL NODAL NUMBER TO
C ------- CONTAINING BOUNDARY NODAL NUMBER.
C
      DO 969 NP=1,NNNP
      NI=NPNB(NP)
C
      DO 965 I=1,NBNP
      NJ=NPBB(I)
      IF(NJ.NE.NI) GO TO 965
      NII=I
      GO TO 967
  965 CONTINUE
C
      WRITE(16,9670) NP
      STOP
  967 NPNB(NP)=NII
C
  969 CONTINUE
C
C --------- PRINT COMPUTED BOUNDARY NODAL NUMBER FOR ALL NEUMANN NODES
C
      LINE=0
      DO 970 I=1,NNNP,6
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,9700)
      NJMN=I
      NJMX=MIN0(I+5,NNNP)
      WRITE(16,9750) (J,NPNB(J),J=NJMN,NJMX)
  970 CONTINUE
C
      END IF
C
      RETURN
C
 1000 FORMAT(A1)
C
 6000 FORMAT('1',5X,' **** RAINFALL-SEEPAGE BOUNDARY CONDITIONS ***'/5X,
     1 'NO. OF VARIABLE BOUNDARY ELEMENT SIDES, NVES . . . ',I5/5X,
     2 'NO. OF VARIABLE BOUNDARY NODAL POINTS, NVNP  . . . ',I5/5X,
     3 'NO. OF RAINFALL PROFILES, NVPR . . . . . . . . . . ',I5/5X,
     4 'NO. OF DATA POINTS ON RAINFALL PROFILES, NVDP  . . ',I5/5X,
     5 'ANLYTICAL RAINFALL INPUT CONTROL . . . . . . . . . ',I5/)
 6100 FORMAT('0'/10X,' --- RAINFALL PROFILE ---')
 6150 FORMAT('0'/5X,' PROFILE NO.',I2/1X,3(4X,'TIME       RAINS  ')/1X,
     > 3(4X,'----      ------  '))
 6155 FORMAT(' ',3(1PD11.3,1PD11.3))
 6300 FORMAT('1','** ERROR READING RAINFALL-SEEPAGE ELEMENT SIDE: STOP')
 6400 FORMAT('0'/10X,' --- INPUTTED VARIABLE SIDE DATA ---'//5X,
     2 2('   MP  GN1  GN2  GN3  GN4 RTYP',5X)/5X,
     3 2('   --  ---  ---  ---  --- ----',5X))
 6450 FORMAT(' ',4X,2(6I5,5X))
 6500 FORMAT('0'/10X,' --- INPUTTED VARIABLE NODE DATA ---'//1X,
     1 2(1X,'    I NPVB     HCON        HMIN   ',1X)/1X,
     2 2(1X,'    - ----     ----        ----   ',1X))
 6550 FORMAT(' ',2(1X,2I5,2D12.4,1X))
 6570 FORMAT('1','*** CANNOT FIND A BOUNDARY SIDE COINCIDING'/1X,
     1 'WITH',I3,'-TH VARIABLE BOUNDARY SIDE: STOP ***')
 6670 FORMAT('1',' *** CANNOT FIND A BOUNDARY NODAL NUMBER FOR'/1X,
     1 I3,'-TH VARIABLE BOUNDARY NODE: STOP ***')
 6700 FORMAT('0'//10X,'COMPUTED BOUNDARY NODAL NUMBER OF ALL VB NODES'
     1 //5X,6('    I NPVB',2X)/5X,6('    - ----',2X))
 6750 FORMAT(' ',4X,6(2I5,2X))
 6751 FORMAT('0',' *** CAN NOT FIND A COMPRESSED RS NODE FOR'/1X,
     1 I2,'-TH POINT OF',I4,'-TH RS SIDE: STOP ***')
 6900 FORMAT('0'/10X,' --- COMPUTED VB SIDE DATA ---'//1X,
     1 2('   MP CNP1 CNP2 CNP3 CNP4  MPB RTYP',1X)/1X,
     2 2('   -- ---- ---- ---- ---- ---- ----',1X))
 6950 FORMAT(' ',2(7I5,1X))
C
 7000 FORMAT('1'/5X,' **** DIRICHLET BOUNDARY CONDITIONS ****'/5X,
     > 'NO. OF DIRICHLET NODES, NDNP . . . . . . . . . . . ',I5/5X,
     > 'NO. OF DIRICHLET PROFILES, NDPR  . . . . . . . . . ',I5/5X,
     > 'NO. OF DATA POINTS ON DIRICHLET PROFILES, NDDP . . ',I5/5X,
     > 'ANALYTICAL DIRICHLET BV INPUT CONTROL  . . . . . . ',I5/)
 7100 FORMAT('0'/5X,' PROFILE NO.',I2/1X,3(4X,'TIME       HEAD   ')/1X,
     > 3(4X,'----      ------  '))
 7155 FORMAT(' ',3(1PD11.3,1PD11.3))
 7200 FORMAT('0'//10X,'GLOBAL NODAL NUMBER AND PROFILE TYPE OF ',
     1 'DIRICHLET BOUNDARY NODES'/1X,4('    I NPDB TYPE',4X)/1X,
     2 4('    - ---- ----',4X))
 7250 FORMAT(' ',4(3I5,4X))
C
 8000 FORMAT('1'/5X,' **** CAUCHY  BOUNDARY CONDITIONS ****'/5X,
     6 'NO. OF CAUCHY BOUNDARY ELEMENT SIDES, NCES . . . . ',I5/5X,
     7 'NO. OF CAUCHY BOUNDARY NODAL POINTS, NCNP  . . . . ',I5/5X,
     8 'NO. OF CAUCHY FLUX PROFILES, NCPR  . . . . . . . . ',I5/5X,
     9 'NO. OF DATA POINTS ON CAUCHY FLUX PROFILES, NCDP . ',I5/5X,
     A 'ANALYTICAL CAUCHY FLUX INPUT CONTROL . . . . . . . ',I5/)
 8100 FORMAT('0'/5X,' PROFILE NO.',I2/1X,3(4X,'TIME       FLUX   ')/1X,
     > 3(4X,'----      ------  '))
 8155 FORMAT(' ',3(1PD11.3,1PD11.3))
 8300 FORMAT('0',10X,'*** ERROR READING CAUCHY BOUNDARY ELEMENT SIDE',
     1 ' SIDE: STOP ***')
 8400 FORMAT('0'/10X,' --- INPUTTED CAUCHY SIDE DATA ---'//5X,
     1 2('   MP  GN1  GN2  GN3  GN4 CTYP',5X)/5X,
     3 2('   --  ---  ---  ---  --- ----',5X))
 8450 FORMAT(' ',4X,2(6I5,5X))
 8500 FORMAT('0'/10X,' --- INPUTTED CAUCHY NODE DATA ---'//5X,
     1 6('    I NPCB',2X)/5X,6('    - ----',2X))
 8550 FORMAT(' ',4X,6(2I5,2X))
 8570 FORMAT('1','*** CANNOT FIND A BOUNDARY SIDE COINCIDING '/1X,
     1 ' WITH',I3,'-TH CAUCY BOUNDARY SIDE:  STOP ***')
 8670 FORMAT('1','*** CANNOT FIND A BOUNDARY NODAL NUMBER FOR'/1X,
     1 I3,'-TH CAUCHY BOUNDARY NODE:  STOP')
 8700 FORMAT('0'/10X,' --- COMPUTED CAUCHY NODE DATA ---'//5X,
     1 6('    I NPCB',2X)/5X,6('    - ----',2X))
 8750 FORMAT(' ',4X,6(2I5,2X))
C
 9000 FORMAT('1'/5X,' **** NEUMANN BOUNDARY CONDITIONS ****'/5X,
     B 'NO. OF NEUMANN BOUNDARY ELEMENT SIDES, NNES  . . . ',I5/5X,
     C 'NO. OF NEUMANN BOUNDARY NODAL POINTS, NNNP . . . . ',I5/5X,
     D 'NO. OF NEUMANN FLUX PROFILES, NNPR . . . . . . . . ',I5/5X,
     E 'NO. OF DATA POINTS ON NEUMANN FLUX PROFILES, NNDP. ',I5/5X,
     F 'ANALYTICAL NEUMANN FLUX INPUT CONTROL  . . . . . . ',I5/)
 9100 FORMAT('0'/5X,' PROFILE NO.',I2/1X,3(4X,'TIME       FLUX   ')/1X,
     > 3(4X,'----      ------  '))
 9155 FORMAT(' ',3(1PD11.3,1PD11.3))
 9300 FORMAT('1',' *** ERROR READING NEUMANN BOUNDARY ELEMENT',
     > ' SIDE: STOP ***')
 9400 FORMAT('0'/10X,' --- INPUTTED NEUMANN SIDE DATA ---'//5X,
     1 2('   MP  GN1  GN2  GN3  GN4 NTYP',5X)/5X,
     2 2('   --  ---  ---  ---  ---  ---',5X))
 9450 FORMAT(' ',4X,2(6I5,5X))
 9500 FORMAT('0'/10X,' --- INPUTTED NEUMANN NODE DATA ---'//5X,
     1 6('    I NPNB',2X)/5X,6('    - ----',2X))
 9550 FORMAT(' ',4X,6(2I5,2X))
 9570 FORMAT('1','*** CANNOT FIND A BOUNDARY SIDE COINCIDING'/1X,
     1 'WITH',I3,'-TH NEUMANN BOUNDARY SIDE:  STOP ***')
 9670 FORMAT('1',' *** CANNOT FIND A BOUNDARY NODAL NUMBER FOR'/1X,
     1 I3,'-TH NEUMANN BOUNDARY NODE: STOP ***')
 9700 FORMAT('0'/10X,' --- COMPUTED NEUMANN NODE DATA ---'//5X,
     1 6('    I NPNB',2X)/5X,6('    - ----',2X))
 9750 FORMAT(' ',4X,6(2I5,2X))
C
      END
C
C
C
C
      SUBROUTINE ESSFCT(PR,TPRF,PRF,T,MXPR,MXDP,NPR,NDP,KANALY)
C
C ------- FIND PROFILE VALUE AT TIME T EITHER BY TABLUR INTERPORATION
C ------- OR BY ANALYTICAL EXPRESSION.  FOR THE LATTER CASE, THE USER
C ------- MUST SUPPLY THE FUNCTION AFTER LINE ESSF 175
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PR(MXPR),TPRF(MXDP,MXPR),PRF(MXDP,MXPR)
C
      IF(KANALY.EQ.1) GO TO 200
C
C ------- THE PROFILE VALUE IS OBTAINED BY INTERPOLATION OF THE INPUT
C ------- TABULAR VALUES
C
      DO 160 I=1,NPR
      DO 140 J=2,NDP
      IF(TPRF(J-1,I).LE.T .AND. T.LE.TPRF(J,I)) GO TO 120
      GO TO 140
  120 RFJM1=PRF(J-1,I)
      TRFJM1=TPRF(J-1,I)
      RFJ=PRF(J,I)
      TRFJ=TPRF(J,I)
      ABC=RFJ-RFJM1
      ABCD=TRFJ-TRFJM1
      PR(I)=RFJM1+(T-TRFJM1)*ABC/ABCD
      GO TO 160
  140 CONTINUE
      PR(I)=0.0
  160 CONTINUE
      RETURN
C
C ------- PROFILE VALUE IS OBTAINED ANALYTICALLY, THE USER MUST SUPPLY
C ------- THE FUNCTION.
C
  200 DO 260 I=1,NPR
      A1=PRF(1,I)
      A2=PRF(2,I)
      A3=PRF(3,I)
      A4=PRF(4,I)
      A5=PRF(5,I)
      PR(I)=A1 + A2*T + A3*T*T + A4*T*T*T + A5*T*T*T*T
  260 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE WSSFCT(PR,TPRF,PRF,T,MXPR,MXDP,NPR,NDP,KANALY)
C
C ------- FIND PROFILE VALUE AT TIME T EITHER BY TABLUR INTERPORATION
C ------- OR BY ANALYTICAL EXPRESSION.  FOR THE LATTER CASE, THE USER
C ------- MUST SUPPLY THE FUNCTION AFTER LINE WSSF 175
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PR(MXPR),TPRF(MXDP,MXPR),PRF(MXDP,MXPR)
C
      IF(KANALY.EQ.1) GO TO 200
C
C ------- THE PROFILE VALUE IS OBTAINED BY INTERPOLATION OF THE INPUT
C ------- TABULAR VALUES
C
      DO 160 I=1,NPR
      DO 140 J=2,NDP
      IF(TPRF(J-1,I).LE.T .AND. T.LE.TPRF(J,I)) GO TO 120
      GO TO 140
  120 RFJM1=PRF(J-1,I)
      TRFJM1=TPRF(J-1,I)
      RFJ=PRF(J,I)
      TRFJ=TPRF(J,I)
      ABC=RFJ-RFJM1
      ABCD=TRFJ-TRFJM1
      PR(I)=RFJM1+(T-TRFJM1)*ABC/ABCD
      GO TO 160
  140 CONTINUE
      PR(I)=0.0
  160 CONTINUE
      RETURN
C
C ------- PROFILE VALUE IS OBTAINED ANALYTICALLY, THE USER MUST SUPPLY
C ------- THE FUNCTION.
C
  200 DO 260 I=1,NPR
      A1=PRF(1,I)
      A2=PRF(2,I)
      A3=PRF(3,I)
      A4=PRF(4,I)
      A5=PRF(5,I)
      PR(I)=A1 + A2*T + A3*T*T + A4*T*T*T + A5*T*T*T*T
  260 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE DBVFCT(PR,TPRF,PRF,T,MXPR,MXDP,NPR,NDP,KANALY)
C
C ------- FIND PROFILE VALUE AT TIME T EITHER BY TABLUR INTERPORATION
C ------- OR BY ANALYTICAL EXPRESSION.  FOR THE LATTER CASE, THE USER
C ------- MUST SUPPLY THE FUNCTION AFTER LINE DBVF 175
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PR(MXPR),TPRF(MXDP,MXPR),PRF(MXDP,MXPR)
C
      IF(KANALY.EQ.1) GO TO 200
C
C ------- THE PROFILE VALUE IS OBTAINED BY INTERPOLATION OF THE INPUT
C ------- TABULAR VALUES
C
      DO 160 I=1,NPR
      DO 140 J=2,NDP
      IF(TPRF(J-1,I).LE.T .AND. T.LE.TPRF(J,I)) GO TO 120
      GO TO 140
  120 RFJM1=PRF(J-1,I)
      TRFJM1=TPRF(J-1,I)
      RFJ=PRF(J,I)
      TRFJ=TPRF(J,I)
      ABC=RFJ-RFJM1
      ABCD=TRFJ-TRFJM1
      PR(I)=RFJM1+(T-TRFJM1)*ABC/ABCD
      GO TO 160
  140 CONTINUE
      PR(I)=0.0
  160 CONTINUE
      RETURN
C
C ------- PROFILE VALUE IS OBTAINED ANALYTICALLY, THE USER MUST SUPPLY
C ------- THE FUNCTION.
C
C
  200 DO 260 I=1,NPR
      A1=PRF(1,I)
      A2=PRF(2,I)
      A3=PRF(3,I)
c     A4=PRF(4,I)
c     A5=PRF(5,I)
c     PR(I)=A1 + A2*T + A3*T*T + A4*T*T*T + A5*T*T*T*T
      pr(i)=a1+a2*DSIN(T*3.141592654/a3)
  260 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE VBVFCT(PR,TPRF,PRF,T,MXPR,MXDP,NPR,NDP,KANALY)
C
C ------- FIND PROFILE VALUE AT TIME T EITHER BY TABLUR INTERPORATION
C ------- OR BY ANALYTICAL EXPRESSION.  FOR THE LATTER CASE, THE USER
C ------- MUST SUPPLY THE FUNCTION AFTER LINE VBVF 175
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PR(MXPR),TPRF(MXDP,MXPR),PRF(MXDP,MXPR)
C
      IF(KANALY.EQ.1) GO TO 200
C
C ------- THE PROFILE VALUE IS OBTAINED BY INTERPOLATION OF THE INPUT
C ------- TABULAR VALUES
C
      DO 160 I=1,NPR
      DO 140 J=2,NDP
      IF(TPRF(J-1,I).LE.T .AND. T.LE.TPRF(J,I)) GO TO 120
      GO TO 140
  120 RFJM1=PRF(J-1,I)
      TRFJM1=TPRF(J-1,I)
      RFJ=PRF(J,I)
      TRFJ=TPRF(J,I)
      ABC=RFJ-RFJM1
      ABCD=TRFJ-TRFJM1
      PR(I)=RFJM1+(T-TRFJM1)*ABC/ABCD
      GO TO 160
  140 CONTINUE
      PR(I)=0.0
  160 CONTINUE
      RETURN
C
C ------- PROFILE VALUE IS OBTAINED ANALYTICALLY, THE USER MUST SUPPLY
C ------- THE FUNCTION.
C
  200 DO 260 I=1,NPR
      A1=PRF(1,I)
      A2=PRF(2,I)
      A3=PRF(3,I)
      A4=PRF(4,I)
      A5=PRF(5,I)
      PR(I)=A1 + A2*T + A3*T*T + A4*T*T*T + A5*T*T*T*T
  260 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE CBVFCT(PR,TPRF,PRF,T,MXPR,MXDP,NPR,NDP,KANALY)
C
C ------- FIND PROFILE VALUE AT TIME T EITHER BY TABLUR INTERPORATION
C ------- OR BY ANALYTICAL EXPRESSION.  FOR THE LATTER CASE, THE USER
C ------- MUST SUPPLY THE FUNCTION AFTER LINE CBVF 175
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PR(MXPR),TPRF(MXDP,MXPR),PRF(MXDP,MXPR)
C
      IF(KANALY.EQ.1) GO TO 200
C
C ------- THE PROFILE VALUE IS OBTAINED BY INTERPOLATION OF THE INPUT
C ------- TABULAR VALUES
C
      DO 160 I=1,NPR
      DO 140 J=2,NDP
      IF(TPRF(J-1,I).LE.T .AND. T.LE.TPRF(J,I)) GO TO 120
      GO TO 140
  120 RFJM1=PRF(J-1,I)
      TRFJM1=TPRF(J-1,I)
      RFJ=PRF(J,I)
      TRFJ=TPRF(J,I)
      ABC=RFJ-RFJM1
      ABCD=TRFJ-TRFJM1
      PR(I)=RFJM1+(T-TRFJM1)*ABC/ABCD
      GO TO 160
  140 CONTINUE
      PR(I)=0.0
  160 CONTINUE
      RETURN
C
C ------- PROFILE VALUE IS OBTAINED ANALYTICALLY, THE USER MUST SUPPLY
C ------- THE FUNCTION.
C
  200 DO 260 I=1,NPR
      A1=PRF(1,I)
      A2=PRF(2,I)
      A3=PRF(3,I)
      A4=PRF(4,I)
      A5=PRF(5,I)
      PR(I)=A1 + A2*T + A3*T*T + A4*T*T*T + A5*T*T*T*T
  260 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE NBVFCT(PR,TPRF,PRF,T,MXPR,MXDP,NPR,NDP,KANALY)
C
C ------- FIND PROFILE VALUE AT TIME T EITHER BY TABLUR INTERPORATION
C ------- OR BY ANALYTICAL EXPRESSION.  FOR THE LATTER CASE, THE USER
C ------- MUST SUPPLY THE FUNCTION AFTER LINE NBVF 175
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PR(MXPR),TPRF(MXDP,MXPR),PRF(MXDP,MXPR)
C
      IF(KANALY.EQ.1) GO TO 200
C
C ------- THE PROFILE VALUE IS OBTAINED BY INTERPOLATION OF THE INPUT
C ------- TABULAR VALUES
C
      DO 160 I=1,NPR
      DO 140 J=2,NDP
      IF(TPRF(J-1,I).LE.T .AND. T.LE.TPRF(J,I)) GO TO 120
      GO TO 140
  120 RFJM1=PRF(J-1,I)
      TRFJM1=TPRF(J-1,I)
      RFJ=PRF(J,I)
      TRFJ=TPRF(J,I)
      ABC=RFJ-RFJM1
      ABCD=TRFJ-TRFJM1
      PR(I)=RFJM1+(T-TRFJM1)*ABC/ABCD
      GO TO 160
  140 CONTINUE
      PR(I)=0.0
  160 CONTINUE
      RETURN
C
C ------- PROFILE VALUE IS OBTAINED ANALYTICALLY, THE USER MUST SUPPLY
C ------- THE FUNCTION.
C
  200 DO 260 I=1,NPR
      A1=PRF(1,I)
      A2=PRF(2,I)
      A3=PRF(3,I)
      A4=PRF(4,I)
      A5=PRF(5,I)
      PR(I)=A1 + A2*T + A3*T*T + A4*T*T*T + A5*T*T*T*T
  260 CONTINUE
C
      RETURN
      END
C
C
C
C
      SUBROUTINE READR(F,MAXNOD,NNP)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO AUTOMATICALLY GENERATE REAL NUMBER INPUT.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION F(MAXNOD)
C
      NODES=0
  150 READ(15,*) NI,NSEQ,NAD,FNI,FAD,FRD
      IF(NI.EQ.0) GO TO 170
      NJ=NI+NSEQ
      DO 160 N=NI,NJ
      NODES=NODES+1
      I=NI+(N-NI)*NAD
      IF(FRD.NE.0.0) GO TO 155
      F(I)=FNI+FAD*  dble(N-NI)
      GO TO 160
  155 I1=I-NAD
      IF(N.EQ.NI) F(I)=FNI
      IF(N.EQ.NI) DINC=1.0D0
      IF(N.GT.NI) DINC=DINC*(1.0D0+FRD)
      IF(N.GT.NI) F(I)=F(I1)+FAD*DINC
  160 CONTINUE
      GO TO 150
  170 IF(NODES.EQ.NNP) GO TO 180
      WRITE(16,1100)
      STOP
  180 IF(NNP.LE.MAXNOD) GO TO 190
      WRITE(16,1200)
      STOP
  190 CONTINUE
C
 1100 FORMAT('1'/'    *** ERROR IN EXECUTING READR SINCE NODES  .NE.',
     > ' NNP:   STOP ***')
 1200 FORMAT('1'/'   *** NNP .GT. MAXNOD  IN EXECUTING READR: STOP')
C
      RETURN
      END
C
C
C
C
      SUBROUTINE READN(INDTYP,MXTYP,NTYPE)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO AUTOMATICALLY GENERATE INTEGER INPUT.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION INDTYP(MXTYP)
C
      NTYPES=0
  110 READ(15,*) NI,NSEQ,NAD,NITYP,NTYPAD
      IF(NI.EQ.0) GO TO 130
      NJ=NI+NSEQ
      DO 120 N=NI,NJ
      I=NI+(N-NI)*NAD
      INDTYP(I)=NITYP + (N-NI)*NTYPAD
      NTYPES=NTYPES+1
  120 CONTINUE
      GO TO 110
  130 IF(NTYPES.EQ.NTYPE) GO TO 140
      WRITE(16,1100)
      STOP
  140 IF(NTYPE.LE.MXTYP) GO TO 150
      WRITE(16,1200)
      STOP
  150 CONTINUE
C
 1100 FORMAT('1'/'   *** ERROR IN EXECUTING READN SINCE NTYPES .NE. ',
     > 'NTYPE: STOP ***')
 1200 FORMAT('1'/'   *** NTYPE .GT. MXTYP IN EXECUTING READN: STOP')
C
      RETURN
      END
C
C
      SUBROUTINE LRL3D
     I                 (IE,NEL,NNP,MAXEL,MAXNP,MXKBD,NO2,
     O                  NLRL,LRL)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,NO2),NLRL(MAXNP),LRL(MXKBD,MAXNP)
C
      DO 10 I=1,NNP
        NLRL(I)=0
        DO K=1,MXKBD
          LRL(K,I)=0
        ENDDO
  10  CONTINUE
      DO 400 M=1,NEL
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,I,I)
        DO 200 I=1,NODE
          IEM=IE(M,I)
          NLRL(IEM)=NLRL(IEM)+1
          NLRLN=NLRL(IEM)
          LRL(NLRLN,IEM)=M
 200    CONTINUE
 400  CONTINUE
      RETURN
      END
      SUBROUTINE LRN3D
     I                 (IE,NEL,NNP,MAXEL,MAXNP,MXADNP,MXJBD,
     O                  NLRN,LRN)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IE(MAXEL,11)
      DIMENSION NLRN(MAXNP),LRN(MXJBD,MXADNP)
C
      DO 10 I=1,NNP
        NLRN(I)=0
        DO J=1,MXJBD
          LRN(J,I)=0
        ENDDO
   10 CONTINUE
C
      DO 400 M=1,NEL
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,I,I)
        DO 200 I=1,NODE
          IEM=IE(M,I)
          DO 100 J=1,NODE
            IEMJ=IE(M,J)
            IF(NLRN(IEM).EQ.0)THEN
              NLRN(IEM)=NLRN(IEM)+1
              LRN(1,IEM)=IEMJ
            ELSE
              DO K=1,NLRN(IEM)
                N=LRN(K,IEM)
                IF(N.EQ.IEMJ)GOTO 100
              ENDDO
              NLRN(IEM)=NLRN(IEM)+1
              NLRNN=NLRN(IEM)
              LRN(NLRNN,IEM)=IEMJ
            ENDIF
 100      CONTINUE
 200    CONTINUE
 400  CONTINUE
      RETURN
      END
      SUBROUTINE SURF(X,IE,LRL,NLRL, DCOSB,ISB,NPBB, IGEOM,LRN,NLRN)
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- TO GENERATE BOUNDARY GEOMETRY.
C
C********1*********2*********3*********4*********5*********6*********7**
C
C ------- INPUT: X(NNP,3), IE(NEL,9), LRL(MXKBD,NNP), NLRL(NNP).
C
C ------- OUTPUT: DCOSB(3,NBES), ISB(6,NBES), NPBB(NBNP), NBES, NBNP.
C
C********1*********2*********3*********4*********5*********6*********7**
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
C
      DIMENSION X(MAXNP,3),IE(MAXEL,11),LRL(MXKBD,MAXNP),NLRL(MAXNP)
      DIMENSION DCOSB(3,MAXBES),ISB(6,MAXBES),NPBB(MAXBNP)
      DIMENSION LRN(MXJBD,MXADNP),NLRN(MAXNP)
C
      DIMENSION IEMI(4)
C
      DIMENSION KGB(4,6,3)
C
      DATA KGB/1,4,8,5, 1,2,6,5, 2,3,7,6, 4,3,7,8, 1,2,3,4, 5,6,7,8,
     >         1,3,6,4, 1,4,5,2, 2,5,6,3, 1,2,3,0, 4,5,6,0, 0,0,0,0,
     >         4,3,2,0, 4,1,3,0, 4,2,1,0, 1,2,3,0, 0,0,0,0, 0,0,0,0/
C
      NBES=0
      NBNP=0
C
C ----- INITIATE NLRN(1..NNP) FOR BEING AS A WORKING ARRAY
C
      DO I=1,NNP
        NLRN(I)=0
      ENDDO
C
      DO 390 MI=1,NEL
        CALL ELENOD
     I      (IE(MI,5),IE(MI,7),
     O       NODE,NSIDE,IK)
        DO 380 LS=1,NSIDE
C
C ------- STORE FOUR GLOBAL NODAL NUMBERS OF LS-TH SIDE OF MI-TH
C ------- ELEMENT IN ARRAY IEMI(4) FOR LATER USE.
C
          DO 120 IQ=1,4
            K=KGB(IQ,LS,IK)
            IF(K.NE.0)THEN
              IEMI(IQ)=IE(MI,K)
            ELSE
              IEMI(IQ)=0
            ENDIF
  120     CONTINUE
C
C ------- CHECK IF THE LS-TH SIDE OF ELEMENT MI IS A BOUNDARY
C ------- SIDE BY LOOPING OVER ALL ELEMENTS CONNECTING TO THE
C ------- FOUR NODES OF THE SIDE TO SEE IF ANY OF THOSE ELEMENTS
C ------- CONTAINS THE SAME FOUR NODES.  IF YES, THEN LS IS NOT
C ------- A BOUNDARY SIDE.  IF NO, THEN LS IS A BOUNDAY SIDE.
C
          NOD1=IEMI(1)
          NOD2=IEMI(2)
          NOD3=IEMI(3)
          DO 220 MJ1=1,NLRL(NOD1)
            MM1=LRL(MJ1,NOD1)
            IF(MM1.EQ.MI)GOTO 220
            DO 210 MJ2=1,NLRL(NOD2)
              MM2=LRL(MJ2,NOD2)
              IF(MM1.NE.MM2)GOTO 210
              DO 200 MJ3=1,NLRL(NOD3)
                MM3=LRL(MJ3,NOD3)
                IF(MM2.EQ.MM3)GOTO 380
  200         CONTINUE
  210       CONTINUE
  220     CONTINUE
C
C ------- AFTER LOOPING OVER ALL ELEMENTS CONNECTED TO THE FOUR NODES
C ------- OF THE LS-TH SIDE OF ELEMENT MI, WE CANNOT FIND ANY OF THOSE
C ------- ELEMENTS CONTAINING THE SAME FOUR NODES OF THE LS-SIDE OF
C ------- OF ELEMENT MI, HENCE THIS SIDE IS A BOUNDARY SIDE.
C
          NBES=NBES+1
          ISB(5,NBES)=LS
          ISB(6,NBES)=MI
C
C ------- COMPUTE DIRECTIONAL COSINES FOR THE NBES-TH SIDE.
C
          NI=IEMI(1)
          NJ=IEMI(2)
          A1=X(NJ,1)-X(NI,1)
          A2=X(NJ,2)-X(NI,2)
          A3=X(NJ,3)-X(NI,3)
          NJ=IEMI(3)
          B1=X(NJ,1)-X(NI,1)
          B2=X(NJ,2)-X(NI,2)
          B3=X(NJ,3)-X(NI,3)
          AB23=A2*B3-A3*B2
          AB31=A3*B1-A1*B3
          AB12=A1*B2-A2*B1
          AREA=DSQRT(AB23*AB23+AB31*AB31+AB12*AB12)
          DCOSB(1,NBES)=AB23/AREA
          DCOSB(2,NBES)=AB31/AREA
          DCOSB(3,NBES)=AB12/AREA
          IF(NODE.EQ.8 .AND. (LS.EQ.2 .OR. LS.EQ.3 .OR. LS.EQ.6))
     >      GOTO 305
          IF(NODE.EQ.6 .AND. LS.EQ.5)GOTO 305
          DCOSB(1,NBES)=-DCOSB(1,NBES)
          DCOSB(2,NBES)=-DCOSB(2,NBES)
          DCOSB(3,NBES)=-DCOSB(3,NBES)
C
  305     CONTINUE
          DO 310 IQ=1,4
              NI=IEMI(IQ)
              IF(NI.EQ.0)GOTO 310
              IF(NLRN(NI).NE.0)GOTO 310
              NBNP=NBNP+1
              NPBB(NBNP)=NI
C ----- NOTE: NLRN IS TO BE USED AS A WORKING ARRAY. IT WILL BE
C             RECOMPUTED AT THE END OF THIS SUBROUTINE
              NLRN(NI)=NBNP
  310     CONTINUE
C
  380   CONTINUE
  390 CONTINUE
C     WRITE(16,*)'NBNP=',NBNP
C     DO NP=1,NNP
C       WRITE(16,*)NP,NLRN(NP)
C     ENDDO
C
C ------ PRINT BOUNDARY NODE INFORMATION
C
      IF(MOD(IGEOM,2).EQ.0)GOTO 396
      LINE=0
      DO 395 I=1,NBNP,6
        LINE=LINE+1
        IF(MOD(LINE-1,50).EQ.0) WRITE(16,3900)
        NJMN=I
        NJMX=MIN0(I+5,NBNP)
        WRITE(16,3950) (J,NPBB(J),J=NJMN,NJMX)
  395 CONTINUE
  396 CONTINUE
C
C ------- COMPUTE THE COMPRESSED BOUNDARY NODE NUMBER FOR EACH OF THE
C ------- FOUR NODES OF A BOUNDARY SIDE
C
      DO 490 MP=1,NBES
        LS=ISB(5,MP)
        M=ISB(6,MP)
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,IQ,IK)
        DO 460 IQ=1,4
          I=KGB(IQ,LS,IK)
          IF(I.NE.0)THEN
            NI=IE(M,I)
          ELSE
            NI=0
          ENDIF
          IF(NI.EQ.0)THEN
            ISB(IQ,MP)=0
            GOTO 460
          ENDIF
          NII=NLRN(NI)
          IF(NII.GT.NBNP)THEN
            WRITE(16,4000) IQ,MP
 4000       FORMAT('0',10X,' CAN NOT FIND A COMPRESSED BOUNDARY NODE',
     1             ' FOR',I2,'-TH POINT OF',I4,'-TH BOUNDARY SIDE',
     2             ' --- STOP')
            STOP
          ENDIF
          ISB(IQ,MP)=NII
  460   CONTINUE
  490 CONTINUE
C
C ------ RECOMPUTE NLRN
C
      DO 500 NP=1,NNP
        DO K=1,MXJBD
          IF(LRN(K,NP).EQ.0)THEN
            NLRN(NP)=K-1
            GOTO 500
          ENDIF
        ENDDO
        NLRN(NP)=MXJBD
  500 CONTINUE
C
C ------- PRINT BOUNDARY SIDE INFORMATION
C
      IF(MOD(IGEOM,2).EQ.0) GO TO 696
      LINE=0
      DO 695 MP=1,NBES
      LINE=LINE+1
      IF(MOD(LINE-1,50).EQ.0) WRITE(16,6900)
      WRITE(16,6950) MP,(DCOSB(I,MP),I=1,3),(ISB(I,MP),I=1,6)
  695 CONTINUE
  696 CONTINUE
C
 3900 FORMAT('1'//10X,' **** COMPUTED BOUNDARY NODE DATA ****'//5X,
     1 6('    I NPBB',2X)/5X,6('    - ----',2X))
 3950 FORMAT(' ',4X,6(2I5,2X))
 6900 FORMAT('1'//10X,' *** COMPUTED BOUNDARY ELEMENT SIDE INFORMATION',
     1 '***'//1X,'   MP     DCOSXB         DCOSYB        DCOSZB   ',
     2 '  BP1  BP2  BP3  BP4   LS    M'/1X,
     3           '   --     ------         ------        ------   ',
     4 '  ---  ---  ---  ---   --    -')
 6950 FORMAT(' ',I5,3D14.6,6I5)
C
      RETURN
      END
C
C
      SUBROUTINE PAGEN(LRN,NLRN,LRL,NLRL,LNOJCN,LMAXDF,NTNPLR,GNLR,IE,
     >                 NNPLR,IGEOM,ND)
C
C $$$$$ TO GENERATE POINTER ARRAYS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*4 GNLR
C
      COMMON /SGEOM/ MAXEL,MAXNP,MXADNP,MAXBES,MXTUBS,MAXBNP,MXJBD,
     >               MXKBD,MXNTI,MXDTC
      COMMON /CGEOM/ NNP,NEL,NBNP,NTUBS,NBES,NTI,NDTCHG,ISHAPE
      COMMON /LGEOM/ LTMXNP,LMXNP,LMXBW,MXREGN,NREGN
      COMMON /NOPTN/ ILUMP,IMID,IWET,IOPTIM,KSORP,LGRN,IQUAR
      COMMON /FINTE/ NCYLF,NITERF,NPITERF,KSP,KGRAV,IPNTSF
      COMMON /TINTE/ NCMt,KVIt,NITERt,NPITERt,IPNTSt,MICONF,IFLUX
C
      DIMENSION LRN(MXJBD,MXADNP),LRL(MXKBD,MAXNP)
      DIMENSION NLRN(MAXNP),NLRL(MAXNP),ND(MAXNP)
      DIMENSION IE(MAXEL,11)
      DIMENSION NTNPLR(MXREGN),NNPLR(MXREGN),LMAXDF(MXREGN)
      DIMENSION GNLR(LTMXNP,MXREGN),LNOJCN(MXJBD,LMXNP,MXREGN)
      DIMENSION KK(100),KJ(100)
C
C ***** GENERATE NLRL(MAXNP), LRL(MXKBD,MAXNP), NLRN(MAXNP), AND
C       LRN(MXJBD,MAXNP) BASED ON IE(1..NEL,1..8)
C
      CALL LRL3D
     I         (IE,NEL,NNP,MAXEL,MAXNP,MXKBD,11,
     O          NLRL,LRL)
      CALL LRN3D
     I         (IE,NEL,NNP,MAXEL,MAXNP,MXADNP,MXJBD,
     O          NLRN,LRN)
C
      IF(IPNTSF.EQ.3 .OR. IPNTST.EQ.3)THEN
C
C ----- REARRANGE LRN IN ACENDING ORDER AND COMPUTE ND IF EITHER
C       IPNTSF=3 OR IPNTST=3
C
        DO 150 NP=1,NNP
          NCOUNT=0
          DO 140 I=1,NLRN(NP)
            NI=LRN(I,NP)
            IF(NCOUNT.EQ.0)THEN
              KK(1)=I
              KJ(1)=NI
              NCOUNT=NCOUNT+1
            ELSE
              DO J=NCOUNT,1,-1
                NJ=KJ(J)
                IF(J.EQ.1)THEN
                  IF(NI.LT.NJ)GOTO 100
                ELSE
                  NK=KJ(J-1)
                  IF(NI.LT.NJ .AND. NI.GT.NK)GOTO 100
                ENDIF
              ENDDO
              NCOUNT=NCOUNT+1
              KK(NCOUNT)=I
              KJ(NCOUNT)=NI
              GOTO 140
  100         CONTINUE
              NCOUNT=NCOUNT+1
              DO JJ=NCOUNT,J+1,-1
                KK(JJ)=KK(JJ-1)
                KJ(JJ)=KJ(JJ-1)
              ENDDO
              KK(J)=I
              KJ(J)=NI
            ENDIF
  140     CONTINUE
C
          DO I=1,NLRN(NP)
            LRN(I,NP)=KJ(I)
            IF(KJ(I).EQ.NP)ND(NP)=I
          ENDDO
C
  150   CONTINUE
      ENDIF
C
      NMAX=0
      NOCUR=0
      DO 195 NP=1,NNP
        NN=NLRL(NP)
        IF(NN.GT.NMAX) THEN
          NMAX=NN
          NOCUR=NP
        ENDIF
  195 CONTINUE
      IF(NMAX.GT.MXKBD) THEN
        WRITE(16,1100) NOCUR,NMAX,MXKBD
        STOP
      ENDIF
C
      NMAX=0
      NOCUR=0
      DO 200 NP=1,NNP
        NN=NLRN(NP)
        IF(NN.GT.NMAX) THEN
          NMAX=NN
          NOCUR=NP
        ENDIF
  200 CONTINUE
      IF(NMAX.GT.MXJBD) THEN
        KONT=NMAX-1
        JBND=MXJBD-1
        WRITE(16,1000) NP,KONT,JBND
        STOP
      ENDIF
C
C ----- PRINT GENERATED ARRAY LRN
C
      IF(MOD(IGEOM,2).NE.0)THEN
        LINE=0
        DO 510 NP=1,NNP
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,5000)
          WRITE(16,5100) NP,(LRN(I,NP),I=1,NLRN(NP))
  510   CONTINUE
      ENDIF
C
C ----- PRINT GENERATED ARRAY LRL
C
      IF(MOD(IGEOM,2).NE.0)THEN
        LINE=0
        DO 530 NP=1,NNP
          LINE=LINE+1
          IF(MOD(LINE-1,50).EQ.0) WRITE(16,5200)
          WRITE(16,5100) NP,(LRL(I,NP),I=1,NLRL(NP))
  530   CONTINUE
      ENDIF
C
C ---------------------------------------------------------------------
C
      IF(IPNTSF.GE.1 .AND. IPNTST.GE.1) RETURN
C
C ***** 1. FILL-UP GNLR(I,MXKR) FROM I=MXNR+1 TO I=MXTNR BASED ON
C          IE(MAXEL,1..8) AND ON GNLR(I,MXKR) FROM I=1 TO I=MXNR,
C ***** 2. GENERATE NTNPLR(MXKR) ALSO BASED ON IE(MAXEL,8) AND ON
C          GNLR(I,MXKR) FROM I=1 TO I=MXNR.
C NOTE: IN THIS BLOCK, NLRL(1..NNP) AND NLRN(1..NNP) ARE USED AS
C       WORKING ARRAYS. THEY WILL BE RECOVERED LATER ON.
C       NLRL(NP)=THE SUBREGION NUMBER THAT GLOBAL NODE NP BELONGS TO.
C       NLRN(NP)=THE SUBREGIONAL NODE NUMBER THAT GLOBAL NODE NP IS
C                RELATED TO.
C
      DO 610 K=1,NREGN
        NTNPLR(K)=NNPLR(K)
        DO I=1,NNPLR(K)
          NP=GNLR(I,K)
          NLRL(NP)=K
          NLRN(NP)=I
        ENDDO
  610 CONTINUE
C
C ----- INITIATE LNOJCN
C
      DO 620 K=1,NREGN
        LNNP=NNPLR(K)
        DO 615 J=1,MXJBD
          DO LI=1,LNNP
            LNOJCN(J,LI,K)=0
          ENDDO
  615   CONTINUE
  620 CONTINUE
C
C ----- START PREPARING POINTER ARRAYS FOR ALL THE SUBREGIONS
C
      DO 690 M=1,NEL
        CALL ELENOD
     I      (IE(M,5),IE(M,7),
     O       NODE,I,I)
C
        DO 680 I=1,NODE
          IEMI=IE(M,I)
          KI=NLRL(IEMI)
          KIJ=NLRN(IEMI)
          DO 670 J=1,NODE
            IEMJ=IE(M,J)
            DO 650 JJ=1,MXJBD
              LNJJ=LNOJCN(JJ,KIJ,KI)
              IF(LNJJ.EQ.0)GOTO 660
              GNJJ=GNLR(LNJJ,KI)
              IF(IEMJ.EQ.GNJJ)GOTO 670
  650       CONTINUE
C
C ----- IT IS NOT POSSIBLE TO FIND NO POINT THE SAME AS GLOBAL NODE
C       IEMJ AFGER WE MXJBD ARE RUN OUT OF ==> ERROR OCCURRED
C
            WRITE(16,*)'ERROR IN PAGEN --- RUN OUT OF MXJBD, CANNOT ',
     >                 'FIND A POINT THE SAME AS IEMJ'
            STOP
C
  660       CONTINUE
            DO NT=1,NTNPLR(KI)
              IF(GNLR(NT,KI).EQ.IEMJ)GOTO 665
            ENDDO
            NTNPLR(KI)=NTNPLR(KI)+1
            NT=NTNPLR(KI)
            GNLR(NT,KI)=IEMJ
C
  665       CONTINUE
            LNOJCN(JJ,KIJ,KI)=NT
  670     CONTINUE
  680   CONTINUE
  690 CONTINUE
C
      IF(IPNTSF.EQ.3 .OR. IPNTST.EQ.3)THEN
C
C ----- REARRANGE LNOJCN IF EITHER IPNTSF=3 OR IPNTST=3
C
        DO 700 K=1,NREGN
          LNNP=NNPLR(K)
          DO 695 I=1,LNNP
            NI=GNLR(I,K)
            NII=NLRN(NI)
            DO J=1,MXJBD
              KJ(J)=LNOJCN(J,NII,K)
            ENDDO
C
            DO 693 J=1,MXJBD
              NJ=LRN(J,NI)
              IF(NJ.EQ.0)GOTO 695
              DO JJ=1,MXJBD
                NJJ=KJ(JJ)
                NJN=GNLR(NJJ,K)
                IF(NJN.EQ.NJ)THEN
                  LNOJCN(J,NII,K)=NJJ
                  GOTO 693
                ENDIF
              ENDDO
  693       CONTINUE
C
  695     CONTINUE
  700   CONTINUE
      ENDIF
C
C ----- RECOVER NLRN AND NLRL
C
      DO 710 NP=1,NNP
        DO K=1,MXJBD
          IF(LRN(K,NP).EQ.0)THEN
            NLRN(NP)=K-1
            GOTO 710
          ENDIF
        ENDDO
        NLRN(NP)=MXJBD
  710 CONTINUE
      DO 720 NP=1,NNP
        DO K=1,MXKBD
          IF(LRL(K,NP).EQ.0)THEN
            NLRL(NP)=K-1
            GOTO 720
          ENDIF
        ENDDO
        NLRL(NP)=MXKBD
  720 CONTINUE
C
C ----- GENERATE LMAXDF(MXKR)
C
      DO 790 K=1,NREGN
        MAXDF=0
        LNNP=NNPLR(K)
        DO 780 LI=1,LNNP
          DO 780 J=1,MXJBD
            LJ=LNOJCN(J,LI,K)
            IF(LJ.GT.LNNP .OR. LJ.EQ.0) GOTO 780
            IDIF=LI-LJ
            IF(IDIF.LT.0) IDIF=LJ-LI
            IF(MAXDF.LT.IDIF) MAXDF=IDIF
  780   CONTINUE
        LMAXDF(K)=MAXDF
  790 CONTINUE
C
C ***** PRINT GENERATED ARRAYS LNOJCN(J,AND,MXNR,MXKR) AND
C        GNLR(MXTNR,MXKR).
C
      IF(MOD(IGEOM,2).EQ.0) GOTO 895
      DO 890 K=1,NREGN
        WRITE(16,8000) K
C
        LNNP=NNPLR(K)
        LNNP1=LNNP+1
        LTNNP=NTNPLR(K)
C
        DO 820 I=1,LNNP
          WRITE(16,8200) I,GNLR(I,K),(LNOJCN(J,I,K),J=1,MXJBD)
  820   CONTINUE
C
        DO 830 I=LNNP1,LTNNP
          WRITE(16,8200) I,GNLR(I,K)
  830   CONTINUE
C
  890 CONTINUE
C
  895 CONTINUE
      WRITE(16,9000) (LMAXDF(K),K=1,NREGN)
C
 1000 FORMAT('0'//5X,' ***',I4,'-TH NODE HAS',I4,' NODES SURROUNDING ',
     1 'IT, WHICH IS MORE THAN MXJBD - 1 =',I5,'  STOP ***')
 1100 FORMAT('0'/5X,'*** NUMBER OF ELEMENTS CONNECTING TO NODE ',I6,
     1 ' IS',I3, ' WHICH IS GREATER THAN MXKBD =',I3,'  STOP')
 5000 FORMAT('1'/1X,'** GENERATED CONNECTING NODES OF ALL NODES *'//1X,
     1 '   NP    1    2    3    4    5    6    7    8    9   10',
     2 '   11   12   13'/1X,9X,'   14   15   16   17   18   19   20',
     3 '   21   22   23   24   25   26   27'/1X,14('   --')/1X,9X,
     4 14('   --')/)
 5100 FORMAT(' ',14I5/1X,9X,14I5)
 5200 FORMAT('1'/1X,' ** GENERATED CONNECTING ELEMENTS OF ALL NODES *'
     1 //1X,'   NP    1    2    3    4    5    6    7    8'/1X,
     4 9('  ---')/)
 8000 FORMAT('1',10X,' *** ARRAY GNLR AND LNOJCN ***'///5X,
     1 ' -- SUBREGION K =',I3,' --'//1X,
     2 ' LNODE GNODE  SURROUNDING AND INCLUDING LOCAL NODES'/1X,
     3 ' ----- -----  -------------------------------------')
 8200 FORMAT(' ',2I6,14I4/17X,13I4)
 9000 FORMAT('0'//5X,' LIST OF HALF BAND WIDTH FOR ALL SUBREGIONS'/5X,
     1 (15I4))
C
      RETURN
      END
