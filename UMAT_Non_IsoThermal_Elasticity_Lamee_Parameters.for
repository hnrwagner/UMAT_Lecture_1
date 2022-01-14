C     ABAQUS Subroutine for Isotropic material
C
C     More Infos can be found at 
C     https://abaqus-docs.mit.edu/2017/English/SIMACAESUBRefMap/simasub-c-umat.htm

C---------------------------------------------------------------------------------------

C Start of Base code, DO NOT change

C---------------------------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
C---------------------------------------------------------------------------------------

C End of Base code

C---------------------------------------------------------------------------------------

C End of Base code

C---------------------------------------------------------------------------------------

C Start of USER code

C---------------------------------------------------------------------------------------

C LOCAL ARRAYS
C ----------------------------------------------------------------
C           EELAS - ELASTIC STRAINS
C           ETHERM - THERMAL STRAINS
C           DTHERM - INCREMENTAL THERMAL STRAINS
C           DELDSE - CHANGE IN STIFFNESS DUE TO TEMPERATURE CHANGE
C ----------------------------------------------------------------
        DIMENSION EELAS(6), ETHERM(6), DTHERM(6), DELDSE(6,6)
C
        PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0)
        
        REAL myE, myNu, myLambda, myMu, myLambda_0, myMu_0
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC THERMO-ELASTICITY WITH LINEARLY VARYING
C MODULI - CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
C           PROPS(1) = Elasticity Modulus - myE - function of Temperatur T_0
C           PROPS(2) = Poisson's Ratio - myNu - function of Temperatur T_0
C           PROPS(3) = Temperatur T_0 
C           PROPS(4) = Elasticity Modulus - myE - function of Temperatur T_1
C           PROPS(5) = Poisson's Ratio - myNu - function of Temperatur T_1
C           PROPS(6) = Temperatur T_1 
C           PROPS(7) = ALPHA - thermal expansion coefficient
C           PROPS(8) = Temperatur T_INITIAL

C           Lamé parameter 1 = myLambda
C           Lamé parameter 2 = Shear Modulus = myMu 

C           ELASTIC PROPERTIES AT START OF INCREMENT

C           DTEMP = Increment of temperature
C           TEMP = Temperature at the start of the increment
C
        FAC1=(TEMP-PROPS(3))/(PROPS(6)-PROPS(3))
        IF (FAC1 .LT. ZERO) FAC1=ZERO
        IF (FAC1 .GT. ONE) FAC1=ONE
        FAC0=ONE-FAC1
        
        myE=FAC0*PROPS(1)+FAC1*PROPS(4)
        myNu=FAC0*PROPS(2)+FAC1*PROPS(5)

        myMu_0= myE / (TWO * (ONE + myNu)) 
        myLambda_0= (myE * myNu) / ((ONE + myNu) * (ONE - TWO * myNu))
C
C           ELASTIC PROPERTIES AT END OF INCREMENT
C
        FAC1=(TEMP+DTEMP-PROPS(3))/(PROPS(6)-PROPS(3))
        IF (FAC1 .LT. ZERO) FAC1=ZERO
        IF (FAC1 .GT. ONE) FAC1=ONE
        FAC0=ONE-FAC1
        
        myE=FAC0*PROPS(1)+FAC1*PROPS(4)
        myNu=FAC0*PROPS(2)+FAC1*PROPS(5)
        
        myMu= myE / (TWO * (ONE + myNu)) 
        myLambda= (myE * myNu) / ((ONE + myNu) * (ONE - TWO * myNu))

C           DDSDDE(NTENS,NTENS)
C           Jacobian matrix of the constitutive model

C           NTENS
C           Size of the stress or strain component array (NDI + NSHR).

C           NDI
C           Number of direct stress components at this point.

C           NSHR
C           Number of engineering shear stress components at this point.

C           DELDSE - CHANGE IN STIFFNESS DUE TO TEMPERATURE CHANGE

C           ELASTIC STIFFNESS AT END OF INCREMENT AND STIFFNESS CHANGE

C           Stiffness Matrix for Isotropic Non-Isothermal Elasticity with Lamé parameters 

        DO K1=1,NDI
            DO K2=1,NDI
                DDSDDE(K2,K1)=myLambda
                DELDSE(K2,K1)=myLambda-myLambda_0
            END DO
            DDSDDE(K1,K1)=TWO * myMu + myLambda
            DELDSE(K1,K1)=TWO * myMu + myLambda - TWO * myMu_0 + myLambda_0
        END DO
        DO K1=NDI+1,NTENS
            DDSDDE(K1,K1)=myMu
            DELDSE(K1,K1)=myMu-myMu_0
        END DO
C
C           CALCULATE THERMAL EXPANSION

C           ETHERM - THERMAL STRAINS
C           DTHERM - INCREMENTAL THERMAL STRAINS

C           PROPS(7) = ALPHA - thermal expansion coefficient
C           PROPS(8) = Temperatur T_INITIAL

C           DTEMP = Increment of temperature
C           TEMP = Temperature at the start of the increment

C
        DO K1=1,NDI
            ETHERM(K1)=PROPS(7)*(TEMP-PROPS(8))
            DTHERM(K1)=PROPS(7)*DTEMP
        END DO

        DO K1=NDI+1,NTENS
            ETHERM(K1)=ZERO
            DTHERM(K1)=ZERO
        END DO
C
C           CALCULATE STRESS, ELASTIC STRAIN AND THERMAL STRAIN
C
C           STRESS(NTENS)
C           This array is passed in as the stress tensor at the beginning of the increment and must be updated
C           in this routine to be the stress tensor at the end of the increment.
C      
C           DSTRAN(NTENS)
C           Array of strain increments.

C           EELAS - ELASTIC STRAINS

C
        DO K1=1, NTENS
            DO K2=1, NTENS
C                STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*(DSTRAN(K1)-DTHERM(K1))+DELDSE(K2,K1)*(STRAN(K1)-ETHERM(K1))
                STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*(DSTRAN(K1)-DTHERM(K1))
            END DO
            ETHERM(K1)=ETHERM(K1)+DTHERM(K1)
            EELAS(K1)=STRAN(K1)+DSTRAN(K1)-ETHERM(K1)
        END DO

C---------------------------------------------------------------------------------------

C End of USER code

C---------------------------------------------------------------------------------------

         RETURN
         END