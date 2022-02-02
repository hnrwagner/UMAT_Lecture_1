C     ABAQUS Subroutine for Plasticity with Kinematic Hardening
C     routine is written for linear hardening because the classical
C			Prager-Ziegler theory is limited to this case
C
C 		More infos at:
C     https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-hardening.htm

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

C---------------------------------------------------------------------------------------

C Start of USER code

C---------------------------------------------------------------------------------------

C LOCAL ARRAYS

C ----------------------------------------------------------------
C EELAS - ELASTIC STRAINS
C EPLAS - PLASTIC STRAINS
C ALPHA - SHIFT TENSOR
C FLOW - PLASTIC FLOW DIRECTIONS
C OLDS - STRESS AT START OF INCREMENT
C OLDPL - PLASTIC STRAINS AT START OF INCREMENT
C
 			DIMENSION EELAS(6), EPLAS(6), ALPHA(6), FLOW(6), OLDS(6), OLDPL(6)
C
 			PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0, ENUMAX=.4999D0, TOLER=1.0D-6)
C
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC ELASTICITY AND MISES PLASTICITY
C WITH KINEMATIC HARDENING - CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
C PROPS(1) - E = Elasticity Modulus 
C PROPS(2) - NU = Poisson's Ratio 
C PROPS(3) - SYIELD = Yield Stress
C PROPS(4) - HARD = Hardening Modulus 
C PROPS(5) - SULTIMATE = Ultimate Stress

C ----------------------------------------------------------------
C
C ELASTIC PROPERTIES
C
		 EMOD=PROPS(1)
		 ENU=MIN(PROPS(2), ENUMAX)
		 EBULK3=EMOD/(ONE-TWO*ENU)
		 EG2=EMOD/(ONE+ENU)
C Lamé parameter 2 - EG - Shear Modulus		  
		 EG=EG2/TWO
		 EG3=THREE*EG
C Lamé parameter 1 - ELAM   		 
		 ELAM=(EBULK3-EG2)/THREE
C
C ELASTIC STIFFNESS
C
C DDSDDE(NTENS,NTENS) - Jacobian matrix of the constitutive model
C NTENS - Size of the stress or strain component array (NDI + NSHR). 
C NDI - Number of direct stress components at this point. 
C

 			DO K1=1, NDI
 				DO K2=1, NDI
 					DDSDDE(K2, K1)=ELAM
 				END DO
 				DDSDDE(K1, K1)=EG2+ELAM
 			END DO
 
 
 			DO K1=NDI+1, NTENS
 				DDSDDE(K1, K1)=EG
 			END DO
 
C
C RECOVER ELASTIC STRAIN, PLASTIC STRAIN AND SHIFT TENSOR AND ROTATE
C NOTE: USE CODE 1 FOR (TENSOR) STRESS, CODE 2 FOR (ENGINEERING) STRAIN
C 
C Calls Function ROTSIG
C
C EELAS - ELASTIC STRAINS
C EPLAS - PLASTIC STRAINS
C ALPHA - SHIFT TENSOR
C NSHR - Number of engineering shear stress components at this point.
C

		 CALL ROTSIG(STATEV( 1), 			  DROT, EELAS, 2, NDI, NSHR)
		 CALL ROTSIG(STATEV( NTENS+1),  DROT, EPLAS, 2, NDI, NSHR)
		 CALL ROTSIG(STATEV(2*NTENS+1), DROT, ALPHA, 1, NDI, NSHR)
		 
C
C SAVE STRESS AND PLASTIC STRAINS AND
C CALCULATE PREDICTOR STRESS AND ELASTIC STRAIN
C
C OLDS - STRESS AT START OF INCREMENT
C OLDPL - PLASTIC STRAINS AT START OF INCREMENT
C DSTRAN(NTENS) - Array of strain increments. 
C STRESS(NTENS) - This array is passed in as the stress tensor at the beginning of the increment and must be updated
C                 in this routine to be the stress tensor at the end of the increment.


		 DO K1=1, NTENS
			 OLDS(K1)=STRESS(K1)
			 OLDPL(K1)=EPLAS(K1)
			 EELAS(K1)=EELAS(K1)+DSTRAN(K1)
		 	 	DO K2=1, NTENS
		 			STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
		 		END DO
		 END DO
 
C
C CALCULATE EQUIVALENT VON MISES STRESS (Elastic Predictor)
C
 		SMISES=(STRESS(1)-ALPHA(1)-STRESS(2)+ALPHA(2))**2+(STRESS(2)-ALPHA(2)-STRESS(3)+ALPHA(3))**2+(STRESS(3)-ALPHA(3)-STRESS(1)+ALPHA(1))**2
 		
 		DO K1=NDI+1,NTENS
 			SMISES=SMISES+SIX*(STRESS(K1)-ALPHA(K1))**2
 		END DO
 		SMISES=SQRT(SMISES/TWO)
C
C GET YIELD STRESS AND HARDENING MODULUS
C
	 SYIELD=PROPS(3)
	 HARD=PROPS(4)
C
C DETERMINE IF ACTIVELY YIELDING
C Plastic flow occurs if the elastic predictor is larger than the yield stress !
C
 		IF(SMISES.GT.(ONE+TOLER)*SYIELD) THEN
C
C ACTIVELY YIELDING

C SEPARATE THE HYDROSTATIC FROM THE DEVIATORIC STRESS
C CALCULATE THE FLOW DIRECTION
C
C FLOW - PLASTIC FLOW DIRECTIONS

 			SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
 			DO K1=1,NDI
 				FLOW(K1)=(STRESS(K1)-ALPHA(K1)-SHYDRO)/SMISES
 			END DO
 
 			DO K1=NDI+1,NTENS
 				FLOW(K1)=(STRESS(K1)-ALPHA(K1))/SMISES
 			END DO
C
C SOLVE FOR EQUIVALENT PLASTIC STRAIN INCREMENT
C
 			DEQPL=(SMISES-SYIELD)/(EG3+HARD)
 
C
C UPDATE SHIFT TENSOR, ELASTIC AND PLASTIC STRAINS AND STRESS
C
C EELAS - ELASTIC STRAINS
C EPLAS - PLASTIC STRAINS
C ALPHA - SHIFT TENSOR
C FLOW - PLASTIC FLOW DIRECTIONS

			 DO K1=1,NDI
				 ALPHA(K1)=ALPHA(K1)+HARD*FLOW(K1)*DEQPL
				 EPLAS(K1)=EPLAS(K1)+THREE/TWO*FLOW(K1)*DEQPL
				 EELAS(K1)=EELAS(K1)-THREE/TWO*FLOW(K1)*DEQPL
				 STRESS(K1)=ALPHA(K1)+FLOW(K1)*SYIELD+SHYDRO
			 END DO
			
			 DO K1=NDI+1,NTENS
				 ALPHA(K1)=ALPHA(K1)+HARD*FLOW(K1)*DEQPL
				 EPLAS(K1)=EPLAS(K1)+THREE*FLOW(K1)*DEQPL
				 EELAS(K1)=EELAS(K1)-THREE*FLOW(K1)*DEQPL
				 STRESS(K1)=ALPHA(K1)+FLOW(K1)*SYIELD
			 END DO
			 
C			 
C CALCULATE PLASTIC DISSIPATION
C
			 SPD=ZERO
			 DO K1=1,NTENS
			 	SPD=SPD+(STRESS(K1)+OLDS(K1))*(EPLAS(K1)-OLDPL(K1))/TWO
			 END DO
C
C FORMULATE THE JACOBIAN (MATERIAL TANGENT)
C FIRST CALCULATE EFFECTIVE MODULI
C
			 EFFG=EG*(SYIELD+HARD*DEQPL)/SMISES
			 EFFG2=TWO*EFFG
			 EFFG3=THREE*EFFG
			 EFFLAM=(EBULK3-EFFG2)/THREE
			 EFFHRD=EG3*HARD/(EG3+HARD)-EFFG3
			 	DO K1=1, NDI
			 		DO K2=1, NDI
			 			DDSDDE(K2, K1)=EFFLAM
			 		END DO
			 		DDSDDE(K1, K1)=EFFG2+EFFLAM
			 	END DO
			 
			 DO K1=NDI+1, NTENS
			 	DDSDDE(K1, K1)=EFFG
			 END DO
			 DO K1=1, NTENS
			 	DO K2=1, NTENS
			 		DDSDDE(K2, K1)=DDSDDE(K2, K1)+EFFHRD*FLOW(K2)*FLOW(K1)
			 	END DO
			 END DO
			 ENDIF
 
C
C STORE ELASTIC STRAINS, PLASTIC STRAINS AND SHIFT TENSOR
C IN STATE VARIABLE ARRAY


C
			 DO K1=1,NTENS
			 	STATEV(K1)=EELAS(K1)
			 	STATEV(K1+NTENS)=EPLAS(K1)
			 	STATEV(K1+2*NTENS)=ALPHA(K1)
			 END DO
			 
C GET Ultimate Stress - SULTIMATE=PROPS(5)
C
C Define STATEV (19) as the STRESS (2) in Loading Direction Y
C IF STATEV (19) is higher than the Ultimate Stress we have element deletion	
C 
C Define STATEV (20) for Element Deletion 
C - STATEV(20) = 0 - Element is deleted
C - STATEV(20) = 1 - Element is active
C		 
       SULTIMATE=PROPS(5)
			 STATEV(19) = STRESS(2)
			 IF (STATEV(19) .GE. SULTIMATE) THEN
			 	IF (STATEV(20) .NE. 0) THEN
			 		STATEV(20) = 0
			 	END IF
			 END IF
			 
C 
			 RETURN
			 END