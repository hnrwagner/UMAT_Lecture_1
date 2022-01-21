C     ABAQUS Subroutine for Neo-Hookean Hyperelasticity
C
C     More Infos can be found at 
C     
C			http://130.149.89.49:2080/v6.8/books/stm/default.htm?startat=ch04s06ath123.html

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
C BBAR - DEVIATORIC CAUCHY-GREEN TENSOR
C BBARP - PRINCIPAL VALUES OF BBAR
C BBARN - PRINCIPAL DIRECTION OF BBAR 
C DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C ----------------------------------------------------------------
C
 					DIMENSION BBAR(6), BBARP(3), BBARN(3, 3), DISTGR(3,3)
C
 					PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0, SIX=6.D0)
C
C ----------------------------------------------------------------
C UMAT FOR COMPRESSIBLE NEO-HOOKEAN HYPERELASTICITY
C CANNOT BE USED FOR PLANE STRESS, 2D SHELL ELEMENTS AND ELEMENTS WITH REDUCED INTEGRATION
C Tested for C3D8 & C3D20
C ----------------------------------------------------------------
C PROPS(1) - C10 Coefficiant  = G/2 ( G = shear modulus)
C PROPS(2) - D1 Coefficiant  = 2/K ( K = bulk modulus)

C ----------------------------------------------------------------
C
C INPUT
C

 			C10=PROPS(1)
 			D1=PROPS(2)
C
C JACOBIAN AND DISTORTION TENSOR
C
 			DET=DFGRD1(1, 1)*DFGRD1(2, 2)*DFGRD1(3, 3)-DFGRD1(1, 2)*DFGRD1(2, 1)*DFGRD1(3, 3)
C	Dimension Check - 3D or 2D elements ? 			
 			IF(NSHR.EQ.3) THEN
 				DET=DET+DFGRD1(1, 2)*DFGRD1(2, 3)*DFGRD1(3, 1)+DFGRD1(1, 3)*DFGRD1(3, 2)*DFGRD1(2, 1)-DFGRD1(1, 3)*DFGRD1(3,1)*DFGRD1(2, 2)-DFGRD1(2, 3)*DFGRD1(3, 2)*DFGRD1(1, 1)
 			END IF
 			SCALE=DET**(-ONE/THREE)
 
 			DO K1=1, 3
 				DO K2=1, 3
 					DISTGR(K2, K1)=SCALE*DFGRD1(K2, K1)
 				END DO
 			END DO

C 			
C CALCULATE DEVIATORIC LEFT CAUCHY-GREEN DEFORMATION TENSOR
C
 
 			BBAR(1)=DISTGR(1, 1)**2+DISTGR(1, 2)**2+DISTGR(1, 3)**2
 			BBAR(2)=DISTGR(2, 1)**2+DISTGR(2, 2)**2+DISTGR(2, 3)**2
 			BBAR(3)=DISTGR(3, 3)**2+DISTGR(3, 1)**2+DISTGR(3, 2)**2
 			BBAR(4)=DISTGR(1, 1)*DISTGR(2, 1)+DISTGR(1, 2)*DISTGR(2, 2)+DISTGR(1, 3)*DISTGR(2, 3)
C	Dimension Check - 3D or 2D elements ? 	
 			IF(NSHR.EQ.3) THEN 
 			BBAR(5)=DISTGR(1, 1)*DISTGR(3, 1)+DISTGR(1, 2)*DISTGR(3, 2)+DISTGR(1, 3)*DISTGR(3, 3)
 			BBAR(6)=DISTGR(2, 1)*DISTGR(3, 1)+DISTGR(2, 2)*DISTGR(3, 2)+DISTGR(2, 3)*DISTGR(3, 3)
 
 			END IF
 			
C
C CALCULATE THE STRESS
C

 			TRBBAR=(BBAR(1)+BBAR(2)+BBAR(3))/THREE
 			EG=TWO*C10/DET
 			EK=TWO/D1*(TWO*DET-ONE)
 			PR=TWO/D1*(DET-ONE)
 			DO K1=1,NDI
 				STRESS(K1)=EG*(BBAR(K1)-TRBBAR)+PR
 			END DO
 
  		DO K1=NDI+1,NDI+NSHR
 				STRESS(K1)=EG*BBAR(K1)
 			END DO
 			
C 			
C CALCULATE THE STIFFNESS
C
 			EG23=EG*TWO/THREE
 			DDSDDE(1, 1)= EG23*(BBAR(1)+TRBBAR)+EK
 			DDSDDE(2, 2)= EG23*(BBAR(2)+TRBBAR)+EK
 			DDSDDE(3, 3)= EG23*(BBAR(3)+TRBBAR)+EK
 			DDSDDE(1, 2)=-EG23*(BBAR(1)+BBAR(2)-TRBBAR)+EK
 			DDSDDE(1, 3)=-EG23*(BBAR(1)+BBAR(3)-TRBBAR)+EK
 			DDSDDE(2, 3)=-EG23*(BBAR(2)+BBAR(3)-TRBBAR)+EK
 			DDSDDE(1, 4)= EG23*BBAR(4)/TWO
 			DDSDDE(2, 4)= EG23*BBAR(4)/TWO
 			DDSDDE(3, 4)=-EG23*BBAR(4)
 			DDSDDE(4, 4)= EG*(BBAR(1)+BBAR(2))/TWO
C	Dimension Check - 3D or 2D elements ? 	
 			IF(NSHR.EQ.3) THEN
	 			DDSDDE(1, 5)= EG23*BBAR(5)/TWO
 				DDSDDE(2, 5)=-EG23*BBAR(5)
 				DDSDDE(3, 5)= EG23*BBAR(5)/TWO
 				DDSDDE(1, 6)=-EG23*BBAR(6)
 				DDSDDE(2, 6)= EG23*BBAR(6)/TWO
 				DDSDDE(3, 6)= EG23*BBAR(6)/TWO
 				DDSDDE(5, 5)= EG*(BBAR(1)+BBAR(3))/TWO
 
 				DDSDDE(6, 6)= EG*(BBAR(2)+BBAR(3))/TWO
 				DDSDDE(4,5)= EG*BBAR(6)/TWO
 				DDSDDE(4,6)= EG*BBAR(5)/TWO
 				DDSDDE(5,6)= EG*BBAR(4)/TWO
 			END IF
 			
 			DO K1=1, NTENS
 				DO K2=1, K1-1
 					DDSDDE(K1, K2)=DDSDDE(K2, K1)
 				END DO
 			END DO

C---------------------------------------------------------------------------------------

C End of USER code

C---------------------------------------------------------------------------------------

 		RETURN
 		END