C     ABAQUS Subroutine for Isotropic material - Plane Stress
C     by Dr.Ing. Ronald Wagner Siemens Mobility GmbH Braunschweig / Germany
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

        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

C---------------------------------------------------------------------------------------

C Start of USER code

C---------------------------------------------------------------------------------------
C           PROPS(NPROPS)
C           User-specified array of material constants associated with this user material.
C
C           NPROPS
C           User-defined number of material constants associated with this user material.
C
C           PROPS(1) = Elasticity Modulus
C           PROPS(2) = Poisson's Ratio     
     
           E=PROPS(1)
           v=PROPS(2)
           
            
    			DO i = 1, NTENS
    			 DO j = 1, NTENS
    			 DDSDDE(i,j)= ZERO
    			 END DO
    			END DO

C           DDSDDE(NTENS,NTENS)
C           Jacobian matrix of the constitutive model

C           NTENS
C           Size of the stress or strain component array (NDI + NSHR).

C           NDI
C           Number of direct stress components at this point.

C           NSHR
C           Number of engineering shear stress components at this point.
				   
C          Stiffness Matrix for Plane Stress
		       			
	     	   DDSDDE(1,1) = E/(ONE-v*v)
    		   DDSDDE(2,2) = E/(ONE-v*v)
    	       DDSDDE(3,3) = E*(ONE-v)/(TWO*(ONE-v*v))
    		   
    		   DDSDDE(1,2) = E*v/(ONE-v*v)
    		   DDSDDE(2,1) = E*v/(ONE-v*v)


C           Calculation of STRESSES 
C
C           STRESS(NTENS)
C           This array is passed in as the stress tensor at the beginning of the increment and must be updated
C           in this routine to be the stress tensor at the end of the increment.
C      
C           DSTRAN(NTENS)
C           Array of strain increments.	           			

               DO i=1, NTENS
                   DO j=1, NTENS
                   STRESS(i)=STRESS(i)+DDSDDE(i, j)*DSTRAN(j)
                   END DO
               END DO

C---------------------------------------------------------------------------------------

C End of USER code

C---------------------------------------------------------------------------------------

      RETURN
      END
