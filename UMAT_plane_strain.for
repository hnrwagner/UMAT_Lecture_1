C     ABAQUS Subroutine for linear elastic material
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
C           G = Shear Modulus 

       E=PROPS(1)
       v=PROPS(2)
       G=E/(2.D0*(1.D0+v))

C           DDSDDE(NTENS,NTENS)
C           Jacobian matrix of the constitutive model

C           NTENS
C           Size of the stress or strain component array (NDI + NSHR).

C           NDI
C           Number of direct stress components at this point.

C           NSHR
C           Number of engineering shear stress components at this point.

C          Stiffness Matrix for Plane Strain
		
       DO i=1, NDI
           DO j=1, NDI
               DDSDDE(j, i)=(E*v)/((ONE+v)*(ONE-TWO*v))
           END DO
           DDSDDE(i, i)=(E*(ONE-v))/((ONE+v)*(ONE-TWO*v))
       END DO
       
       DO i=NDI+1, NTENS
           DDSDDE(i ,i)=G
       END DO

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
           STRESS(j)=STRESS(j)+DDSDDE(j, i)*DSTRAN(i)
       END DO
       END DO

C---------------------------------------------------------------------------------------

C End of USER code

C---------------------------------------------------------------------------------------


      RETURN
      END