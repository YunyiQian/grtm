	subroutine bessj012d(x,j0,j1,j2,j0d,j1d,j2d)

!c	--------------------------------------------	
!c	called subroutines: bessj0,bessj1,bessj	
	implicit none
	real*8      j0,j1,j2,j3,j0d,j1d,j2d,x  
	real*8   bessj0,bessj1,bessj
	
	j0  = bessj0(x)
	j1  = bessj1(x)
	j2  = bessj(2,x)
      j3  = bessj(3,x)
!c	if(j3.ne.0.0) then
!c	print*,'j3 is not equal to 0,and j3=',j3
!c	pause
!c	end if
	j0d = -j1
	j1d = (j0-j2)/2.0
      j2d = (j1-j3)/2.0

	return
	end 


	subroutine bessj01d(x,j0,j1,j0d,j1d)
	implicit none

!c	--------------------------------------------	
!c	called subroutines: bessj0,bessj1,bessj	
	real*8      j0,j1,j2,j0d,j1d ,x
	real*8  bessj0,bessj1,bessj
	
	j0  = bessj0(x)
	j1  = bessj1(x)
	j2  = bessj(2,x)
	j0d =-j1
	j1d = (j0-j2)/2.0

	return
	end 

      subroutine bessj0d(x,j0,j0d)
!c       --------------------------------------------
!c       called subroutines: bessj0,bessj1
	implicit none
      real*8      j0,j0d ,x
      real*8  bessj0,bessj1
 
      j0  = bessj0(x)
      j0d =-bessj1(x)

      return
      end


!c-----------------------------------------------------------
!c	BESSEL FUNCTION OF FIRST KIND OF ORDER N(N.GT.1)   C
!c	------------------------------------------------   C

	real*8 FUNCTION BESSJ(N,X)
	implicit none
	real*8 x,bigno,bigni,ax,tox, bjm,bj,bjp,jsum,sum 
	integer n,iacc,j,m
	real*8  bessj0,bessj1
	PARAMETER (IACC=40,BIGNO=1.E10,BIGNI=1.E-10)

	IF(N.EQ.0)THEN
	BESSJ=BESSJ0(X)
	RETURN
	ENDIF

	IF(N.EQ.1)THEN
	BESSJ=BESSJ1(X)
	RETURN
	ENDIF
	
	AX=dabs(X)
	IF (AX.EQ.0.0) THEN
	   BESSJ=0.0
	ELSE IF(AX.GT.FLOAT(N)) THEN	
!c #Use upwards recurrence from Jo(x) and J1(x)
	   TOX=2.0/AX		        
           BJM=BESSJ0(AX)
           BJ=BESSJ1(AX)
	   DO J=1,N-1
	      BJP=J*TOX*BJ-BJM
	      BJM=BJ
	      BJ=BJP
	   ENDDO
	   BESSJ=BJ
	ELSE				
!c #Use downwards recurrence from an even value M here computed.
	   TOX=2.0/AX	     	        
	   M=2*((N+INT(sqrt(FLOAT(IACC*N))))/2)
!c  Make IACC large to increase accuracy.
	   BESSJ=0.0
           JSUM=0
	   SUM=0.
	   BJP=0.
	   BJ=1.
	   DO J=M,1,-1		
!c # downwards recurrence.
	      BJM=J*TOX*BJ-BJP
	      BJP=BJ
	      BJ=BJM
	      IF (dabs(BJ).GT.BIGNO) THEN  
!c # renormalize to prevent overflows.
		 BJ=BJ*BIGNI
		 BJP=BJP*BIGNI
		 BESSJ=BESSJ*BIGNI
		 SUM=SUM*BIGNI
	      ENDIF
	      IF (JSUM.NE.0) SUM=SUM+BJ     
!c # accumulate the sum.
	      JSUM=1-JSUM		  
!c change 0 to 1 or vice-versa.
	      IF (J.EQ.N) BESSJ=BJP	  
!c Save the unnormalized answer.
	   ENDDO	
	   SUM=2.0*SUM-BJ		 
!c compute (5.4.6)
	   BESSJ=BESSJ/SUM		  
!c and use it to normalize the answer.
	ENDIF
	IF (X.LT.0.0.and.MOD(N,2).EQ.1) BESSJ=-BESSJ
	RETURN
	END

!c-------------------------------------------------------------------
!c 	SUBROUTINES OF BESSJ0 and BESSJ1
!c	-------------------------------------

	real*8 FUNCTION BESSJ0(X)
	implicit none
	real*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     &	       S1,S2,S3,S4,S5,S6,x,ax,z,xx
	DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     &    -.2073370639D-5,.2093887211D-6/,Q1,Q2,Q3,Q4,Q5/-.1562499995D-1
     &   ,.1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
     	DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,
     &    651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0/,
     &    S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,9494680.718D0,
     &    59272.64853D0,267.8532712D0,1.D0/

	IF (dabs(X).LT.8.0) THEN	
!c direct rational function fit.
	   Y=X**2
	   BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/
     &            (S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))) 
	ELSE
	   AX=dabs(X)
	   Z=8./AX
	   Y=Z**2
	   XX=AX-.785398164
	 BESSJ0=dsqrt(.636619772/AX)*(dcos(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)
     &            )))-Z*dsin(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
	ENDIF
	RETURN
	END

!c----------------------------------------------------------------------

	real*8 FUNCTION BESSJ1(X)
	implicit none
	real*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     &	       S1,S2,S3,S4,S5,S6,x,ax,z,xx
	DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,
     &	  .2457520174D-5,-.240337019D-6/,Q1,Q2,Q3,Q4,Q5/.04687499995D0
     &   ,-.2002690873D-3,.8449199096D-5,-.88228987D-6,.105787412D-6/
     	DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,
     &    242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0/,
     &	 S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,18583304.74D0
     &	  ,99447.43394D0,376.9991397D0,1.D0/

	IF (dabs(X).LT.8.0) THEN	
!c direct rational function fit.
	   Y=X**2
	   BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))/
     &            (S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))) 
	ELSE
	   AX=dabs(X)
	   Z=8./AX
	   Y=Z**2
	   XX=AX-2.356194491
	 BESSJ1=dsqrt(.636619772/AX)*(dcos(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)
     &        )))-Z*dsin(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))*SIGN(P1,X)
	ENDIF
	RETURN
	END













