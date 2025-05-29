	real*8 function stability(mm1,mm2,mm3,ein,eout,inc)

*       Three-body stability test (Mardling & Aarseth 1999)
*       
        implicit real*8 (a-h,m,o-z)
        REAL*8 inc
	
*       Employ the semi-analytical stability criterion (MA 1999).
	Q = mm3/(mm1 + mm2)
	IF (EOUT.LT.1.0) THEN
	   XFAC = (1.0 + Q)*(1.0 + EOUT)/SQRT(1.0 - EOUT)
	ELSE
	   XFAC = 40.0*(1.0 + Q)
	END IF
	PCRIT = 2.8*XFAC**0.4
*       
*       Include the inclination fudge factor.
	YFAC = 1.0 - 0.3*INC/3.1415
	PCRIT = YFAC*PCRIT
*       
        stability=PCRIT

	end
