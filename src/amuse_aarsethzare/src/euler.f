	subroutine euler(angles,v)

*	Rotates vector given wrt (i,j,k) through Euler angles (i,w,Omega).

	implicit real*8 (a-h,o-z)
	real*8 v(3),angles(3),B(3,3),temp(3)

	cosi=cos(angles(1))
	sini=sin(angles(1))
	cosw=cos(angles(2))
	sinw=sin(angles(2))
	cosOm=cos(angles(3))
	sinOm=sin(angles(3))

        B(1,1) = cosw*cosOm - sinw*cosi*sinOm
        B(1,2) = cosw*sinOm + sinw*cosi*cosOm
        B(1,3) = sinw*sini
        B(2,1) = -sinw*cosOm - cosw*cosi*sinOm
        B(2,2) = -sinw*sinOm + cosw*cosi*cosOm
        B(2,3) = cosw*sini
        B(3,1) = sini*sinOm
        B(3,2) = -sini*cosOm
        B(3,3) = cosi

	do i=1,3
	   sum=0
	   sum2=0
	   do j=1,3
	      sum=sum + B(j,i)*v(j)
	   enddo
	   temp(i)=sum
	enddo
	do i=1,3
	   v(i)=temp(i)
	enddo

	end

