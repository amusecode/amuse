        subroutine  readmm
	
	real r,drhalo,dr,fixedpot,radhalo(10000),rhohalo(10000),
     &	 xmhalo(10000),upothalo(10000),f2halo(10000),fhalo(10000),masshalo
        integer i,in,in1,nhalo
	
	character*16 halofile
	logical firstcall
	
	common/fixed/nhalo,masshalo,drhalo,radhalo,rhohalo,xmhalo,
     &   upothalo,fhalo,f2halo
	
		 halofile='MODEL.MASS'
	 OPEN(19,FILE=halofile,
     &	    STATUS='OLD',iostat=ioerror)
	 IF(ioerror.NE.0) THEN
	  print*,' stop -- error reading halofile:',halofile
	  stop
	 ENDIF 
	 READ(19,*) nhalo
	 
	 IF(nhalo.GT.10000) THEN
	  print*,' nhalo overflow'
	  stop
	 ENDIF

         DO i=1,nhalo
          READ(19,*) radhalo(i),rhohalo(i),xmhalo(i),upothalo(i),fhalo(i),
     &	  f2halo(i)
	 ENDDO
	 
         CLOSE(19)
	 
	 masshalo=xmhalo(nhalo)
         drhalo=radhalo(5)-radhalo(4)
		
	end
	
	function fixedpot(r)
	
	real r,drhalo,dr,fixedpot,radhalo(10000),rhohalo(10000),
     &	 xmhalo(10000),upothalo(10000),f2halo(10000),fhalo(10000),masshalo
        integer i,in,in1,nhalo
	
	character*16 halofile
	logical firstcall
	
	common/fixed/nhalo,masshalo,drhalo,radhalo,rhohalo,xmhalo,
     &   upothalo,fhalo,f2halo
	common/firstt/ firstcall
	data firstcall/.TRUE./
	
	
	if(firstcall) then
	 firstcall = .FALSE.
         call readmm	
	endif


       if(r.GE.radhalo(1)) then
        if(r.LE.radhalo(nhalo)) then
         in=(r-radhalo(1))/drhalo
         in1=in+1
         dr=(r-radhalo(1))/drhalo-in
         fixedpot=(1-dr)*upothalo(in1)+dr*upothalo(in1+1)
        else
         fixedpot=upothalo(nhalo)
        endif
       else
        fixedpot=upothalo(1)
       endif
       
       end

        function halofixeddens(r,z) result(x)
         real r,z,x
         x=fixeddens(sqrt(r*r+z*z))
        end function

        function fixeddens(r)

	real r,drhalo,dr,fixeddens,radhalo(10000),rhohalo(10000),
     &	 xmhalo(10000),upothalo(10000),f2halo(10000),fhalo(10000),masshalo
        integer i,in,in1,nhalo
	
	character*16 halofile
	logical firstcall
	
	common/fixed/nhalo,masshalo,drhalo,radhalo,rhohalo,xmhalo,
     &   upothalo,fhalo,f2halo
	common/firstt/ firstcall
			
	if(firstcall) then
	 firstcall = .FALSE.
         call readmm		
	endif


       if(r.GE.radhalo(1)) then
        if(r.LE.radhalo(nhalo)) then
         in=(r-radhalo(1))/drhalo
         in1=in+1
         dr=(r-radhalo(1))/drhalo-in
         fixeddens=(1-dr)*rhohalo(in1)+dr*rhohalo(in1+1)
        else
         fixeddens=rhohalo(nhalo)
        endif
       else
        fixeddens=rhohalo(1)
       endif
       
       end

        function fixedforce(r)

	real r,drhalo,dr,fixedforce,radhalo(10000),rhohalo(10000),
     &	 xmhalo(10000),upothalo(10000),f2halo(10000),fhalo(10000),masshalo
        integer i,in,in1,nhalo
	
	character*16 halofile
	logical firstcall
	
	common/fixed/nhalo,masshalo,drhalo,radhalo,rhohalo,xmhalo,
     &   upothalo,fhalo,f2halo
	common/firstt/ firstcall
		
	if(firstcall) then
	 firstcall = .FALSE.
         call readmm	
	endif

       if(r.GE.radhalo(1)) then
        if(r.LE.radhalo(nhalo)) then
         in=(r-radhalo(1))/drhalo
         in1=in+1
         dr=(r-radhalo(1))/drhalo-in
         fixedforce=(1-dr)*fhalo(in1)+dr*fhalo(in1+1)
        else
         fixedforce=fhalo(nhalo)
        endif
       else
        fixedforce=fhalo(1)
       endif
       
       end

        function fixedf2(r)

	real r,drhalo,dr,fixedf2,radhalo(10000),rhohalo(10000),
     &	 xmhalo(10000),upothalo(10000),f2halo(10000),fhalo(10000),masshalo
        integer i,in,in1,nhalo
	
	character*16 halofile
	logical firstcall
	
	common/fixed/nhalo,masshalo,drhalo,radhalo,rhohalo,xmhalo,
     &   upothalo,fhalo,f2halo
	common/firstt/ firstcall
			
	if(firstcall) then
	 firstcall = .FALSE.
         call readmm		
	endif


       if(r.GE.radhalo(1)) then
        if(r.LE.radhalo(nhalo)) then
         in=(r-radhalo(1))/drhalo
         in1=in+1
         dr=(r-radhalo(1))/drhalo-in
         fixedf2=(1-dr)*f2halo(in1)+dr*f2halo(in1+1)
        else
         fixedf2=f2halo(nhalo)
        endif
       else
        fixedf2=f2halo(1)
       endif
       
       end


       
!       program test
!       real r,p1,p2,p3,p4
!       integer i
!10     read*,r
!       if(r.EQ.0) stop
!       p2=fixeddens(r)
!       p1=fixedpot(r)
!       p3=fixedforce(r)
!       p4=fixedf2(r)
!       print*,r,p1,p2,p3,p4
!       goto 10
!        do i=0,1000
!	 r=.01*i
!	 p1=fixedf2(r)
!	 print*,r,p1
!	enddo 
!       end
              
