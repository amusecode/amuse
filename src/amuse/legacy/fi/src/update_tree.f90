! code to reduce an existing tree structure parallel
! counts the number of cells
recursive subroutine update_reduction(start,tthread,option,ind)
  include 'globals.h'
 integer, intent(in) :: start,tthread
 character (len=4), intent(in) :: option
 integer,intent(out) :: ind(*)
 integer :: facell, lacell, nacell, k, i

 ind(1) = start
 facell = 1;lacell = 1;nacell = 1
 do while(lacell-facell+1.GT.0)  
       do k = 1, nsubcell 
        do i = facell,lacell 
         if (subp(ind(i),k).GE.nbods1) then 
          nacell = nacell + 1 
          ind(nacell) = subp(ind(i),k) 
         endif 
        enddo 
       enddo
       facell = lacell + 1 
       lacell = nacell 
       if(tthread.GT.1.AND.lacell-facell+1.GT.4*tthread) then
!$omp parallel do private(i,k)       
        do i=facell,lacell
         k=ind(i)
         call update_reduction(k,1,option,srlist)
        enddo
	exit
       endif		
 enddo
 call tree_reduction(start,facell-1,option)
 end subroutine
