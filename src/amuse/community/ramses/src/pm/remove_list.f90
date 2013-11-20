!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_list(ind_part,list1,ok,np)
  use amr_commons
  use pm_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part,list1
  logical,dimension(1:nvector)::ok
  !----------------------------------------------------
  ! Remove particles from their original linked lists
  !----------------------------------------------------
  integer::j
  do j=1,np
     if(ok(j))then
        if(prevp(ind_part(j)) .ne. 0) then
           if( nextp(ind_part(j)) .ne. 0 )then
              nextp(prevp(ind_part(j)))=nextp(ind_part(j))
              prevp(nextp(ind_part(j)))=prevp(ind_part(j))
           else
              nextp(prevp(ind_part(j)))=0
              tailp(list1(j))=prevp(ind_part(j))
           end if
        else
           if(nextp(ind_part(j)) .ne. 0)then
              prevp(nextp(ind_part(j)))=0
              headp(list1(j))=nextp(ind_part(j))
           else
              headp(list1(j))=0
              tailp(list1(j))=0
           end if
        end if
        numbp(list1(j))=numbp(list1(j))-1
     end if
  end do
end subroutine remove_list
!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_free(ind_part,np)
  use amr_commons
  use pm_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  !-----------------------------------------------
  ! Get np particle from free memory linked list
  !-----------------------------------------------
  integer::j,ipart
  do j=1,np
     ipart=headp_free
     ind_part(j)=ipart
     numbp_free=numbp_free-1
     if(numbp_free<0)then
        write(*,*)'No more free memory'
        write(*,*)'in PE ',myid
        write(*,*)'Increase npartmax'
        call clean_stop
     end if
     headp_free=nextp(headp_free)
  end do
  npart=npartmax-numbp_free
end subroutine remove_free
