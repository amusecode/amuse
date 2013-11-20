!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list(ind_part,list2,ok,np)
  use amr_commons
  use pm_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part,list2
  logical,dimension(1:nvector)::ok
  !
  ! Add particles to their new linked lists
  !
  integer::j

  do j=1,np
     if(ok(j))then
        if(numbp(list2(j))>0)then
           ! Add particle at the tail of its linked list
           nextp(tailp(list2(j)))=ind_part(j)
           prevp(ind_part(j))=tailp(list2(j))
           nextp(ind_part(j))=0
           tailp(list2(j))=ind_part(j)
           numbp(list2(j))=numbp(list2(j))+1
        else
           ! Initialise linked list
           headp(list2(j))=ind_part(j)
           tailp(list2(j))=ind_part(j)
           prevp(ind_part(j))=0
           nextp(ind_part(j))=0
           numbp(list2(j))=1
        end if
     end if
  end do

end subroutine add_list
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free(ind_part,np)
  use amr_commons
  use pm_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  !
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=0.0
        vp(ind_part(j),idim)=0.0
     end do
  end do
  do j=1,np
     mp(ind_part(j))=0.0
     idp(ind_part(j))=0
     levelp(ind_part(j))=0
  end do
  if(star.or.sink)then
     do j=1,np
        tp(ind_part(j))=0.0
     end do
     if(metal)then
        do j=1,np
           zp(ind_part(j))=0.0
        end do
     end if
  end if

  do j=1,np
     if(numbp_free>0)then
        ! Add particle at the tail of its linked list
        nextp(tailp_free)=ind_part(j)
        prevp(ind_part(j))=tailp_free
        nextp(ind_part(j))=0
        tailp_free=ind_part(j)
        numbp_free=numbp_free+1
     else
        ! Initialise linked list
        headp_free=ind_part(j)
        tailp_free=ind_part(j)
        prevp(ind_part(j))=0
        nextp(ind_part(j))=0
        numbp_free=1
     end if
  end do
  npart=npartmax-numbp_free

end subroutine add_free
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free_cond(ind_part,ok,np)
  use amr_commons
  use pm_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  logical,dimension(1:nvector)::ok
  !
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xp(ind_part(j),idim)=0.0
           vp(ind_part(j),idim)=0.0
        endif
     end do
  end do
  do j=1,np
     if(ok(j))then
        mp(ind_part(j))=0.0
        idp(ind_part(j))=0
        levelp(ind_part(j))=0
     endif
  end do
  if(star.or.sink)then
     do j=1,np
        if(ok(j))then
           tp(ind_part(j))=0.0
        endif
     end do
     if(metal)then
        do j=1,np
           if(ok(j))then
              zp(ind_part(j))=0.0
           endif
        end do
     end if
  end if

  do j=1,np
     if(ok(j))then
        if(numbp_free>0)then
           ! Add particle at the tail of its linked list
           nextp(tailp_free)=ind_part(j)
           prevp(ind_part(j))=tailp_free
           nextp(ind_part(j))=0
           tailp_free=ind_part(j)
           numbp_free=numbp_free+1
        else
           ! Initialise linked list
           headp_free=ind_part(j)
           tailp_free=ind_part(j)
           prevp(ind_part(j))=0
           nextp(ind_part(j))=0
           numbp_free=1
        end if
     endif
  end do
  npart=npartmax-numbp_free

end subroutine add_free_cond
