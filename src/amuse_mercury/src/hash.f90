module hashMod
! simple integer hash

private
public :: initHash, find, endHash

 integer, save, allocatable :: nextkey(:)
 integer, save :: nkey

 contains
 
function hashfunc(dat)
 integer :: dat,hashfunc
 hashfunc=ABS(mod(dat,nkey))+1
end function
 
subroutine initHash(nk,n,dat)
 integer :: nk,n,dat(n)
 nkey=nk
 if(allocated(nextkey)) call endHash() 
 allocate(nextkey(nkey+n))
 call insert(n,dat)
end subroutine

subroutine insert(n,dat)
 integer :: i,n,next,dat(n)
 nextkey(1:nkey)=0
 do i=1,n
   next=hashfunc(dat(i))
   do while(nextkey(next).NE.0)
      next=nextkey(next)
   enddo
   nextkey(next)=nkey+i
   nextkey(nkey+i)=0
 enddo
end subroutine

function find(item,dat) result(next)
 integer :: item,dat(*)
 integer :: next
 next=hashfunc(item)
 next=nextkey(next)
 do
   if(next.EQ.0) exit
   if(dat(next-nkey).EQ.item) exit
   next=nextkey(next)
 enddo
 if(next.NE.0) next=next-nkey  
end function

subroutine endHash()
  deallocate(nextkey)
end subroutine

end module


subroutine test
 use hashMod
 integer :: n,i,dat(13)
 
 n=13
 dat=(/12,2,3,4,5,5,6,12,13,16,-6,0,-8383/)
 
 call inithash(10,6,dat(8:13))
 
10 read*,i
  if(i.eq.-999) stop
  i=find(i,dat(8:13))
  if(i.NE.0) print*,i,dat(i+7)  
  goto 10
 
 call endHash

end subroutine
