module hashMod
! simple integer hash

private
public :: initHash, find, endHash, hash_type

type hash_type
  integer,  allocatable :: nextkey(:)
  integer :: nkey
end type

 contains
 
function hashfunc(nkey,dat)
 integer :: dat,nkey,hashfunc
 hashfunc=ABS(mod(dat,nkey))+1
end function
 
subroutine initHash(nk,n,dat,hash)
 integer :: nk,n,dat(n)
 type(hash_type) :: hash
 hash%nkey=nk
 if(allocated(hash%nextkey)) call endHash(hash) 
 allocate(hash%nextkey(hash%nkey+n))
 call insert(n,dat,hash)
end subroutine

subroutine insert(n,dat,hash)
 integer :: i,n,next,dat(n)
 type(hash_type) :: hash
 hash%nextkey(1:hash%nkey)=0
 do i=1,n
   next=hashfunc(hash%nkey,dat(i))
   do while(hash%nextkey(next).NE.0)
      next=hash%nextkey(next)
   enddo
   hash%nextkey(next)=hash%nkey+i
   hash%nextkey(hash%nkey+i)=0
 enddo
end subroutine

function find(item,dat,hash) result(next)
 integer :: item,dat(*)
 integer :: next
 type(hash_type) :: hash
 next=hashfunc(hash%nkey,item)
 next=hash%nextkey(next)
 do
   if(next.EQ.0) exit
   if(dat(next-hash%nkey).EQ.item) exit
   next=hash%nextkey(next)
 enddo
 if(next.NE.0) next=next-hash%nkey  
end function

subroutine endHash(hash)
  type(hash_type) :: hash
  deallocate(hash%nextkey)
end subroutine

end module


subroutine test
 use hashMod
 type(hash_type) :: hash
 integer :: n,i,dat(13)
 
 n=13
 dat=(/12,2,3,4,5,5,6,12,13,16,-6,0,-8383/)
 
 call inithash(5,13,dat(1:n),hash)
 
10 read*,i
  if(i.eq.-999) stop
  i=find(i,dat(1:n),hash)
  if(i.NE.0) print*,i,dat(i)  
  goto 10
 
 call endHash(hash)

end subroutine
