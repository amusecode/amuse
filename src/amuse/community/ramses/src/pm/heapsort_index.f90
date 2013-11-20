subroutine swap(ind,size,x,y)
  integer::x,y,z,size
  integer,dimension(1:size)::ind
  z = ind(x)
  ind(x) = ind(y)
  ind(y) = z
end subroutine swap
 
subroutine sift_down(a,ind,end,i,size) 
  integer::i,size,ln,j,end
  real(kind=8),dimension(1:size)::a
  integer,dimension(1:size)::ind
 
  j = i
  do while (j .le. (end/2))
     ln = 2*j
     if (ln .lt. end) then
        if ((a(ind(ln+1))) .gt. (a(ind(ln)))) then 
           ln = ln+1
        end if
     end if
     if (a(ind(ln)) .gt. a(ind(j))) then
        call swap(ind,size,j,ln)
        j = ln
     else
        j = end
     end if
  end do
 
  return
end subroutine sift_down
 
subroutine heapify(a,ind,size)
  integer::size,i
  real(kind=8),dimension(1:size)::a
  integer,dimension(1:size)::ind
 
  do i = size/2,1,-1
     call sift_down(a,ind,size,i,size)
  end do

end subroutine heapify

subroutine heapsort_index(a,ind,size) 
  integer::size,i
  real(kind=8),dimension(1:size)::a
  integer,dimension(1:size)::ind

  do i=1,size 
     ind(i)=i
  end do

  call heapify(a,ind,size)

  do i=size,1,-1
     call swap(ind,size,1,i)
     call sift_down(a,ind,i-1,1,size)
  end do


end subroutine heapsort_index
