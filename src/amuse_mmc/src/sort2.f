      subroutine sort2(n,ra,rb)
*
*
*       heapsort method (press p. 231).
*       -------------------------------
*
      real*8 ra,rra
*
      integer rb,rrb,l,n,ir,i,j
*
      dimension  ra(n),rb(n)
*
*
      l = n/2+1
      ir=n
   10 continue
*
      if(l.gt.1)then
	  l=l-1
	  rra=ra(l)
	  rrb=rb(l)
      else
          rra=ra(ir)
	  rrb=rb(ir)
	  ra(ir)=ra(1)
	  rb(ir)=rb(1)
          ir=ir-1
*
	  if(ir.eq.1)then
	      ra(1)=rra
	      rb(1)=rrb
	      return
          endif
*
      endif
*
      i=l
      j=l+l
*
   20 if(j.le.ir)then
*
	  if(j.lt.ir)then
	      if(ra(j).lt.ra(j+1))j=j+1
          endif
*
	  if(rra.lt.ra(j))then
	       ra(i)=ra(j)
	       rb(i)=rb(j)
	       i=j
	       j=j+j
           else
	       j=ir+1
           endif
*
           go to 20
      endif
*
      ra(i)=rra
      rb(i)=rrb
*
      go to 10
*
      end
*
*
*
*
