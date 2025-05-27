      subroutine sort3(n,ra,rb,rc)
*
*
*       heapsort method (press p. 231).
*       -------------------------------
*
      real*8 ra,rra
*
      integer rb,rrb,rc,rrc,l,n,ir,i,j
*
      dimension  ra(n),rb(n),rc(n)
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
          rrc=rc(l)
      else
          rra=ra(ir)
	  rrb=rb(ir)
          rrc=rc(ir)
	  ra(ir)=ra(1)
	  rb(ir)=rb(1)
          rc(ir)=rc(1)
          ir=ir-1
*
	  if(ir.eq.1)then
	      ra(1)=rra
	      rb(1)=rrb
              rc(1)=rrc
	      return
          endif
*
      end if
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
               rc(i)=rc(j)
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
      rc(i)=rrc
*
      go to 10
*
      end
*
*
*
*
