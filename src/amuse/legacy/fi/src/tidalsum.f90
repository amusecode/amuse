subroutine evaluate_tidalfield(peps,ppos,ptide)
 include 'globals.h'
 real,intent(in) :: ppos(3),peps
 real,intent(inout) :: ptide(6)
 
 ptide(1:6)=0
 if(selfgrav) call system_tidalfield(peps,ppos,ptide)
 if(fixthalo) call halotidalfield(ppos,ptide)

end subroutine

subroutine system_tidalfield(peps,ppos,ptide)
 include 'globals.h'
 real,intent(in) :: ppos(3),peps
 real,intent(inout) :: ptide(6)
 integer i,nterms
! periodic - 
  if(.NOT.directsum) then
   nterms=0
   call treewalk(root,ppos,peps,0.,0.,nterms)
  else
   nterms=nbodies
   do i=1,nbodies
    bodlist(i)=i
   enddo
  endif 
  call tidalsum(peps,ppos,ptide,nterms)
end subroutine

subroutine tidalsum(peps,ppos,ptide,nterms)
 include 'globals.h'
 integer, intent(in):: nterms
 real,intent(in) :: peps,ppos(3)
 real,intent(inout):: ptide(6)
 real :: dx,dy,dz,dist2,dist,epseff,epseffi,drdeldrg, &
  &  rinveff,r3inveff,r5inveff,drsm,rhinv
 integer :: i,q,smindex
 do i=1,nterms
   q=bodlist(i)
   dx=ppos(1)-pos(q,1)
   dy=ppos(2)-pos(q,2)
   dz=ppos(3)-pos(q,3)
   dist2=dx*dx+dy*dy+dz*dz
   dist=SQRT(dist2)     
   epseff=peps+epsgrav(q)  
   if(dist.GE.epseff) then
    rinveff=1./dist
    r3inveff=rinveff**3
    r5inveff=rinveff**5
   else
    epseffi=1./epseff
    rhinv=dist*epseffi
    if(rhinv.LT.0.5) then
      r3inveff=4./3.+rhinv**2*(4.*rhinv-4.8)
      r5inveff=12*rhinv-9.6
    else
      r3inveff=8./3.-6.*rhinv+4.8*rhinv**2- &
            4./3.*rhinv**3-1./120./rhinv**3
      r5inveff=9.6-4.*rhinv-6./rhinv+1./40./rhinv**5
    endif   
    r3inveff=r3inveff*8*epseffi**3
    r5inveff=-r5inveff*8./3.*epseffi**5
   endif

   ptide(1)=ptide(1)+mass(q)*(-r3inveff+3*r5inveff*dx**2)              
   ptide(4)=ptide(4)+mass(q)*(-r3inveff+3*r5inveff*dy**2)              
   ptide(6)=ptide(6)+mass(q)*(-r3inveff+3*r5inveff*dz**2)              
   ptide(2)=ptide(2)+mass(q)*3*r5inveff*dx*dy              
   ptide(3)=ptide(3)+mass(q)*3*r5inveff*dx*dz              
   ptide(5)=ptide(5)+mass(q)*3*r5inveff*dy*dz              

   if(usequad.AND.q.GE.nbods1) then
    call terror("tide+quad tbd")
   endif

 enddo
end subroutine tidalsum


