pro cmp_memory,ndim=ndim,nvar=nvar $
               ,ngridmax=ngridmax,npart=npart,single=single

if not keyword_set(ndim)     then ndim=3L
if not keyword_set(nvar)     then nvar=0L
if not keyword_set(ngridmax) then ngridmax=0L
if not keyword_set(npart)    then npart=0L

ndim=double(ndim)
nvar=double(nvar)
ngridmax=double(ngridmax)
npart=double(npart)

INTEGER=4.d0
if keyword_set(single) then REAL=4.d0 else REAL=8.0d0

memory=   $
;particles (x,v,m,a, level,prevp,nextp)
  (2.0d0*ndim+2.0d0)*npart*REAL+3.0d0*npart*INTEGER + $ 
;particles linked list (headp, tailp, numbp)
  3.0d0*ngridmax*INTEGER*npart/(npart+1.0d0) + $
;amr (xg,son,father,nbor,next,prev
  ndim*ngridmax*REAL+(2.0d0^ndim+1.0d0+2.0d0*ndim+2.0d0)*ngridmax*INTEGER+ $
;amr(flag1, flag2, cpumap1, cpumap2)
  4.0d0*2.0d0^ndim*ngridmax*INTEGER + $
;Hilbert key (always of type double)
  2.0d0^ndim*ngridmax*8.0d0 + $
;hydro (uold, unew)
  2.0d0*nvar*2.0d0^ndim*ngridmax*REAL + $
;poisson (rho, phi, f)
  (2.0d0+ndim)*2.0d0^ndim*ngridmax*REAL*npart/(npart+1.0d0) 

print,'Total memory=',double(memory)/(1024.^3.),' Gb'


end
