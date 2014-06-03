

      subroutine contourden(n,d,dens,asp,rscale)
      external dens
c contour function dens(r,z) with nxn bins of size d, and aspect ratio asp.
c assumes that a PGPLOT viewport has been defined before the call; inside
c this space a plot with properly scaled axes is made.
      dimension tr(6),ddens(101,101),clev(25)
      dmax=-1e31
      do i=0,n
         r=d*i
         do j=0,n
            z=d*j*asp
            ddens(i+1,j+1)=dens(r/rscale,z/rscale)
            dmax=max(dmax,ddens(i+1,j+1))
         enddo
      enddo
      tr(1)=-d
      tr(2)=d
      tr(3)=0
      tr(4)=-d*asp
      tr(5)=0
      tr(6)=d*asp
      clev(1)=dmax/2.
      do ic=2,25
         clev(ic)=clev(ic-1)/2.
      enddo
      call pgwnad(0.,n*d,0.,n*d*asp)
      call pgbox('BCNST',0.,0,'BCNST',0.,0)
      call pgcons(ddens,101,101,1,n+1,1,n+1,clev,25,tr)
      return

      end
