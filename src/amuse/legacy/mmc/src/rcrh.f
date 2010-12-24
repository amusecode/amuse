      real*8 t,tf,r005,r01,r02,r03,r04,r05,r1,r2,r3,r4,r5,r6,r7,r8,
     &       r9,r10,smt,rh,rt,smc,rc,vc,roc,x,xx,xxx,t1,t2
      integer is,nt,nb,nc,i,ii,k,kk,kkk,l,ll
      dimension x(14),xx(4),xxx(4),i(4),ii(13)
*
      open(10,file='lagrangi.dat')
      open(11,file='system.dat')
      open(12,file='core.dat')
      open(13,file='rcrh.dat')
*
 10   continue
*
      read(10,100,end=20) is,t,tf,r005,r01,r02,r03,r04,r05,r1,r2,r3,r4,
     &                    r5,r6,r7,r8,r9,r10
      read(11,110) is,t1,tf,smt,(x(k),k=1,14),rh,rt,nt,(i(l),l=1,4),nb,
     &             (ii(ll),ll=1,13),(xx(kk),kk=1,4)
      read(12,120) is,t2,tf,smc,rc,vc,roc,(xxx(kkk),kkk=1,4),nc
      write(13,130) is,t,tf,r005,r01,r02,r03,r04,r05,r1,r2,r3,r4,
     &              r5,r6,r7,r8,r9,r10,rc,rh,rt,smc,vc,roc,nc,nt,nb
*
      if(t.ne.t1) then
        print*,'t,t1,t2 = ', t,t1,t2
        stop
      endif
 100  format(1x,i5,1p18e12.4)
 110  format(1x,i5,1p7e12.4,1pe15.7,1p11e12.4,19i10,1p4e12.4)
 120  format(1x,i5,1p10e12.4,i6)
 130  format(1x,i5,1p24e12.4,3i7)
*
      go to 10
*
20    continue
*
      close(10)
      close(11)
      close(12)
      close(13)  
*
      stop      
*
      end      