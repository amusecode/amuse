      implicit double precision (a-h,o-z)
      double precision m1,m2,time
      common /zset/ zini
      read (5,*) m1,m2,a,e,z
      in = 1
      id = 1
      if (m2.gt.0.d0) then
         call init_binaries(in,id,a,e,m1,m2)
         call evbinary(id,tnow)
 10      continue
         call bsupdatetime(id,time)
         call evbinary(id,time)
         if (time.lt.5e3) go to 10
      else
      stop
      end
