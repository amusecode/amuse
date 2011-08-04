      function erfcen(x)
c      write(*,*) 'erfcen: x=',x
      if (x.gt.5.) then 
         erfcen=2.5066283
      elseif (x.lt.-5.) then
         erfcen=-2.5066283
      else
         erfcen=erf(x/1.4142136)*2.5066283
         endif
      return
      end
