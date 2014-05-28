! wat als idiskflag=2 etc?
      function dens(r,z)
      parameter (pi=3.1415926535)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
     
      if( idiskflag .eq. 1 .or. idiskflag.eq.3) then
	 addens = appdiskdens(r,z)
         dens=totdens(r,z)
         dens = dens - addens
!         if(dens.LT.0) dens=0
      else
         dens=totdens(r,z)
      endif
      return
      end
