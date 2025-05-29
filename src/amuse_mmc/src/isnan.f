      function isnan(x)
      implicit double precision (a-h,o-z)
      logical isnan
      if (.not.x.ge.0.d0.and..not.x.lt.0.d0) then
         isnan = .true.
      else
         isnan = .false.
      endif
      return
      end

      
