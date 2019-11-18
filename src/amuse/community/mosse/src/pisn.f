***
      real*8 FUNCTION pisn(kw,mche,mt)
*
*     Author : N. Giacobbo
*     Date :   27th March 2017
*
*
*     Pari Instability and Pulsation Pair Instability during the 
*     supernova explosions.
*     Ref: Woosley 2016 and SEVN implementation (Spera & Mapelli 2017)
*     ----------------------------------------------------------
*     pisn = pair instability correction
*
*
      implicit none
      real*8 hefrac,kpar,val,mche,mt
      integer kw
*
* No correction
*
      pisn = 1.d0
*     
*     Pulsation Pair Instability PPISN
*     
*     use Table 2 fit from Woosley 2006
*
      if(kw.lt.7)then

         if(mche.gt.32.d0 .and. mche.lt.64.d0)then
*
*     I suppose the total he_tot is well described by the
*     the core mass
*
            hefrac = mche/mt
            kpar = 0.67d0*hefrac + 0.1d0
            if(mche.le.37.d0)then
               pisn = (kpar - 1.d0)/5.d0*mche + 
     &				(37.d0 - 32.d0*kpar) / 5.d0
            elseif(mche.gt.37.d0 .and. mche.le.60.d0)then
               pisn = kpar
            else
               pisn = - (kpar/4.d0)*mche + 16.d0*kpar
            endif
*
*     Pair-instability supernova: disintegrated
*
         elseif(mche.ge.64.d0 .and. mche.lt.135.d0)then
            pisn = 0.d0
*
*     end of PISN... again standard direct collapse
*
         else
            pisn = 1.d0
         endif


*
*     use WR table 1 fit from Woosley 2016
*
      elseif(kw.ge.7 .and. kw.le.9)then
*
*     I suppose the is composed only by he and co
*
         if(mt.gt.32.d0 .and. mt.lt.64.d0)then
            
            hefrac = 1.d0
            if(mt.le.37.d0)then
                pisn = (0.5226d0*hefrac - 0.52974d0)*(mt - 32.d0) +
     &              1.d0
            elseif(mt.gt.37.d0 .and. mt.le.56.d0)then
                val = (0.5226d0*hefrac - 0.52974d0)*5.d0 + 1.d0
                if(val.lt.0.82916d0)then
                   pisn = val
                else
                   pisn = (- 0.1381d0*hefrac + 0.1309d0)*(mt - 56.d0)
     &                  + 0.82916d0
                endif
            else
                pisn = - 0.103645d0*mt + 6.63328d0
            endif
*
*     Pair-instability supernova: disintegrated
*
         elseif(mt.ge.64.d0 .and. mt.lt.135.d0)then
            pisn = 0.d0
*
*     end of PISN... again standard direct collapse
*
         else
            pisn = 1.d0
         endif

      else
         write(*,*) 'Which are those cases?'

      endif
*
      return
**
      end
