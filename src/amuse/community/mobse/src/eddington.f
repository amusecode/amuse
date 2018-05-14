***
      real*8 FUNCTION eddington(mt,lum,kw)
*
*     Author : N. Giacobbo
*     Date :   19th March 2017
*
      integer kw
      real*8 mt,Xh,lum
*
*       Computes the Eddigton factor.
*       Ref: Graefener et al. 2011
*
* WR stars
      if(kw.ge.7 .and. kw.lt.10)then
         Xh = 0.d0
      else
         if(kw.eq.0 .or. kw.eq.1)then
            Xh = 0.7d0
         elseif(kw.eq.2)then
            Xh = 0.6d0
         elseif(kw.eq.3)then
            Xh = 0.5d0
         elseif(kw.eq.4)then
            Xh = 0.4d0
         elseif(kw.eq.5 .or. kw.eq.6)then
            Xh = 0.2d0
        endif
      endif
      
* equation (8) in Graefener+ 2011 but used in both Graenefer+ and Vink+
      eddington = 10.d0**(-4.813d0 + Log10(1.d0 + Xh) + Log10(lum) -
     &               Log10(mt))
*
      return
*
      end
***      
