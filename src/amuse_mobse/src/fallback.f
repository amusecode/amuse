***
      real*8 FUNCTION fallback(sne,mc,mt)
*
*     Author : N. Giacobbo
*     Date :   25th March 2017
*
*     Fallback of matter after the SNe. 
*     Ref: Fryer et al. 2012 and Spera et al. 2015
*     --------------------------------------------------
*
*
      implicit none
      real*8 mt,mc
      real*8 mproto,alpha_R,beta_R,alpha_D,beta_D
      integer sne

      if(sne.eq.1)then
* StarTrack      
         if(mc.lt.5.d0)then
	    fallback = 0.d0
	 elseif(5.d0.le.mc .and. mc.lt.7.6d0)then
	    fallback = 0.378d0*mc - 1.889d0
	 elseif(7.6d0.le.mc)then
	    fallback = 1.d0
	 endif
* 
      elseif(sne.eq.2)then
* Rapid
	 mproto = 1.d0
	 alpha_R = 0.25d0 - 1.275d0/(mt - mproto)
	 beta_R = 1.d0 - 11.d0*alpha_R
*	 
	 if(mc.lt.2.5d0)then
	    fallback = 0.2d0/(mt - mproto)
	 elseif(2.5d0.le.mc .and. mc.lt.6.d0)then
	    fallback = (0.286d0*mc - 0.514d0)/(mt - mproto)
	 elseif(6.d0.le.mc .and. mc.lt.7.d0)then
	    fallback = 1.d0
	 elseif(7.d0.le.mc .and. mc.lt.11.d0)then
	    fallback = alpha_R*mc + beta_R
	 elseif(11.d0.le.mc)then
	    fallback = 1.d0
	 endif
*
      elseif(sne.eq.3)then
* Delayed
	 if(mc.lt.3.5d0)then
	    mproto = 1.2d0
	 elseif(3.5d0.le.mc .and. mc.lt.6.d0)then
	    mproto = 1.3d0
	 elseif(6.d0.le.mc .and. mc.lt.11.d0)then
	    mproto = 1.4d0
	 elseif(11.d0.le.mc)then
	    mproto = 1.6d0
	 endif
*
	 alpha_D = 0.133d0 - 0.093d0/(mt - mproto)
	 beta_D = 1.d0 - 11.d0*alpha_D
*	 
	 if(mc.lt.2.5d0)then
	    fallback = 0.2d0/(mt - mproto)
	 elseif(2.5d0.le.mc .and. mc.lt.3.5d0)then
	    fallback = (0.5d0*mc - 1.05d0)/(mt - mproto)
	 elseif(3.5d0.le.mc .and. mc.lt.11.d0)then
	    fallback = alpha_D*mc + beta_D
	 elseif(11.d0.le.mc)then
	    fallback = 1.d0
	 endif
*
      elseif(sne.eq.4)then
* Belczynski 2002
	 fallback = 0.d0
      elseif(sne.eq.5)then
* Only remnants
	 fallback = 0.d0
      else
         write(*,*) "Invalid SNe mechanism"
      endif
*
      return
**
      end