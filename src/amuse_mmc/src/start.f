      subroutine start
*
*
*       initialization of initial model
*       -------------------------------
*
      include 'common.h'
*
      real*4 runtim,cpp      
*
      integer n
*
*
*       inicialize cpu time
*
      cpp = runtim(0.0)
      cpu = cpp
      itime = 0
      timeold = 0.0
*
      n = nt
      ntnew = nt
*
*       initialize global scalars,arrays, counters & useful constants
*
      call zero
*
*       set initial conditions: body(i), r(i), vr(i), vt(i)
*
      print*,'calling data'
      call flush(6)
      call data
      print*,'called data'
*
*       scale initial conditions to N-body units
*
      call scale0
*
*       define mean mass in scaled units 
*
      bodym = smt/float(n)
*
*
      return
*
      end
*
*
*
*
