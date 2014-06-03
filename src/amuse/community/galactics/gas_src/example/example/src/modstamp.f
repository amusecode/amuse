      subroutine modstamp
      character*60 string
      open(20,file='modname',err=23,status='old')
      read(20,'(a60)') string
      close(20)
      call pgqvp(0,x1,x2,y1,y2)
      call pgqwin(xw1,xw2,yw1,yw2)
      call pgsvp(0.,1.,0.,1.)
      call pgswin(0.,1.,0.,1.)
      call pgtext(0.01,0.03,string)
      call pgsvp(x1,x2,y1,y2)
      call pgswin(xw1,xw2,yw1,yw2)
 23   return
      end
