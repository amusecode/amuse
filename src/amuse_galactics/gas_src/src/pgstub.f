! nicked and adapted from IFEFFIT, the  Interactive XAFS Analysis Library
       SUBROUTINE PGVPORT (XLEFT, XRIGHT, YBOT, YTOP)
       REAL XLEFT, XRIGHT, YBOT, YTOP
       END

       SUBROUTINE PGWNAD (X1, X2, Y1, Y2)
       REAL X1, X2, Y1, Y2
       END
       
       SUBROUTINE PGCONS (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR)
       INTEGER IDIM, JDIM, I1, I2, J1, J2, NC
       REAL    A(IDIM,JDIM), C(*), TR(6)
       END
       
       INTEGER FUNCTION PGBEGIN (UNIT, FILE, NXSUB, NYSUB)
       INTEGER       UNIT
       CHARACTER*(*) FILE
       INTEGER       NXSUB, NYSUB
       PGBEGIN=1
       END
       
       SUBROUTINE PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
       REAL XLEFT, XRIGHT, YBOT, YTOP
       END
       
       SUBROUTINE PGQWIN (X1, X2, Y1, Y2)
       REAL X1, X2, Y1, Y2
       END
       
       SUBROUTINE PGQVP (UNITS, X1, X2, Y1, Y2)
       INTEGER UNITS
       REAL    X1, X2, Y1, Y2
       END

       SUBROUTINE PGLAB (XLBL, YLBL, TOPLBL)
       CHARACTER*(*) XLBL, YLBL, TOPLBL
       END

       SUBROUTINE PGLABEL (XLBL, YLBL, TOPLBL)
       CHARACTER*(*) XLBL, YLBL, TOPLBL
       END
       
       SUBROUTINE PGENV (XMIN, XMAX, YMIN, YMAX, JUST, AXIS)
       REAL XMIN, XMAX, YMIN, YMAX
       INTEGER JUST, AXIS
       END

       subroutine pgend 
       return
       end

       subroutine pgclos 
       return
       end

       integer function pgopen(i)
       integer i
       pgopen = 0
       return
       end

       subroutine  pgqid(i)
       integer i
       i = 0
       return
       end
c
       subroutine pgslct(i)
       integer i
       return
       end
c
       subroutine  pgbeg(i, file, j, k)
       integer i,j, k
       character*(*) file
       print*, ' pgplot not installed'
       return
       end
c
       subroutine pgqinf(type, arg, i)
       integer i
       character*(*) type, arg
       i = 0
       return
       end
c
       subroutine pgpage
       return
       end
c
       subroutine  pgbbuf 
       return
       end
c
       subroutine pgask(flag)
       logical flag
       return
       end
c
       subroutine pgeras
       return
       end
c
       subroutine pgsls(i)
       integer i
       return
       end
c
       subroutine pgsch(x)
       real x
       return
       end
c
       subroutine pgscf(i)
       integer i
       return
       end
c
       subroutine pgslw(i)
       integer i
       return
       end
c
       subroutine pgvstd
       return
       end
c
       subroutine pgpt1(x1,x2,i)
       real x1,x2
       integer i
       return
       end
c
       subroutine pgrnge(x1,x2,x3,x4)
       real x1,x2,x3,x4
       return
       end
c
       subroutine pgsci(i)
       integer i
       return
       end
c
       subroutine pgswin(x1,x2,x3,x4)
       real x1,x2,x3,x4
       return
       end
c
       subroutine pgbox(s1,x1,i1,s2,x2,i2)
       integer i1, i2
       character*(*) s1, s2
       real x1, x2
       return
       end
c
       subroutine pgmtxt(s1, x1, x2, x3, s2)
       character*(*) s1, s2
       real x1, x2, x3
       return
       end
c
       subroutine pgline(i1,x1,x2)
       integer i1
       real x1(*), x2(*)
       return
       end
c
       subroutine pgpt(i1,x1,x2,i2)
       integer i1, i2
       real x1(*), x2(*)
       return
       end
c
       subroutine pgtext(x1, x2, s1)
       character*(*) s1
       real x1, x2
       return
       end
c
       subroutine pgebuf
       return
       end
c
       subroutine pgscrn(i1,s1,i2)
       character*(*) s1
       integer i1, i2
       return
       end
c
       subroutine pgscr(i,r,g,b)
       integer i
       real r, g, b
       return
       end
c
       integer function pgcurs(x,y,c)
       real x,y
       character*(*) c
       pgcurs = 1
       x = 0
       y = 0
       c = 'a'
       return
       end
c      
       subroutine  pgsah(i,x,y)
       real x,y
       integer i
       return
       end
       subroutine  pgarro(x,y,u,v)
       real x,y,u,v
       return
       end

       integer function pgband(m,n,x1,y1,x,y,c)
       real x1,y1,x,y
       character*(*) c
       integer m,n
       pgband = 1
       x = x1
       y = y1
       c = 'a'
       m = 0
       return
       end
       
       subroutine pgqndt(i)
       integer i
       return
       end

       subroutine pgqdt(i,s1,i2,s2,i3,i4)
       integer i, i2, i3, i4
       character s1*(*), s2*(*)
       return
       end

       subroutine pgerry(n,x,y1,y2,t)
       integer i
       real x(*), y1(*), y2(*), t
       return
       end
       subroutine pgerrx(n,x,y1,y2,t)
       integer i
       real x(*), y1(*), y2(*), t
       return
       end
