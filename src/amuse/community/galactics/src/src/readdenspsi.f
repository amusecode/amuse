      subroutine readdenspsibulge

      include 'commonblocks'

      open(file='denspsibulge.dat',unit=50,status='old')

      do i=1,npsi
         read(50,*) tableE(i),denspsibulge(i)
      enddo
      close(50)

      open(file='dfsersic.dat',unit=50,status='old')

      do i=1,npsi
         read(50,*) tableE(i),dfsersic(i)
      enddo

      close(50)
      return
      end

      subroutine readdenspsihalo

      include 'commonblocks'

      open(file='denspsihalo.dat',unit=50,status='old')

      do i=1,npsi
         read(50,*) tableE(i),denspsihalo(i)
      enddo
      close(50)

      open(file='dfnfw.dat',unit=50,status='old')

      do i=1,npsi
         read(50,*) tableE(i),dfnfw(i)
      enddo

      close(50)
      return
      end
