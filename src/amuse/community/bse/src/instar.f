***
      SUBROUTINE instar
*
*
*       Initialization of collision matrix.
*       ------------------------
*
      implicit none
      integer i,j,ktype(0:14,0:14)
      common /TYPES/ ktype
*
*       Initialize stellar collision matrix.
*
      ktype(0,0) = 1
      do 10 , j = 1,6
         ktype(0,j) = j
         ktype(1,j) = j
 10   continue
      ktype(0,7) = 4
      ktype(1,7) = 4
      do 15 , j = 8,12
         if(j.ne.10)then
            ktype(0,j) = 6
         else
            ktype(0,j) = 3
         endif
         ktype(1,j) = ktype(0,j)
 15   continue
      ktype(2,2) = 3
      do 20 , i = 3,14
         ktype(i,i) = i
 20   continue
      ktype(5,5) = 4
      ktype(7,7) = 1
      ktype(10,10) = 15
      ktype(13,13) = 14
      do 25 , i = 2,5
         do 30 j = i+1,12
            ktype(i,j) = 4
 30      continue
 25   continue
      ktype(2,3) = 3
      ktype(2,6) = 5
      ktype(2,10) = 3
      ktype(2,11) = 5
      ktype(2,12) = 5
      ktype(3,6) = 5
      ktype(3,10) = 3
      ktype(3,11) = 5
      ktype(3,12) = 5
      ktype(6,7) = 4
      ktype(6,8) = 6
      ktype(6,9) = 6
      ktype(6,10) = 5 
      ktype(6,11) = 6
      ktype(6,12) = 6
      ktype(7,8) = 8
      ktype(7,9) = 9
      ktype(7,10) = 7
      ktype(7,11) = 9
      ktype(7,12) = 9
      ktype(8,9) = 9
      ktype(8,10) = 7
      ktype(8,11) = 9
      ktype(8,12) = 9
      ktype(9,10) = 7
      ktype(9,11) = 9
      ktype(9,12) = 9
      ktype(10,11) = 9
      ktype(10,12) = 9
      ktype(11,12) = 12
      do 35 , i = 0,12
         ktype(i,13) = 13
         ktype(i,14) = 14
 35   continue
      ktype(13,14) = 14
*
* Increase common-envelope cases by 100.
      do 40 , i = 0,9
         do 45 , j = i,14
            if(i.le.1.or.i.eq.7)then
               if(j.ge.2.and.j.le.9.and.j.ne.7)then
                  ktype(i,j) = ktype(i,j) + 100
               endif
            else
               ktype(i,j) = ktype(i,j) + 100
            endif
 45      continue
 40   continue
*
*       Assign the remaining values by symmetry.
      do 50 , i = 1,14
         do 55 , j = 0,i-1
            ktype(i,j) = ktype(j,i)
 55      continue
 50   continue
*
      return
      end
***
