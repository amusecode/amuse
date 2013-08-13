      subroutine sort3(arr, ind, n)
c
c -------------------------------------------------------------------------
c     SORT3: Sort array of indices ind(1:n) using array of values,
c     so that for all i <= j, arr(ind(i)) <= arr(ind(j)).
c
c     arr        real*8 array (*)   in    keys for sorting
c     ind        INT  array (*)     out   indices sorted according to arr
c     n          INT                in    number of data
c--------------------------------------------------------------------------
c
      real*8 arr(*)
      integer ind(*), n
      integer LSTK, sp, i, j, k, l, ij, m
      parameter(LSTK = 32)
      integer iu(LSTK), il(LSTK)
      real*8 tij
c
      sp = 1
      i = 1
      j = n
 5    continue
      if (i .lt. j) then
         ij = (j + i)/2
         if (arr(ind(i)) .gt. arr(ind(ij))) then
            m = ind(ij)
            ind(ij) = ind(i)
            ind(i) = m
         endif
         if (arr(ind(j)) .lt. arr(ind(ij))) then
            m = ind(ij)
            ind(ij) = ind(j)
            ind(j) = m
            if (arr(ind(i)) .gt. arr(ind(ij))) then
               m = ind(ij)
               ind(ij) = ind(i)
               ind(i) = m
            endif
         endif
         tij = arr(ind(ij))
         k = i
         l = j
 40      continue
         l = l - 1
         if (arr(ind(l)) .gt. tij) go to 40
 50      continue
         k = k + 1
         if (arr(ind(k)) .lt. tij) go to 50
         if (k .le. l) then
            m = ind(l)
            ind(l) = ind(k)
            ind(k) = m
            go to 40
         endif
         if (l - i .gt. j - k) then
            il(sp) = i
            iu(sp) = l
            i = k
         else
            il(sp) = k
            iu(sp) = j
            j = l
         endif
         sp = sp + 1
      else
         sp = sp - 1
         if (sp .eq. 0) return
         i = il(sp)
         j = iu(sp)
      endif
      go to 5
      end
