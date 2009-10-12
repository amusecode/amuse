! Solver for the TWIN stellar evolution code.
! Written by Evert Glebbeek (2009)
!
! This solver uses the same API as the original solver package in order
! that it can be plugged in without modifying other parts of the code.
!
! Canonically, the jacobian that needs to be inverted during the solution
! process (also called "Henyey matrix" in astrophysical literature) has a
! block diagonal structure. In TWIN the matrix is more complicated for two
! reasons. First of all by the presence of second order equations, secondly
! by the presence of "eigenvalues". The matrix has the following schematic 
! form (columns label variables, rows label equations and the matrix
! elements are the derivatives of the equations to the variables):
!
! ccc             e     ?
! ccc             e     ?
! ... ...         e     ?
! xxx xxx         e     ?
! xxx xxx         e     ?
! ... ... ...     e     ?
!     xxx xxx     e     ?
!     xxx xxx     e     ?
!     ... ... ... e     ?
!         xxx xxx e     ?
!         xxx xxx e     ?
!         ... ... e     ?
!             sss e     ?
!             sss e     ?
! Here, 'ccc' represents the block of central boundary conditions,
! 'sss' represents the block of surface boundary conditions, the 'xxx'
! blocks refer to stellar structure (1st order diff. eqn.) blocks and the
! '...' blocks refer to composition (2nd order diff eqn.) blocks.
! Gridpoints are seperated by a space. The points labelled 'e' involve the
! eigenvalues (the location of this column within the matrix is essentially
! but not practically arbitrary). The RHS of the equation is indicated by '?'.
! To proceed, we divide the matrix in four parts: the square block diagonal
! part without the eigenvalues, the column with the eigenvalue
! derivatives, the bottom rows below the block diagonal that needs to be
! removed to make the block diagonal part square and the lower right corner
! where the removed rows and columns join. To solve the system of equations
! we will perform the following steps in order: inversion of the square
! band-diagonal part of the matrix, Gaussian elimination of the rows below
! the block matrix, inversion of the small matrix in the lower right corner
! and backsubstitution along the column.
! While calculating the inverse of the block matrix we need to keep track
! changes this causes in the far right column. In Gauss-Jordan elimination
! this is done by performing the same row operations as on the rest of the
! matrix, but this is equivalent to treating the column as a right hand
! side vector. During the first step we therefore pass both the real RHS
! and the right-most columns of the matrix to the routine that solves the
! block diagonal part.
! The next step is to complete the Gaussian elimination of the bottom rows
! below the block diagonal; this is straightforward because the rows above
! it have been cleared to unity.
! The small NEV x NEV matrix in the corner needs to be inverted next. Its
! size varies between 1 (single non-rotating star) and 9 (TWIN binary). We
! should use a regular matrix solver for this (we may need pivoting) but
! for now we don't.
! Finally the column is cleared by backsubstitution from the lower
! righthand block.
      MODULE JACOBIAN
      DOUBLE PRECISION, ALLOCATABLE :: BANDED_JACOBIAN(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: RHS(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EV_BLOCK(:,:)
      INTEGER, ALLOCATABLE :: PIVOT(:)
      INTEGER :: NSUP, NSUB, LDA, NBAND, NJAC
      INTEGER :: NSBC, NRHS

      ! Block jacobian: NROWS x 3xNVAR + NEV + RHS
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK_JACOBIAN(:,:)
      ! Number of blank columns at the start of each row and number of
      ! non-zero columns in each blok.
      INTEGER, ALLOCATABLE :: ROWSTART(:)
      INTEGER, ALLOCATABLE :: ROWWIDTH(:)
      ! First column where we may search for a pivot in this row
      INTEGER, ALLOCATABLE :: ROWPIVOT(:)
      ! Ditto, with rows and columns interchanged:
      ! Number of blank rows at the start of each column and number of
      ! non-zero rows in each blok.
      INTEGER, ALLOCATABLE :: COLSTART(:)
      INTEGER, ALLOCATABLE :: COLHEIGHT(:)
      ! Row and column pivoting
      INTEGER, ALLOCATABLE :: RP(:)
      INTEGER, ALLOCATABLE :: RHP(:)
      INTEGER, ALLOCATABLE :: CP(:)
      ! Maximum row to search for a suitable pivot from this row
      INTEGER, ALLOCATABLE :: MAXPIVOT(:)
      ! Shape of the block Jacobian: eigenvalue and RHS columns
      INTEGER :: BLOCKW, CEV, CRHS, NBLOCK
      contains

      subroutine set_jacobian(irow, jcol, x)
      implicit none
      integer, intent(in) :: irow, jcol
      double precision, intent(in) :: x
      integer :: jj

      if (irow < 1 .or. jcol < 1 .or. irow>NJAC .or. jcol > NJAC) return

      !print *, irow, jcol, x
      if (jcol > NBAND) then               ! Right-hand column
         RHS(irow, jcol - NBAND) = x
      else if (irow > NBAND) then          ! Bottom rows
         EV_BLOCK(jcol, irow-NBAND) = x
      else if (max(1, jcol-NSUP) <= irow .and. irow<=min(NBAND,jcol+NSUB)) then
         ! Block diagonal
         jj = LDA-NSUB+irow-jcol
         BANDED_JACOBIAN(jj, jcol) = x
      else if (x/=0.0d0) then
         !print *, 'Should never get here!', irow, jcol, x
         !pause
      end if

      ! Block jacobian
      call set_block_jacobian(irow, jcol, x)
      end subroutine

      function get_jacobian(irow, jcol)
      implicit none
      double precision :: get_jacobian
      integer, intent(in) :: irow, jcol
      double precision :: x1, x2
      integer :: jj
      get_jacobian = 0.0d0

      x1 = 0.0d0
      if (jcol > NBAND) then               ! Right-hand column
         x1 = RHS(irow, jcol - NBAND)
      else if (irow > NBAND) then          ! Bottom rows
         x1 = EV_BLOCK(jcol, irow-NBAND)
      else if (max(1, jcol-NSUP) <= irow .and. irow<=min(NBAND,jcol+NSUB)) then
         ! Block diagonal
         jj = LDA-NSUB+irow-jcol
         x1 = BANDED_JACOBIAN(jj, jcol)
      else
!         print *, 'Attempting to read matrix element outside band jacobian',
!     &            irow, jcol
      end if
      get_jacobian = x1
      end function

      ! At each row, set elements that are < 1.e-16*max element in row to
      ! since they are essentially 0 to numerical accuracy anyway
      subroutine regularise_jacobian
      implicit none
      integer :: irow, jcol
      double precision :: x, y(NJAC)
      double precision, parameter :: eps = 1.0d-16

      do irow = 1, NJAC
         ! Find maximum element in row
         x = 0.0d0
         do jcol = 1, NJAC
            y(jcol) = dabs(get_jacobian(irow, jcol))
            x = max(x, y(jcol))
         end do
         ! Set small elements to 0
         do jcol = 1, NJAC
            if (y(jcol)>0.0d0 .and. y(jcol)<x*eps) then
               !print *, jcol, y(jcol), x
               call set_jacobian(irow, jcol, 0.0d0)
            end if
         end do
      end do
      end subroutine

      END MODULE

 
! Compute (and store) unperturbed functions for all meshpoints from k1
! to k2 in steps of dk. In effect, computes the RHS of the matrix
! equation
      subroutine compute_functions (k1, k2, dk)
      use mesh
      use jacobian
      implicit none
      integer, intent(in) :: k1, k2, dk
      integer :: jk
      integer :: keq, kj2, kj5, kj6, kj10, kj12, ki4, kee, ia, ke1, ke2,
     &        ke3, kbc, kev, kfn, kl, jh1, jh2, jh3, kd, kq, kvb, kvc
      double precision :: C, S, ER
      common /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :         KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC, 
     :         KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC

C Evaluate functions at first meshpoint
      JK = k1 + KL
      CALL DIFRNS ( JK, KEQ - KBC+1, KEQ, 80 + KEV, 1 ) 
      JK = JK + dk
      IF ( k1 == k2 ) return
C Ditto second
      CALL DIFRNS ( JK, 1, KEQ, 40, 1 )
      DO JK = k1 + KL + 2*dk, k2 + KL, dk
C Ditto remaining meshpoints except last
         CALL DIFRNS ( JK, 1, KEQ, 40, 1)
      END DO
C Ditto last meshpoint
      CALL DIFRNS ( JK, 1, KEQ-KBC + KEV, 80, 1) 
      end subroutine

! Compute the Jacobian for all meshpoints from k1 to k2 in steps of dk.
! In effect, computes the RHS of the matrix equation.
      subroutine compute_jacobian(k1, k2, dk)
      use mesh
      use jacobian
      implicit none
      integer, intent(in) :: k1, k2, dk
      integer :: jk
      integer :: keq, kj2, kj5, kj6, kj10, kj12, ki4, kee, ia, ke1, ke2,
     &        ke3, kbc, kev, kfn, kl, jh1, jh2, jh3, kd, kq, kvb, kvc
      double precision :: C, S, ER
      common /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :         KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC, 
     :         KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC

C Evaluate functions at first meshpoint
      JK = k1 + KL
      CALL DIFRNS ( JK, KEQ - KBC+1, KEQ, 80 + KEV, 0 ) 
      JK = JK + dk
      IF ( k1 == k2 ) return
C Ditto second
      CALL DIFRNS ( JK, 1, KEQ, 40, 0 )
      DO JK = k1 + KL + 2*dk, k2 + KL, dk
C Ditto remaining meshpoints except last
         CALL DIFRNS ( JK, 1, KEQ, 40, 0)
      END DO
C Ditto last meshpoint
      CALL DIFRNS ( JK, 1, KEQ-KBC + KEV, 80, 0) 
      end subroutine

! get_block_jacobian:
!  Read a matrix element from the block Jacobian. Returns 0 if the requested
!  element is outside the Jacobian.
      function get_block_jacobian(iirow, jjcol)
      use jacobian
      implicit none
      double precision :: get_block_jacobian
      integer, intent(in) :: iirow, jjcol
      integer :: irow, jcol
      double precision :: x

      irow = rp(iirow)
      jcol = cp(jjcol)

      x = 0.0d0
      if (jcol > NBLOCK) then     ! EV column
         x = BLOCK_JACOBIAN(irow, CEV + jcol - NBLOCK - 1)
      else if (jcol<=rowstart(irow) .or. jcol>rowstart(irow)+rowwidth(irow)) then
!         print *, 'Attempting to read matrix element outside block jacobian',
!     &            irow, jcol, rowstart(irow), rowwidth(irow)
      else
         !print *, 'get>>>', irow, jcol - rowstart(irow)
!      else if (jcol>rowstart(irow).and.jcol<=rowstart(irow)+rowwidth(irow)) then
         x = BLOCK_JACOBIAN(irow, jcol - rowstart(irow))
      end if
      get_block_jacobian = x
      end function

! sel_block_jacobian:
!  Set the value of an element in the block Jacobian. Does nothing if the
!  element to be set is outside the Jacobian.
      subroutine set_block_jacobian(iirow, jjcol, x)
      use jacobian
      implicit none
      integer, intent(in) :: iirow, jjcol
      integer :: irow, jcol
      double precision, intent(in) :: x

      irow = rp(iirow)
      jcol = cp(jjcol)

      if (jcol > NBLOCK) then     ! EV column
         BLOCK_JACOBIAN(irow, CEV + jcol - NBLOCK - 1) = x
      else if (jcol<=rowstart(irow) .or. jcol>rowstart(irow)+rowwidth(irow)) then
!         print *, 'Attempting to read matrix element outside block jacobian',
!     &            irow, jcol, rowstart(irow), rowwidth(irow)
      else
         !print *, 'get>>>', irow, jcol - rowstart(irow)
!      else if (jcol>rowstart(irow).and.jcol<=rowstart(irow)+rowwidth(irow)) then
         BLOCK_JACOBIAN(irow, jcol - rowstart(irow)) = x
      end if
      end subroutine

! swap_rows:
!  Swaps the elements in the row pivot array
      subroutine swap_rows(row1, row2)
      use jacobian
      implicit none
      integer, intent(in) :: row1, row2
      integer :: row
      if (row1 /= row2) then
         row = rp(row1)
         rp(row1) = rp(row2)
         rp(row2) = row
      end if
      end subroutine


      subroutine swap_columns(column1, column2)
      use jacobian
      implicit none
      integer, intent(in) :: column1, column2
      integer :: column
      double precision :: x
      if (column1 /= column2) then
         column = cp(column1)
         cp(column1) = cp(column2)
         cp(column2) = column

         ! Swap elements in RHS vector
         !x = block_jacobian(cp(column1), CRHS)
         !block_jacobian(cp(column1), CRHS) = block_jacobian(cp(column2), CRHS)
         !block_jacobian(cp(column2), CRHS) = x
      end if
      end subroutine

      subroutine locate_pivot_element(pivot_row, pivot_col)
      use jacobian
      implicit none
      integer, intent(in) :: pivot_row
      integer, intent(out) :: pivot_col
      double precision, external :: get_block_jacobian
      integer :: i, j, row, col, new_pivot_row, new_pivot_col, first, last
      integer :: imax(pivot_row:maxpivot(pivot_row))
      double precision :: x, xmax, sx, xx
      ! we perform partial pivoting, so we may do row interchanges while
      ! looking for the pivot.
      ! We keep track of this by permuting the "rp" array, which is
      ! accessed by all functions acting directly on the data of the
      ! jacobian.
      ! TODO: use an actually sane pivoting algorithm...
      new_pivot_row = pivot_row
      new_pivot_col = pivot_row
      pivot_col = pivot_row
      if (.false.) then    ! "Full" pivoting, doesn't quite work (yet)
! Locate the most significant (remaining) row. `Significance' of a row is 
! the ratio of largest |element| = VM to sum of remaining |elements| = VS 
      xmax = 0.0d0
      imax(:) = 0
      print *, 'Search pivot in rows', pivot_row, maxpivot(pivot_row)
      do i = pivot_row, maxpivot(pivot_row)
         ! Sanity check
         if (rowstart(rp(i)) /= rowstart(rp(pivot_row)))
     &      print *, 'Error in matrix'
         ! Determine max element in row and sum of remaining elements in
         ! this row
         row = rp(i)    ! Find row in original (unscrambled) matrix

         first = rowpivot(row)
         last = rowpivot(row)+blockw-1
         ! Something like this should work, but it blatently doesn't
         ! with ifort on OS X at least... so use a regular 80s-style
         ! search...
         !ii(1) = maxloc(dabs(block_jacobian(row, first:last)))
         !x  =   dabs(block_jacobian(row, ii(1)))
         sx =    sum(dabs(block_jacobian(row, first:last)))
         col = 0
         x = 0.0d0
         !print *, 'Search block columns', first, last
         do j=first, last
            if (dabs(block_jacobian(row, j)) > x) then
               !print *, j, x, dabs(block_jacobian(row, j)) 
               x = dabs(block_jacobian(row, j)) 
               col = j
            end if
         end do
         if (col == 0) then
            print *, 'Encountered singular matrix at row', i
            print '("sss",1P, I5, 35E13.5)', row,
     &      block_jacobian(row, :)
            stop
         end if
         imax(i) = col
         sx = sx - x

         if (x > 0.0d0 .and. sx == 0.0d0) then
            x = 2.0d0
         else
            x = x / (x+sx)
         end if
         if (x > xmax) then
            xmax = x
            new_pivot_row = i
         end if
      end do

      new_pivot_col = rowstart(rp(new_pivot_row)) + imax(new_pivot_row)
      ! We now know the new pivot column in the original, unscrambled
      ! matrix - but this is not good enough: we need to know the column
      ! in the current, scrambled matrix. To find this, we essentially need
      ! the inverse of the cp() array - but there is a clever way to look
      ! for this that we use here for the time being: because the column
      ! interchanges are, ultimately, cyclic, we can follow the permutations
      ! until we find the one that will return us to the starting point.
      ! This is the current column number of the column in the unscrambled
      ! matrix we found above.
      col = new_pivot_col
      do while (cp(col) /= new_pivot_col)
         col = cp(col)
      end do
      new_pivot_col = col
      !print *, new_pivot_row, rp(new_pivot_row), imax(new_pivot_row)
      print *, 'Found pivot in row, col', new_pivot_row, new_pivot_col

      ! Swap rows
      call swap_rows(pivot_row, new_pivot_row)

      ! Now swap columns such that the pivot element is on the diagonal
      call swap_columns(pivot_col, new_pivot_col)
      return
      end if

      if (get_block_jacobian(pivot_row, pivot_col) == 0.0d0) then
         do i = pivot_row, maxpivot(pivot_row)
            if (get_block_jacobian(i, pivot_col) /= 0.0d0) then
               call swap_rows(pivot_row, i)
               return
            end if
         end do
      end if
      end subroutine

      ! Divide a row by the value of the pivot element located at
      ! pivot_row, pivot_col
      subroutine divide_row(pivot_row, pivot_col)
      use jacobian
      implicit none
      integer, intent(in) :: pivot_row, pivot_col
      double precision :: x
      double precision, external :: get_block_jacobian

      ! Divide row by pivot
      x = get_block_jacobian(pivot_row, pivot_col)
      if (x /= 0.0d0) then
         x = 1.0d0/x
         block_jacobian(rp(pivot_row), :) = x*block_jacobian(rp(pivot_row), :)
      else
         print *, 'Encountered 0 pivot element - singlar matrix?'
         print *, pivot_row, pivot_col
         print *, rp(pivot_row), cp(pivot_col)
         stop
      end if

      ! The pivot element should now be 1.0, but set this explicitly
      ! (because we don't trust floating point calculations and don't
      ! want to accumulate rounding errors)
      call set_block_jacobian(pivot_row, pivot_col, 1.0d0)
      end subroutine

      subroutine eliminate_row(irow, pivot_row, pivot_col)
      use jacobian
      implicit none
      integer, intent(in) :: irow, pivot_row, pivot_col
      integer :: col1, col2, width, j
      double precision :: x
      double precision, external :: get_block_jacobian

      ! Subtract overlapping parts of the block matrix 
      ! We only need to worry about columns to the right in the scrambled
      ! matrix, those to the left should all have been set to 0 already
      col1 = cp(pivot_col) - rowstart(rp(irow))
      col2 = cp(pivot_col) - rowstart(rp(pivot_row)) 

      ! If the unscrambled pivot column is to the left of the first
      ! element in the current row, there is nothing to do, so just return.
      ! NB: col2<1 indicates and error in the way the matrix is set up
      ! because it means that the pivot is not in the row that contains it.
      if (col1 < 1) return
      if (col2 < 1) then   ! Sanity check
         print *, 'Pivot is not on the row', pivot_row, pivot_col
         print *, "I'm confused now."
         stop
      end if

      ! Width of the block to work with:
      !  ..... .p...
      !        ..... .....
      !        <--->
      ! We extend to the left of the pivot because, due to column
      ! interchanges, the column to the left in the unscrambled matrix may
      ! not have been cleared yet.
      width = min(rowstart(rp(irow))+rowwidth(rp(irow)), 
     &               rowstart(rp(pivot_row))+rowwidth(rp(pivot_row)))
      width = width - max(rowstart(rp(irow)), rowstart(rp(pivot_row)))
      col2 = rowstart(rp(irow)) - rowstart(rp(pivot_row)) + 1
      !print *, 'Block columns', col1, col2, width
!      if (col2 + width > CEV) then
!         print *, 'width exceeds width of matrix... bug?'
!         stop
!      end if
!      if (width < 1) then
!         print *, '*******'
!         print *, 'width = 0 for row, unscrambled', irow, rp(irow)
!         print *, 'row start, width:', rowstart(rp(irow)), rowwidth(rp(irow))
!         print *, 'pivot start, width:', rowstart(rp(pivot_row)), rowwidth(rp(pivot_row))
!         print *, 'pivot, unscrambled:', pivot_row, rp(pivot_row)
!         !return
!      end if

      x = get_block_jacobian(irow, pivot_col)
      if (x /= 0.0d0) then
         ! Subtract pivot from row
         if (width >= 1) block_jacobian(rp(irow), 1:width) = 
     &      block_jacobian(rp(irow), 1:width) - x *
     &      block_jacobian(rp(pivot_row), col2:col2+width-1)

         ! Subtract EV column on the right and operate on RHS
         block_jacobian(rp(irow), CEV:CRHS) = 
     &      block_jacobian(rp(irow), CEV:CRHS) - x *
     &      block_jacobian(rp(pivot_row), CEV:CRHS)

         ! Leading element should now be 0, but we don't trust floating-
         ! point operations to be perfectly accurate and we are afraid of
         ! rounding errors, so force it to 0 identically.
         call set_block_jacobian(irow, pivot_col, 0.0d0)
      end if
      end subroutine

      subroutine back_substitution(jcol)
      use jacobian
      implicit none
      integer, intent(in) :: jcol
      integer :: irow, i, minrow
      double precision :: x, y
      double precision, external :: get_block_jacobian

      ! Find pivot row, which is the same as the pivot column
      ! Remap the row in the scrambled matrix onto the location in the
      ! original matrix, which is where it is still stored in memory
      irow = cp(jcol)
      !irow = rp(jcol)

      ! RHS element in the pivot row
      x = block_jacobian(irow, CRHS)

      if (cp(jcol) > NBLOCK) then    ! EV column
         minrow = 1
      else
         minrow = max(1, jcol-2*NSUP)
      end if
      !minrow = colstart(irow)

      ! Back substitution, eliminate super rows (starting at jcol-1)
      do i = jcol-1, minrow, -1
         y = get_block_jacobian(i, jcol)
         if (y /= 0.0d0) then
            block_jacobian(rp(i), CRHS) = block_jacobian(rp(i), CRHS) - x*y
            call set_block_jacobian(i, jcol, 0.0d0)
         end if
      end do
      end subroutine

      subroutine invert_block_jacobian
      use jacobian
      implicit none
      integer :: irow, jcol
      integer :: iirow, ii, row_swap
      integer :: jcol_pivot, j, i, col1, col2, width, firstrow
      integer :: unscc(NJAC)
      double precision :: x
      double precision, external :: get_block_jacobian

      ! Do the matrix inversion, loop over all rows
      do irow = 1, NJAC
         ! Find the pivot element for this row. This may involve some
         ! row interchanges, but since we don't want to move rows and
         ! columns in the original matrix around, we keep track of this
         ! by using tag lists "rp" and "cp" (for row and column permutation)
         ! (only row interchanges are implemented at the moment)
         call locate_pivot_element(irow, jcol)
         rhp(irow) = cp(jcol)
         ! The pivot element is at irow, jcol
         call divide_row(irow, jcol)
         ! Now eliminate rows below the pivot
         do i = irow+1, colheight(cp(jcol))
            if (i == irow) cycle
            call eliminate_row(i, irow, jcol)

            ! Sanity check
            do j=rowstart(rp(i))+1, jcol - 1
               if (.false..and.j /= i .and. get_block_jacobian(i, j) /= 0.0d0) then
                  print *, 'Non-zero elements to left of pivot in row, col',i,j
                  print *, 'unscrambled row, col', rp(i), cp(j)
                  print *, 'pivot (scrambled)', irow, jcol
                  print *, 'pivot (unscrambled)', rp(irow), cp(jcol)
                  print *, 'element, pivot, start', cp(j) - rowstart(rp(i)), cp(jcol) - rowstart(rp(i)), rowstart(rp(i))
                  print *, get_block_jacobian(i, j)
                  print '("***",1P, I5, 35E13.5)', rp(i),
     &                  (block_jacobian(rp(i), ii), ii = 1, CEV-1),
     &                  block_jacobian(rp(i), CEV:)
                  stop
               end if
            end do
         end do
      end do

      ! Back substitution. This is easy, since all elements to the left of
      ! the pivot have been eliminated (by Gauss-Jordan elimination) and
      ! all elements to the right of the pivot have been eliminated by
      ! back substitution
      ! In fact, we only need to worry about the solution vector.
      do jcol=NJAC, 1, -1
         ! FIXME: back substiution should be done in the order of the
         ! original matrix: so jcol here refers to the column number in the
         ! original matrix, not in the scrambled matrix
         j = cp(jcol)
         do while (cp(j) /= jcol)
            j = cp(j)
         end do
         ! FIXME: Checking the values in the matrix after the row
         ! elimination, it appears that the full pivoting does not
         ! give the correct results at the moment.
         ! This problem occurs already in the first block of first
         ! order equations before any back substitution is done
         call back_substitution(j)
      end do
      end subroutine

      SUBROUTINE SOLVER ( ITER, IG, KT5, JO )
      USE MESH
      USE SETTINGS
      USE JACOBIAN
      use compare_floats
c     Explicitly type all variables (Steve, 5/08):
      IMPLICIT NONE

c     Arguments:
      INTEGER :: iter, ig, kt5, jo

c     Local variables:
      INTEGER :: i, ii, ij, ik, ikfst, iklst, ikm, jk, k, ke4, ki1,
     &        ki2, ki3, kj1, kj11, kj3, kj4, kj7, kj8, kj9, kk, info
      DOUBLE PRECISION :: df, err, errpr, ert, es, fac, frr, vx, vz

      DIMENSION IKM(NVAR), ERT(NVAR), ES(NEQ), IG(130)
      DOUBLE PRECISION, ALLOCATABLE :: AFB(:,:), RSCALE(:), CSCALE(:)
      DOUBLE PRECISION, ALLOCATABLE :: X(:,:), FERR(:), BERR(:), WORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: EV_MATRIX(:,:), B(:)
      INTEGER, ALLOCATABLE :: IWORK(:)
      DOUBLE PRECISION :: RCOND
      INTEGER :: EQUED, LDAFB

c     Common variables:
      INTEGER :: kh, ktw, kw
      INTEGER :: jmod, jb, jnn, jter, joc, jkh
      INTEGER :: keq, kj2, kj5, kj6, kj10, kj12, ki4, kee, ia, ke1, ke2,
     &        ke3, kbc, kev, kfn, kl, jh1, jh2, jh3, kd, kq, kvb, kvc

      DOUBLE PRECISION :: h, dh, eps, del, dh0
      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: c, s, er
      DOUBLE PRECISION :: errs, errppr, besterr, besth, bestdh
      DOUBLE PRECISION :: fval, fjac, xd
      LOGICAL :: RECOMPUTE

      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :         KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC, 
     :         KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC
      COMMON /SOLMON/ ERRS, ERRPPR, BESTERR, BESTH(NVAR,NM), BESTDH(NVAR,NM)
      COMMON /JAC   / FVAL(NM, NFUNC), FJAC(NM,NVAR,NFUNC), XD(NM), RECOMPUTE

      CHARACTER*7 VARNAMES(1:NVAR)
      DATA VARNAMES/'   ln f','   ln T','    X16','      M','     X1',
     &              '      C','   ln r','      L','     X4','    X12',
     &              '    X20','      I','   Prot','    phi','  phi_s',
     &              '    X14','  H_orb','      e','     xi','    M_B',
     &              ' *1:X24','  22   ','  23   ','  24   ','*2:ln f',
     &              '*2:ln T','*2: X16','*2:   M','*2:  X1','*2:   C',
     &              '*2:ln r','*2:   L','*2:  X4','*2: X12','*2: X20',
     &              '*2:   I','*2:Prot','*2: phi','*2:phis','*2: X14',
     &              '  41   ','  42   ','  43   ','  44   ','  45   ',
     &              '  46   ','  47   ','  48   ','  49   ','  50   ',
     &              '  51   ','  52   ','  53   ','  54   ','  55   ',
     &              '  56   ','  57   ','  58   ','  59   ','  60   ',
     &              '  61   ','  62   ','  63   ','  64   ','  65   ',
     &              '  66   ','  67   ','  68   ','  69   ','  70   ',
     &              '  71   ','  72   ','  73   ','  74   ','  75   ',
     &              '  76   ','  77   ','  78   ','  79   ','  80   ',
     &              '  81   ','  82   ','  83   ','  84   ','  85   ',
     &              '  86   ','  87   ','  88   ','  89   ','  90   ',
     &              '  91   ','  92   ','  93   ','  94   ','  95   ',
     &              '  96   ','  97   ','  98   ','  99   ',' 100   ',
     &              ' 101   ',' 102   ',' 103   ',' 104   ',' 105   ',
     &              ' 106   ',' 107   ',' 108   ',' 109   ',' 110   ',
     &              ' 111   ',' 112   ',' 113   ',' 114   ',' 115   ',
     &              ' 116   ',' 117   ',' 118   ',' 119   ',' 120   '/

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      ii = 0
      ij = 0
      ik = 0
      ikfst = 0
      iklst = 0
      ikm = 0
      jk = 0
      k = 0
      ke4 = 0
      ki1 =0
      ki2 = 0
      ki3 = 0
      kj1 = 0
      kj11 = 0
      kj3 = 0
      kj4 = 0
      kj7 = 0
      kj8 = 0
      kj9 = 0
      kk = 0
      df = 0
      err = 0
      errpr = 0
      ert = 0
      es = 0
      fac = 0
      frr = 0
      vx = 0
      vz = 0
c-----------------------------------------------------------------------

      ! Set number of equations/variables.
      ! Used to be a block copy, but does not work reliably when using
      ! aggressive optimisations in the compiler.
      KE1 = IG(1)
      KE2 = IG(2)
      KE3 = IG(3)
      KBC = IG(4)
      KEV = IG(5)
      KFN = IG(6)
      KL  = IG(7)
      JH1 = IG(8)
      JH2 = IG(9)
      JH3 = IG(10)
      KD(1:120) = IG(11:130)
      KQ = 1 - 2*KL
      KEQ = KE1 + KE2
      ERR = 0.0D0 
      IF ( KEQ.EQ.0 ) RETURN        ! NO EQUNS TO SOLVE
      KVB = KEQ + KEV
      KVC = KVB + KVB
      KE4 = KE3 + KE2
C LOCATE VERTICAL (KJ) AND HORIZONTAL (KI) PARTITIONS IN MATRIX
      KJ1 = KBC + 1
      KJ2 = KEQ + 1
      KJ3 = KEQ + KBC
      KJ4 = KJ3 + 1
      KJ5 = KEQ + KEQ
      KJ6 = KJ5 + 1
      KJ7 = KJ5 + KBC
      KJ8 = KJ7 + 1
      KJ9 = KJ5 + KEQ
      KJ10 = KJ9 + 1
      KJ11 = KJ9 + KEV
      KJ12 = KJ10 + KEV
      KI1 = KEQ - KBC
      KI2 = KI1 + 1
      KI3 = KI1 + KEV
      KI4 = KI3 + 1
C FIRST MESHPOINT IS IKFST, LAST IS IKLST
      IKFST = 1 + KL*(KH - 1)
      IKLST = KH + 1 - IKFST
C DETERMINE 'TYPICAL', USUALLY MAXIMAL, VALUES FOR EACH VARIABLE
      DO IJ = 1, KVB
         I = KD(IJ)
         IF ( JNN.GT.1 ) ES(I) = ER(I) 
         ER(I) = 1.0D0
         DO K = 1, KH
            ER(I) = DMAX1(ER(I), DABS(H(I, K)))
         END DO
         IF ( feq(ER(I), 0.0D0) ) ER(I) = 1.0D0
         IF (JNN.GT.1 .AND. nfeq(ES(I), 1.0D0)) ER(I) = 0.5D0*(ER(I) + ES(I))
      END DO
! If there are no second order equations, then we only need to take
! the current and next meshpoints into account when calculating the
! Jacobian (all other elements are 0), otherwise we also need to take
! into account the previous meshpoint.
      IF (KE4 > 0) THEN
         KEE = 1
      ELSE
         KEE = 2
      END IF

c Allocate space for storage of banded jacobian
c FIXME: this should probably only be done once when the program starts
      NSBC = KE1 - KBC + KEV           ! Number of linear equations at surface
      IF (KEV < KBC) NSBC = KEV
      NSUP = KEQ + KE2 + KBC + 1       ! Superdiagonals
      NSUB = KEQ + KE2 + KBC + 1       ! Subdiagonals
      LDA = NSUB+NSUP+1                ! Leading dimension of jacobian
      LDAFB = LDA+NSUB                 ! Leading dimension for workspace
      NJAC = KH*KEQ + KEV              ! Total number of equations to solve
      NBAND = NJAC - NSBC              ! Number of equations in band matrix
      NRHS = 1 + NSBC                  ! Number of columns to the right
      ALLOCATE(BANDED_JACOBIAN(LDA,NBAND))
      ALLOCATE(RHS(NJAC, NRHS))
      ALLOCATE(EV_BLOCK(NBAND, NSBC))
      ALLOCATE(PIVOT(NBAND))

      ALLOCATE(AFB(LDAFB,NBAND))
      ALLOCATE(RSCALE(NBAND))
      ALLOCATE(CSCALE(NBAND))
      ALLOCATE(X(NBAND,NRHS))
      ALLOCATE(FERR(NRHS))
      ALLOCATE(BERR(NRHS))
      ALLOCATE(WORK(3*NBAND))
      ALLOCATE(IWORK(NBAND))
      IF (NSBC > 0) THEN
         ALLOCATE(EV_MATRIX(NSBC, NSBC))
         ALLOCATE(B(NSBC))
      END IF

      BANDED_JACOBIAN(:, :) = 0.0D0
      RHS(:,:) = 0.0D0
      EV_BLOCK(:,:) = 0.0D0
      C(:,:,:) = 0.0

! Allocate space for the block jacobian
      ALLOCATE(BLOCK_JACOBIAN(NJAC, 3*KEQ + KEV + 1))
      ALLOCATE(ROWSTART(NJAC))
      ALLOCATE(ROWWIDTH(NJAC))
      ALLOCATE(COLSTART(NJAC))
      ALLOCATE(COLHEIGHT(NJAC))
      ALLOCATE(RP(NJAC))
      ALLOCATE(RHP(NJAC))
      ALLOCATE(CP(NJAC))
      ALLOCATE(MAXPIVOT(NJAC))
      ALLOCATE(ROWPIVOT(NJAC))

! Store shape of the block Jacobian, that is to say, for each row, store
! the number of empty columns at the beginning and the length of the
! stretch of non-empty elements
      CEV = 3*KEQ + 1
      CRHS = 3*KEQ + KEV + 1
      NBLOCK = KH*KEQ
      BLOCKW = KEQ
! Central boundary condition block: offset 0, KBC x KEQ wide
      DO IJ = 1, KBC
         ROWSTART(IJ) = 0
         ROWWIDTH(IJ) = KEQ
         MAXPIVOT(IJ) = KBC
         ROWPIVOT(IJ) = 1
      END DO
! Second order equations at central meshpoint
      DO IJ = 1, KE2
         ROWSTART(IJ + KBC) = 0
         ROWWIDTH(IJ + KBC) = 2*KEQ
         MAXPIVOT(IJ + KBC) = KBC + KE2
         ROWPIVOT(IJ + KBC) = 1 + KBC
      END DO
! First order equations at central meshpoint
      DO IJ = 1, KE1
         ROWSTART(IJ + KE2 + KBC) = 0
         ROWWIDTH(IJ + KE2 + KBC) = 2*KEQ
         MAXPIVOT(IJ + KE2 + KBC) = KBC + KE2 + KE1
         ROWPIVOT(IJ + KE2 + KBC) = 1 + KBC
      END DO
! Subsequent meshpoints, except last
      DO IK = 2, KH-1
! Second order equations at meshpoint IK
         DO IJ = 1, KE2
            ROWSTART(IJ + KBC + (IK - 1)*KEQ) = (IK - 1)*KEQ - KEQ
            ROWWIDTH(IJ + KBC + (IK - 1)*KEQ) = 3*KEQ
            MAXPIVOT(IJ + KBC + (IK - 1)*KEQ) = KBC + (IK - 1)*KEQ + KE2
            ROWPIVOT(IJ + KBC + (IK - 1)*KEQ) = 1 + KEQ + KBC
         END DO
! First order equations at meshpoint IK
         DO IJ = 1, KE1
            !ROWSTART(IJ + KBC + KE2 + (IK - 1)*KEQ) = (IK - 1)*KEQ
            !ROWWIDTH(IJ + KBC + KE2 + (IK - 1)*KEQ) = 2*KEQ
            !ROWPIVOT(IJ + KBC + KE2 + (IK - 1)*KEQ) = 1 + KBC
            ROWSTART(IJ + KBC + KE2 + (IK - 1)*KEQ) = (IK - 1)*KEQ - KEQ
            ROWWIDTH(IJ + KBC + KE2 + (IK - 1)*KEQ) = 3*KEQ
            MAXPIVOT(IJ + KBC + KE2 + (IK - 1)*KEQ) = KBC + (IK-1)*KEQ+KE2+KE1
            ROWPIVOT(IJ + KBC + KE2 + (IK - 1)*KEQ) = 1 + KEQ + KBC
         END DO
      END DO
! Second order equations, last meshpoint (surface)
      DO IJ = 1, KE2
         ROWSTART(IJ + KBC + (KH - 1)*KEQ) = (KH - 1)*KEQ - KEQ
         ROWWIDTH(IJ + KBC + (KH - 1)*KEQ) = 2*KEQ
         ROWPIVOT(IJ + KBC + (KH - 1)*KEQ) = 1 + KBC
         !ROWSTART(IJ + KBC + (KH - 1)*KEQ) = (KH - 2)*KEQ - KEQ
         !ROWWIDTH(IJ + KBC + (KH - 1)*KEQ) = 3*KEQ
         MAXPIVOT(IJ + KBC + (KH - 1)*KEQ) = KBC + (KH - 1) * KEQ + KE2
         ROWPIVOT(IJ + KBC + (KH - 1)*KEQ) = 1 + KEQ + KBC
         !ROWPIVOT(IJ + KBC + (KH - 1)*KEQ) = 1 + KEQ + KBC
      END DO
! Surface boundary conditions and EV
! This block is special in the sense that we let it start at the same
! location as the below-surface meshpoint. This is so that the block will
! actually by physically adjacent to the EV column on the right. This way,
! we do not need to write special-case code to invert this last block.
! In effect, we make the SBC block "right aligned" within the stored
! sub-matrix. An alternative is to do this for all blocks except the CBC
! block, which logically aligns first and second order blocks. See how this
! works out first.
      DO IJ = 1, KE1 - KBC + KEV
         !ROWSTART(IJ + KBC + KE2 + (KH - 1)*KEQ) = (KH - 1)*KEQ
         !ROWWIDTH(IJ + KBC + KE2 + (KH - 1)*KEQ) = KEQ
         !ROWPIVOT(IJ + KBC + KE2 + (KH - 1)*KEQ) = 1 + KBC
         ROWSTART(IJ + KBC + KE2 + (KH - 1)*KEQ) = (KH - 1)*KEQ - 2*KEQ
         ROWWIDTH(IJ + KBC + KE2 + (KH - 1)*KEQ) = 3*KEQ
         MAXPIVOT(IJ + KBC + KE2 + (KH - 1)*KEQ) = NJAC
         ROWPIVOT(IJ + KBC + KE2 + (KH - 1)*KEQ) = 1 + CEV - KEQ
      END DO
! Column shapes: first and last non-zero row in each column, for the purpose
! of finding out which rows to loop over when doing Gaussian elimination.
! Because of column interchanges, these need to extend beyond the row where
! each column has data in the unscrambled matrix, because an interchange
! can put non-zero elements where they usually wouldn't be.
      COLSTART(1:KEQ) = 1
      COLHEIGHT(1:KEQ) = KBC + KE2 + KEQ + KEQ
      DO IJ = KEQ+1, NJAC-KEV
         COLSTART(IJ) = KBC + (IJ / KEQ - 1) * KEQ
         COLHEIGHT(IJ) = MIN(NJAC, COLSTART(IJ) + KE2 + 2 * KEQ + KEQ)
      END DO
      IF (KEV > 0) THEN
         COLSTART(NJAC-KEV:NJAC) = 1
         COLHEIGHT(NJAC-KEV:NJAC) = NJAC
      END IF
      do ij=1, 50
         !print *, ij, rowstart(ij), colstart(ij), colheight(ij)
!         print '(6I6)', ij, rowstart(ij), colstart(ij), rowstart(ij)
!     &    + rowwidth(ij), colheight(ij), max(rowstart(ij) + rowwidth(ij), colheight(ij))
      end do

C BEGIN ITERATIVE LOOP
      FAC = 1.0D0
      BESTERR = 1.0D0
      RECOMPUTE = .TRUE.
      DO JTER = 1, ITER
      ERR = 0.0D0
      JKH = 0
      IF ( JNN.EQ.JH1 .AND. JTER.EQ.JH2 ) JKH = 1
      XD(1:NM) = 0.0D0

! Initialise block Jacobian
      BLOCK_JACOBIAN(:,:) = 0.0D0
      DO IK=1, NJAC
         RP(IK) = IK
         RHP(IK) = IK
         CP(IK) = IK
      END DO

      call compute_functions(IKFST, IKLST, KQ)
      call compute_jacobian(IKFST, IKLST, KQ)
      !call regularise_jacobian

      ! Gauss-Jordan elimination of the block Jacobian
      call invert_block_jacobian

      if (.true.) then
      ! Solve the banded part of the Jacobian
      call dgbsvx('E', 'N',  NBAND, NSUB, NSUP, NRHS, 
     &               BANDED_JACOBIAN, LDA, AFB, LDAFB,
     &               PIVOT, EQUED, RSCALE, CSCALE,
     &               RHS, NJAC, X,  NBAND,
     &               RCOND,  FERR,  BERR,  WORK,  IWORK, INFO)
      rhs(1:NBAND,:) = x(1:NBAND,:)
      print *, 'info = ',  info
      print *, 'Condition number: ', 1/RCOND
      if (info /= 0) exit           ! Singular matrix
      ! Now invert the remaining part of the matrix
      ! First complete the Gauss-Jordan elimination of the block below the
      ! block diagonal. Since this entire block is below the main diagonal
      ! this means all the elements in this block miust be reduced to 0.
      ! Since we don't need to use this for anything else, we only perform
      ! the row operations on the columns in RHS.
      do ij=1, NSBC
         do ik = NBAND - 2*KEQ, NBAND
            RHS(NBAND + ij, 1:NRHS) = 
     &         RHS(NBAND + ij, 1:NRHS) - EV_BLOCK(ik, ij)*RHS(ik, 1:NRHS)
         end do
      end do

      ! The NSBC x NSBC block
      if (NSBC > 0) then
         EV_MATRIX = TRANSPOSE(RHS(NBAND+1:NJAC, 1:NSBC))
         B(1:NSBC) = RHS(NBAND+1:NJAC, NRHS)
         call dgesv ( NSBC, 1, EV_MATRIX, NSBC, PIVOT, B, NSBC, INFO)
         RHS(NBAND+1:NJAC, NRHS) = B(1:NSBC)
         print *, 'info = ',  info
         if (info /= 0) exit        ! Singular matrix

      ! We will need the inverse of this matrix for the back substitution
         RHS(NBAND+1:NJAC, 1:NSBC) = 0.0
         do ij = 1, NSBC
            rhs(NBAND+ij, ij) = 1.0
         end do
      end if
      ! Back substitution: for each EV, eliminate the columns above it.
      do ij = 1, NSBC
         do ik = 1, NJAC - ij
            RHS(ik, 1:NRHS) = RHS(ik, 1:NRHS) - RHS(ik, ij)*RHS(NBAND+ij,1:NRHS)
         end do
      end do
      ! Unscramble the RHS vector
      do ij = 1, NJAC
         ik = rp(rhp(ij))
         print *, RHS(ij, NRHS), BLOCK_JACOBIAN(ik, CRHS), 1.0d0 - dabs(RHS(ij, NRHS)/BLOCK_JACOBIAN(ik, CRHS))
         !RHS(ij, NRHS) = BLOCK_JACOBIAN(rp(ij), CRHS)
      end do

      stop
      ! Test the solution
      ij=1
      do ik = 1, KH
         do I = 1, KEQ
            C(kh+1-ik, I, 1) = RHS(ij, NRHS)
            ij = ij+1
         end do
      end do
      do I = 1, KEV
         C(1:KH, KEQ + I, 1) = RHS(KH*KEQ + I, NRHS)
      end do
      else
      ! Test the solution
      ij=1
      do ik = 1, KH
         do I = 1, KEQ
            C(kh+1-ik, I, 1) = BLOCK_JACOBIAN(rp(ij), CRHS)
            ij = ij+1
         end do
      end do
      do I = 1, KEV
         C(1:KH, KEQ + I, 1) = BLOCK_JACOBIAN(rp(KH*KEQ + I), CRHS)
      end do
      end if
      DO IJ = 1, KVB
C ESTIMATE ACCURACY OF ITERATION
         VX = 0.0D0
         DO IK = 1, KH
            VZ = DABS(C(IK, IJ, 1))
            IF ( VZ > VX ) VX = VZ
            IF ( VX.EQ.VZ ) IKM(IJ) = IK
            ERR = ERR + VZ
         END DO
         ERT(IJ) = VX
      END DO
      ERR = ERR/(KVB*KH)
      If (JTER > 3 .AND. ERR > ERRPR) THEN
         FAC = 0.8D0*FAC;
      ELSE
         FAC = MIN(DEL/DMAX1(ERR, DEL), 1.1D0*FAC)
      ENDIF
      FRR = DLOG10(ERR)
      ERT(1:KVB) = DLOG10(ERT(1:KVB) + 1.3D-10)
      IF ( JTER == KT5+1 )
     & WRITE (JB, '(A4,2X,A3,2X,A3,1X, 17(1X,A7),/,3(15X,17(1X,A7),/))') 
     &        'Iter', 'Err', 'Fac', (VARNAMES(KD(IJ)), IJ=1,KVB)
      IF ( JTER.GT.KT5 )
     : WRITE (JB, 992) JTER, FRR, FAC, (IKM(IJ), ERT(IJ), IJ = 1, KVB) 
  992 FORMAT (1X, I2, 2F6.2, 17(I4, F4.1),/, 3(15X, 17(I4, F4.1),/))
      CALL FLUSH ( JB )
C APPLY CORRECTIONS, SCALED DOWN IF TOO LARGE
      DO I = 1, KVB
         IJ = KD(I)
         DH(IJ, 1:KH) = DH(IJ, 1:KH) - MIN(C(1:KH, I, 1), CLIMIT)*FAC*ER(IJ)
      END DO
      CALL CHECKS
      IF ( ERR < BESTERR ) THEN
         BESTH(1:NVAR,1:KH) = H(1:NVAR,1:KH)
         BESTDH(1:NVAR,1:KH) = DH(1:NVAR,1:KH)
         BESTERR = ERR
      END IF
      ERRS = ERR
      IF ( ERR.GT.1.0D1*DEL ) JO = 1
      IF ( (ERR.LT.EPS .AND. FAC>0.1D0) .OR. JO.EQ.1 ) EXIT
      ERRPPR = ERRPR
      ERRPR = ERR
C CONTINUE ITERATING IF NOT YET ACCURATE ENOUGH
      END DO
! Free memory
      DEALLOCATE(BANDED_JACOBIAN)
      DEALLOCATE(RHS)
      DEALLOCATE(EV_BLOCK)
      DEALLOCATE(PIVOT)

      DEALLOCATE(AFB)
      DEALLOCATE(RSCALE)
      DEALLOCATE(CSCALE)
      DEALLOCATE(X)
      DEALLOCATE(FERR)
      DEALLOCATE(BERR)
      DEALLOCATE(WORK)
      DEALLOCATE(IWORK)
      IF (NSBC > 0) THEN
         DEALLOCATE(EV_MATRIX)
         DEALLOCATE(B)
      END IF

      DEALLOCATE(BLOCK_JACOBIAN)
      DEALLOCATE(ROWSTART)
      DEALLOCATE(ROWWIDTH)
      DEALLOCATE(COLSTART)
      DEALLOCATE(COLHEIGHT)
      DEALLOCATE(RP)
      DEALLOCATE(RHP)
      DEALLOCATE(CP)
      DEALLOCATE(MAXPIVOT)
      DEALLOCATE(ROWPIVOT)

      IF ( info == 0 .and. ERR<=EPS ) RETURN
      JO = 1
      RETURN
      END

      ! Compute FUNCS and derivatives
      SUBROUTINE COMPFUNCS ( K )
      USE MESH
      USE JACOBIAN
      IMPLICIT NONE
c     Parameter:
      INTEGER :: k

c     Local variables:
      INTEGER :: i, ivx, ji
      DOUBLE PRECISION :: dx, dvx, vx

c     Common variables:
      INTEGER :: kh, ktw, kw
      INTEGER :: jmod, jb, jnn, jter, joc, jkh
      INTEGER :: keq, kj2, kj5, kj6, kj10, kj12, ki4, kee, ia, ke1, ke2,
     &        ke3, kbc, kev, kfn, kl, jh1, jh2, jh3, kd, kq, kvb, kvc

      DOUBLE PRECISION :: h, dh, eps, del, dh0
      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: c, s, er
      DOUBLE PRECISION :: var, dvar, fn1, dfn1
      DOUBLE PRECISION :: fn2, dfn2, equ, dequ
      DOUBLE PRECISION :: fval, fjac, xd
      LOGICAL :: RECOMPUTE

      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :             KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC,  
     :       KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), FN1(NFUNC), DFN1(NVAR,NFUNC)
      COMMON /INE   / FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)
      COMMON /JAC   / FVAL(NM, NFUNC), FJAC(NM,NVAR,NFUNC), XD(NM), RECOMPUTE

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      ivx = 0
      ji = 0
      dx = 0
      dvx = 0
      vx = 0
c-----------------------------------------------------------------------

      DVAR(1:NVAR) = DH(1:NVAR, K - KL)
       VAR(1:NVAR) =  H(1:NVAR, K - KL) + DVAR(1:NVAR)
      IF ( JOC.EQ.1 ) THEN
         XD(1:KVC) = XD(1 + KVB:KVC + KVB)
C EVALUATE AND STORE THE REQUIRED FUNCS OF THE INDEP. VBLES
         CALL FUNCS1 ( K - KL, 0 )
         FVAL(K - KL, 1:KFN) = FN1(1:KFN)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The effect of the above code is the following: After the call to FUNCS1, the
!  array FN1() contains the values for independent variables and boundary cond.
! This is mapped onto the array FN2 for three subsequent meshpoints. The end
!  result is that the order of arrays in NAMEIN in EQUNS1 corresponds to the
!  order of variables in NAMEOUT in FUNCS1.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C VARYING INDEP. VBLES IN TURN, EVALUATE NUMERIC DERIVATIVES OF FUNCS
         DO I = 1, KVB
            JI = KD(I)
            VX = VAR(JI)
            DVX = DVAR(JI)
            IVX = 1
            IF ( DVX.LT.0.0D0 ) IVX = -1
            DX = IVX*DH0*MAX(ABS(VX), 1.0D0)  
            XD(I + KVC) = ER(JI)/DX
            VAR(JI) = VX + DX
            DVAR(JI) = DVX + DX
            CALL FUNCS1 ( K - KL, JI )
            FJAC(K - KL, I, 1:KFN) = FN1(1:KFN)
            DVAR(JI) = DVX
            VAR(JI) = VX
         END DO
      ELSE
C ALTERNATIVE TO BE USED IF DERIVS COMPUTED ANALYTICALLY
         CALL FUNCS2 ( K - KL, 1, KH )
         FVAL(K - KL, 1:KFN) = FN1(1:KFN)
         FJAC(K - KL, 1:KVB, 1:KFN) = DFN1(1:KVB, 1:KFN)
      END IF
      END SUBROUTINE COMPFUNCS

      ! JZ is the `start' index in the list of permutations
      SUBROUTINE DIFRNS ( K, JX, JY, JZ, J00 )
      USE MESH
      USE JACOBIAN
      IMPLICIT NONE

c     Parameters:
      INTEGER :: k, jx, jy, jz, j00

c     Local variables:
      INTEGER :: i, iee, ieq, ii, ivb, j, jeq, jj, jvb3, jvb4
      INTEGER :: i1, j1
      DOUBLE PRECISION :: ds

c     Common variables:
      INTEGER :: kh, ktw, kw
      INTEGER :: jmod, jb, jnn, jter, joc, jkh
      INTEGER :: keq, kj2, kj5, kj6, kj10, kj12, ki4, kee, ia, ke1, ke2,
     &        ke3, kbc, kev, kfn, kl, jh1, jh2, jh3, kd, kq, kvb, kvc

      DOUBLE PRECISION :: h, dh, eps, del, dh0
      DOUBLE PRECISION :: ml, ql, xl, uc
      DOUBLE PRECISION :: c, s, er
      DOUBLE PRECISION :: var, dvar, fn1, dfn1
      DOUBLE PRECISION :: fn2, dfn2, equ, dequ
      DOUBLE PRECISION :: fval, fjac, xd
      LOGICAL :: RECOMPUTE

      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, KJ5, 
     :             KJ6, KJ10, KJ12, KI4, KEE, IA(1), KE1, KE2, KE3, KBC,  
     :       KEV, KFN, KL, JH1, JH2, JH3, KD(120), KQ, KVB, KVC
      COMMON /INF   / VAR(NVAR), DVAR(NVAR), FN1(NFUNC), DFN1(NVAR,NFUNC)
      COMMON /INE   / FN2(3,NFUNC), DFN2(3,NVAR,NFUNC), EQU(NEQ), DEQU(NVAR,3,NEQ)
      COMMON /JAC   / FVAL(NM, NFUNC), FJAC(NM,NVAR,NFUNC), XD(NM), RECOMPUTE

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      iee = 0
      ieq = 0
      ii = 0
      ivb = 0
      j = 0
      jeq = 0
      jj = 0
      jvb3 = 0
      jvb4 = 0
      ds = 0
c-----------------------------------------------------------------------

      IF ( JKH.GE.1 .AND. ( K.LT.4.OR.K.GT.KH-2 )) THEN
         JKH = JH3 + 1
      ELSE
         JKH = 1
      END IF
C REDEFINE PREVIOUS CURRENT MESHPOINT AS CURRENT PREVIOUS MESHPOINT
      IF ( (KL.EQ.0.AND.K.LE.KH) .OR. (KL.EQ.1.AND.K.GE.2) ) THEN
          FN2(KEE:2, 1:KFN)        =  FN2(KEE+1:3, 1:KFN)
         DFN2(KEE:2, 1:KVB, 1:KFN) = DFN2(KEE+1:3, 1:KVB, 1:KFN)

C Call function to perform the following tasks:
C * evaluate argument list of vbles and increments at current meshpoint
C * evaluate and store the required funcs of the indep. vbles
C * varying indep. vbles in turn, evaluate numeric derivatives of funcs
C * Or: Store result of derivs computed analytically
         CALL COMPFUNCS( K )
         FN1(1:KFN) = FVAL(K - KL, 1:KFN)
         DFN1(1:KVB, 1:KFN) = FJAC(K - KL, 1:KVB, 1:KFN)
         FN2(3, 1:KFN) = FN1(1:KFN)
         DFN2(3, 1:KVB, 1:KFN) = DFN1(1:KVB, 1:KFN)
      END IF
      IF ( JOC.EQ.1 ) THEN
C EVALUATE AND STORE THE DIFFERENCE EQUNS WHICH ARE TO BE SATISFIED
         IF (J00 == 1) THEN
            CALL EQUNS1 ( K, KL, KQ )
            DO IEQ = JX, JY
               JEQ = KD(IEQ + JZ)
               S(IEQ, KJ12) = EQU(JEQ)
               RHS(KBC + (KH-K)*KEQ + IEQ, NRHS) = EQU(JEQ)
               BLOCK_JACOBIAN(KBC + (KH-K)*KEQ + IEQ, CRHS) = EQU(JEQ)
            END DO
            RETURN
         END IF
C VARYING INDEP. VBLES IN TURN, EVALUATE NUMERIC DERIVATIVES OF EQUNS
         JVB3 = KEE*KVB - KVB
         DO IEE = KEE, 3
            DO IVB = 1, KVB
               JVB3 = JVB3 + 1
               FN1(1:KFN) = FN2(IEE, 1:KFN)
               FN2(IEE, 1:KFN) = DFN2(IEE, IVB, 1:KFN)
               CALL EQUNS1 ( K, KL, KQ )
               FN2(IEE, 1:KFN) = FN1(1:KFN)
               DO IEQ = JX, JY
                  JEQ = KD(IEQ + JZ)
                  DS = (EQU(JEQ) - RHS(KBC + (KH-K)*KEQ + IEQ, NRHS))*XD(JVB3)
                  ! Indices for banded Jacobian
                  j1 = KBC + (KH - K)*KEQ + IEQ          ! Equation (row)
                  if (IVB <= KEQ) then       ! Not an eigenvalue, in band
                     i1 = KH*KEQ - (K + (2 - IEE))*KEQ + IVB    ! Variable (col)
                     ! Surface boundary condition block needs to be shifted
                     ! one block to the left, so that it is below the
                     ! block for the next-to-surface meshpoint.
                     if (k == 1) i1 = i1 - KEQ
                  else                       ! Eigenvalue, outside band
                     i1 = KH*KEQ + IVB - KEQ
                  end if
                  call set_jacobian(j1, i1, DS)
               END DO
            END DO
         END DO
      ELSE
C ALTERNATIVE TO BE USED IF DERIVS COMPUTED ANALYTICALLY
         CALL EQUNS2 ( K, 1, KH, KEQ )
         DO J = JX, JY
            JJ = KD(J + JZ)
            DO I = 1, KEQ
               II = KD(I)
               S(J, I) = DEQU(I, 1, JJ)*ER(II)
               S(J, I + KEQ) = DEQU(I, 2, JJ)*ER(II)
               S(J, I + KJ5) = DEQU(I, 3, JJ)*ER(II)
            END DO
            IF ( KEV.NE.0 ) THEN
               DO I = KEQ+1, KVB
                  S(J, I + KJ5) = (DEQU(I, 1, JJ) + DEQU(I, 2, JJ) 
     :                         + DEQU(I, 3, JJ)) * ER(KD(I))
               END DO
            END IF
            S(J, KJ12) = EQU(JJ)
         END DO
      END IF
      RETURN
      END

      SUBROUTINE DIVIDE ( ID1, JD1, ID2, K, JO )
      USE MESH
      use compare_floats;
      IMPLICIT NONE

c     Parameters:
      INTEGER :: id1, jd1, id2, k, jo

c     Local variables:
      INTEGER :: i, ii, im, j, jc, jd2, jd3, jj, jl, jmi
      DOUBLE PRECISION :: vm, vs, vt, vx

c     Common variables:
      INTEGER :: keq, jw, kj10, kj12, ki4, kee, jx
      INTEGER :: itype, jmax, idone
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: c, s, er
      DOUBLE PRECISION :: ml, ql, xl, uc

      COMMON /SOLV  / C(NM+1, NVAR, NVAR+1), S(NEQ, 121), ER(NVAR), KEQ, JW(3), 
     :                KJ10, KJ12, KI4, KEE, JX(134) 
      COMMON /DIV   / ITYPE(NVAR), JMAX(NVAR+1), IDONE(NVAR)
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      ii = 0
      im = 0
      j = 0
      jc = 0
      jd2 = 0
      jd3 = 0
      jj = 0
      jl = 0
      jmi = 0
      vm = 0
      vs = 0
      vt = 0
      vx = 0
c-----------------------------------------------------------------------

      JC = 1
      JD2 = JD1 + ID2 - ID1
      IF ( ID1.GT.ID2 .OR. JD2 + 1.GT.KJ12 ) RETURN
      !IF ( JKH.GE.2 ) CALL PRINTS ( K, JKH, ID1, ID2, 0 )
      IDONE(ID1:ID2) = 1
      ITYPE(ID1:ID2) = 5
   2  JC = JC + 1
      DO I = ID1, ID2
         IF ( ITYPE(I).LT.3 ) ITYPE(I) = ITYPE(I) + 2
         IF ( ITYPE(I).LT.5 ) IDONE(JMAX(I)) = 0
      END DO
      VT = 0.0
c Locate the most significant (remaining) row. `Significance' of a row is 
c the ratio of largest |element| = VM to sum of remaining |elements| = VS 
      DO I = ID1, ID2
         IF ( ITYPE(I).GE.5 ) THEN
            VM = 0.0
            DO J = JD1, JD2
               JJ = J - JD1 + ID1
               IF ( IDONE(JJ).NE.0 ) THEN
                  VX = DABS(S(I, J))
                  IF ( VX.GE.VM ) THEN
                     VM = VX
                     JL = JJ
                  END IF
               END IF
               VS = 0.0D0
            END DO
            IF ( JL.LT.1 .OR. JL.GT.NVAR+1 ) GO TO 8
            JMAX(I) = JL
            DO J = JD1, JD2
               IF ( J - JD1 + ID1.NE.JL ) VS = VS + DABS(S(I, J))
            END DO
            IF ( FEQ(VS, 0.0D0) ) ITYPE(I) = 1
            IF ( FEQ(VM, 0.0D0) ) ITYPE(I) = 6
            IF ( FEQ(VS, 0.0D0) .AND. NFEQ(VM, 0.0D0) ) IDONE(JL) = 0
            IF ( FEQ(VM, 0.0D0) ) VX = 0.0D0
            IF ( FEQ(VS, 0.0D0) .AND. VM.GT.0.0D0 ) VX = 2.0D0
            IF ( VS.GT.0.0D0 ) VX = VM/(VM + VS)
            IF ( VX.GE.VT ) THEN
               VT = VX
               IM = I
            END IF
         END IF
      END DO
      IF ( IM.LT.1 .OR. IM.GT.NVAR ) GO TO 8
      IF ( ITYPE(IM).EQ.5 ) ITYPE(IM) = 2
c Largest element moduluswise of most significant row is the leading 
c pivot; eliminate elements above and below it
      DO I = ID1, ID2
         IM = ITYPE(I)
         IF ( IM.LT.3 ) THEN
            JMI = JMAX(I) + JD1 - ID1
            JD3 = JD1
            IF ( IM.EQ.1 ) JD3 = JD2 + 1
            VX = 1.0D0/S(I, JMI)
            S(I, JD3:KJ12) = VX*S(I, JD3:KJ12)
            S(I, JMI) = 1.0D0
            DO II = ID1, ID2
               IF ( ITYPE(II).GT.3 ) THEN
                  VX = S(II, JMI)
                  S(II, JD3:KJ12) = S(II, JD3:KJ12) - VX*S(I, JD3:KJ12)
                  S(II, JMI) = 0.0D0
               END IF
            END DO
         END IF
         IDONE(I) = 1
      END DO
      DO I = ID1, ID2
         IF ( (ITYPE(I).EQ.5 .OR. ITYPE(I).LE.2) .AND. JC.LT.5000 ) 
     :      GO TO 2
      END DO
      DO I = ID1, ID2
         DO J = JD2 + 1, KJ12
            C(K, JMAX(I), J - KJ12 + KI4) = S(I, J)
         END DO
      END DO
      IF ( JKH.GE.2 ) GO TO 9
      IF ( JC.LT.5000 ) RETURN
   8  JO = 1
c......../........./........./........./........./........./........./..
c Some emergency printout
   9  WRITE (10, 991) K, JNN, JTER, JKH, JL, IM, JC, ITYPE, JMAX, IDONE
  991 FORMAT (7I5, /, 40I3, /, 41I3, /, 40I3)
      !CALL PRINTC ( K, JKH, ID1, ID2, JD2 + 1 - KJ12 + KI4, KI4 )
      RETURN
      END

      SUBROUTINE ELIMN8 ( IS1, JS1, IS2, JS2, IS3, JS4, JS5, K, JCL )
      USE MESH
      IMPLICIT NONE

c     Parameters:
      INTEGER :: is1, js1, is2, js2, is3, js4, js5, k, jcl

c     Local variables:
      INTEGER :: i, jj, jjj
      DOUBLE PRECISION :: vx

c     Common variables:
      INTEGER :: keq, jw, kj10, kj12, ki4, kee, jx
      INTEGER :: jmod, jb, jnn, jter, joc, jkh

      DOUBLE PRECISION :: c, s, er
      DOUBLE PRECISION :: ml, ql, xl, uc

      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, JW(3), 
     :                KJ10, KJ12, KI4, KEE, JX(134) 
      COMMON /QUERY / ML, QL, XL, UC(21), JMOD, JB, JNN, JTER, JOC, JKH

c-----------------------------------------------------------------------
c     Explicitly initialize all local variables (Steve, 5/08).

      i = 0
      jj = 0
      jjj = 0
      vx = 0
c-----------------------------------------------------------------------

      !IF ( JKH.GE.2 ) CALL PRINTS ( K, JKH, IS1, IS2, JCL )
      IF ( JS1.GT.JS2 .OR. 1.GT.IS2 ) RETURN
      DO JJ = JS1, JS2
         JJJ = JJ - JS1 + IS3
         DO I = IS1, IS2
            VX = S(I, JJ)
            IF ( JS4 <= JS5 ) THEN
               S(I, JS4:JS5) = S(I, JS4:JS5) - VX*C(K, JJJ, 1:JS5 - JS4 + 1)
            END IF

            S(I, KJ10:KJ12) = S(I, KJ10:KJ12) - VX*C(K, JJJ, KJ10 - KJ12 + KI4:KI4)
         END DO
      END DO
      RETURN
      END

      SUBROUTINE PRINTS ( JK, JKH, IPS1, IPS2, JCL )
      USE MESH
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ, KJ2, 
     :                KJ5, KJ6, KJ10, KJ12, KI4, KEE, JX(134)
      WRITE (10,*) 'printS', JCL, JK, JKH, KEE, KEQ, KJ2, KJ5, KJ6, KJ12
      IF ( KEE.NE.2 ) THEN
         DO I = IPS1, IPS2
            WRITE (10, 991) JK, (S(I, J), J = 1, KEQ)
         END DO
         WRITE (10, 991) JK
      END IF
      DO I = IPS1, IPS2
         WRITE (10, 991) JK, (S(I, J), J = KJ2, KJ5)
      END DO
      WRITE (10, 991) JK
      DO I = IPS1, IPS2
         WRITE (10, 991) JK, (S(I, J), J = KJ6, KJ12)
      END DO
      WRITE (10, 991) JK
  991 FORMAT (I5, 1P, 17D9.1, /, 17D9.1)
      RETURN
      END

      SUBROUTINE PRINTC ( JK, JKH, IPC1, IPC2, JPC1, JPC2 )
      USE MESH
      IMPLICIT REAL*8 (A-H, L-Z)
      COMMON H(NVAR,NM), DH(NVAR,NM), EPS, DEL, DH0, KH, KTW, KW(260)
      COMMON /SOLV  / C(NM+1,NEQ,NVAR+1), S(NEQ,121), ER(NVAR), KEQ(142)
      IF ( JKH.LT.2 .OR. (JK.GE.4 .AND. JK.LE.KH - 2) ) RETURN
      WRITE (10,*) 'printC', JK, JKH, JPC1, JPC2
      DO I = IPC1, IPC2
         WRITE (10, 99002) (C(JK, I, J), J = JPC1, JPC2)
      END DO
99002 FORMAT (1P, 15D10.2)
      RETURN
      END

