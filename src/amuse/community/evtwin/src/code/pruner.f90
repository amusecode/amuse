subroutine pruner ( jn, jt, isb )
   use real_kind
   use ltm2ubv
   ! Reads in a long output file (eg fort.1 or 2) and prunes it down
   ! to 5 lines per timestep
   
   implicit none
   integer :: jn, jt, isb
   integer :: i,j,k
   
   integer, allocatable :: im(:), ib(:)
   real(double), allocatable :: y(:), z(:,:)
   character :: char(150)

   allocate(im(20000))
   allocate(ib(20000))
   allocate(y(10000))
   allocate(z(20000,85))
   
98 format (3(1p, e16.9, 5e10.3, 0p, f8.3, 7f8.5, 2f8.4, f9.5, /), &
        (1p, e16.9, 5e10.3, 0p, f8.3, 7f8.3, 2f8.4, f9.5, /), &
        1p, e16.9, 5e10.3, 0p, 8f8.3, 24x, i6, i2)
   k=1
   do
      read (jn, '(150a1)', end=7) char
      if ( char(8).ne.'M' ) cycle
      read (jn, '(150a1)') char
      read (jn, 98, end=7) (z(k, j), j = 1, 82), im(k), ib(k)
      k = k + 1
   end do
   
   !  1-M*;    2-Porb; 3-xi ;  4-tn;   5-LH;    6-P(cr); 7-McHe;  8-CXH;   9-CXHe;
   !  10-CXC; 11-CXN; 12-CXO; 13-CXNe; 14-CXMg; 15-Cpsi; 16-Crho; 17-CT; 
   !  
   ! 18-dty;  19-Prot; 20-zet; 21-tKh; 22-LHe;  23-RCZ;  24-McCO; 25-TXH; 26-TXHe;
   ! 27-TXC;  28-TXN;  29-TXO; 30-TXNe; 31-TXMg; 32-Tpsi; 33-Trho; 34-TT; 
   !
   ! 35-age;  36-ecc;  37-mdt; 38-tET; 39-LCO;  40-DRCZ; 41-McNe; 42-SXH; 43-SXHe;
   ! 44-SXC;  45-SXN;  46-SXO; 47-SXNe; 48-SXMg; 49-Spsi; 50-Srho; 51-ST;
   ! 
   ! 52-cM*;  53-RLF1; 54-RLF2; 55-DLT;  56-Lnu; 57-RA/R; 58-MH; 
   ! 59 to 66-conv. bdries;     67-logR; 68-logL
   !
   ! 69-Horb; 70-F1; 71-DF21; 72-BE;  73-Lth;  74-Bp;  75-MHe;  76 to 81-semiconv.
   ! bdries;  82-k**2; 83-MV; 84- B-V; 85- U-B; then JMOD
7  k = k - 1
   y(1:k) = z(1:k, 35)
   
   ! age = Y(J). Look for *small* backward jumps in age: due to convergence 
   ! failure and backtracking with timestep reduced.
   do j = 1, k - 1
      if ( dabs(y(j)).lt.1.0d-8 ) y(j) = -1.0d0
      if ( dabs(y(j + 1)).lt.1.0d-8 ) y(j + 1) = -1.0d0
      if ( y(j).gt.y(j + 1) ) im(j) = -1
   end do
   
99 format (3(1p, e16.9, 5e10.3, 0p, f8.3, 7f8.5, 2f8.4, f9.5, /), &
        (1p, e16.9, 5e10.3, 0p, f8.3, 7f8.3, 2f8.4, f9.5, /), &
        1p, e16.9, 5e10.3, 0p, 11f8.3, i6, i2)
   
   do j = 1, k
      if ( j /= k ) then
         if ( y(j).lt.0.0d0 .and. y(j + 1).lt.0.0d0) cycle
      end if
      ! add In Bol. Corrn, colours./........./........./......../........./..
      call lt2ubv (z(j,68), z(j,51), z(j,1), z(j,83), z(j,84), z(j,85))
      if ( im(j).gt.0 .or. y(j).lt.0.0d0 ) write (jt, 99) &
           (z(j, i), i = 1, 34), y(j), (z(j, i), i = 36, 85), im(j), ib(k)
   end do
   
   if ( .not. ((jn.eq.1 .and. isb.eq.2) .or. jt.ne.9 ) ) then
      z(1, 1:85) = 0.0d0
      write (9,99) (z(1, i), i = 1, 85), 0, isb
   end if

   deallocate(im)
   deallocate(ib)
   deallocate(y)
   deallocate(z)
   
end subroutine pruner

