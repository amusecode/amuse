module create_new_model
   ! Create a new particle from a user supplied model (non-ZAMS, e.g. merger product)
   
   use real_kind
   
   ! variables for storing the user supplied model
   implicit none
   private
   
   real(double), pointer :: new_model(:,:)
   logical :: new_model_defined = .false.
   integer :: n_shells_new_model
   
   public :: new_stellar_model
   public :: finalize_stellar_model
   
   contains
   
   ! Each shell has: Mass coordinate [Msun], Radius [Rsun], density [cgs],
   !  pressure [cgs], XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE
   ! Mespoints should be passed in from *surface* to *centre*
   function new_stellar_model(mass, radius, rho, pressure, &
         XH, XHE, XC, XN, XO, XNE, XMG, XSI, XFE, n)
      implicit none
      integer :: new_stellar_model
      integer, intent(in) :: n
      double precision, intent(in) :: mass(n), radius(n), rho(n), pressure(n), &
         XH(n), XHE(n), XC(n), XN(n), XO(n), XNE(n), XMG(n), XSI(n), XFE(n)
      double precision :: logrho(n), logpressure(n)
      
      if (new_model_defined) then
         new_stellar_model = -30
         return
      endif
      
      new_stellar_model = new_empty_model(n)
      if (new_stellar_model /= 0) then
         return
      endif
      
      new_model(1, 1:n) = mass(1:n)
      new_model(2, 1:n) = radius(1:n)
      logrho(1:n) = log(rho(1:n))
      logpressure(1:n) = log(pressure(1:n))
      new_model(3, 1:n) = rho(1:n)
      new_model(4, 1:n) = pressure(1:n)
      new_model(5, 1:n) = XH
      new_model(6, 1:n) = XHE
      new_model(7, 1:n) = XC
      new_model(8, 1:n) = XN
      new_model(9, 1:n) = XO
      new_model(10, 1:n) = XNE
      new_model(11, 1:n) = XMG
      new_model(12, 1:n) = XSI
      new_model(13, 1:n) = XFE
      new_model_defined = .true.
      n_shells_new_model = n
      
   end function
   
   ! Start with an empty model with the correct size
   function new_empty_model(n_shell)
      use mesh, only: NM
      implicit none
      integer :: new_empty_model
      integer, intent(in) :: n_shell
      
      select case (n_shell)
         case ( : 2)
            new_empty_model = -31   ! Too few shells in new empty model
         case (NM+1 : )
            new_empty_model = -32   ! Too many shells in new empty model
         case default
            allocate(new_model(13, n_shell))
            new_empty_model = 0    ! A new empty model was defined
      end select
   end function
   
   function finalize_stellar_model(star_id, age_tag)
      use import
      implicit none
      integer :: finalize_stellar_model
      integer, intent(out) :: star_id
      double precision, intent(in) :: age_tag
      
      if (.not. new_model_defined) then
         finalize_stellar_model = -35
         return
      endif
      star_id = import_stellar_merger(n_shells_new_model, 13, new_model, age_tag)
      call flush()
      
      if (star_id .eq. 0) then
         finalize_stellar_model = -36
      else
         deallocate(new_model)
         new_model_defined = .false.
         finalize_stellar_model = 0
      end if
      
   end function

end module create_new_model
