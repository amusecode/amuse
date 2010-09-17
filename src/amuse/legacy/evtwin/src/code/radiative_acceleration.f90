! FORTRAN 90 module to calculate radiative accelerations.
! Original code based on the Opacity Project (OP) code "ax.f"
! It should not be hard to replace the backend with a different
! implementation, for instance the "single value approximation" by
! Alecian & LeBlanc.
! The main performance hit is due to loading of the input files,
! which happens on-the-fly. Loaded files are cached to relieve the
! pain as much as possible (the cache size should perhaps be increased
! from that used here for practical use).
! There are some experimental OpenMP directives in the code, but all they
! seem to do right now is increase the CPU load of the code; the runtime
! performance is the same.
!
! Evert Glebbeek 2010
module radiative_acceleration
   use real_kind
   use constants

   character, private :: op_dir*132                   ! Should point to the OPCD mono/ directory. Initialise once.
   integer, private :: op_dir_length                  ! Length of the character string. Used in the original OP code.

   ! Data type for caching an entry that was read from disk so that we can reuse it
   ! We store separate entries for all elements.
   type, private :: radacc_cache_t
      sequence
      integer :: id, it, ine, iz
      real(single), pointer :: ff(:)         ! Typically 10000 elements
      real(single), pointer :: ta(:)         ! Typically 10000 elements
   end type radacc_cache_t
   type(radacc_cache_t), allocatable, private :: cache(:)
   integer, parameter :: cache_size = 17*256
   integer :: last_cache_entry

   ! Data for the frequency mesh, for the integration of the cross sections
   real(single), allocatable, private :: umesh(:)
   real(single), private :: dv
   integer, private :: ntot

   ! List of nuclear charges Z for which input files are available
   integer, parameter, private :: num_acc_elem = 17
   integer, parameter, private :: acceleration_z(num_acc_elem) = (/ 1, 2, 6, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20, 24, 25, 26, 28 /)
   character, parameter, private :: zlabel(num_acc_elem)*2 = (/'01','02','06','07','08','10','11','12','13',&
                                                               '14','16','18','20','24','25','26','28'/)

   ! Temperature "indices", for temperatures in the range logT=3.5..logT=8.0, in steps of 2
   integer, parameter, private :: min_temperature_index = 140
   integer, parameter, private :: max_temperature_index = 320

   ! Range of "density indices" for a given temperature, for all isotopes
   integer, parameter, private :: min_density_index(140:320) = (/       &
      14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 18, 19, 22, 23,   &
      26, 27, 30, 31, 34, 34, 34, 34, 36, 36, 36, 36, 36, 36, 38, 38,   &
      38, 38, 38, 39, 40, 40, 40, 40, 40, 41, 42, 42, 42, 42, 42, 43,   &
      44, 44, 44, 44, 44, 45, 46, 46, 46, 46, 46, 47, 48, 48, 48, 48,   &
      48, 49, 50, 50, 50, 50, 52, 52, 52, 52, 52, 52, 54, 54, 54, 54,   &
      54, 54, 56, 56, 56, 56, 56, 56, 58, 58, 58, 58, 58, 58, 60, 60,   &
      60, 60, 60, 60, 62, 62, 62, 62, 62, 63, 64, 64, 64, 64, 64, 65,   &
      66, 66, 66, 66, 66, 67, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70,   &
      70, 71, 72, 72, 72, 72, 72, 73, 74, 74, 74, 74, 76, 76, 76, 76,   &
      76, 76, 78, 78, 78, 78, 78, 78, 80, 80, 80, 80, 80, 80, 82, 82,   &
      82, 82, 82, 82, 84, 84, 84, 84, 84, 84, 86, 86, 86, 86, 86, 87,   &
      88, 88, 88, 88, 88 /)

   integer, parameter, private :: max_density_index(140:320) = (/       &
      52, 56, 56, 58, 58, 60, 60, 60, 60, 62, 62, 64, 64, 66, 66, 68,   &
      68, 70, 70, 72, 72, 74, 74, 74, 74, 76, 76, 78, 78, 80, 80, 80,   &
      80, 81, 82, 82, 82, 82, 82, 83, 84, 84, 84, 84, 84, 86, 86, 86,   &
      86, 86, 86, 88, 88, 88, 88, 89, 90, 90, 90, 88, 88, 88, 88, 92,   &
      92, 92, 92, 94, 94, 94, 94, 94, 94, 89, 90, 90, 90, 90, 90, 91,   &
      92, 92, 92, 92, 92, 94, 94, 90, 90, 92, 92, 92, 92, 92, 92, 93,   &
      94, 94, 94, 94, 94, 94, 94, 96, 96, 96, 96, 96, 96, 98, 98, 98,   &
      98, 98, 98, 99,100,100,100,100,100,100,100,102,102,102,102,102,   &
     102,103,104,104,104,104,104,104,104,106,106,106,106,106,106,107,   &
     108,108,108,108,108,108,108,110,110,110,110,110,110,111,112,112,   &
     112,112,112,112,112,114,114,114,114,114,114,115,116,116,116,116,   &
     116,116,116,118,118 /)

   ! Temperature labels for use in constructing file names
   character, parameter, private :: tlabel(140:320)*3 = (/                                                                 &
         '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155',   &
         '156', '157', '158', '159', '160', '161', '162', '163', '164', '165', '166', '167', '168', '169', '170', '171',   &
         '172', '173', '174', '175', '176', '177', '178', '179', '180', '181', '182', '183', '184', '185', '186', '187',   &
         '188', '189', '190', '191', '192', '193', '194', '195', '196', '197', '198', '199', '200', '201', '202', '203',   &
         '204', '205', '206', '207', '208', '209', '210', '211', '212', '213', '214', '215', '216', '217', '218', '219',   &
         '220', '221', '222', '223', '224', '225', '226', '227', '228', '229', '230', '231', '232', '233', '234', '235',   &
         '236', '237', '238', '239', '240', '241', '242', '243', '244', '245', '246', '247', '248', '249', '250', '251',   &
         '252', '253', '254', '255', '256', '257', '258', '259', '260', '261', '262', '263', '264', '265', '266', '267',   &
         '268', '269', '270', '271', '272', '273', '274', '275', '276', '277', '278', '279', '280', '281', '282', '283',   &
         '284', '285', '286', '287', '288', '289', '290', '291', '292', '293', '294', '295', '296', '297', '298', '299',   &
         '300', '301', '302', '303', '304', '305', '306', '307', '308', '309', '310', '311', '312', '313', '314', '315',   &
         '316', '317', '318', '319', '320' /)

   ! Actual range for different isotopes (temperature index is packed, (320-140)/2+1=91)
   integer, private :: density_start_index(num_acc_elem, 91)
   integer, private :: density_end_index(num_acc_elem, 91)

   ! Data for calculating partial ionisation
   real(single), allocatable :: bound_electrons(:,:,:)

   ! Private functions and subroutines
   private :: find_electron_density
   private :: get_bound_electrons
   private :: read_crosssections
   private :: geteta
   private :: interpolate
   private :: fint
   private :: solve_cubic

contains

! ------------------------------------------------------------------------------
! get_free_file_unit:
!  Helper function to get an available file handle for the next file we want
!  to open.
!  Beats having to remember which file numbers are in use by any part of the
!  code. Not to mention that this makes the code more easily portable.
! ------------------------------------------------------------------------------
!  Return value:
!   A free file handle
! ------------------------------------------------------------------------------
   function get_free_file_unit()
   implicit none
   integer :: get_free_file_unit
   integer :: n
   logical :: in_use
   n = 1
   do
      inquire(unit = n, opened = in_use)
      if (.not. in_use) exit
      n = n+1
   end do
   get_free_file_unit = n
   end function get_free_file_unit



! ------------------------------------------------------------------------------
! SET_OP_DIRECTORY:
! Set the location of the OP 'mono/' directory.
! How to figure this out is another problem that we let the caller worry about. Possibilities include:
!  - /usr/local/OPCD/
!  - ~/OPCD/
!  - $(OPCD)
! ------------------------------------------------------------------------------
   subroutine set_op_directory(dir)
   implicit none
   character(len=*), intent(in) :: dir

   op_dir_length = len(trim(adjustl(dir)))
   op_dir(1:op_dir_length) = trim(adjustl(dir))
   end subroutine set_op_directory



! ------------------------------------------------------------------------------
! INIT_RADIATIVE_ACCELERATION:
! Initialise global data for the calculation of radiative accelerations.
! Call once at startup.
! ------------------------------------------------------------------------------
   function init_radiative_acceleration()
   implicit none
   integer :: n, k, i, j, jn1, jn2, ite
   integer :: stat
   real(single) :: epa        ! bound electrons per atom
   logical :: init_radiative_acceleration

   init_radiative_acceleration = .false.

! Determine the frequence mesh for calculating the integrals over the cross sections
   n = get_free_file_unit()
   open(unit=n, file=op_dir(1:op_dir_length)//'m01.mesh', form='unformatted', status='old', iostat=stat)
   if (stat /= 0) return      ! Failed
   read(n) dv, ntot
   allocate(umesh(1:ntot))
   ! FORTRAN file I/O is retarded in that you cannot seek within a file unless you open it for "direct" access, in which case
   ! all "records" must be the same length and you cannot read more than one at a time ("sequential access"). *sigh*
   ! So we have to seek back to the beginning of the file and re-read the data there (which we just read).
   ! Am sorely tempted to write a C wrapper for this...
   rewind(n)
   read(n) dv, ntot, (umesh(k), k = 1,ntot)
   close(n)

! Read in data for calculating the electron density (ie, the ionisation states for each element)
   allocate(bound_electrons(num_acc_elem,                                              &
                               maxval(max_density_index - min_density_index)/2+1,      &
                               (max_temperature_index-min_temperature_index)/2+1))

   do k = 1, num_acc_elem
      n = get_free_file_unit()
      open(unit = n, file = op_dir(1:op_dir_length)//'m'//zlabel(k)//'.smry', status='old', iostat=stat)
      if (stat /= 0) return      ! Failed
      ! Skip first two lines
      read(n, *)
      read(n, *)
      do i = min_temperature_index, max_temperature_index, 2   ! Loop over temperature
         read(n, *) ite, jn1, jn2                              ! Read temperature index (redundant) and density range for element
         density_start_index(k, (i-min_temperature_index)/2+1) = jn1
         density_end_index(k, (i-min_temperature_index)/2+1) = jn2
         do j = jn1, jn2, 2
            read(n, *) ite, epa                                ! Partial ionisation
            bound_electrons(k, (j - min_density_index(i))/2 + 1, (i-min_temperature_index)/2+1) = epa
         end do
      end do
      close(n)
   end do

   ! Initialise the cache
   if (cache_size > 0 .and. .not. allocated(cache)) then
      allocate(cache(1:cache_size))
      last_cache_entry = 0
      do n=1, cache_size
         allocate(cache(n)%ff(1:ntot))
         allocate(cache(n)%ta(1:ntot))
      end do
   end if

   init_radiative_acceleration = .true.
   end function init_radiative_acceleration



   subroutine store_cache(it, ine, iz, ff, ta)
   implicit none
   integer, intent(in) :: it, ine, iz
   real(single), intent(in) :: ff(ntot), ta(ntot)
   type(radacc_cache_t) :: tmp
   integer :: id
   integer :: n

   if (cache_size == 0) return

!$omp critical
   ! Construct a unique ID, which is easier to compare with
   id = iz + 100*ine + 100000*it

   ! Are we already in the list?
   do n = 1, last_cache_entry
      if (cache(n)%id == id) then   ! Yes, swap with first item
         if (n == 1) goto 1         ! YUCK!!!! - replace the goto with return?   - No, because we need to close the omp critical, EG
         tmp = cache(1)
         cache(1) = cache(n)
         cache(n) = tmp
         goto 1                     ! YUCK!!!! - replace the goto with return?   - No, because we need to close the omp critical, EG
      end if
   end do

   ! Not yet in the list
   if (last_cache_entry < cache_size) last_cache_entry = last_cache_entry+1
   tmp = cache(last_cache_entry)
   do n = last_cache_entry, 2, -1
      cache(n) = cache(n-1)
   end do
   cache(1) = tmp

   cache(1)%id = id
   cache(1)%it = it
   cache(1)%ine = ine
   cache(1)%iz = iz
   cache(1)%ff(1:ntot) = ff(1:ntot)
   cache(1)%ta(1:ntot) = ta(1:ntot)

1 continue                       ! YUCK!!!!
!$omp end critical
   end subroutine store_cache



   function retrieve_cache(it, ine, iz, ff, ta)
   implicit none
   integer, intent(in) :: it, ine, iz
   real(single), intent(out) :: ff(ntot), ta(ntot)
   logical :: retrieve_cache
   integer :: id
   integer :: n

   retrieve_cache = .false.
   if (cache_size == 0) return
!$omp critical

   ! Construct a unique ID, which is easier to compare with
   id = iz + 100*ine + 100000*it

   ! Search through the list. Ultimately, we should use a hash table for this
   do n = 1, last_cache_entry
      if (cache(n)%id == id) then
         ff = cache(n)%ff(1:ntot)
         ta = cache(n)%ta(1:ntot)
         retrieve_cache = .true.
         exit
      end if
   end do
!$omp end critical
   end function retrieve_cache



! ------------------------------------------------------------------------------
! GET_BOUND_ELECTRONS
! Helper function to get the number of free electrons for given temperature and
! density indices and a particular particle species.
! ------------------------------------------------------------------------------
! Input:
!  it    - Temperature index (140..320)
!  jn    - Density index (range depends on it)
!  n     - Particle species (1..17)
! Returns:
!  The degree of ionisation
! ------------------------------------------------------------------------------
   pure function get_bound_electrons(it, jn, n)
   implicit none
   integer, intent(in) :: it, jn, n
   real(single) :: get_bound_electrons
   get_bound_electrons = bound_electrons(n, (jn - min_density_index(it))/2 + 1, (it-min_temperature_index)/2+1)
   end function



! ------------------------------------------------------------------------------
! GET_RADIATIVE_ACCELERATIONS:
! Calculate the radiative accelerations for all elements for given
! temperature and density, up to a factor 1/A (1/atomic weight in AMU).
! This factor is different for different isotopes of the same species, whereas
! the accelerations calculated here otherwise don't depend on the particular
! isotope considered (ignoring fine structure).
!
! An alternative approach is to pass in isotopic abundances and figure out
! elemental abundances from those, but that is a bit annoying to do and
! probably easier for the caller (with knowledge of which isotope belongs
! to which element).
!
! Does some file I/O to read the required input files. There should *maybe* be
! an option to hold all of this data in memory, but the full data is rather
! enormous, making this impractical in general. Assuming there is enough
! free memory to do so to begin with, we can perhaps rely on the OS caching
! the files in memory.
!
! NB: the OP code uses SINGLE PRECISION. Our function interface uses DOUBLE
! PRECISION.
!
! FIXME: we need to return the "status", for instance, "temperature out of range" to the caller.
! ------------------------------------------------------------------------------
!  Input:
!   logT       - log10 T (in K)
!   logrho     - log10 rho (in g/cm^3)
!   mub        - Mean molecular weight per baryon [AMU]
!   F          - Local radiative flux (Lrad/r^2, erg/cm^2/s)
!   nelem      - Number of elements
!   na(nelem)  - Abundance fraction (by number) for all elements
!   z(nelem)   - Nuclear charge for all isotopes (in units of e)
!  Output:
!   grad(nelem)- radiative acceleration for all elements (without the factor 1/A)
!   logkappa   - log Rosseland average opacity (in cm^2/g)
! ------------------------------------------------------------------------------
   subroutine get_radiative_accelerations(logT, logrho, mub, F, nelem, na, z, grad, logkappa)
   implicit none
   real(double),  intent(in) :: logT, logrho, mub, F
   integer,       intent(in) :: nelem
   integer,       intent(in) :: z(nelem)
   real(double),  intent(in) :: na(nelem)
   real(double), intent(out) :: grad(nelem), logkappa
   real(single) :: x
   real(single) :: xi         ! Temperature interpolation factor
   real(single) :: eta        ! Density interpolation factor
   real(single) :: abund(nelem)  ! Abundance for calculation
   real(single) :: epa(4, 4)     ! Electrons per atom
   real(single) :: gam(nelem), ross
   real(single) :: flne, fmu, flmu, slogrho
   real(single), allocatable :: ss(:,:,:)
   real(single), allocatable :: ff(:,:,:,:)
   real(single), allocatable :: ta(:,:,:,:)
   real(single), allocatable :: slogkappa(:,:)
   real(single), allocatable :: logg(:,:,:)
   integer :: perm(nelem), acc_perm(nelem)
   integer :: nae
   integer :: ih(4), jh(4)    ! Temperature and density indices
   integer :: i, j, n

   grad = 0.0d0
   logkappa = 0.0d0

!  Map the abundances that have been passed in onto the 17 elements for which we have data.
!  We ignore the rest (TODO: can we do better by some sort of "averaging" over neighbouring Z values?)
   nae = 0
   perm = 0
   do i = 1, nelem
      do j = 1, num_acc_elem
         if (acceleration_z(j) == z(i)) then
            nae = nae + 1
            abund(nae) = na(i)         ! Abundance list for calculation of acceleration
            perm(i) = nae              ! Input element -> index in computational array
            acc_perm(nae) = j          ! index in computational array -> element ID for calculation
         end if
      end do
   end do

   if (nae == 0) return                      ! No known isotopes

!  Mean molecular weight and log mean molecular weight, in g
   fmu = mub * AMU
   flmu = log10(fmu)
   slogrho = logrho

!  Determine temperature indices for interpolation
   if (logT < 3.5 .or. logT > 8.0) return    ! Temperature out of range
   x = 40.0 * logT / 2.0

   ih(2) = min(max(int(x), 140/2+1), 320/2-2)
   do i = 1, 4
      ih(i) = ih(2)+i-2
   enddo
   xi=2.*(x-ih(2))-1

!  Find electron density and density indices.
!  FIXME: wouldn't it be better to pass this data in as an argument, since we know what it is in the evolution code?
!  One problem with that is that we also need to know the ionisation states for neighbouring points in the interpolation, which
!  would be quite expensive to calculate. Doing it this way is also the easiest way to compare the results with the original
!  OP code because it uses the exact same data.
!  UPDATE: should be fine to pass it in, actually, but keep this for now for easier comparison. Once we remove this, one-for-one
!  comparison becomes much harder

!  Determine electron density indices
   if(.not. find_electron_density(ih, nae, abund(1:nae), flmu, slogrho, xi, jh, flne) ) then
      return   ! Out of range, return 0
   end if

!  Find interpolation factors for density interpolation
   epa = 0
   do n=1,nae
      forall(i=1:4, j=1:4)
         epa(i, j) = epa(i, j) + get_bound_electrons(2*ih(i), 2*jh(j), n)*abund(n)
      end forall
   end do
   eta = geteta(slogrho, xi, flmu, epa, jh)

!  Read atomic data for current parameters and all isotopes
!  NB: we can (and do) re-use the calculation for isotopes with the same Z
   allocate(ss(ntot,4,4))
   allocate(ff(ntot,nae,4,4))
   allocate(ta(ntot,nae,4,4))
   call read_crosssections(ih, jh, nae, acc_perm(1:nae), ff, ta)

!  Calculate integrals over the composition, ie, abundance-weighted cross-sections
   forall (n = 1:ntot, i=1:4, j=1:4)
      ss(n, i, j) = dot_product(ff(n,1:nae,i,j), abund(1:nae))
   end forall

!  Calculate Rosseland opacities and accelerations on grid points
!  Unlike the original code, we don't calculate the derivatives
   allocate(slogkappa(4,4))
   allocate(logg(4,4,nelem))
   call get_rosseland_opacity(nae, flmu, ss, ta, slogkappa, logg(1:4,1:4,1:nae))

!  Interpolate grid points to get value for current parameters
   call interpolate(nae, slogkappa, logg(1:4,1:4,1:nae), xi, eta, ross, gam(1:nae))
   logkappa = ross
   do i = 1, nelem
      if (perm(i) > 0) then
         if (gam(perm(i)) > -30) then  ! FORTRAN doesn't short-circuit boolean logic...
            grad(i) = mub*F*10**( gam(perm(i)) + logkappa ) / CL
         end if
      end if
   end do

!  Clean up
   deallocate(ss)
   deallocate(ff)
   deallocate(ta)
   deallocate(slogkappa)
   deallocate(logg)
   end subroutine get_radiative_accelerations



! ------------------------------------------------------------------------------
! FIND_ELECTRON_DENSITY:
! Find the number density of electrons in the medium for a given composition,
! temperature and overall density of atoms.
! ------------------------------------------------------------------------------
!  Input:
!   ih(4)      - Array of temperature points, to be used for interpolation
!   nelem      - Number of chemical elements (species) to consider
!   na(nelem)  - Abundance of each species
!   flmu       - log10 mu_B, mean molecular weight per baryon [g]
!   logrho     - log density [g/cm^3]
!   xi         - Temperature interpolation factor
!  Output:
!   flne       - log10 ne, log electron density [number / cm^3]
!   jh(4)      - Array of density indices, for interpolation
!  Return value:
!   TRUE or FALSE, depending on whether the function succeeded or not.
! ------------------------------------------------------------------------------
   function find_electron_density(ih, nelem, na, flmu, logrho, xi, jh, flne)
   implicit none
!  Input
   integer, intent(in) :: ih(4), nelem
   real(single), intent(in) :: na(nelem), flmu, logrho, xi
!  Output
   real(single), intent(out) :: flne
   integer, intent(out) :: jh(4)
   logical :: find_electron_density
!  Local variable
   real(single) :: efa(4, 7:118)
   real(single) :: flrho(4, 7:118)
   integer :: n, i, j, jn1, jn2, it
   integer :: jhmin, jhmax
   real(single) :: logrhomin, logrhomax
   real(single) :: u(4), flnei(4), y, zeta, dum

!  Determine minimum and maximum density indices
!  Find density range for selected temperature range. Not what we want yet, but it's a start
   jhmin=max(min_density_index(ih(1)*2), min_density_index(ih(2)*2), min_density_index(ih(3)*2), min_density_index(ih(4)*2))/2
   jhmax=min(max_density_index(ih(1)*2), max_density_index(ih(2)*2), max_density_index(ih(3)*2), max_density_index(ih(4)*2))/2

   efa = 0.0
!$omp parallel do private(n,it,jn1,jn2,j)
   do i = 1, 4
      do n = 1, nelem
         it = 2*(ih(1) + i - 1)
         jn1 = density_start_index(n, (it-min_temperature_index)/2+1)
         jn2 = density_end_index(n, (it-min_temperature_index)/2+1)
         forall (j = jn1:jn2:2)
            efa(i, j/2)=efa(i, j/2) + na(n)*bound_electrons(n, (j - min_density_index(it))/2 + 1, (it-min_temperature_index)/2+1)
         end forall
      end do
   end do

!  Get range where efa > 0 (ie, not fully ionised)
   do j = jhmin, jhmax
      do i = 1,4
         if (efa(i,j) <= 0.) then
            jhmax = MIN(jhmax, j - 1)
            exit
         end if
      end do
   end do

!  Get flrho
   forall(j=jhmin:jhmax)
      flrho(1:4,j)=flmu+0.5*j-log10(efa(1:4,j))
   end forall

!  Find minimum and maximum values for the density, abort if the density is outside this range
   logrhomin = maxval(flrho(1:4, jhmin))
   logrhomax = minval(flrho(1:4, jhmax))
   if (logrho < logrhomin .or. logrho > logrhomax) then
      find_electron_density = .false.
      return
   end if

!  Find density interval
   do j = jhmin-1, jhmax-1
      if (flrho(2, j+1) > logrho) exit
   end do
   if (j == jhmax) then
      find_electron_density = .false.
      return
   end if
   jn1 = min(max(j, jhmin+1), jhmax-2)

   do i=1,4
      forall (j=1:4)
         u(j) = flrho(i, jn1+j-2)
      end forall
      if (.not. solve_cubic(u,logrho,zeta,dum) ) then
         find_electron_density = .false.
         return
      end if
      y=jn1+0.5*(zeta+1)
      flnei(i)=.25*2.0*y
   enddo

!  Interpolate in termperature to find electron density
   flne=fint(flnei,xi)
   find_electron_density = .true.

!  Find density indices
   j = 4.0 * flne / 2.0
   j = min(max(j, jhmin+1), jhmax-2)
   forall (i=1:4)
      jh(i) = j+i-2
   end forall
   end function find_electron_density



! ------------------------------------------------------------------------------
! READ_CROSSSECTIONS
! Read data for the frequency-dependent cross sections for the selected isotopes
! for the given temperature and density indices
! ------------------------------------------------------------------------------
! Input:
!  ih(4)       - Temperature indices
!  jh(4)       - Density indices
!  nelem       - Number of elements
! Output:
!  ff(ntot,nelem,4,4)   - Frequency dependent cross sections
!  ta(ntot,nelem,4,4)   - Correction term for cross section (Seaton 2005, eq. 10)
! ------------------------------------------------------------------------------
   subroutine read_crosssections(ih, jh, nelem, perm, ff, ta)
   implicit none
   integer, intent(in) :: ih(4), jh(4), nelem, perm(nelem)
   real(single), intent(out) :: ff(ntot,nelem,4,4)
   real(single), intent(out) :: ta(ntot,nelem,4,4)
   integer :: i, n

!  Initialise
   ff = 0.0
   ta = 0.0

! NB: test whether parallelising this loop really gives a benefit or not... it's not clear at the moment.
! We're doing a lot of I/O in the inner loops, if we can't really do that in parallel (may depend on the OS) then trying to
! parallise this loop will not give any benefit.
!$omp parallel do private(n)
   do i = 1, 4             ! Loop over temperature indices
      do n = 1, nelem      ! Loop over elements
         call read_freq_dep_cross_section(ih(i), i, 1, 4, jh, n, perm(n), nelem, ff, ta)
      end do
   end do

   end subroutine read_crosssections




! ------------------------------------------------------------------------------
! READ_FREQ_DEP_CROSS_SECTION
! Helper function to read the cross section for one particular element for one
! particular temperature and density.
! It is possible to read data for one specific choise of T and Ne by setting
! the density indices appropriately.
! ------------------------------------------------------------------------------
! Input:
!  it          - Temperature index/2 (70..160)
!  ih          - Temperature point in cross section list
!  j1          - First density point to read in (normally 1)
!  j2          - Last density point to read in (normally 4)
!  jh(4)       - List of density indices
!  elem        - Element number for which to read data (in the list of abundances)
!  pelem       - perm(elem), the index of elem in the list of known elements (1..17)
!  nelem       - Total number of elements in the list
! Output:
!  ff(ntot,nelem,4,4)   - Frequency dependent cross section; set (ih,jh,*,elem)
!  ta(ntot,nelem,4,4)   - Correction term for cross section; set (ih,jh,*,elem)
! ------------------------------------------------------------------------------
   subroutine read_freq_dep_cross_section(it, ih, j1, j2, jh, elem, pelem, nelem, ff, ta)
   implicit none
   integer, intent(in) :: it, ih, j1, j2, jh(4), elem, pelem, nelem
   real(single), intent(out) :: ff(ntot,nelem,4,4)
   real(single), intent(out) :: ta(ntot,nelem,4,4)
   real(single) :: yx(ntot)
   integer :: nx(ntot)
   integer :: f, stat
   integer :: m, j, k
   integer :: izp,ite,ncrse,ntott,jn1,jn2,jn3, np
   real(single) :: amamu,umin,umax,dpack, d
   real(double) :: se, u
   logical :: cached(j1:j2)
   integer :: cache_count

!  Check if data is available from the cache
   cached = .false.
   cache_count = 0
   do j = j1, j2
      if (retrieve_cache(it, jh(j), elem, ff(1:ntot,elem,ih,j), ta(1:ntot,elem,ih,j))) then
         cached(j) = .true.
         cache_count = cache_count+1
      end if
   end do
   if (cache_count == (j2-j1+1)) return      ! All entries found in the data cache

!  Read mono opacity
!$omp critical
   f = get_free_file_unit()
   open(unit=f, file=op_dir(1:op_dir_length)//'m'//zlabel(pelem)//'.'//tlabel(it*2), form='unformatted', status='old', iostat=stat)
!$omp end critical

   read(f) izp,ite,amamu,umin,umax,ncrse,ntott,dpack,jn1,jn2,jn3
!  Skip to the density we want
   do m = jn1, 2*(jh(j1)-1), 2
      read(f)
      read(f)
      read(f)
   end do

   do j=j1, j2
      if (cached(j)) then
         read(f)
         read(f)
         read(f)
         cycle
      end if

      read(f)
      read(f) np           ! Do we need to interpolate the data?
      if(np == 0)then      ! No
         read(f) ff(1:ntot,elem,ih,j)
      else                 ! Yes
         read(f) (nx(k),yx(k),k=1,np)
!  Interpolate yx(k) to get ff(k,n,i,j)
!  This seems very wasteful because it means we spend more time than we need when calculating Rosseland opacities. The calculation
!  for accelerations could probably also be accelerated by doing this differently, but it's a bit more tricky (obviously, the
!  easiest case is if all integrants are discretised on the same frequency grid)
         ff(1,elem,ih,j)=yx(1)
!$omp parallel do private(d,k)
         do m=2,np
            d=(yx(m)-yx(m-1))/float(nx(m)-nx(m-1))
            forall (k = nx(m-1)+1:nx(m)-1)
               ff(k,elem,ih,j)=yx(m-1) + (k-nx(m-1))*d
            end forall
            ff(nx(m),elem,ih,j)=yx(m)
         end do
      end if
   end do
   close(f)

!  Read opta file
!  FIXME: skip elements that are fully ionised, like H and He at higher T (the acceleration always comes back as 0 anyway).
!$omp critical
   f = get_free_file_unit()
   open(unit=f, file=op_dir(1:op_dir_length)//'a'//zlabel(pelem)//'.'//tlabel(it*2), form='unformatted', status='old', iostat=stat)
!$omp end critical
   if (stat == 0) then
      read(f) jn1,jn2,jn3
      do m=jn1, 2*(jh(1)-1),2
        read(f)
        read(f)
      enddo
      do j=j1,j2
         if (cached(j)) then
            read(f)
            read(f)
            cycle
         end if

         read(f) k, np
         read(f) (nx(k),yx(k),k=1,np)
!  Interpolate yx(k) to get ta(j,k,i)
         ta(1,elem,ih,j)=yx(1)
!$omp parallel do private(d,k)
         do m=2,np
            d=(yx(m)-yx(m-1))/float(nx(m)-nx(m-1))
            forall (k = nx(m-1)+1:nx(m)-1)
               ta(k,elem,ih,j)=yx(m-1) + (k-nx(m-1))*d
            end forall
            ta(nx(m),elem,ih,j)=yx(m)
         end do
!  Complete calculation of ta(j,k,it)
         do k=1,ntot
            u=umesh(k)
            se=1.-exp(-u)
            ta(k,elem,ih,j) = se*ff(k,elem,ih,j) - ta(k,elem,ih,j)
         end do
      end do
   end if
   close(f)

!  Store data in the cache so we can retrieve it from there if needed
   do j = j1, j2
      if (.not. cached(j)) call store_cache(it, jh(j), elem, ff(1:ntot,elem,ih,j), ta(1:ntot,elem,ih,j))
   end do

   end subroutine



! ------------------------------------------------------------------------------
! GET_ROSSELAND_OPACITY
! Calculate Rosseland opacities and radiative accelerations on the interpolation
! points of the grid.
! ------------------------------------------------------------------------------
! Input:
!  nelem                - Number of elements
!  lfmu                 - Log10 mean molecular weight [g]
!  s(ntot,4,4)          - Total frequency-dependent cross section
!  ta(ntot,nelem,4,4)   - Modified frequency-dependent cross section for elements
! Output:
!  logkappa(4,4)        - Rosseland opacity
!  logg(4,4,nelem)      - Gravitational acceleration
! ------------------------------------------------------------------------------
   subroutine get_rosseland_opacity(nelem, lfmu, s, ta, logkappa, logg)
   implicit none
   integer, intent(in) :: nelem
   real(single), intent(in) :: lfmu
   real(single), intent(in) :: s(ntot,4,4)
   real(single), intent(in) :: ta(ntot,nelem,4,4)
   real(single), intent(out) :: logkappa(4,4)
   real(single), intent(out) :: logg(4,4,nelem)
   integer :: i, j, n
   real(double) :: drs, dgm(nelem)
   real(single) :: invs(ntot,4,4)

   invs = 1.0 / s
!$omp parallel do private(i,j,n,drs, dgm)
   do j=1,4
      do i=1,4
         drs = sum(invs(1:ntot,i,j))
         forall (n=1:nelem)
            dgm(n) = dot_product(ta(1:ntot,n,i,j), invs(1:ntot,i, j))
         end forall
         logkappa(i,j)=-log10(drs*dv)-16.55280-lfmu
         where (dgm>0)
            dgm = log10(dgm * dv)
         elsewhere
            dgm = -30.0
         end where
         logg(i,j,1:nelem)=dgm(1:nelem)
      end do
   end do
   end subroutine get_rosseland_opacity



   function geteta(slogrho, xi, flmu, epa, jh)
   implicit none
   real(single), intent(in) :: slogrho, xi, flmu, epa(4,4)
   integer, intent(in) :: jh(4)
   real(single) :: geteta
   integer :: i, j
   real(single) :: u(4), zeta(4), dummy
   logical :: res

   do i =1, 4
      forall (j=1:4)
         u(j) = (0.25*2.*jh(j)) + flmu - log10(epa(i, j))
      end forall
      res = solve_cubic(u, slogrho, zeta(i), dummy)
   end do
   geteta = fint(zeta, xi)
   end function geteta



! Cubic fit through 4 data points
   pure function fint(u,r)
   implicit none
   real(single), intent(in) :: u(4), r
   real(single) :: fint

!  If  P(R) =   u(1)  u(2)  u(3)  u(4)
!  for   R  =    -3    -1     1     3
!  then a cubic fit is:
   fint=(                                 &
      27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(   &
      27*(u(3)-u(2))-(u(4)-u(1))   +R*(   &
      -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(   &
      -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
   end function fint




! ------------------------------------------------------------------------------
! INTERPOLATE
! Bi-cubic interpolation for the opacity and the accelerations
! ------------------------------------------------------------------------------
! Input:
!  nelem             - Number of elements
!  slogkappa(4,4)    - Opacity values
!  logg(4,4,nelem)   - Acceleration values
!  xi                - Temperature interpolation factor
!  eta               - Density interpolation factor
! Output:
!  ross              - Rosseland opacity
!  loggam(nelem)     - Dimensionless accelerations for all elements ('gamma')
! ------------------------------------------------------------------------------
   pure subroutine interpolate(nelem, slogkappa, logg, xi, eta, ross, loggam)
   implicit none
   integer, intent(in) :: nelem
   real(single), intent(in) :: slogkappa(4,4), logg(4,4,nelem)
   real(single), intent(in) :: xi, eta
   real(single), intent(out) :: ross, loggam(nelem)
   integer :: i, j, n
   real(single) :: u(4), v(4)

!  Interpolation of opacity
   forall (i=1:4)
      v(i) = fint(slogkappa(i, 1:4), eta)
   end forall
   ross = fint(v, xi)

!  Ditto accelerations
   loggam = 0.0
!$omp parallel do private(i,j,u,v)
   do n = 1, nelem
outer:do i=1,4
         do j=1,4
            if (logg(i,j,n)  <= -30.)then
               loggam(n) = -30.
               exit outer
            end if
            u(j)=logg(i,j,n)
         end do
         v(i)=fint(u,eta)
      end do outer
      if (loggam(n) /= -30.) loggam(n)=fint(v,xi)
   end do

   end subroutine interpolate




! Solve an equation where the function to be inverted is represented by a cubic fit using Newton-Raphson iterations
   function solve_cubic(u,v,z,uz)
      implicit none
      real(single), intent(in) :: u(4), v
      real(single), intent(out) :: z, uz
      logical :: solve_cubic
      integer :: k
      real(single) :: d

!  Find value of z giving P(R)=v
!  First estimate
      z=(2.*v-u(3)-u(2))/(u(3)-u(2))
!  Newton-Raphson iterations
      solve_cubic = .true.
      do k=1,10
         uz=pp(z)
         d=(v-p(z))/uz
         z=z+d
         if (abs(d) < 1.0e-4) return
      enddo
      solve_cubic = .false.
      contains

!  If  P(R) =   u(1)  u(2)  u(3)  u(4)
!  for   R  =    -3    -1    1     3
!  then a cubic fit is:
        pure function P(R)
        implicit none
        real(single), intent(in) :: r
        real(single) :: p
        p = (&
           27*(u(3)+u(2))-3*(u(1)+u(4)) +R*(&
           27*(u(3)-u(2))-(u(4)-u(1))   +R*(&
           -3*(u(2)+u(3))+3*(u(4)+u(1)) +R*(&
           -3*(u(3)-u(2))+(u(4)-u(1)) ))))/48.
        end function

!  First derivative is:
        pure function PP(R)
        implicit none
        real(single), intent(in) :: r
        real(single) :: pp
        pp = (&
           27*(u(3)-u(2))-(u(4)-u(1))+ 2*R*(&
           -3*(u(2)+u(3))+3*(u(4)+u(1)) +3*R*(&
           -3*(u(3)-u(2))+(u(4)-u(1)) )))/48.
        end function

      end function solve_cubic
end module
