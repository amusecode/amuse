! Storage module for stars in a stellar evolution code.
! Only used for storing and retrieving, has no checks/calculations.
! Units used:
! - time: julian years
! - mass: solar masses
! - radius: solar radii
! - luminosity: solar luminosities
! - temperature: Kelvin
!
! - Steven Rieder

module store_stars
  use iso_c_binding
  implicit none

  type, public:: star
    private
    integer:: id  ! never-changing identifier
    integer:: track_id  ! track of the star in metisse-may change if a star is removed/a track is deallocated??
    real(c_double):: age
    real(c_double):: CO_core_mass
    real(c_double):: core_mass
    real(c_double):: core_radius
    real(c_double):: convective_envelope_mass
    real(c_double):: convective_envelope_radius
    real(c_double):: epoch
    real(c_double):: initial_mass
    real(c_double):: luminosity
    real(c_double):: main_sequence_lifetime
    real(c_double):: mass
    real(c_double):: metallicity
    real(c_double):: radius
    real(c_double):: spin
    real(c_double):: temperature
    real(c_double):: time_step
    integer:: stellar_type
  end type star

  type, public:: stars
    private
    type(star), allocatable:: star_array(:)
    integer:: num_stars = 0  ! number of stars in the system
    integer:: next_star_id = 1  ! the id of the next star, should only ever increase
  contains
    procedure, public:: initialize
    procedure, public:: new_star
    procedure, public:: remove_star
    procedure, private:: resize
    procedure, private:: lookup_star_id

    procedure, private:: get_property_double
    procedure, private:: get_property_int
    procedure, private:: set_property_double
    procedure, private:: set_property_int
    
    procedure, public:: get_number_of_stars

    ! Every property has a public getter and a setter, listed alphabetically here.
    ! 'id' is only used internally, so it is not exposed.
    procedure, public:: get_age
    procedure, public:: get_CO_core_mass
    procedure, public:: get_core_mass
    procedure, public:: get_core_radius
    procedure, public:: get_convective_envelope_mass
    procedure, public:: get_convective_envelope_radius
    procedure, public:: get_epoch
    procedure, public:: get_initial_mass
    procedure, public:: get_luminosity
    procedure, public:: get_main_sequence_lifetime
    procedure, public:: get_mass
    procedure, public:: get_metallicity
    procedure, public:: get_radius
    procedure, public:: get_spin
    procedure, public:: get_stellar_type
    procedure, public:: get_temperature
    procedure, public:: get_time_step

    procedure, public:: set_age
    procedure, public:: set_CO_core_mass
    procedure, public:: set_core_mass
    procedure, public:: set_core_radius
    procedure, public:: set_convective_envelope_mass
    procedure, public:: set_convective_envelope_radius
    procedure, public:: set_epoch
    procedure, public:: set_initial_mass
    procedure, public:: set_luminosity
    procedure, public:: set_main_sequence_lifetime
    procedure, public:: set_mass
    procedure, public:: set_metallicity
    procedure, public:: set_radius
    procedure, public:: set_spin
    procedure, public:: set_stellar_type
    procedure, public:: set_temperature
    procedure, public:: set_time_step
  end type stars

contains

  subroutine initialize(self)
    class(stars), intent(inout):: self
    allocate(self%star_array(0))    
    self%num_stars = 0
    self%next_star_id = 1
  end subroutine

  subroutine get_number_of_stars(self, number_of_stars)
    class(stars), intent(in):: self
    integer, intent(out):: number_of_stars
    number_of_stars = self%num_stars
  end subroutine

  function new_star(self, initial_mass) result(new_id)
    class(stars), intent(inout):: self
    real(c_double), intent(in):: initial_mass
    integer:: new_id
    integer:: i

    self%num_stars = self%num_stars+1
    !write(*,*) "adding new star  ! so resizing to ", self%num_stars
    !call flush(6)
    call self%resize(self%num_stars)
    i = self%num_stars
    new_id = self%next_star_id

    self%star_array(i)%id = new_id

    self%star_array(i)%age = 0.0_c_double
    self%star_array(i)%CO_core_mass = 0.0_c_double
    self%star_array(i)%core_mass = 0.0_c_double
    self%star_array(i)%core_radius = 0.0_c_double
    self%star_array(i)%convective_envelope_mass = 0.0_c_double
    self%star_array(i)%convective_envelope_radius = 0.0_c_double
    self%star_array(i)%epoch = 0.0_c_double
    self%star_array(i)%initial_mass = initial_mass
    self%star_array(i)%luminosity = 0.0_c_double
    self%star_array(i)%main_sequence_lifetime = 0.0_c_double
    self%star_array(i)%mass = initial_mass
    self%star_array(i)%metallicity = 0.0_c_double
    self%star_array(i)%radius = 0.0_c_double
    self%star_array(i)%spin = 0.0_c_double
    self%star_array(i)%stellar_type = 0_c_int
    self%star_array(i)%time_step = 1.0_c_double
    self%star_array(i)%temperature = 0.0_c_double

    self%next_star_id = new_id+1

  end function new_star

  subroutine remove_star(self, id)
    class(stars), intent(inout):: self
    integer, intent(in):: id

    integer:: i

    do i = 1, self%num_stars
      if (self%star_array(i)%id == id) then
        if (i /= self%num_stars) then
          self%star_array(i:self%num_stars-1) = self%star_array(i+1:self%num_stars)
        end if
        self%num_stars = self%num_stars-1
        !write(*,*) "resizing to ", self%num_stars
        !call flush(6)
        call self%resize(self%num_stars)
        exit
      end if
    end do

  end subroutine remove_star

  subroutine resize(self, required_size)
    class(stars), intent(inout):: self
    integer, intent(in):: required_size
    type(star), allocatable:: temp(:)
    integer:: current_size, new_capacity
 
    if (required_size <= 0) then
      self%num_stars = 0
      if (allocated(self%star_array)) deallocate(self%star_array)
      return
    end if

    if (allocated(self%star_array)) then
      current_size = size(self%star_array)
    else
      current_size = 0
    end if
    if (required_size .lt. current_size) return
  
    new_capacity = current_size
    do while (required_size .gt. new_capacity)
      new_capacity = max(100, int(new_capacity*1.1))
    end do

    if (required_size .lt. current_size) then
      ! Reduce the size of the array
      if (allocated(self%star_array)) then
        deallocate(self%star_array)
      endif
      allocate(self%star_array(required_size))
    else
      ! Increase the size of the array
      if (.not. allocated(self%star_array)) then
        allocate(self%star_array(required_size))
      else
        ! Allocate a temporary array to hold the new data
        !type(star), allocatable:: temp(:)
        allocate(temp(required_size))
        ! Copy the old data to the temporary array
        temp(1:current_size) = self%star_array
        ! Deallocate the old memory
        deallocate(self%star_array)
        ! Allocate new memory for the array
        allocate(self%star_array(required_size))
        ! Copy the data from the temporary array to the new array
        self%star_array = temp
        ! Deallocate the temporary array
        deallocate(temp)
      endif
    endif
 
    !if (.not. allocated(self%star_array) .or. new_capacity > current_size) then
    !  write(*,*) "not allocated OR resizing needed"
    !  call flush(6)
    !  if (allocated(self%star_array)) then
    !    write(*,*) "allocated but resizing needed"
    !    call flush(6)
    !    allocate(temp(current_size))
    !    temp = self%star_array
    !    deallocate(self%star_array)
    !  end if
    !  write(*,*) "allocating array of needed size (", new_capacity, ")"
    !  call flush(6)
    !  write(*,*) "allocated? ", allocated(self%star_array)
    !  call flush(6)
    !  if (.not. allocated(self%star_array)) then
    !      allocate(self%star_array(new_capacity))
    !  else
    !      print *, "Error: Memory already allocated"
    !      stop
    !  endif
    !  !allocate(self%star_array(new_capacity))
    !  if (allocated(temp)) then
    !    write(*,*) "copying data from temp"
    !    call flush(6)
    !    self%star_array(1:current_size) = temp
    !    deallocate(temp)
    !  end if
    !end if
  
  end subroutine resize

  function lookup_star_id(self, id) result(index_of_the_star)
    class(stars), intent(in):: self
    integer, intent(in):: id
    integer:: index_of_the_star
    integer:: i

    do i = 1, self%num_stars
      if (self%star_array(i)%id == id) then
        index_of_the_star = i
        return
      end if
    end do

    index_of_the_star = 0
  end function

  ! Getters for all the stellar properties
  subroutine get_property_double(self, id, property_name, value, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    character(len=*), intent(in):: property_name
    real(c_double), intent(out):: value
    integer:: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      value = 0.0_c_double
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('age')
        value = self%star_array(i)%age
      case ('CO_core_mass')
        value = self%star_array(i)%CO_core_mass
      case ('core_mass')
        value = self%star_array(i)%core_mass
      case ('core_radius')
        value = self%star_array(i)%core_radius
      case ('convective_envelope_mass')
        value = self%star_array(i)%convective_envelope_mass
      case ('convective_envelope_radius')
        value = self%star_array(i)%convective_envelope_radius
      case ('epoch')
        value = self%star_array(i)%epoch
      case ('initial_mass')
        value = self%star_array(i)%initial_mass
      case ('luminosity')
        value = self%star_array(i)%luminosity
      case ('main_sequence_lifetime')
        value = self%star_array(i)%main_sequence_lifetime
      case ('mass')
        value = self%star_array(i)%mass
      case ('metallicity')
        value = self%star_array(i)%metallicity
      case ('radius')
        value = self%star_array(i)%radius
      case ('spin')
        value = self%star_array(i)%spin
      case ('temperature')
        value = self%star_array(i)%temperature
      case ('time_step')
        value = self%star_array(i)%time_step
      case default
        value = 0.0_c_double
        error = -2  ! property not found
    end select
    error = 0
  end subroutine get_property_double

  subroutine get_property_int(self, id, property_name, value, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    character(len=*), intent(in):: property_name
    integer, intent(out):: value
    integer:: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      value = 0.0_c_double
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('stellar_type')
        value = self%star_array(i)%stellar_type
      case default
        value = 0
        error = -2  ! property not found
        return
    end select
    error = 0
  end subroutine get_property_int

  ! Setters for all the stellar properties
  subroutine set_property_double(self, id, property_name, value, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    character(len=*), intent(in):: property_name
    real(c_double), intent(in):: value
    integer:: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('age')
        self%star_array(i)%age = value
      case ('CO_core_mass')
        self%star_array(i)%CO_core_mass = value
      case ('core_mass')
        self%star_array(i)%core_mass = value
      case ('core_radius')
        self%star_array(i)%core_radius = value
      case ('convective_envelope_mass')
        self%star_array(i)%convective_envelope_mass = value
      case ('convective_envelope_radius')
        self%star_array(i)%convective_envelope_radius = value
      case ('epoch')
        self%star_array(i)%epoch = value
      case ('initial_mass')
        self%star_array(i)%initial_mass = value
      case ('luminosity')
        self%star_array(i)%luminosity = value
      case ('main_sequence_lifetime')
        self%star_array(i)%main_sequence_lifetime = value
      case ('mass')
        self%star_array(i)%mass = value
      case ('metallicity')
        self%star_array(i)%metallicity = value
      case ('radius')
        self%star_array(i)%radius = value
      case ('spin')
        self%star_array(i)%spin = value
      case ('temperature')
        self%star_array(i)%temperature = value
      case ('time_step')
        self%star_array(i)%time_step = value
      case default
        error = -2  ! property not found
        return
    end select
    error = 0
  end subroutine

  subroutine set_property_int(self, id, property_name, value, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    character(len=*), intent(in):: property_name
    integer, intent(in):: value
    integer:: i, error
  
    i = lookup_star_id(self, id)
    if (i == 0) then
      error = -1  ! star not found
      return
    end if
  
    select case (trim(property_name))
      case ('stellar_type')
        self%star_array(i)%stellar_type = value
      case default
        error = -2  ! property not found
        return
    end select
    error = 0
  end subroutine


  ! getters for all the stellar properties
  subroutine get_age(self, id, age, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: age
    integer:: error
    call get_property_double(self, id, 'age', age, error)
  end subroutine

  subroutine get_CO_core_mass(self, id, CO_core_mass, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: CO_core_mass
    integer:: error
    call get_property_double(self, id, 'CO_core_mass', CO_core_mass, error)
  end subroutine

  subroutine get_core_mass(self, id, core_mass, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: core_mass
    integer:: error
    call get_property_double(self, id, 'core_mass', core_mass, error)
  end subroutine

  subroutine get_core_radius(self, id, core_radius, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: core_radius
    integer:: error
    call get_property_double(self, id, 'core_radius', core_radius, error)
  end subroutine

  subroutine get_convective_envelope_mass(self, id, convective_envelope_mass, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: convective_envelope_mass
    integer:: error
    call get_property_double(self, id, 'convective_envelope_mass', convective_envelope_mass, error)
  end subroutine

  subroutine get_convective_envelope_radius(self, id, convective_envelope_radius, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: convective_envelope_radius
    integer:: error
    call get_property_double(self, id, 'convective_envelope_radius', convective_envelope_radius, error)
  end subroutine

  subroutine get_epoch(self, id, epoch, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: epoch
    integer:: error
    call get_property_double(self, id, 'epoch', epoch, error)
  end subroutine

  subroutine get_initial_mass(self, id, initial_mass, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: initial_mass
    integer:: error
    call get_property_double(self, id, 'initial_mass', initial_mass, error)
  end subroutine

  subroutine get_luminosity(self, id, luminosity, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: luminosity
    integer:: error
    call get_property_double(self, id, 'luminosity', luminosity, error)
  end subroutine

  subroutine get_main_sequence_lifetime(self, id, main_sequence_lifetime, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: main_sequence_lifetime
    integer:: error
    call get_property_double(self, id, 'main_sequence_lifetime', main_sequence_lifetime, error)
  end subroutine

  subroutine get_mass(self, id, mass, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: mass
    integer:: error
    call get_property_double(self, id, 'mass', mass, error)
  end subroutine

  subroutine get_metallicity(self, id, metallicity, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: metallicity
    integer:: error
    call get_property_double(self, id, 'metallicity', metallicity, error)
  end subroutine

  subroutine get_radius(self, id, radius, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: radius
    integer:: error
    call get_property_double(self, id, 'radius', radius, error)
  end subroutine

  subroutine get_spin(self, id, spin, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: spin
    integer:: error
    call get_property_double(self, id, 'spin', spin, error)
  end subroutine

  subroutine get_stellar_type(self, id, stellar_type, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    integer, intent(out):: stellar_type
    integer:: error
    call get_property_int(self, id, 'stellar_type', stellar_type, error)
  end subroutine

  subroutine get_temperature(self, id, temperature, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: temperature
    integer:: error
    call get_property_double(self, id, 'temperature', temperature, error)
  end subroutine

  subroutine get_time_step(self, id, time_step, error)
    class(stars), intent(in):: self
    integer, intent(in):: id
    real(c_double), intent(out):: time_step
    integer:: error
    call get_property_double(self, id, 'time_step', time_step, error)
  end subroutine


  ! setters for all the stellar properties (in the same order as the getters)

  subroutine set_age(self, id, age, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: age
    integer:: error
    call set_property_double(self, id, 'age', age, error)
  end subroutine

  subroutine set_CO_core_mass(self, id, CO_core_mass, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: CO_core_mass
    integer:: error
    call set_property_double(self, id, 'CO_core_mass', CO_core_mass, error)
  end subroutine

  subroutine set_core_mass(self, id, core_mass, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: core_mass
    integer:: error
    call set_property_double(self, id, 'core_mass', core_mass, error)
  end subroutine

  subroutine set_core_radius(self, id, core_radius, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: core_radius
    integer:: error
    call set_property_double(self, id, 'core_radius', core_radius, error)
  end subroutine

  subroutine set_convective_envelope_mass(self, id, convective_envelope_mass, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: convective_envelope_mass
    integer:: error
    call set_property_double(self, id, 'convective_envelope_mass', convective_envelope_mass, error)
  end subroutine

  subroutine set_convective_envelope_radius(self, id, convective_envelope_radius, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: convective_envelope_radius
    integer:: error
    call set_property_double(self, id, 'convective_envelope_radius', convective_envelope_radius, error)
  end subroutine

  subroutine set_epoch(self, id, epoch, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: epoch
    integer:: error
    call set_property_double(self, id, 'epoch', epoch, error)
  end subroutine

  subroutine set_initial_mass(self, id, initial_mass, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: initial_mass
    integer:: error
    call set_property_double(self, id, 'initial_mass', initial_mass, error)
  end subroutine

  subroutine set_luminosity(self, id, luminosity, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: luminosity
    integer:: error
    call set_property_double(self, id, 'luminosity', luminosity, error)
  end subroutine

  subroutine set_main_sequence_lifetime(self, id, main_sequence_lifetime, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: main_sequence_lifetime
    integer:: error
    call set_property_double(self, id, 'main_sequence_lifetime', main_sequence_lifetime, error)
  end subroutine

  subroutine set_mass(self, id, mass, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: mass
    integer:: error
    call set_property_double(self, id, 'mass', mass, error)
  end subroutine

  subroutine set_metallicity(self, id, metallicity, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: metallicity
    integer:: error
    call set_property_double(self, id, 'metallicity', metallicity, error)
  end subroutine

  subroutine set_radius(self, id, radius, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: radius
    integer:: error
    call set_property_double(self, id, 'radius', radius, error)
  end subroutine

  subroutine set_spin(self, id, spin, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: spin
    integer:: error
    call set_property_double(self, id, 'spin', spin, error)
  end subroutine

  subroutine set_stellar_type(self, id, stellar_type, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    integer, intent(in):: stellar_type
    integer:: error
    call set_property_int(self, id, 'stellar_type', stellar_type, error)
  end subroutine

  subroutine set_temperature(self, id, temperature, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: temperature
    integer:: error
    call set_property_double(self, id, 'temperature', temperature, error)
  end subroutine

  subroutine set_time_step(self, id, time_step, error)
    class(stars), intent(inout):: self
    integer, intent(in):: id
    real(c_double), intent(in):: time_step
    integer:: error
    call set_property_double(self, id, 'time_step', time_step, error)
  end subroutine

end module store_stars
