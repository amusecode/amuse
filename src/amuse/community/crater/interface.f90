MODULE CraterRadius

  double precision :: target_density
  double precision :: target_gravity
  double precision :: projectile_diameter, projectile_density 
  double precision :: crater_diameter, formation_time
  double precision :: impact_angle, impact_velocity
  integer :: crater_type, target_type, projectile_type

CONTAINS

  FUNCTION set_projectile_type(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_projectile_type
    INTEGER, intent(in) :: inputvalue
    projectile_type = inputvalue
    set_projectile_type=0
  END FUNCTION 

  FUNCTION get_projectile_type(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_projectile_type
    INTEGER, intent(out) :: outputvalue
    outputvalue = projectile_type  
    get_projectile_type=0
  END FUNCTION 
  
  FUNCTION set_projectile_density(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_projectile_density
    DOUBLE PRECISION, intent(in) :: inputvalue
    projectile_density = inputvalue
    set_projectile_density=0
  END FUNCTION 

  FUNCTION get_projectile_density(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_projectile_density
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = projectile_density  
    get_projectile_density=0
  END FUNCTION 

  
  FUNCTION get_target_density(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_target_density
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = target_density
    get_target_density=0
  END FUNCTION 

  
  FUNCTION set_projectile_diameter(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_projectile_diameter
    DOUBLE PRECISION, intent(in) :: inputvalue
    projectile_diameter = inputvalue
    set_projectile_diameter=0
  END FUNCTION 

  FUNCTION get_projectile_diameter(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_projectile_diameter
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = projectile_diameter  
    get_projectile_diameter=0
  END FUNCTION 
  
  FUNCTION set_impact_angle(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_impact_angle
    DOUBLE PRECISION, intent(in) :: inputvalue
    impact_angle = inputvalue
    set_impact_angle=0
  END FUNCTION 

  FUNCTION get_impact_angle(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_impact_angle
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = impact_angle  
    get_impact_angle=0
  END FUNCTION 
  
  FUNCTION set_impact_velocity(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_impact_velocity
    DOUBLE PRECISION, intent(in) :: inputvalue
    impact_velocity = inputvalue
    set_impact_velocity=0
  END FUNCTION 

  FUNCTION get_impact_velocity(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_impact_velocity
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = impact_velocity  
    get_impact_velocity=0
  END FUNCTION 
  
  FUNCTION set_target_density(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_target_density
    DOUBLE PRECISION, intent(in) :: inputvalue
    target_density = inputvalue
    set_target_density=0
  END FUNCTION 
  
  FUNCTION set_target_type(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_target_type
    INTEGER, intent(in) :: inputvalue
    target_type = inputvalue
    set_target_type=0
  END FUNCTION 

  FUNCTION get_target_type(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_target_type
    INTEGER, intent(out) :: outputvalue
    outputvalue = target_type  
    get_target_type=0
  END FUNCTION 

  
  FUNCTION set_target_gravity(inputvalue)
    IMPLICIT NONE
    INTEGER :: set_target_gravity
    DOUBLE PRECISION, intent(in) :: inputvalue
    target_gravity = inputvalue
    set_target_gravity=0
  END FUNCTION 
  
  FUNCTION get_target_gravity(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_target_gravity
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = target_gravity  
    get_target_gravity=0
  END FUNCTION 

  FUNCTION get_crater_diameter(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_crater_diameter
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = crater_diameter  
    get_crater_diameter=0
  END FUNCTION 

  FUNCTION get_formation_time(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_formation_time
    DOUBLE PRECISION, intent(out) :: outputvalue
    outputvalue = formation_time
    get_formation_time=0
  END FUNCTION 

  FUNCTION get_crater_type(outputvalue)
    IMPLICIT NONE
    INTEGER :: get_crater_type
    INTEGER, intent(out) :: outputvalue
    outputvalue = crater_type  
    get_crater_type=0
  END FUNCTION 

  FUNCTION crater_radius(rhoproj, L, v, theta, &
       rhotarget, g, targtype, cratertype, Dyield, Dgault, Tform, &
       Dfinal, Lyield, Lgault)
    implicit none
    integer :: crater_radius
    integer,intent(in) :: targtype
    DOUBLE PRECISION,intent(in) :: L,rhoproj,v,theta,rhotarget,g
    DOUBLE PRECISION,intent(out) :: Dyield,Dgault,Tform
    DOUBLE PRECISION,intent(out) :: Dfinal, Lyield,Lgault
    integer,intent(out) :: cratertype

    target_density = rhotarget
    target_type = targtype
    projectile_density = rhoproj
    projectile_diameter = L
    impact_angle = theta
    impact_velocity = v
    call crater(rhoproj, &
         L, v, theta, rhotarget, g, targtype, &
         crater_type, Dyield,Dgault, formation_time,&
         crater_diameter, Lyield,Lgault)
    crater_radius = 0
  END FUNCTION 

!  FUNCTION evolve_model(end_time)
!    integer :: evolve_model
!    integer,intent(in) :: end_time
!    
!    write (*,*) "evolve model."
!    evolve_model = 1
  !  END FUNCTION evolve_model
  
  FUNCTION evolve_model()
    
    write (*,*) "evolve model."
    evolve_model = 1
  END FUNCTION evolve_model
  
END MODULE
