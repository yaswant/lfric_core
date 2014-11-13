!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Defines various constants.

!> @details Various physical and geometrical constants are defined in this module.
!! Their values are also set here.
module constants_mod
use, intrinsic :: iso_fortran_env, only : real64, int32
implicit none

!Define default application-defined kinds for all intrinsic data types
!Reals
!< default real for the application
integer,       parameter    :: r_def     = real64

real, private               :: r_val
double precision, private   :: dp_val

!< native kind for real 
integer,       parameter    :: r_native  = kind(r_val)
!< native kind for double precision
integer,       parameter    :: dp_native = kind(dp_val)

!Complex

!Integers
!< default integer for the application
integer,       parameter    :: i_def    = int32

integer, private            :: i_val
!< native kind for integer 
integer,       parameter    :: i_native = kind(i_val)

!Logical
logical, private            :: l_val
integer,       parameter    :: l_def    = kind(l_val)
integer,       parameter    :: l_native = kind(l_val)

!Character
character, private          :: c_val
integer,       parameter    :: c_def    = kind(c_val)
integer,       parameter    :: c_native = kind(c_val)

!String Length
integer,       parameter    :: str_def          = 128
integer,       parameter    :: str_long         = 255
integer,       parameter    :: str_max_filename = 255

! Default numbers
real(kind=r_def), parameter :: large_real = huge(0.0_r_def)

!Numerical constants
real(kind=r_def), parameter :: pi  = 3.141592654_r_def    !< pi value
real(kind=r_def), parameter :: eps = 3.0e-15_r_def        !< relative precision

! Linear solver constants
integer (kind=i_def), parameter :: max_iter = 99 !< maximum iteration number for solver
integer (kind=i_def), parameter :: cg_solver     = 1, &
                                   bicg_solver   = 2, &
                                   jacobi_solver = 3, &
                                   gmres_solver  = 4, &
                                   gcr_solver    = 5
integer (kind=i_def), parameter :: solver_option = bicg_solver
integer (kind=i_def), parameter :: no_pre_cond       = -1, &
                                   diagonal_pre_cond = 1
real(kind=r_def),     parameter :: solver_tol = 1.0e-4_r_def
integer (kind=i_def), parameter :: gcrk  = 4

! physical constants
real(kind=r_def), parameter :: gravity = 9.80665_r_def
real(kind=r_def), parameter :: rd = 287.05_r_def 
real(kind=r_def), parameter :: cp = 1005.0_r_def 
real(kind=r_def), parameter :: kappa = 287.05_r_def/1005.0_r_def 
real(kind=r_def), parameter :: p_zero = 100000.0_r_def
real(kind=r_def), parameter :: n_sq = 0.0001_r_def
real(kind=r_def), parameter :: omega_unscaled = 7.292116E-5_r_def
real(kind=r_def), parameter :: earth_radius_unscaled = 6371229.0_r_def

! Small earth scalings
real(kind=r_def), parameter :: earth_scaling = 1.0_r_def
real(kind=r_def)            :: omega = omega_unscaled*earth_scaling
real(kind=r_def)            :: earth_radius = earth_radius_unscaled*earth_scaling

end module constants_mod

