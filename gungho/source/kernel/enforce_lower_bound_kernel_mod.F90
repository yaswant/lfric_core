!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Enforces (non-conservatively) a lower bound on a field, such that
!>        if the value of the field is below this value, then the field is set
!>        to zero.  The rational behind this is that if the threshold is set
!>        to be O(SPACING(field)) then any small floating point errors will 
!>        be removed. 
module enforce_lower_bound_kernel_mod

  use argument_mod,  only : arg_type, func_type,       &
                            GH_FIELD, GH_READ, GH_INC, &
                            GH_BASIS, GH_REAL,         &
                            CELLS,ANY_SPACE_1
  use constants_mod, only : i_def, r_def
  use kernel_mod,    only : kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  ! The type declaration for the kernel. Contains the metadata needed by the
  ! Psy layer.
  !
  type, public, extends(kernel_type) :: enforce_lower_bound_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/    &
        arg_type(GH_FIELD,   GH_INC,  ANY_SPACE_1), &
        arg_type(GH_REAL,    GH_READ ) &
        /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass ::enforce_lower_bound_code
  end type

  !---------------------------------------------------------------------------
  ! Constructors
  !---------------------------------------------------------------------------

  ! overload the default structure constructor
  interface enforce_lower_bound_kernel_type
    module procedure enforce_lower_bound_kernel_constructor
  end interface

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public enforce_lower_bound_code

contains

type(enforce_lower_bound_kernel_type) &
function enforce_lower_bound_kernel_constructor() result(self)
implicit none
return
end function enforce_lower_bound_kernel_constructor

!> @brief Limits the field dofs by some measure of CFL limit
!! @param[in] nlayers Number of layers
!! @param[inout] field Field
!! @param[in] lower_bound The lower bound
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column
subroutine enforce_lower_bound_code(nlayers, field, lower_bound, &
                                    ndf, undf, map)

  implicit none

  !Arguments
  integer, intent(in) :: nlayers, ndf, undf
  integer, dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: field
  real(kind=r_def), intent(in)    :: lower_bound

  !Internal variables
  integer          :: df, k

  do k = 0, nlayers-1
    do df = 1, ndf
        
      ! Clip field
      if (field(map(df)+k) < lower_bound) field(map(df)+k) = 0.0_r_def
        
    end do
  end do

end subroutine enforce_lower_bound_code

end module enforce_lower_bound_kernel_mod
