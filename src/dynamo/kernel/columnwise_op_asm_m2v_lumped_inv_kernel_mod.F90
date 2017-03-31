!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which constructs a CMA representation of the diagonally
!> lumped inverse of the vertical velocity mass matrix
!> @details Extract the diagonal of a locally assembled matrix
!> (LMA) for the vertical velocity mass matrix. The inverse of this diagonal is
!> then assembled into a CMA. Note that this CMA is diagonal, i.e. it has
!> parameters \f$\alpha=\beta=1\f$ and \f$\gamma_-=\gamma_+=0\f$.
!> 

module columnwise_op_asm_m2v_lumped_inv_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,                    &
                                    GH_OPERATOR, GH_COLUMNWISE_OPERATOR,    &
                                    GH_READ, GH_WRITE,                      &
                                    ANY_SPACE_1,                            &
                                    GH_COLUMN_BANDED_DOFMAP,                &
                                    CELLS 

use constants_mod,           only : r_def, i_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_asm_m2v_lumped_inv_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                       &
       arg_type(GH_OPERATOR,           GH_READ,  ANY_SPACE_1, ANY_SPACE_1), &
       ! NOT CURRENTLY SUPPORTED BY PSY
       arg_type(GH_COLUMNWISE_OPERATOR, GH_WRITE, ANY_SPACE_1, ANY_SPACE_1) &
       /)
  type(func_type) :: meta_funcs(1) =  (/                                    &
       func_type(ANY_SPACE_1, GH_COLUMN_BANDED_DOFMAP)                      &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: columnwise_op_asm_m2v_lumped_inv_kernel_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface columnwise_op_asm_m2v_lumped_inv_kernel_type
   module procedure columnwise_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public columnwise_op_asm_m2v_lumped_inv_kernel_code
contains
  
type(columnwise_op_asm_m2v_lumped_inv_kernel_type) function columnwise_constructor() result(self)
  implicit none
  return
end function columnwise_constructor

!> @brief The subroutine which is called directly from the PSY layer and
!> assembles the LMA into a CMA
!> @details Given an LMA representation of local operator for the 
!> vertical velocity mass matrix, assemble the columnwise matrix
!> which represents the inverse lumped mass matrix
!>
!> @param [in]  cell Horizontal cell index
!> @param [in]  nlayers Number of vertical layers
!> @param [in]  ncell_3d Total number of cells
!> @param [in]  ncell_2d Number of cells in 2d grid
!> @param [in]  local_stencil Locally assembled matrix
!> @param [out] columnwise_matrix Banded matrix to assemble into
!> @param [in]  nrow Number of rows in the banded matrix
!> @param [in]  ncol Number of columns in the banded matrix
!> @param [in]  bandwidth Bandwidth of the banded matrix
!> @param [in]  alpha Banded matrix parameter \f$\alpha=1\f$
!> @param [in]  beta Banded matrix parameter \f$\beta=1\f$
!> @param [in]  gamma_m Banded matrix parameter \f$\gamma_-=0\f$
!> @param [in]  gamma_p Banded matrix parameter \f$\gamma_+=0\f$
!> @param [in]  ndf_to Number of degrees of freedom per cell for the to-space
!> @param [in]  ndf_from Number of degrees of freedom per cell for the from-sp
!> @param [in]  column_banded_dofmap_to List of offsets for to-space
!> @param [in]  column_banded_dofmap_from List of offsets for from-space
subroutine columnwise_op_asm_m2v_lumped_inv_kernel_code(cell,                    &
                                                        nlayers,                 &
                                                        ncell_3d,                &
                                                        ncell_2d,                &
                                                        local_stencil,           &
                                                        columnwise_matrix,       &
                                                        nrow,                    &
                                                        ncol,                    &
                                                        bandwidth,               &
                                                        alpha,                   &
                                                        beta,                    &
                                                        gamma_m,                 &
                                                        gamma_p,                 &
                                                        ndf_to,                  &
                                                        ndf_from,                &
                                                        column_banded_dofmap_to, &
                                                        column_banded_dofmap_from)

  implicit none
  
  ! Arguments
  integer(kind=i_def),                                      intent(in)  :: cell
  integer(kind=i_def),                                      intent(in)  :: nlayers
  integer(kind=i_def),                                      intent(in)  :: ncell_3d
  integer(kind=i_def),                                      intent(in)  :: ncell_2d
  integer(kind=i_def),                                      intent(in)  :: nrow
  integer(kind=i_def),                                      intent(in)  :: ncol
  integer(kind=i_def),                                      intent(in)  :: bandwidth
  integer(kind=i_def),                                      intent(in)  :: ndf_to
  integer(kind=i_def),                                      intent(in)  :: ndf_from
  integer(kind=i_def),                                      intent(in)  :: alpha
  integer(kind=i_def),                                      intent(in)  :: beta
  integer(kind=i_def),                                      intent(in)  :: gamma_m
  integer(kind=i_def),                                      intent(in)  :: gamma_p
  integer(kind=i_def), dimension(ndf_to,nlayers),           intent(in)  :: column_banded_dofmap_to
  integer(kind=i_def), dimension(ndf_from,nlayers),         intent(in)  :: column_banded_dofmap_from
  real   (kind=r_def), dimension(ndf_to,ndf_from,ncell_3d), intent(in)  :: local_stencil
  real   (kind=r_def), dimension(bandwidth,nrow,ncell_2d),  intent(out) :: columnwise_matrix

  ! Internal parameters
  integer(kind=i_def) :: df1    ! Loop index for dofs
  integer(kind=i_def) :: i      ! Row and column index
  integer(kind=i_def) :: ik     ! ncell3d counter
  integer(kind=i_def) :: k      ! nlayers  counter

  k = alpha + beta + gamma_m + gamma_p + ncol

  ! Initialise matrix to zero
  columnwise_matrix( :, :, cell ) = 0.0_r_def
  ! Loop over all vertical layers add add up diagonal entries
  do k = 1, nlayers
    ik = (cell-1)*nlayers + k ! cell index in 3d
    do df1 = 1, ndf_to
      i = column_banded_dofmap_to( df1, k )
      columnwise_matrix( 1, i, cell ) = columnwise_matrix( 1, i, cell ) &
                                      + local_stencil( df1 ,df1, ik )
    end do
  end do

  ! Calculate inverse of diagonal entries
  columnwise_matrix(:,:,cell) = 1.0_r_def/columnwise_matrix(:,:,cell)

end subroutine columnwise_op_asm_m2v_lumped_inv_kernel_code

end module columnwise_op_asm_m2v_lumped_inv_kernel_mod
