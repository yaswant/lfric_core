!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which calculates the diagonal contribution to the term D_h*M_{2v,lumped,inv}*D_h.
!> @details Takes the operator D_h in LMA representation and the field representation of
!> the diagonally lumped horizontal velocity mass matrix. Based on this, the kernel
!> assembles a CMA which contains the diagonal couplings of the term
!> D_h*M_{2v,lumped,inv}*D_h^T

module columnwise_op_asm_diag_hmht_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,       &
                                    GH_OPERATOR, GH_FIELD,     &
                                    GH_COLUMNWISE_OPERATOR,    &
                                    GH_READ, GH_WRITE,         &
                                    ANY_SPACE_1, ANY_SPACE_2,  &
                                    GH_COLUMN_BANDED_DOFMAP,   &
                                    CELLS 

use constants_mod,           only : r_def, i_def

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: columnwise_op_asm_diag_hmht_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                       &
       arg_type(GH_OPERATOR,         GH_READ,  ANY_SPACE_1, ANY_SPACE_2),   &
       arg_type(GH_FIELD,            GH_READ,  ANY_SPACE_2),                &
       arg_type(GH_COLUMNWISE_OPERATOR, GH_WRITE, ANY_SPACE_1, ANY_SPACE_1) &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass :: columnwise_op_asm_diag_hmht_kernel_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! Overload the default structure constructor for function space
interface columnwise_op_asm_diag_hmht_kernel_type
   module procedure columnwise_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public columnwise_op_asm_diag_hmht_kernel_code
contains
  
type(columnwise_op_asm_diag_hmht_kernel_type) function columnwise_constructor() result(self)
  implicit none
  return
end function columnwise_constructor

!> @brief The subroutine which is called directly from the PSY layer and
!> assembles the LMA into a CMA
!> @details Given an LMA representation of the operator mapping between two
!> horizontally discontinuous spaces, assemble the columnwise matrix
!> representation of the operator.
!>
!> @param [in]  cell Horizontal cell index
!> @param [in]  nlayers Number of vertical layers
!> @param [in]  ncell_2d Number of cells in 2d grid
!> @param [in]  ncell_3d Total number of cells
!> @param [in]  local_stencil Locally assembled matrix
!> @param [in]  field_data Values of the field representing the
!>              diagonally lumped W2h velocity mass matrix
!> @param [out] columnwise_matrix Banded matrix to assemble into
!> @param [in]  nrow Number of rows in the banded matrix
!> @param [in]  ncol Number of columns in the banded matrix
!> @param [in]  bandwidth Bandwidth of the banded matrix
!> @param [in]  alpha banded Matrix parameter \f$\alpha\f$
!> @param [in]  beta banded Matrix parameter \f$\beta\f$
!> @param [in]  gamma_m Banded matrix parameter \f$\gamma_-\f$
!> @param [in]  gamma_p Banded matrix parameter \f$\gamma_+\f$
!> @param [in]  ndf_w3 Number of dofs per cell for the W_3 space
!> @param [in]  ndf_w2h Number of dofs per cell for the W_{2h} space
!> @param [in]  undf_w2h Number of unique dofs per column for the W_{2h} space
!> @param [in]  dofmap Indirection map for W2h space
!> @param [in]  column_banded_dofmap List of offsets for W3-space
subroutine columnwise_op_asm_diag_hmht_kernel_code(cell,              &
                                                   nlayers,           &
                                                   ncell_2d,          &
                                                   ncell_3d,          &
                                                   local_stencil,     &
                                                   field_data,        &
                                                   columnwise_matrix, &
                                                   nrow,              &
                                                   ncol,              &
                                                   bandwidth,         &
                                                   alpha,             &
                                                   beta,              &
                                                   gamma_m,           &
                                                   gamma_p,           &
                                                   ndf_w3,            &
                                                   ndf_w2h,           &
                                                   undf_w2h,          &
                                                   dofmap,            &
                                                   column_banded_dofmap)

  implicit none
    
  ! Arguments
  integer(kind=i_def),                                     intent(in)  :: cell
  integer(kind=i_def),                                     intent(in)  :: nlayers
  integer(kind=i_def),                                     intent(in)  :: ncell_3d
  integer(kind=i_def),                                     intent(in)  :: ncell_2d
  integer(kind=i_def),                                     intent(in)  :: alpha
  integer(kind=i_def),                                     intent(in)  :: beta
  integer(kind=i_def),                                     intent(in)  :: gamma_m
  integer(kind=i_def),                                     intent(in)  :: gamma_p
  integer(kind=i_def),                                     intent(in)  :: undf_w2h
  integer(kind=i_def),                                     intent(in)  :: nrow
  integer(kind=i_def),                                     intent(in)  :: ncol
  integer(kind=i_def),                                     intent(in)  :: bandwidth
  integer(kind=i_def),                                     intent(in)  :: ndf_w3
  integer(kind=i_def),                                     intent(in)  :: ndf_w2h
  integer(kind=i_def), dimension(ndf_w2h),                 intent(in)  :: dofmap
  integer(kind=i_def), dimension(ndf_w3,nlayers),          intent(in)  :: column_banded_dofmap
  real   (kind=r_def), dimension(ndf_w3,ndf_w2h,ncell_3d), intent(in)  :: local_stencil
  real   (kind=r_def), dimension(undf_w2h),                intent(in)  :: field_data
  real   (kind=r_def), dimension(bandwidth,nrow,ncell_2d), intent(out) :: columnwise_matrix


  ! Internal parameters
  integer(kind=i_def) :: df1, df2, df3  ! Loop indices for dofs
  integer(kind=i_def) :: i,j            ! Row and column index index
  integer(kind=i_def) :: j_minus        ! First column in a row
  integer(kind=i_def) :: ik             ! ncell3d counter
  integer(kind=i_def) :: k              ! nlayers  counter
  real   (kind=r_def) :: tmp            ! Local contribution

  k = gamma_m + ncol

  ! Initialise matrix to zero
  columnwise_matrix( :, :, cell ) = 0.0_r_def
  ! Loop over all vertical layers
  do k = 1, nlayers
    ik = (cell-1)*nlayers + k ! cell index in 3d
    do df1 = 1, ndf_w3
      i = column_banded_dofmap( df1, k )
      j_minus = ceiling((alpha*i-gamma_p)/(1.0_8*beta),i_def)
      do df2 = 1, ndf_w3
        tmp = 0.0_r_def
        do df3 = 1, ndf_w2h
          tmp = tmp                           &
              + local_stencil( df1 ,df3, ik ) &
              * local_stencil( df2 ,df3, ik ) &
              / field_data(dofmap(df3) + k-1)
        end do
        j = column_banded_dofmap( df2, k )
        columnwise_matrix( j-j_minus+1, i, cell )      &
           = columnwise_matrix( j-j_minus+1, i, cell ) + tmp
      end do
    end do
  end do

end subroutine columnwise_op_asm_diag_hmht_kernel_code

end module columnwise_op_asm_diag_hmht_kernel_mod
