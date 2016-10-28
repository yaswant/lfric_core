!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes boundary part of the weighted divergence operator

!> @details The kernel computes the boundary part of the weighted divergence operator
!>         This consists of \f[<\sigma,\theta*\mathbf{v}\cdot\mathbf{n}> \f]
!>          where sigma is the W3 test function, v is the W2 trial function,
!>          theta is the potential temperature, and \mathbf{n} is the outward pointing normal.

module weighted_div_bd_kernel_mod
  use kernel_mod,              only : kernel_type
  use argument_mod,            only : arg_type, func_type,                 &
                                      GH_OPERATOR, GH_FIELD,               &
                                      GH_READ, GH_WRITE,                   &
                                      W2, W3, Wtheta, GH_BASIS,            &
                                      CELLS
  use constants_mod,           only : r_def, i_def
  use reference_element_mod,    only: nfaces_h, out_face_normal


  implicit none

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: weighted_div_bd_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                               &
       arg_type(GH_OPERATOR, GH_WRITE, W3, W2),                       &
       arg_type(GH_FIELD,    GH_READ, Wtheta)                         &
      /)
    type(func_type) :: meta_funcs(3) = (/                             &
       func_type(W3, GH_BASIS),                                       &
       func_type(W2, GH_BASIS),                                       &
       func_type(Wtheta, GH_BASIS)                               &
      /)
    integer :: iterates_over = CELLS
  contains
    procedure, nopass :: weighted_div_bd_code
  end type

  !-------------------------------------------------------------------------------
  ! Constructors
  !-------------------------------------------------------------------------------

  ! overload the default structure constructor for function space
  interface weighted_div_bd_kernel_type
    module procedure weighted_div_bd_kernel_constructor
  end interface

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public weighted_div_bd_code
contains

  type(weighted_div_bd_kernel_type) function weighted_div_bd_kernel_constructor() result(self)
    return
  end function weighted_div_bd_kernel_constructor

  !> @brief Compute the boundary terms in the weighted divergence for the Helmholtz lhs
  !! @param[in] cell Cell number
  !! @param[in] nlayers Number of layers
  !! @param[in] ncell_3d ncell*ndf
  !! @param[in] div Local stencil of the div operator
  !! @param[in] theta Potential temperature
  !! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
  !! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
  !! @param[in] ndf_wtheta Number of degrees of freedom per cell for wtheta
  !! @param[in] undf_wtheta Number unique of degrees of freedom  for wtheta
  !! @param[in] map_wtheta Dofmap for wtheta
  !! @param[in] nqp_v Number of quadrature points in the vertical
  !! @param[in] nqp_h_1d Number of quadrature points in a single horizontal direction
  !! @param[in] wqp_v Vertical quadrature weights
  !! @param[in] w2_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] w3_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces
  !! @param[in] wtheta_basis_face Basis functions evaluated at gaussian quadrature points on horizontal faces

  subroutine weighted_div_bd_code(cell, nlayers, ncell_3d,                                    &
                                  div, theta,                                                 &
                                  ndf_w2, ndf_w3,                                             &
                                  ndf_wtheta, undf_wtheta,                                    &
                                  map_wtheta,                                                 &
                                  nqp_v, nqp_h_1d, wqp_v, w2_basis_face, w3_basis_face,       &
                                  wtheta_basis_face)


    ! Arguments
    integer(kind=i_def), intent(in) :: cell, nlayers, nqp_v, nqp_h_1d, ncell_3d
    integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta

    integer(kind=i_def), dimension(ndf_wtheta),  intent(in) :: map_wtheta

    real(kind=r_def), dimension(4,3,ndf_w2,nqp_h_1d,nqp_v), intent(in)     :: w2_basis_face
    real(kind=r_def), dimension(4,1,ndf_w3,nqp_h_1d,nqp_v), intent(in)     :: w3_basis_face
    real(kind=r_def), dimension(4,1,ndf_wtheta,nqp_h_1d,nqp_v), intent(in) :: wtheta_basis_face

    real(kind=r_def), dimension(ndf_w3,ndf_w2,ncell_3d), intent(inout) :: div
    real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta

    real(kind=r_def), dimension(nqp_v), intent(in)      ::  wqp_v

    !Internal variables
    integer(kind=i_def)              :: df, df2, df3, k, ik, face
    integer(kind=i_def)              :: qp1, qp2

    real(kind=r_def), dimension(ndf_wtheta)      :: theta_e

    real(kind=r_def) :: v(3), integrand
    real(kind=r_def) :: theta_at_fquad, bd_term

    do k = 0, nlayers - 1
      ik = k + 1 + (cell-1)*nlayers
      do face = 1, nfaces_h

        bd_term = 0.0_r_def

        do df = 1,ndf_wtheta
          theta_e(df) = theta(map_wtheta(df) + k)
        end do
        do df2 = 1, ndf_w2
          do df3 = 1, ndf_w3
            div(df3,df2,ik) = 0.0_r_def
            do qp2 = 1, nqp_v
              do qp1 = 1, nqp_h_1d
                theta_at_fquad = 0.0_r_def
                do df = 1, ndf_wtheta
                  theta_at_fquad   = theta_at_fquad + theta_e(df)*wtheta_basis_face(face,1,df,qp1,qp2)
                end do
                v  = w2_basis_face(face,:,df2,qp1,qp2)
                bd_term = dot_product(v, out_face_normal(:, face))*theta_at_fquad
                integrand = wqp_v(qp1)*wqp_v(qp2)*w3_basis_face(face,1,df3,qp1,qp2) &
                             * bd_term
                div(df3,df2,ik) = div(df3,df2,ik) + integrand
              end do ! qp1
            end do ! qp2
          end do ! df3
        end do ! df2
      end do ! faces
    end do ! layers

  end subroutine weighted_div_bd_code

end module weighted_div_bd_kernel_mod
