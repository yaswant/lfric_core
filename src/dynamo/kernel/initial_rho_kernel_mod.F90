!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes LHS of Galerkin projection and solves equation in W3 space

module initial_rho_kernel_mod

use argument_mod,               only : arg_type, func_type,            &
                                       GH_FIELD, GH_READ, GH_WRITE,    &
                                       W0, W3,                         &
                                       GH_BASIS, GH_DIFF_BASIS,        &
                                       CELLS
use base_mesh_config_mod,       only : geometry, &
                                       base_mesh_geometry_spherical
use constants_mod,              only : r_def, PI
use coord_transform_mod,        only : xyz2llr, central_angle
use idealised_config_mod,       only : test,                         &
                                       idealised_test_cold_bubble,   &
                                       idealised_test_gaussian_hill, &
                                       idealised_test_cosine_hill,   &
                                       idealised_test_slotted_cylinder
use initial_density_config_mod, only : r1, x1, y1, r2, x2, y2,     &
                                       tracer_max, tracer_background
use kernel_mod,                 only : kernel_type
use planet_config_mod,          only : p_zero, Rd, kappa
use reference_profile_mod,      only : reference_profile

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: initial_rho_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE,  W3),                            &
       arg_type(GH_FIELD*3, GH_READ, W0)                               &
       /)
  type(func_type) :: meta_funcs(2) = (/                                &
       func_type(W3, GH_BASIS),                                        &
       func_type(W0, GH_BASIS, GH_DIFF_BASIS)                          &
       /)
  integer :: iterates_over = CELLS
contains
  procedure, nopass ::initial_rho_code
end type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
interface initial_rho_kernel_type
   module procedure initial_rho_kernel_constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public initial_rho_code
contains

type(initial_rho_kernel_type) function initial_rho_kernel_constructor() result(self)
  return
end function initial_rho_kernel_constructor

!> @brief The subroutine which is called directly by the Psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in] ndf_w3 The number of degrees of freedom per cell
!! @param[in] map_w3 Integer array holding the dofmap for the cell at the base of the column
!! @param[in] w3_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points 
!! @param[inout] rho Real array the data 
!! @param[inout] gq The gaussian quadrature rule 
!! @param[in] ndf_w0 The number of degrees of freedom per cell
!! @param[in] map_w0 Integer array holding the dofmap for the cell at the base of the column
!! @param[in] w0_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[in] w0_diff_basis Real 5-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[inout] chi_1 Real array, the x component of the w0 coordinate field
!! @param[inout] chi_2 Real array, the y component of the w0 coordinate field
!! @param[inout] chi_3 Real array, the z component of the w0 coordinate field
subroutine initial_rho_code(nlayers, rho, chi_1, chi_2, chi_3, &
                            ndf_w3, undf_w3, map_w3, w3_basis, &
                            ndf_w0, undf_w0, map_w0, w0_basis, w0_diff_basis, &
                            nqp_h, nqp_v, wqp_h, wqp_v )

   use matrix_invert_mod,       only : matrix_invert
   use coordinate_jacobian_mod, only : coordinate_jacobian

  ! needs to compute the integral of rho_df * P
  ! P_analytic over a single column

  !Arguments
  integer, intent(in) :: nlayers, ndf_w3, ndf_w0, undf_w3, undf_w0, nqp_h, nqp_v
  integer, dimension(ndf_w3), intent(in) :: map_w3
  integer, dimension(ndf_w0), intent(in) :: map_w0
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v), intent(in)    :: w3_basis
  real(kind=r_def), dimension(undf_w3),              intent(inout) :: rho
  real(kind=r_def), dimension(undf_w0),              intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(3,ndf_w0,nqp_h,nqp_v), intent(in)    :: w0_diff_basis
  real(kind=r_def), dimension(1,ndf_w0,nqp_h,nqp_v), intent(in)    :: w0_basis
  real(kind=r_def), dimension(nqp_h),                intent(in)    :: wqp_h
  real(kind=r_def), dimension(nqp_v),                intent(in)    :: wqp_v

  !Internal variables
  integer               :: df1, df2, k
  integer               :: qp1, qp2

  real(kind=r_def), dimension(ndf_w3)          :: rho_e, rhs_e
  real(kind=r_def), dimension(ndf_w3,ndf_w3)   :: mass_matrix_w3, inv_mass_matrix_w3
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w0)          :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                             :: exner_ref, rho_ref, theta_ref, &
                                                  integrand
  real(kind=r_def)                             :: x(3)
  real(kind=r_def)            :: l, dt
  real(kind=r_def), parameter :: XC = 0.0_r_def, &
                                 XR = 4000.0_r_def, &
                                 ZC_cold = 3000.0_r_def, &
                                 ZR = 2000.0_r_def
  real(kind=r_def)            :: long, lat, radius
  real(kind=r_def)            :: l1, l2
  real(kind=r_def)            :: h1, h2

  ! compute the RHS & LHS integrated over one cell and solve
  do k = 0, nlayers-1
    do df1 = 1, ndf_w0
      chi_1_e(df1) = chi_1( map_w0(df1) + k)
      chi_2_e(df1) = chi_2( map_w0(df1) + k)
      chi_3_e(df1) = chi_3( map_w0(df1) + k)
    end do
    call coordinate_jacobian(ndf_w0, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e, &
                             w0_diff_basis, jac, dj)
! Compute RHS
    do df1 = 1, ndf_w3
      rhs_e(df1) = 0.0_r_def
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h
          x(:) = 0.0_r_def
          do df2 = 1, ndf_w0
            x(1) = x(1) + chi_1_e(df2)*w0_basis(1,df2,qp1,qp2)
            x(2) = x(2) + chi_2_e(df2)*w0_basis(1,df2,qp1,qp2)
            x(3) = x(3) + chi_3_e(df2)*w0_basis(1,df2,qp1,qp2)
          end do
          call reference_profile(exner_ref, rho_ref, theta_ref, x, test)

          if ( geometry == base_mesh_geometry_spherical ) then
            call xyz2llr(x(1),x(2),x(3),long,lat,radius)
            call central_angle(long,lat,x1,y1,l1)
            call central_angle(long,lat,x2,y2,l2)
          else
            long = x(1)
            lat = x(2)
            l1 = sqrt((long-x1)**2 + (lat-y1)**2)
            l2 = sqrt((long-x2)**2 + (lat-y2)**2)
          end if

          select case( test )
            case( idealised_test_cold_bubble )
              l = sqrt( ((x(1)-XC)/XR)**2 + ((x(3)-ZC_cold)/ZR)**2 )
              if ( l <= 1.0_r_def ) then
                dt =  15.0_r_def/2.0_r_def*(cos(PI*l)+1.0_r_def)
                theta_ref = theta_ref - dt/exner_ref
                rho_ref = p_zero/(Rd*theta_ref) * exner_ref**( (1.0_r_def - kappa )/ kappa )
              end if
            case( idealised_test_GAUSSIAN_HILL )
              h1 = tracer_max*exp( -(l1/r1)**2 )
              h2 = tracer_max*exp( -(l2/r2)**2 )
              rho_ref = h1 +h2
            case( idealised_test_cosine_hill )
              if ( l1 < r1 ) then
                h1 = (tracer_max/2.0_r_def)*(1.0_r_def+cos((l1/r1)*PI))
              else
                h1 = tracer_background
              end if
              if (l2 < r2) then
                h2 = (tracer_max/2.0_r_def)*(1.0_r_def+cos((l2/r2)*PI))
              else
                h2 = tracer_background
              end if
              rho_ref = h1+h2
            case( idealised_test_slotted_cylinder )
              ! Cylinder 1
              if ( l1 < r1 ) then
                if (abs(long-x1) > r1/6.0_r_def) then
                  h1 = tracer_max
                else
                  if (lat < y1-r1*5.0_r_def/12.0_r_def) then
                    h1 = tracer_max
                  else
                    h1 = tracer_background
                  end if
                end if
              else
                h1 = tracer_background
              end if
              ! Cylinder 2
              if ( l2 < r2 ) then
                if (abs(long-x2) > r2/6.0_r_def) then
                  h2 = tracer_max
                else
                  if (lat > y2+r2*5.0_r_def/12.0_r_def) then
                    h2 = tracer_max
                  else
                    h2 = tracer_background
                  end if
                end if
              else
                h2 = tracer_background
              end if
              rho_ref = h1 + h2
          end select

          integrand =  w3_basis(1,df1,qp1,qp2) * rho_ref * dj(qp1,qp2)
          rhs_e(df1) = rhs_e(df1) + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
! Commpute LHS
    do df1 = 1, ndf_w3
       do df2 = 1, ndf_w3
          mass_matrix_w3(df1,df2) = 0.0_r_def
          do qp2 = 1, nqp_v
             do qp1 = 1, nqp_h
                integrand =  w3_basis(1,df1,qp1,qp2) * &
                             w3_basis(1,df2,qp1,qp2) * dj(qp1,qp2)
                 mass_matrix_w3(df1,df2) = mass_matrix_w3(df1,df2) &
                                         + wqp_h(qp1)*wqp_v(qp2)*integrand
             end do
          end do
       end do
    end do
! Solve
    call matrix_invert(mass_matrix_w3,inv_mass_matrix_w3,ndf_w3)
    rho_e = matmul(inv_mass_matrix_w3,rhs_e)
    do df1 = 1,ndf_w3
      rho(map_w3(df1)+k) = rho_e(df1)
    end do
  end do

end subroutine initial_rho_code

end module initial_rho_kernel_mod
