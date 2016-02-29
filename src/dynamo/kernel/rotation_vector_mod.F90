!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel computes the omega (rotation) for the u-equation.

!> @detail The kernel computes the element rotation (omega) vector in W2 space for 
!> the RHS of momentum equation on both F-PLANE and the SPHERE

module rotation_vector_mod

use constants_mod,     only: r_def
use planet_config_mod, only: scaled_omega

implicit none

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the element rotation (omega) vector in W2 space for the
!! RHS of momentum equation on the F-PLANE. 
!! @param[in] ngp_h          Integer. The number of quadrature points in horizontal direction
!! @param[in] ngp_v          Integer. The number of quadrature points in vertical direction
!! @param[in] omegaf         Real scalar. Scaled value of Earth rotation.
!! @param[in] latitude       Real scalar. Pre-set value for the f-plane.
!! @param[out] rotation_vec  Real 3-dim array. Holds the values of the rotation vector on quadrature points
subroutine rotation_vector_fplane(ngp_h, ngp_v, omegaf, latitude, rotation_vec)
!-------------------------------------------------------------------------------
! Compute the rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) on quadrature points 
!-------------------------------------------------------------------------------

integer,          intent(in)  :: ngp_h, ngp_v
real(kind=r_def), intent(in)  :: omegaf, latitude
real(kind=r_def), intent(out) :: rotation_vec(3,ngp_h,ngp_v)

integer :: i, j

rotation_vec = 0.0_r_def

do j = 1, ngp_v
  do i = 1, ngp_h 
     rotation_vec(1,i,j) = 0.0_r_def
     rotation_vec(2,i,j) = 2.0_r_def*omegaf*cos(latitude)
     rotation_vec(3,i,j) = 2.0_r_def*omegaf*sin(latitude)
  end do
end do

end subroutine rotation_vector_fplane

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> Subroutine Computes the element rotation (omega) vector in W2 space for the
!! RHS of momentum equation on the SPHERE. The coordinate field is in W0.
!! @param[in] ndf            Integer. The size of the chi arrays
!! @param[in] ngp_h          Integer. The number of quadrature points in horizontal direction
!! @param[in] ngp_v          Integer. The number of quadrature points in vertical direction
!! @param[in] chi_1          Real 1-dim array. Holds the chi_1 coordinate field
!! @param[in] chi_2          Real 1-dim array. Holds the chi_2 coordinate field
!! @param[in] chi_3          Real 1-dim array. Holds the chi_3 coordinate field
!! @param[in] v0_basis       Real 5-dim array. Holds the W0 basis functions
!! @param[out] rotation_vec  Real 3-dim array. Holds the values of the rotation vector on quadrature points
subroutine rotation_vector_sphere(ndf, ngp_h, ngp_v, chi_1, chi_2, chi_3, v0_basis, rotation_vec)
!-------------------------------------------------------------------------------
! Compute the rotation vector Omega = (0, 2*cos(lat), 2*sin(lat)) on quadrature points 
!-------------------------------------------------------------------------------

use coord_transform_mod, only: xyz2llr, sphere2cart_vector

integer,          intent(in)  :: ndf, ngp_h, ngp_v
real(kind=r_def), intent(in)  :: chi_1(ndf), chi_2(ndf), chi_3(ndf)
real(kind=r_def), intent(in), dimension(1,ndf,ngp_h,ngp_v) :: v0_basis 
real(kind=r_def), intent(out) :: rotation_vec(3,ngp_h,ngp_v)

integer :: i, j, df
real(kind=r_def) :: x, y, z, lat, long, r

lat = 0.0_r_def
long = 0.0_r_def
rotation_vec = 0.0_r_def

do j = 1, ngp_v
  do i = 1, ngp_h  
    x = 0.0_r_def  
    y = 0.0_r_def 
    z = 0.0_r_def   
    do df = 1, ndf
        x = x + chi_1(df)*v0_basis(1,df,i,j)
        y = y + chi_2(df)*v0_basis(1,df,i,j)
        z = z + chi_3(df)*v0_basis(1,df,i,j)        
    end do 
    call xyz2llr(x,y,z,long,lat,r)
    rotation_vec(1,i,j) = 0.0_r_def
    rotation_vec(2,i,j) = 2.0_r_def*scaled_omega*cos(lat)
    rotation_vec(3,i,j) = 2.0_r_def*scaled_omega*sin(lat)

    rotation_vec(:,i,j) = sphere2cart_vector( rotation_vec(:,i,j),(/long, lat, r/) )
  end do
end do

end subroutine rotation_vector_sphere

end module rotation_vector_mod
