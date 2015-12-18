!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief Algorithm to process and dump fields to file
module output_alg_mod
  
  use constants_mod,                     only: r_def, str_max_filename
  use mesh_mod,                          only: mesh_type
  use field_mod,                         only: field_type
  use function_space_mod,                only: function_space_type, W0, W3
  use galerkin_projection_algorithm_mod, only: galerkin_projection_algorithm
  use driver_layer,                      only: interpolated_output
  use quadrature_mod,                    only: quadrature_type, GAUSSIAN
  use operator_mod,                      only: operator_type
  use restart_control_mod,               only: restart_type
  use configuration_mod,                 only: element_order
 
  use psykal_lite_mod,                   only: invoke_set_field_scalar
  implicit none

  private
  public :: output_alg

contains

!> @brief Algorithm to process and dump fields to file
!> @details Projects all fields (or components of vector fields) into a choosen 
!>          scalar space, then samples fields at a given number of points on a 
!>          regular grid and writes them to .m formated files indexed by a 
!>          timestep stamp
!> @param[in] n integer giving the time step index
!> @param[in] rs checkpoint/restart type with timestepping information
!> @param[inout] theta the potential temperature field
!> @param[inout] u the vector wind field
!> @param[inout] rho the density field
!> @param[inout] chi the fem coordinate field array
!> @param[in] mesh  The mesh all fields are on
!> @param[inout] mm_w0 The mass matrix operator for the field to be projected to
  subroutine output_alg(n, rs, theta, u, rho, chi, mesh, mm_w0)

    implicit none
 
    integer,             intent(in)    :: n
    type(field_type),    intent(inout) :: theta, u, rho, chi(3)
    type(mesh_type),     intent(in)    :: mesh
    type(operator_type), intent(inout) :: mm_w0
    type(restart_type),  intent(in)    :: rs

    ! output variables
    integer :: dir
    integer, parameter :: VECTOR_FIELD = 3, &
                          SCALAR_FIELD = 1
    type( field_type ) :: W0_projected_field(3)
    type( field_type ) :: W3_projected_field(1)
    type( quadrature_type )          :: qr
    type( function_space_type )      :: fs
    character(len=str_max_filename)  :: fname

    qr = quadrature_type(element_order+3, GAUSSIAN)

    ! Create fields needed for output (these can be in CG or DG space)
    do dir = 1,3
      W0_projected_field(dir) = field_type( vector_space = fs%get_instance(mesh, W0) )
    end do
    W3_projected_field(1) = field_type( vector_space = fs%get_instance(mesh, W3) )

    call galerkin_projection_algorithm(W0_projected_field(1), theta, mesh, chi, &
                                       SCALAR_FIELD, qr, mm=mm_w0)
    fname=trim(rs%ts_fname("diag_theta",n))//".m"
    call interpolated_output(SCALAR_FIELD, W0_projected_field(1), mesh, chi, &
                             fname)
    call invoke_set_field_scalar(0.0_r_def, W3_projected_field(1)) 
    call galerkin_projection_algorithm(W3_projected_field(1), rho, mesh, chi, &
                                       SCALAR_FIELD, qr)
    fname=trim(rs%ts_fname("diag_rho",n))//".m"
    call interpolated_output(SCALAR_FIELD, W3_projected_field(1), mesh, chi, &
                             fname)
    call galerkin_projection_algorithm(W0_projected_field(:), u, mesh, chi, &
                                       VECTOR_FIELD, qr, mm=mm_w0)
    fname=trim(rs%ts_fname("diag_u",n))//".m"
    call interpolated_output(VECTOR_FIELD, W0_projected_field(:), mesh, chi, &
                             fname)

  end subroutine output_alg

end module output_alg_mod

