!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> An algorithm for computing the galerkin projection 
!> As a first step fields are projected into a continous space

module galerkin_projection_algorithm_mod

  use log_mod,                 only: log_event, log_scratch_space, LOG_LEVEL_INFO
  use solver_mod,              only: solver_algorithm
  use constants_mod,           only: r_def, solver_option  
  use psy,                     only: invoke_copy_field_data, invoke_set_field_scalar, &
                                     invoke_gp_rhs, invoke_gp_vector_rhs
  use field_mod,               only: field_type
  use gaussian_quadrature_mod, only: gaussian_quadrature_type
  use function_space_mod,      only: function_space_type
  use argument_mod,            only: w0

                
  implicit none

  private 
  public :: galerkin_projection_algorithm

contains
!> @brief An algorithm to compute the galerkin projection of a field
!> @details Computes the Galrkin projection of a field f_in into the space
!>          of field f_out. Solves M*f_out = rhs where rhs = int(gamma*f_in)
!>          and gamma is a test function from the same space as f_out and M is 
!>          the mass matrix for the space of f_out.
!>          If f_in is in a vector space then it is decomposed into orthogonal components 
!>          and the galerkin projection of each component is computed.
!> @param[out] f_out A field to project to
!> @param[in]  f_in  The field to project
!> @param[in]  chi   A 3D coordinate field
!> @param[in]  space_dimension The dimension of the space f_in (scalar or vector)
  subroutine galerkin_projection_algorithm(f_out, f_in, chi, space_dimension) 
    
    implicit none

! dimension of space to project
    integer,            intent(in)    :: space_dimension
! Field to output
    type( field_type ), intent(inout) :: f_out(space_dimension)
! Field to intput
    type( field_type ), intent(inout) :: f_in  
! Coodinate fields
    type( field_type ), intent(inout) :: chi(3)  

    integer                          :: out_fs, out_gq
    type( field_type )               :: rhs(space_dimension)
    integer                          :: dir 
    type(function_space_type)        :: fs
    type( gaussian_quadrature_type ) :: gq

! Create continuous fields to project data into
    out_fs = f_out(1)%which_function_space()
    out_gq = f_out(1)%which_gaussian_quadrature()
    do dir = 1,space_dimension
      rhs(dir) = field_type( vector_space = fs%get_instance(out_fs), &
                             gq = gq%get_instance(out_gq) )
      call invoke_set_field_scalar(0.0_r_def, rhs(dir))
    end do
    
! Project field into continuous space
    if ( f_in%which_function_space() == out_fs ) then
      call invoke_copy_field_data(f_in, f_out(1))
    else    
      write( log_scratch_space, '(A)' ) 'Computing Galerkin projection... '
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      if ( space_dimension == 1 ) then
        call invoke_gp_rhs( rhs(1), f_in, chi )
      else
        call invoke_gp_vector_rhs( rhs, f_in, chi ) 
      end if
      do dir = 1,space_dimension
        call solver_algorithm( f_out(dir), rhs(dir), chi, w0, solver_option)
      end do
    end if

  end subroutine galerkin_projection_algorithm  

end module galerkin_projection_algorithm_mod
