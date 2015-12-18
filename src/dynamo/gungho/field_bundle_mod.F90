!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!>@brief Contains wrapper routines for working on arrays of fields
!!       Effectively these are generally wrapper functions to pointwise kernels
module field_bundle_mod
  use field_mod,     only: field_type
  use constants_mod, only: r_def
  implicit none

  contains

!> Create a bundle y of fields on the same function space as bundel x
 subroutine clone_bundle(x, y, mesh, bundle_size)
    use function_space_mod, only: function_space_type
    use mesh_mod,           only: mesh_type
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size)
    type(mesh_type),  intent(in)    :: mesh
    type(field_type), intent(inout) :: y(bundle_size)
    type(function_space_type)       :: fs
    integer :: i

    do i = 1,bundle_size   
      y(i) = field_type( vector_space = &
         fs%get_instance(mesh, x(i)%which_function_space()) )
    end do
  end subroutine clone_bundle
!=============================================================================!

!> Set all fields in a bundle x to the scalar value a
!> @param [in] a The scalar
!> @param [inout] x The field bundle
!> @param [in] bundle_size the number of fields in the bundle x
  subroutine set_bundle_scalar(a, x, bundle_size)
    use psykal_lite_mod, only: invoke_set_field_scalar
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: x(bundle_size)
    real(kind=r_def), intent(in)    :: a
    integer ::i   
 
    do i = 1,bundle_size
      call invoke_set_field_scalar(a, x(i))
    end do
  end subroutine set_bundle_scalar
!=============================================================================!

!> Compute the inner product of two field bundles
!> @param [result] a The inner product
!> @param [in] x The first field bundle
!> @param [in] y The second field bundle
!> @param [in] bundle_size the number of fields in the bundle x
  function bundle_inner_product(x, y, bundle_size) result(a)
    use psykal_lite_mod, only: invoke_inner_prod
    implicit none
    integer,          intent(in) :: bundle_size
    type(field_type), intent(in) :: x(bundle_size), y(bundle_size)
    real(kind=r_def) :: a
    real(kind=r_def) :: b
    integer :: i

    a = 0.0_r_def
    do i = 1,bundle_size
      call invoke_inner_prod(x(i), y(i), b)
      a = a + b
    end do
  end function bundle_inner_product
!=============================================================================!

!> Compute z = a*x + y for bundles x, y, z and scalar a
!> @param [in] a The scalar
!> @param [in] x The first field bundle
!> @param [in] y The second field bundle
!> @param [inout] z The result field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine bundle_axpy(a, x, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_axpy
    implicit none
    integer,          intent(in)    :: bundle_size
    real(kind=r_def), intent(in)    :: a
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_axpy(a, x(i), y(i), z(i))
    end do
  end subroutine bundle_axpy
!=============================================================================!

!> Copy the data from bundle x to bundle y (y = x)
!> @param [in] x The first field bundle
!> @param [inout] y The second field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine copy_bundle(x, y, bundle_size)
    use psykal_lite_mod, only: invoke_copy_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size)
    type(field_type), intent(inout) :: y(bundle_size)  
    integer :: i
    
    do i = 1,bundle_size
      call invoke_copy_field_data(x(i), y(i))
    end do
  end subroutine copy_bundle
!=============================================================================!
!> Compute z = x - y for field bundles x, y and z
!> @param [in] x The first field bundle
!> @param [in] y The second field bundle
!> @param [inout] z The result field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine minus_bundle(x, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_minus_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_minus_field_data(x(i), y(i), z(i))
    end do
  end subroutine minus_bundle
!=============================================================================!
!> Compute y = a*x for bundles x and y and scaler a
!> @param [in] a The scalar
!> @param [in] x The first field bundle
!> @param [inout] y The second field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine bundle_ax(a, x, y, bundle_size)
    use psykal_lite_mod, only: invoke_copy_scaled_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: y(bundle_size)
    type(field_type), intent(in)    :: x(bundle_size)
    real(kind=r_def), intent(in)    :: a
    integer :: i
   
    do i = 1,bundle_size
      call invoke_copy_scaled_field_data( a, x(i), y(i) )
    end do
  end subroutine bundle_ax
!=============================================================================!
!> Divide the dofs in  bundle x by those in bundle y
!> @param [inout] x The first field bundle
!> @param [inout] y The second field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine bundle_divide(x, y, bundle_size)
    use psykal_lite_mod, only: invoke_divide_field
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: x(bundle_size), y(bundle_size)
    integer :: i    

    do i = 1,bundle_size
      call invoke_divide_field(x(i),y(i),x(i))
    end do
  end subroutine bundle_divide
!=============================================================================!
!> Write the min and max values of a field bundle to the log file
!> @param [in] x The first field bundle
!> @param [in] mesh The mesh all the fields in x are defined on
!> @param [in] bundle_size the number of fields in the bundle
  subroutine bundle_minmax(x, mesh, bundle_size)
    use function_space_mod, only: function_space_type
    use mesh_mod,           only: mesh_type
    use psykal_lite_mod,                only: invoke_copy_field_data
    use log_mod,            only: lOG_LEVEL_INFO   
    implicit none
    integer,          intent(in)    :: bundle_size
    type(mesh_type),  intent(in)    :: mesh
    type(function_space_type) :: fs
    type(field_type) :: x(bundle_size)
    type(field_type) :: y
    integer :: i
! This has strange syntax as psykal_lite_modclone doesnt like calls to type bound
! procedures of field arrays
    do i = 1,bundle_size
      y = field_type( vector_space = &
           fs%get_instance(mesh, x(1)%which_function_space()) )
      call invoke_copy_field_data( x(1), y ) 
      call y%log_minmax(LOG_LEVEL_INFO, 'field')
    end do
  end subroutine bundle_minmax

!=============================================================================!
!> Compute z = a*x + b*y for bundles x, y, z and scalars a and b
!> @param [in] a The scalar
!> @param [in] b The scalar
!> @param [in] x The first field bundle
!> @param [in] y The second field bundle
!> @param [inout] z The result field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine bundle_axpby(a, x, b, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_axpby
    implicit none
    integer,          intent(in)    :: bundle_size
    real(kind=r_def), intent(in)    :: a, b
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_axpby(a, x(i), b, y(i), z(i))
    end do
  end subroutine bundle_axpby
!=============================================================================!
!> Compute z = x + y for field bundles x, y and z
!> @param [in] x The first field bundle
!> @param [in] y The second field bundle
!> @param [inout] z The result field bundle
!> @param [in] bundle_size the number of fields in the bundle
  subroutine add_bundle(x, y, z, bundle_size)
    use psykal_lite_mod, only: invoke_axpy
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(in)    :: x(bundle_size), y(bundle_size)
    type(field_type), intent(inout) :: z(bundle_size)
    integer :: i
    
    do i = 1,bundle_size
      call invoke_axpy(1.0_r_def, x(i), y(i), z(i))
    end do
  end subroutine add_bundle
!=============================================================================!
 
end module field_bundle_mod

