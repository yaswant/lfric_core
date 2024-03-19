!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
module psykal_builtin_light_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64
  use constants_mod, only : i_def, i_long, r_def
  use field_mod,     only : field_type, field_proxy_type

  implicit none

  public

contains

  !----------------------------------------------------------------------------
  subroutine invoke_convert_cart2sphere_vector( field, coords)
    use coord_transform_mod, only: cart2sphere_scalar
    implicit none
    type(field_type), intent(inout) :: field(3)
    type(field_type), intent(in)    :: coords(3)

    type(field_proxy_type) :: f_p(3), x_p(3)

    integer :: i, df, undf

    do i = 1,3
      f_p(i) = field(i)%get_proxy()
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = f_p(1)%vspace%get_last_dof_annexed()

!Please see PSyclone issues #1351 regarding this implementation
!$omp parallel default(none)                                                   &
!$omp private(df)                                                              &
!$omp shared(undf,f_p,x_p)
!$omp do schedule(static)
    do df = 1, undf
        call cart2sphere_scalar(                                               &
           x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df) ,                 &
           f_p(1)%data(df), f_p(2)%data(df), f_p(3)%data(df) )
    end do
!$omp end do
!$omp end parallel

    call f_p(1)%set_dirty()
    call f_p(2)%set_dirty()
    call f_p(3)%set_dirty()

  end subroutine invoke_convert_cart2sphere_vector

  !----------------------------------------------------------------------------
  subroutine invoke_pointwise_convert_xyz2llr( coords)
    use coord_transform_mod, only: xyz2llr
    implicit none
    type(field_type), intent(inout) :: coords(3)

    type(field_proxy_type) :: x_p(3)

    integer :: i, df, undf
    real(kind=r_def) :: llr(3)

    do i = 1,3
      x_p(i) = coords(i)%get_proxy()
    end do

    undf = x_p(1)%vspace%get_last_dof_annexed()

    do df = 1, undf
      call xyz2llr(x_p(1)%data(df), x_p(2)%data(df), x_p(3)%data(df), &
                   llr(1), llr(2), llr(3))
      x_p(1)%data(df) = llr(1)
      x_p(2)%data(df) = llr(2)
      x_p(3)%data(df) = llr(3)
    end do

    call x_p(1)%set_dirty()
    call x_p(2)%set_dirty()
    call x_p(3)%set_dirty()

  end subroutine invoke_pointwise_convert_xyz2llr

  subroutine invoke_r32_field_min_max(field_min_norm, &
                                  field_max_norm, &
                                  r32_field)

    use scalar_r32_mod,     only: scalar_r32_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_r32_mod,      only: field_r32_type, field_r32_proxy_type

    implicit none

    real(kind=real32),              intent(out)  :: field_min_norm
    real(kind=real32),              intent(out)  :: field_max_norm
    type(field_r32_type),            intent(in)  :: r32_field
    type(scalar_r32_type)                        :: global_min, global_max
    integer(kind=i_def)                          :: df
    real(kind=real32), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real32), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real32)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start, loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_r32_proxy_type)                   :: field_proxy
    integer(kind=i_def)                          :: max_halo_depth_mesh
    type(mesh_type), pointer                     :: mesh => null()
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = r32_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Create a mesh object
    !
    mesh => field_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    ALLOCATE (l_field_min_norm(nthreads))
    ALLOCATE (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    DO df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    END DO
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    DO th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    END DO
    DEALLOCATE (l_field_min_norm, l_field_max_norm)
    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_r32_field_min_max

  !-------------------------------------------------------------------------------
  subroutine invoke_r64_field_min_max(field_min_norm, &
                                  field_max_norm, &
                                  r64_field)

    use scalar_r64_mod,     only: scalar_r64_type
    use omp_lib,            only: omp_get_thread_num
    use omp_lib,            only: omp_get_max_threads
    use mesh_mod,           only: mesh_type
    use field_r64_mod,      only: field_r64_type, field_r64_proxy_type

    implicit none

    real(kind=real64),               intent(out) :: field_min_norm
    real(kind=real64),               intent(out) :: field_max_norm
    type(field_r64_type),             intent(in) :: r64_field
    type(scalar_r64_type)                        :: global_min, global_max
    integer(kind=i_def)                          :: df
    real(kind=real64), allocatable, dimension(:) :: l_field_min_norm
    real(kind=real64), allocatable, dimension(:) :: l_field_max_norm
    real(kind=real64)                            :: minv, maxv
    integer(kind=i_def)                          :: th_idx
    integer(kind=i_def)                          :: loop0_start, loop0_stop
    integer(kind=i_def)                          :: nthreads
    type(field_r64_proxy_type)                   :: field_proxy
    integer(kind=i_def)                          :: max_halo_depth_mesh
    type(mesh_type), pointer                     :: mesh => null()
    !
    ! Determine the number of OpenMP threads
    !
    nthreads = omp_get_max_threads()
    !
    ! Initialise field and/or operator proxies
    !
    field_proxy = r64_field%get_proxy()
    maxv = huge(maxv)
    minv = -huge(minv)
    !
    ! Create a mesh object
    !
    mesh => field_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()
    !
    ! Set-up all of the loop bounds
    !
    loop0_start = 1
    loop0_stop = field_proxy%vspace%get_last_dof_owned()
    !
    ! Call kernels and communication routines
    !
    ALLOCATE (l_field_min_norm(nthreads))
    ALLOCATE (l_field_max_norm(nthreads))
    !
    l_field_min_norm(:) = maxv
    l_field_max_norm(:) = minv
    !
    !$omp parallel default(shared), private(df,th_idx)
    th_idx = omp_get_thread_num()+1
    !$omp do schedule(static)
    DO df=loop0_start,loop0_stop
      l_field_min_norm(th_idx) = min(l_field_min_norm(th_idx), &
                                 field_proxy%data(df))
      l_field_max_norm(th_idx) = max(l_field_max_norm(th_idx), &
                                 field_proxy%data(df))
    END DO
    !$omp end do
    !$omp end parallel
    !
    ! Find minimum in the partial results sequentially
    !
    field_min_norm = l_field_min_norm(1)
    field_max_norm = l_field_max_norm(1)
    DO th_idx=2,nthreads
      field_min_norm = min(field_min_norm, l_field_min_norm(th_idx))
      field_max_norm = max(field_max_norm, l_field_max_norm(th_idx))
    END DO
    DEALLOCATE (l_field_min_norm, l_field_max_norm)
    global_min%value = field_min_norm
    global_max%value = field_max_norm
    field_min_norm = global_min%get_min()
    field_max_norm = global_max%get_max()
    !
  end subroutine invoke_r64_field_min_max

  !--------------------------------------------------------------------------

end module psykal_builtin_light_mod
