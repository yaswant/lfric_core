!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Provides an implementation of the Psy layer

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_mod

  use field_mod,      only : field_type, field_proxy_type 
  use operator_mod,   only : operator_type, operator_proxy_type
  use quadrature_mod, only : quadrature_type
  use constants_mod,  only : r_def

  implicit none
  public

contains

!-------------------------------------------------------------------------------    
!> invoke_inner_prod: Calculate inner product of x and y
  subroutine invoke_inner_prod(x,y,inner_prod)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x,y
    real(kind=r_def),    intent(out) :: inner_prod

    type( field_proxy_type)          ::  x_p,y_p
    integer                          :: i,undf

    x_p = x%get_proxy()
    y_p = y%get_proxy()

    undf = x_p%vspace%get_undf()
    !sanity check
    if(undf /= y_p%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:inner_prod:x and y live on different w-spaces",LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    inner_prod = 0.0_r_def
    do i = 1,undf
      inner_prod = inner_prod + ( x_p%data(i) * y_p%data(i) )
    end do

  end subroutine invoke_inner_prod
  
!-------------------------------------------------------------------------------   
!> invoke_axpy:  (a * x + y) ; a-scalar, x,y-vector     
  subroutine invoke_axpy(scalar,field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpy:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) + field2_proxy%data(i)
    end do
  end subroutine invoke_axpy
  
!-------------------------------------------------------------------------------   
!> invoke_axmy:  (a * x - y) ; a-scalar, x,y-vector
  subroutine invoke_axmy(scalar,field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axmy:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axmy:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = (scalar * field1_proxy%data(i)) - field2_proxy%data(i)
    end do
  end subroutine invoke_axmy
  
!-------------------------------------------------------------------------------   
!> invoke_copy_field_data: Copy the data from one field to another ( a = b )
  subroutine invoke_copy_field_data(field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:copy_field_data:field1 and field_res live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)
    end do
  end subroutine invoke_copy_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_minus_field_data: Subtract values of field2 from values of field 
!> ( c = a - b )
  subroutine invoke_minus_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:minus_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:minus_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) - field2_proxy%data(i)
    end do
  end subroutine invoke_minus_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_plus_field_data:  Add values of field2 to values of field1
!> ( c = a + b )
  subroutine invoke_plus_field_data(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:plus_field_data:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:plus_field_data:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i) + field2_proxy%data(i)
    end do
  end subroutine invoke_plus_field_data
  
!-------------------------------------------------------------------------------   
!> invoke_set_field_scalar: Set all values in a field to a single value
  subroutine invoke_set_field_scalar(scalar, field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(inout ) :: field_res
    real(kind=r_def),   intent(in )    :: scalar
    type( field_proxy_type)            :: field_res_proxy
    integer                            :: i,undf

    field_res_proxy = field_res%get_proxy()

    undf = field_res_proxy%vspace%get_undf()

    do i = 1,undf
      field_res_proxy%data(i) = scalar
    end do
  end subroutine invoke_set_field_scalar

!-------------------------------------------------------------------------------   
!> invoke_divide_field: Divide the values of field1 by field2 and put result in
!>field_res
!> c = a/b
  subroutine invoke_divide_field(field1,field2,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: field1,field2
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy,field2_proxy      &
                                        , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field2_proxy = field2%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field2_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:divide_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:divide_field:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = field1_proxy%data(i)/field2_proxy%data(i)
    end do
  end subroutine invoke_divide_field

!-------------------------------------------------------------------------------   
!> invoke_copy_scaled_field_data: Copy the scaled data from one field to another ( a = scaler*b )
  subroutine invoke_copy_scaled_field_data(scaler,field1,field_res)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    real( kind=r_def ), intent(in)     :: scaler
    type( field_type ), intent(in )    :: field1
    type( field_type ), intent(inout ) :: field_res
    type( field_proxy_type)            :: field1_proxy , field_res_proxy
    integer                            :: i,undf

    field1_proxy = field1%get_proxy()
    field_res_proxy = field_res%get_proxy()

    !sanity check
    undf = field1_proxy%vspace%get_undf()
    if(undf /= field_res_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:copy_scaled_field_data:field1 and field_res live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      field_res_proxy%data(i) = scaler*field1_proxy%data(i)
    end do
  end subroutine invoke_copy_scaled_field_data


!-------------------------------------------------------------------------------   
!> invoke_sum_field: Sum all values of a field x
  subroutine invoke_sum_field( x, field_sum )
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ),  intent(in ) :: x
    real(kind=r_def),    intent(out) :: field_sum

    type( field_proxy_type)          :: x_p
    integer                          :: df, undf

    x_p = x%get_proxy()   

    undf = x_p%vspace%get_undf()
    
    field_sum = 0.0_r_def
    do df = 1,undf
      field_sum = field_sum + x_p%data(df)
    end do

  end subroutine invoke_sum_field

!-------------------------------------------------------------------------------   
!> invoke_axpby:  z = (a * x + b * y) ; a,b-scalar, x,y-vector     
  subroutine invoke_axpby(a,x,b,y,z)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in )    :: x, y
    type( field_type ), intent(inout ) :: z
    real(kind=r_def),   intent(in )    :: a, b
    type( field_proxy_type)            :: x_proxy,y_proxy      &
                                        , z_proxy
    integer                            :: i,undf

    x_proxy = x%get_proxy()
    y_proxy = y%get_proxy()
    z_proxy = z%get_proxy()

    !sanity check
    undf = x_proxy%vspace%get_undf()
    if(undf /= y_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpby:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif
    if(undf /= z_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:axpby:field1 and result_field live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    endif

    do i = 1,undf
      z_proxy%data(i) = (a * x_proxy%data(i)) + (b * y_proxy%data(i))
    end do
  end subroutine invoke_axpby

!-------------------------------------------------------------------------------   
!> invoke_multiply_field: Compute y = a*x for scalar a and fields y and x
  subroutine invoke_multiply_field(a, x, y)
    use log_mod, only : log_event, LOG_LEVEL_ERROR
    implicit none
    type( field_type ), intent(in)    :: x
    type( field_type ), intent(inout) :: y
    real(kind=r_def),   intent(in)    :: a
    type( field_proxy_type)           :: x_proxy, y_proxy
    integer                           :: i,undf

    x_proxy = x%get_proxy()
    y_proxy = y%get_proxy()

    undf = x_proxy%vspace%get_undf()
    if(undf /= y_proxy%vspace%get_undf() ) then
      ! they are not on the same function space
      call log_event("Psy:multiply_field:field1 and field2 live on different w-spaces" &
                    , LOG_LEVEL_ERROR)
      !abort
      stop
    end if

    do i = 1,undf
      y_proxy%data(i) = a*x_proxy%data(i)
    end do
  end subroutine invoke_multiply_field

!-------------------------------------------------------------------------------   
!> invoke_field_delta: Compute delta, a small perturbation to a field
  subroutine invoke_compute_delta(delta, norm, x)
    implicit none
    type( field_type ), intent(in)    :: x
    real(kind=r_def),   intent(in)    :: norm
    real(kind=r_def),   intent(out)   :: delta
    type( field_proxy_type)           :: x_proxy
    integer                           :: i,undf
    real(kind=r_def), parameter       :: delta0 = 1.0e-6_r_def

    x_proxy = x%get_proxy()

    undf = x_proxy%vspace%get_undf()
    
    delta = 0.0_r_def
    do i = 1,undf
      delta = delta + delta0*abs(x_proxy%data(i)) + delta0
    end do
    delta = delta/(real(undf)*norm)
  end subroutine invoke_compute_delta

!-------------------------------------------------------------------------------
!> Non pointwise Kernels


!------------------------------------------------------------------------------- 
!> invoke_sample_flux_kernel: Retrieve values from flux kernel  
  subroutine invoke_sample_flux_kernel(flux, u, multiplicity, q)
    use sample_flux_kernel_mod, only: sample_flux_code
    implicit none

    type(field_type), intent(in)    :: u, multiplicity, q
    type(field_type), intent(inout) :: flux

    type(field_proxy_type)          :: u_p, m_p, q_p, flux_p
    
    integer                 :: cell, nlayers
    integer                 :: ndf_f, ndf_q
    integer                 :: undf_f, undf_q
    integer                 :: dim_q
    integer, pointer        :: map_f(:), map_q(:) => null()

    real(kind=r_def), allocatable  :: basis_q(:,:,:)

    real(kind=r_def), pointer :: nodes(:,:) => null()

    u_p    = u%get_proxy()
    q_p    = q%get_proxy()
    m_p    = multiplicity%get_proxy()
    flux_p = flux%get_proxy()

    nlayers = flux_p%vspace%get_nlayers()

    ndf_f  = flux_p%vspace%get_ndf( )
    undf_f = flux_p%vspace%get_undf()

    ndf_q  = q_p%vspace%get_ndf( )
    dim_q  = q_p%vspace%get_dim_space( )
    undf_q = q_p%vspace%get_undf()
    allocate(basis_q(dim_q, ndf_q, ndf_f))

    nodes => flux_p%vspace%get_nodes( )
    call q_p%vspace%compute_nodal_basis_function(basis_q, ndf_q, ndf_f, nodes)    

    do cell = 1, flux_p%vspace%get_ncell()
       map_f => flux_p%vspace%get_cell_dofmap( cell )
       map_q => q_p%vspace%get_cell_dofmap( cell )
       call sample_flux_code(nlayers, & 
                             flux_p%data, & 
                             m_p%data, &
                             u_p%data, &
                             q_p%data, &
                             ndf_f, & 
                             undf_f, &
                             map_f, &
                             ndf_q, &
                             undf_q, &
                             map_q, &
                             basis_q &
                            )
    end do

  end subroutine invoke_sample_flux_kernel
!-------------------------------------------------------------------------------  
!> invoke_flux_rhs_kernel: Invoke the RHS of the flux equation Flux = u*f
  subroutine invoke_flux_rhs( rhs, u, f, chi, qr )

    use flux_rhs_kernel_mod, only : flux_rhs_code

    type( field_type ),     intent( in ) :: rhs, f, u
    type( field_type ),     intent( in ) :: chi(3) 
    type( quadrature_type), intent( in ) :: qr

    integer          :: cell
    integer          :: ndf_f, undf_f, dim_f, &
                        ndf_u, undf_u, dim_u, &
                        ndf_chi, undf_chi, dim_diff_chi, &
                        nqp_h, nqp_v

    integer, pointer :: map_f(:) => null(), &
                        map_u(:) => null(), & 
                        map_chi(:) => null(), &
                        boundary_dofs(:,:) => null()

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:)   => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    type( field_proxy_type )        :: rhs_proxy, f_proxy, u_proxy
    type( field_proxy_type )        :: chi_proxy(3) 
    
    real(kind=r_def), dimension(:,:,:,:), allocatable :: &
                                  basis_u, &
                                  basis_f, &
                                  diff_basis_chi

    rhs_proxy    = rhs%get_proxy()
    f_proxy      = f%get_proxy()
    u_proxy      = u%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()
    
    boundary_dofs => rhs_proxy%vspace%get_boundary_dofs()

    ndf_u  = rhs_proxy%vspace%get_ndf( )
    undf_u = rhs_proxy%vspace%get_undf( )
    dim_u  = rhs_proxy%vspace%get_dim_space( ) 

    ndf_f  = f_proxy%vspace%get_ndf( )
    undf_f = f_proxy%vspace%get_undf( )
    dim_f  = f_proxy%vspace%get_dim_space( ) 

    ndf_chi  = chi_proxy(1)%vspace%get_ndf( )
    undf_chi = chi_proxy(1)%vspace%get_undf( )
    dim_diff_chi = chi_proxy(1)%vspace%get_dim_space_diff( ) 

    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    allocate( basis_u(dim_u, ndf_u, nqp_h, nqp_v),           &
              basis_f(dim_f, ndf_f, nqp_h, nqp_v),           &
              diff_basis_chi(dim_diff_chi, ndf_chi, nqp_h, nqp_v) )         
    
    call rhs_proxy%vspace%compute_basis_function( &
         basis_u, ndf_u, nqp_h, nqp_v, xp, zp)  
    call f_proxy%vspace%compute_basis_function( &
         basis_f, ndf_f, nqp_h, nqp_v, xp, zp)  
    call chi_proxy(1)%vspace%compute_diff_basis_function( &
         diff_basis_chi, ndf_chi, nqp_h, nqp_v, xp, zp)  
    
    do cell = 1, rhs_proxy%vspace%get_ncell()
       map_f   => f_proxy%vspace%get_cell_dofmap( cell )
       map_u   => rhs_proxy%vspace%get_cell_dofmap( cell )
       map_chi => chi_proxy(1)%vspace%get_cell_dofmap( cell )

      call flux_rhs_code( rhs_proxy%vspace%get_nlayers(), &
                           ndf_u, &
                           undf_u, &
                           map_u, &
                           basis_u, &
                           boundary_dofs, &
                           rhs_proxy%data, &
                           u_proxy%data, &
                           ndf_f, &
                           undf_f, &
                           map_f, &
                           basis_f, &                             
                           f_proxy%data, &
                           ndf_chi, &
                           undf_chi, &
                           map_chi, &
                           diff_basis_chi, &   
                           chi_proxy(1)%data, &
                           chi_proxy(2)%data, &
                           chi_proxy(3)%data, &
                           nqp_h, &
                           nqp_v, &
                           wh, &
                           wv &
                           )           
    end do 
  end subroutine invoke_flux_rhs
 
!-------------------------------------------------------------------------------  
!> invoke_ru_kernel: Invoke the RHS of the u equation
  subroutine invoke_linear_ru_kernel( r_u, u, rho, theta, phi, chi, qr )

    use linear_ru_kernel_mod, only : linear_ru_code

    type( field_type ), intent( in ) :: r_u, u, rho, theta, phi
    type( field_type ), intent( in ) :: chi(3)
    type( quadrature_type), intent( in ) :: qr

    integer                 :: cell, nlayers, nqp_h, nqp_v
    integer                 :: ndf_w0, ndf_w2, ndf_w3
    integer                 :: undf_w0, undf_w2, undf_w3
    integer                 :: dim_w0, diff_dim_w0, dim_w2, diff_dim_w2,dim_w3
    integer, pointer        :: map_w3(:), map_w2(:), map_w0(:) => null()
    integer, pointer        :: boundary_dofs(:,:) => null()

    type( field_proxy_type )        :: r_u_proxy, u_proxy, rho_proxy, theta_proxy, phi_proxy
    type( field_proxy_type )        :: chi_proxy(3)
    
    real(kind=r_def), allocatable  :: basis_w3(:,:,:,:), &
                                      basis_w2(:,:,:,:), &
                                      basis_w0(:,:,:,:), &
                                      diff_basis_w0(:,:,:,:), &
                                      diff_basis_w2(:,:,:,:) 

    real(kind=r_def), pointer :: xp(:,:) => null()
    real(kind=r_def), pointer :: zp(:) => null()
    real(kind=r_def), pointer :: wh(:), wv(:) => null()

    r_u_proxy    = r_u%get_proxy()
    u_proxy      = u%get_proxy()
    rho_proxy    = rho%get_proxy()
    theta_proxy  = theta%get_proxy()
    phi_proxy    = phi%get_proxy()
    chi_proxy(1) = chi(1)%get_proxy()
    chi_proxy(2) = chi(2)%get_proxy()
    chi_proxy(3) = chi(3)%get_proxy()

    boundary_dofs => r_u_proxy%vspace%get_boundary_dofs()

    nlayers = rho_proxy%vspace%get_nlayers()
    nqp_h=qr%get_nqp_h()
    nqp_v=qr%get_nqp_v()
    zp=>qr%get_xqp_v()
    xp=>qr%get_xqp_h()
    wh=>qr%get_wqp_h()
    wv=>qr%get_wqp_v()

    ndf_w3  = rho_proxy%vspace%get_ndf( )
    dim_w3  = rho_proxy%vspace%get_dim_space( )
    undf_w3 = rho_proxy%vspace%get_undf()
    allocate(basis_w3(dim_w3,ndf_w3,nqp_h,nqp_v))

    ndf_w2      = r_u_proxy%vspace%get_ndf( )
    dim_w2      = r_u_proxy%vspace%get_dim_space( )
    diff_dim_w2 = r_u_proxy%vspace%get_dim_space_diff( )
    undf_w2     = r_u_proxy%vspace%get_undf()
    allocate(basis_w2(dim_w2,ndf_w2,nqp_h,nqp_v))
    allocate(diff_basis_w2(diff_dim_w2,ndf_w2,nqp_h,nqp_v))

    ndf_w0      = theta_proxy%vspace%get_ndf( )
    dim_w0      = theta_proxy%vspace%get_dim_space( )
    diff_dim_w0 = theta_proxy%vspace%get_dim_space_diff( )
    undf_w0     = theta_proxy%vspace%get_undf()
    allocate(basis_w0(dim_w0,ndf_w0,nqp_h,nqp_v))
    allocate(diff_basis_w0(diff_dim_w0,ndf_w0,nqp_h,nqp_v))

    call rho_proxy%vspace%compute_basis_function(basis_w3, ndf_w3,         & 
                                                   nqp_h, nqp_v, xp, zp)    

    call r_u_proxy%vspace%compute_basis_function(basis_w2, ndf_w2,         & 
                                                   nqp_h, nqp_v, xp, zp)    

    call r_u_proxy%vspace%compute_diff_basis_function(                     &
         diff_basis_w2, ndf_w2, nqp_h, nqp_v, xp, zp)

    call theta_proxy%vspace%compute_basis_function(basis_w0, ndf_w0,      & 
                                                   nqp_h, nqp_v, xp, zp)    

    call theta_proxy%vspace%compute_diff_basis_function(                  &
         diff_basis_w0, ndf_w0, nqp_h, nqp_v, xp, zp)


    
    do cell = 1, r_u_proxy%vspace%get_ncell()

       map_w3 => rho_proxy%vspace%get_cell_dofmap( cell )
       map_w2 => r_u_proxy%vspace%get_cell_dofmap( cell )
       map_w0 => theta_proxy%vspace%get_cell_dofmap( cell )


       call linear_ru_code( nlayers,                                      &
                            ndf_w2, undf_w2,                              &
                            map_w2, basis_w2, diff_basis_w2,              &
                            boundary_dofs,                                &
                            r_u_proxy%data,                               &
                            u_proxy%data,                                 &
                            ndf_w3, undf_w3,                              &
                            map_w3, basis_w3,                             &
                            rho_proxy%data,                               &
                            ndf_w0, undf_w0,                              &
                            map_w0, basis_w0, diff_basis_w0,              &   
                            theta_proxy%data,                             &
                            phi_proxy%data,                               &
                            chi_proxy(1)%data,                            &
                            chi_proxy(2)%data,                            &
                            chi_proxy(3)%data,                            &
                            nqp_h, nqp_v, wh, wv                          &
                            )           
    end do

    deallocate(basis_w3, basis_w2, diff_basis_w2, basis_w0, diff_basis_w0)
    
  end subroutine invoke_linear_ru_kernel
 

end module psykal_lite_mod
