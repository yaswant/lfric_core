!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!>@brief Routines for solving the semi-implicit equation set
module si_solver_alg_mod

  use constants_mod,           only: r_def, str_def
  use field_bundle_mod,        only: clone_bundle, &
                                     set_bundle_scalar, &
                                     bundle_axpy, &
                                     copy_bundle, &
                                     minus_bundle, &
                                     bundle_ax, &
                                     bundle_divide, &
                                     bundle_minmax, &
                                     bundle_inner_product
  use field_mod,               only: field_type
  use formulation_config_mod,  only: newton_krylov
  use lhs_alg_mod,             only: lhs_alg
  use log_mod,                 only: log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_ERROR,   &
                                     LOG_LEVEL_DEBUG,   &
                                     LOG_LEVEL_TRACE,   &
                                     lOG_LEVEL_INFO
  use operator_mod,            only: operator_type
  use rhs_alg_mod,             only: rhs_alg
  use runtime_constants_mod,   only: runtime_constants_type
  use solver_config_mod,       only: maximum_iterations, &
                                     tolerance, &
                                     preconditioner, &
                                     solver_preconditioner_none, &
                                     solver_preconditioner_diagonal, &
                                     gcrk
  use timestepping_config_mod, only: dt

  implicit none
  integer, parameter :: bundle_size = 3
  private
  public :: si_solver_alg
 
contains
!>@brief Setup for the semi-implicit solver, extracts mass matrix diagonals and 
!!         sets up terms for the Newton-Krylov method if needed
!>@details Control routine to solve the system L*dx = rhs using a Krlyov method
!!         sets up terms needed either for a prescribed L or to compute it
!!         using a Newton-Krylov method
!>@param[inout] x0 The state array to solve for the increment of 
!>@param[in]    rhs0 Fixed rhs forcing for the solver
!>@param[in]    x_ref A reference state used for computing a proscribed L
!>@param[in]    runtime_constants Container for various constant objects
  subroutine si_solver_alg(x0, rhs0, x_ref, runtime_constants)
    use rhs_alg_mod,            only: rhs_alg
    use psykal_lite_mod,        only: invoke_compute_delta, &
                                      invoke_set_field_scalar
    use mm_diagonal_kernel_mod, only: mm_diagonal_kernel_type
    implicit none

    type(field_type), intent(inout)          :: x0(bundle_size)
    type(field_type), intent(in)             :: rhs0(bundle_size), x_ref(bundle_size)
    type(runtime_constants_type), intent(in) :: runtime_constants
    
    type(field_type)             :: rhs(bundle_size)
    real(kind=r_def)             :: tau_dt    ! tau_dt would eventually be set globally 
                                              ! (probably the same place as alpha)
    real(kind=r_def)             :: delta_i(bundle_size), delta
    integer                      :: i

    ! Set up tau_dt: to be used here and in subsequent algorithms
    tau_dt = -0.5_r_def*dt

    call clone_bundle(x0, rhs, bundle_size)

    if ( newton_krylov ) then
      call rhs_alg(rhs, tau_dt, x0, runtime_constants, .true.)    
      do i = 1, bundle_size
        ! 1.0 should be norm(dx), but this would be 0, so give a 1/0!
        call invoke_compute_delta(delta_i(i), 1.0_r_def, x0(i))       
      end do
      delta = sum(delta_i)/real(bundle_size, r_def)
    end if

    call mixed_gmres_alg(x0, rhs0, rhs, x_ref, delta, tau_dt, runtime_constants)

  end subroutine si_solver_alg

!=============================================================================!
!>@brief GMRES solver adapted for solving the semi-implicit equations
!>@details Standard GMRES algortihm from "Iterative methods for sparse linear
!! systems" by Y Saad, SIAM 2003
!>@param[inout] x0 State to increment 
!>@param[in]    rhs0 Fixed rhs so solve for
!>@param[in]    rhs forcing computed at x0, used to compute the Jacobian for NK
!!                  method
!>@param[in]    x_ref Reference state
!>@param[in]    delta Small distance used to compute the Jacobian for the NK
!!              method
!>@param[in]    tau_dt The offcentering parameter times the timestep
!>@param[in]    runtime_constants Container for various constant objects
!>@param[in]    mm_diagonal fields containing a diagonal approxiamtion to the
!!                          mass matrices
  subroutine mixed_gmres_alg(x0, rhs0, rhs, x_ref, delta, tau_dt, runtime_constants)
    use psykal_lite_mod, only: invoke_inner_prod
    implicit none

    type(field_type),             intent(inout) :: x0(bundle_size),   &
                                                   rhs(bundle_size)
    type(field_type),             intent(in)    :: rhs0(bundle_size), &
                                                   x_ref(bundle_size)
    type(runtime_constants_type), intent(in)    :: runtime_constants
    real(kind=r_def),             intent(in)    :: delta, tau_dt


    ! The temporary fields
    type(field_type)         :: mm_diagonal(bundle_size)
    type(field_type)         :: dx(bundle_size), Ax(bundle_size), &
                                residual(bundle_size), s(bundle_size), &
                                w(bundle_size), v(bundle_size, gcrk)

    ! the scalars
    real(kind=r_def)         :: h(gcrk+1, gcrk), u(gcrk), g(gcrk+1)
    real(kind=r_def)         :: beta, si, ci, nrm, h1, h2, p, q
    ! others
    real(kind=r_def)         :: err, sc_err, init_err
    integer                  :: iter, i, j, k, m
    integer, parameter       :: MAX_GMRES_ITER = 20

    integer                  :: precon = solver_preconditioner_none
    integer                  :: postcon = solver_preconditioner_diagonal


    mm_diagonal(1) = runtime_constants%get_mass_matrix_diagonal(2)
    mm_diagonal(2) = runtime_constants%get_mass_matrix_diagonal(0)
    mm_diagonal(3) = runtime_constants%get_mass_matrix_diagonal(3)

    call clone_bundle(x0, dx, bundle_size)
    call clone_bundle(x0, Ax, bundle_size)
    call clone_bundle(x0, s, bundle_size)
    call clone_bundle(x0, w, bundle_size)
    call clone_bundle(x0, residual, bundle_size)
    do iter = 1,gcrk
      call clone_bundle(x0, v(:,iter), bundle_size)
    end do

    err = bundle_inner_product(rhs0, rhs0, bundle_size)
    sc_err = max( sqrt(err), 1.0e-5_r_def )
    init_err = sc_err

    if (err < tolerance) then
      write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
           "gmres solver_algorithm:converged in ", 0,           &
           " iters, init=", init_err,                           &
           " final=", err
      call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
      return
    else
      write( log_scratch_space, '(A,I2,A, 2E15.8)' ) "solver_algorithm[", 0, &
                                                    "]: residual = ", init_err
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    end if

    ! Initial guess
    call set_bundle_scalar(0.0_r_def, dx, bundle_size)

    call set_bundle_scalar(0.0_r_def, Ax, bundle_size)

    call minus_bundle( rhs0, Ax, residual, bundle_size )

    call bundle_preconditioner(s, residual, precon, mm_diagonal, bundle_size )

    beta = sqrt(bundle_inner_product(s, s, bundle_size)) 

    call bundle_ax( 1.0_r_def/beta, s, v(:,1), bundle_size )

    h(:,:) = 0.0_r_def
    g(:)   = 0.0_r_def
    g(1)   = beta

    do iter = 1, MAX_GMRES_ITER

      do j = 1, gcrk

        call bundle_preconditioner(w, v(:,j), postcon, mm_diagonal, bundle_size)
        call apply_lhs(s, w, rhs, x0, x_ref, delta, runtime_constants, tau_dt, bundle_size)
        call bundle_preconditioner(w, s, precon, mm_diagonal, bundle_size )
        do k = 1, j
          h(k,j) =  bundle_inner_product( v(:,k), w, bundle_size )
          call bundle_axpy( -h(k,j), v(:,k), w, w, bundle_size )
        end do        
        h(j+1,j) = sqrt( bundle_inner_product( w, w, bundle_size ))
        if( j < gcrk ) then
          call bundle_ax(1.0_r_def/h(j+1,j), w, v(:,j+1), bundle_size)
        end if
      end do

      ! Solve (7.2bundle_size) of Wesseling (see Saad's book)
      do m = 1, gcrk
        nrm    = sqrt( h(m,m)*h(m,m) + h(m+1,m)*h(m+1,m) )
        si     = h(m+1,m)/nrm
        ci     = h(m,m)/nrm
        p      = ci*g(m) + si*g(m+1)
        q      = -si*g(m) + ci*g(m+1)
        g(m)   = p
        g(m+1) = q
        do j = m, gcrk
          h1       = ci*h(m,j)   + si*h(m+1,j)
          h2       =-si*h(m,j)   + ci*h(m+1,j)
          h(m,j)   = h1
          h(m+1,j) = h2
        end do
      end do

      u(gcrk) = g(gcrk)/h(gcrk,gcrk)
      do i = gcrk-1, 1, -1
        u(i) = g(i)
        do j = i+1, gcrk
          u(i) = u(i) - h(i,j)*u(j)
        end do
        u(i) = u(i)/h(i,i)
      end do

      do i = 1, gcrk
        call bundle_preconditioner(s, v(:,i), postcon, mm_diagonal, bundle_size)
        call bundle_axpy( u(i), s, dx, dx, bundle_size )
      end do

      ! Check for convergence
      call apply_lhs(Ax, dx, rhs, x0, x_ref, delta, runtime_constants, tau_dt, bundle_size)

      call minus_bundle( rhs0, Ax, residual, bundle_size )

      beta = sqrt(bundle_inner_product(residual, residual, bundle_size))

      err = beta/sc_err
      write( log_scratch_space, '(A,I2,A, E15.8)' ) "solver_algorithm[", iter, &
                                                    "]: residual = ", err
      call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

      if( err <  tolerance ) then
        write( log_scratch_space, '(A, I2, A, E12.4, A, E15.8)' ) &
             "GMRES solver_algorithm:converged in ", iter,        &
             " iters, init=", init_err,                           &
             " final=", err
        call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
        exit
      end if

      call bundle_preconditioner(s, residual, precon, mm_diagonal, bundle_size)
      call bundle_ax(1.0_r_def/beta, s, v(:,1), bundle_size)

      g(:) = 0.0_r_def
      g(1) = beta

    end do

    if( (iter >= MAX_GMRES_ITER .and. err >  tolerance) .or. isnan(err) ) then
      write( log_scratch_space, '(A, I3, A, E15.8)') &
           "GMRES solver_algorithm: NOT converged in", MAX_GMRES_ITER, " iters, Res=", err
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    ! Add increments to field
    call bundle_axpy(1.0_r_def, dx, x0, x0, bundle_size)  

  end subroutine mixed_gmres_alg
!=============================================================================!
!>@brief Computes the action of the lhs matrix L on a field x
!>@details Using either the Newton-Krylov method or a proscribed lhs computes
!!         the action of a matrix L on a field x without explicitly computing the entries
!!         of L
!>@param[inout] ax The product of L*x
!>@param[in]    x The state array to compute the action on
!>@param[in]    rhs The rhs array used to compute L*x with the Newton-Krylov
!!                  method
!>@param[in]    x0  State to perturb around for the Newton-Krylov computation
!>@param[in]    x_ref Reference state used to linearise the equations for the
!!                    proscribed L
!>@param[in]    delta Small distance used to compute the Jacobian for the NK
!!              method
!>@param[in]    runtime_constants Container for various constant objects
!>@param[in]    tau_dt relaxation factor multiplied by the timestep
!>@param[in]    mm_diagonal fields containing a diagonal approxiamtion to the
!!                          mass matrices
!>@param[in]    bundle_size Number of fields the state arrays
  subroutine apply_lhs(ax, x, rhs, x0, x_ref, delta, runtime_constants, &
                       tau_dt, bundle_size)
    implicit none
    integer,                      intent(in)    :: bundle_size
    type(field_type),             intent(inout) :: ax(bundle_size)
    type(field_type),             intent(inout) :: x(bundle_size), rhs(bundle_size), x0(bundle_size)
    type(field_type),             intent(in)    :: x_ref(bundle_size)
    real(kind=r_def),             intent(in)    :: delta
    type(runtime_constants_type), intent(in)    :: runtime_constants
    real(kind=r_def),             intent(in)    :: tau_dt
    type(field_type)                            :: x_new(bundle_size), rhs_new(bundle_size)

    if ( newton_krylov ) then
      call clone_bundle(x0, x_new, bundle_size)
      call clone_bundle(x0, rhs_new, bundle_size)
      call bundle_axpy(delta, x, x0, x_new, bundle_size)
      call rhs_alg(rhs_new, tau_dt, x_new, runtime_constants, .true.)
      call minus_bundle(rhs_new, rhs, ax, bundle_size)
      call bundle_ax(1.0_r_def/delta, ax, ax, bundle_size)
    else
      call lhs_alg(Ax, tau_dt, x, x_ref, runtime_constants)
    end if
  end subroutine apply_lhs
!=============================================================================!
!>@brief Applies a choosen preconditioner to a state x to produce state y
!>@param[inout] y Preconditioned state
!>@param[in]    x Original state
!>@param[in]    option choice of which preconditioner to use
!>@param[in]    mm Arrays containing diagonal approximation to mass matrices
!>@param[in]    bundle_size Number of fields the state arrays
  subroutine bundle_preconditioner(y, x, option, mm, bundle_size)
    use psykal_lite_mod, only: invoke_copy_field_data
    implicit none
    integer,          intent(in)    :: bundle_size
    type(field_type), intent(inout) :: y(bundle_size)
    type(field_type), intent(in)    :: x(bundle_size)
    type(field_type), optional :: mm(bundle_size)
    integer, intent(in) :: option
    integer :: i

    i = option
    if ( option == solver_preconditioner_none ) then
      do i = 1,bundle_size
        call invoke_copy_field_data( x(i), y(i) )
      end do
    elseif ( option == solver_preconditioner_diagonal ) then
      do i = 1,bundle_size
        call invoke_copy_field_data( x(i), y(i) )
      end do
      call bundle_divide(y, mm, bundle_size)
    end if
  end subroutine bundle_preconditioner
!=============================================================================!
end module si_solver_alg_mod
