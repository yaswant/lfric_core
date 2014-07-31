!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Solver algorithm - Krylov subspace iterative solver for A x = b
!! @details Only BiCGstab so far but can contain other solvers so not typing at the momoent
!! Simply contains a single subroutine which implements BiCGStab for a hard-code matrix-vector
!! (the operation Ax) routine
!! @parameter lhs field_type the solution, x
!! @parameter rhs field_type the right hand side, b

module solver_mod
  use constants_mod,           only : r_def, str_def, max_iter
  use log_mod,                 only : log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
                                      LOG_LEVEL_DEBUG
  use field_mod,               only : field_type
  use function_space_mod,      only : function_space_type
  use gaussian_quadrature_mod, only : gaussian_quadrature_type

  use psy,             only : invoke_inner_prod, invoke_matrix_vector,         &
                              invoke_axpy, invoke_minus_field_data,            &
                              invoke_copy_field_data, invoke_set_field_scalar

  implicit none
  private

  public :: solver_algorithm

contains
!> @brief BiCGStab solver with no preconditioning. 
!! @details solves A.x = b where the operation A.x is encoded in a kernel
!! @param rhs_field_p The input b
!! @param lhs_field_p The answser, x

  subroutine solver_algorithm(lhs, rhs)

    type(field_type), intent(inout)    :: lhs
    type(field_type), intent(in)       :: rhs

    character(len=str_def)             :: cmessage
    ! The temporary fields
    type(field_type)                   :: res, cr, p, v, s, t

    ! the scalars 
    ! the greeks - standard BiCGStab
    real(kind=r_def)                   :: rho,alpha,omega,beta, norm
    real(kind=r_def)                   :: ts,tt
    ! others
    real(kind=r_def)                   :: err,sc_err, tol, init_err
    integer                            :: iter
    integer                            :: rhs_fs
    type(function_space_type)          :: fs
    type( gaussian_quadrature_type )   :: gq 

    tol = 1.0e-8_r_def
    ! compute the residual this is a global sum to the PSy ---
    !PSY call invoke ( inner_prod(rhs,rhs,sc_err))
    call invoke_inner_prod(rhs,rhs,sc_err)
    sc_err = sqrt(sc_err)
    write(cmessage,'("solver_algorithm: starting ... ||b|| = ",E15.8)') sc_err
    call log_event(trim(cmessage), LOG_LEVEL_INFO)
    !PSY call invoke ( set_field_scalar(0.0_r_def, lhs))
    call invoke_set_field_scalar(0.0_r_def, lhs)

    rhs_fs = rhs%which_function_space()
    v = field_type(fs%get_instance(rhs_fs),                                &
         gq%get_instance() )
    !PSY call invoke ( set_field_scalar(0.0_r_def, v))
    call invoke_set_field_scalar(0.0_r_def, v)

    call invoke_matrix_vector(v,lhs)
    !PSY call invoke ( inner_prod(v,v,err))
    call invoke_inner_prod(v,v,err)

    res = field_type(fs%get_instance(rhs_fs) )    
    !PSY call invoke ( minus_field_data(rhs,v,res))
    call invoke_minus_field_data(rhs,v,res)

    !PSY call invoke ( inner_prod(res,res,err))
    call invoke_inner_prod(res,res,err)
    err = sqrt(err)/sc_err
    init_err=err

    alpha  = 1.0_r_def
    omega  = 1.0_r_def
    norm   = 1.0_r_def

    cr = field_type(fs%get_instance(rhs_fs) )
    !PSY call invoke ( copy_field_data(res,cr))
    call invoke_copy_field_data(res,cr)

    p = field_type(fs%get_instance(rhs_fs) )
    !PSY call invoke ( set_field_scalar(0.0_r_def, p))
    call invoke_set_field_scalar(0.0_r_def, p)

    t = field_type(fs%get_instance(rhs_fs),                                &
         gq%get_instance() )
    s = field_type(fs%get_instance(rhs_fs) )
    !PSY call invoke ( set_field_scalar(0.0_r_def, v))
    call invoke_set_field_scalar(0.0_r_def, v)

    do iter = 1, max_iter

      !PSY call invoke ( inner_prod(cr,res,rho))
      call invoke_inner_prod(cr,res,rho)
      beta = (rho/norm) * (alpha/omega)
      !PSY call invoke ( axpy((-beta*omega),v,res,t))
      call invoke_axpy((-beta*omega),v,res,t)

      !      ! this is where the preconitioner would go
      !PSY call invoke ( copy_field_data(t,s))
      call invoke_copy_field_data(t,s)
      !PSY call invoke ( axpy(beta,p,s,p))
      call invoke_axpy(beta,p,s,p)
      !PSY call invoke ( set_field_scalar(0.0_r_def, v))
      call invoke_set_field_scalar(0.0_r_def, v)
      call invoke_matrix_vector(v,p)

      !PSY call invoke ( inner_prod(cr,v,norm))
      call invoke_inner_prod(cr,v,norm)
      alpha = rho/norm
      !PSY call invoke ( axpy(-alpha,v,res,s))
      call invoke_axpy(-alpha,v,res,s)

      !precon cs, s - again no preconditioner
      ! either use a cs or zero t first as its an inc field!
      !PSY call invoke ( set_field_scalar(0.0_r_def, t))
      call invoke_set_field_scalar(0.0_r_def, t)
      call invoke_matrix_vector(t,s)

      !PSY call invoke ( inner_prod(t,t,tt), &
      !PSY               inner_prod(t,s,ts))
      call invoke_inner_prod(t,t,tt)
      call invoke_inner_prod(t,s,ts)

      omega = ts/tt

      !      lhs = lhs + omega * s + alpha * p
      !PSY call invoke ( axpy(omega,s,lhs,lhs)
      !PSY             , axpy(alpha,p,lhs,lhs)
      !PSY             , axpy(-omega,t,s,res) )
      call invoke_axpy(omega,s,lhs,lhs)
      call invoke_axpy(alpha,p,lhs,lhs)
      call invoke_axpy(-omega,t,s,res)
      norm = rho

      ! check for convergence
      !PSY call invoke ( inner_prod(res,res,err))
      call invoke_inner_prod(res,res,err)
      err = sqrt(err)/sc_err

      write(cmessage,'("solver_algorithm[",I2,"]: res = ", E15.8)')        &
            iter, err
      call log_event(trim(cmessage), LOG_LEVEL_DEBUG)

      if (err < tol) then 
        write(cmessage,'("solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') iter,init_err,err
        call log_event(trim(cmessage),LOG_LEVEL_INFO)
        exit
      end if
    end do

    if(iter.ge.max_iter) then
      write(cmessage,'("solver_algortihm: NOT converged in", I3," iters, Res=",E15.8)') &
            iter, err
      call log_event(trim(cmessage),LOG_LEVEL_ERROR)
      call log_event(" ... time to flee, bye.",LOG_LEVEL_ERROR)
      stop
    end if

  end subroutine solver_algorithm

end module solver_mod
