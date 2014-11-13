!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Contains methods and algorithms for solving a system A.x = b for known
!! input field b and matrix A and returns field x
!!
!! @details Contains a selction of solvers for inverting the matrix vector
!! system A.x = b to return x = A^{-1}.b Depending upom the type of system to
!! solve a number of iterative solver algorithms are possible or for
!! discontinuous systems an exact solver can be used 
module solver_mod
  use constants_mod,           only : r_def, str_def, max_iter, solver_tol, &
                                      cg_solver, bicg_solver, jacobi_solver, &
                                      gmres_solver, gcr_solver, no_pre_cond
  use log_mod,                 only : log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
                                      LOG_LEVEL_DEBUG, log_scratch_space
  use field_mod,               only : field_type
  use function_space_mod,      only : function_space_type
  use gaussian_quadrature_mod, only : gaussian_quadrature_type, GQ3

  use psy,             only : invoke_inner_prod,                               &
                              invoke_axpy, invoke_minus_field_data,            &
                              invoke_copy_field_data, invoke_set_field_scalar, &
                              invoke_matrix_vector_w0,                         &
                              invoke_matrix_vector_w1,                         &
                              invoke_matrix_vector_w2,                         &
                              invoke_w3_solver_kernel,                         &
                              invoke_divide_field,                             &
                              invoke_copy_scaled_field_data
  use argument_mod,    only : w0, w1, w2, w3

  implicit none
  private

  public :: solver_algorithm

contains

!> @brief Wrapper for specific solver routines for solving system A.x = b
!! @details solves A.x = b for using a choice of solver where A is a mass
!! matrix for a given space and x and b are fields belonging to that space.
!! For a discontinous space the element mass matrix is exactly inverted, for 
!! continuous spaces an iterative solver is used. 
!! The current choices of iteratives solver are:
!! cg: Conjugate gradient method without preconditioning
!! bicgstab: bi-conjugate gradient, stabilised without preconditioning
!! jacobi: a fixed number of jacobi iterations
!> @param[inout] lhs The field to be solved for (x)
!> @param[inout] rhs The right hand side field (b)
!> @param[in]    chi The coordinate array fields
!> @param[in]    space The function space that lhs and rhs are defined on
!> @param[in]    solver_type (optional) The type of iterative solver to use for 
!>               continuous systems
  subroutine solver_algorithm(lhs, rhs, chi, space, solver_type)
    implicit none
    type(field_type), intent(inout)    :: lhs
    type(field_type), intent(inout)    :: rhs
    type(field_type), intent(in)       :: chi(3)
    integer, intent(in)                :: space
    integer, intent(in)                :: solver_type
    integer, parameter                 :: num_jacobi_iters = 5

    select case ( space )
      case ( w3 ) 
        call invoke_w3_solver_kernel(lhs, rhs, chi)
      case ( w0, w1, w2 )
        select case ( solver_type )
          case ( cg_solver )
            call cg_solver_algorithm(lhs, rhs, space)
          case ( bicg_solver ) 
            call bicg_solver_algorithm(lhs, rhs, space)
          case ( jacobi_solver ) 
            call jacobi_solver_algorithm(lhs, rhs, space, num_jacobi_iters)
          case ( gmres_solver )
            call gmres_solver_algorithm(lhs, rhs, space) 
          case ( gcr_solver )
            call gcr_solver_algorithm(lhs, rhs, space)        
          case default
            write( log_scratch_space, '(A)' )  'Invalid linear solver choice, stopping'
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
      case default
       write( log_scratch_space, '(A)' )  'Invalid space for linear solver, stopping'
       call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end subroutine solver_algorithm

!> @brief BiCGStab solver with no preconditioning. 
!! @details solves A.x = b where the operation A.x is encoded in a kernel using
!! the stabilised bi-conjugate gradient method. The choice of matrix is 
!! encoded in the matrix vector kernel that is called
!! @param[in]  rhs_field The input b
!! @param[inout] lhs_field The answer, x
!! @param[in] space The function space lhs and rhs are defined on
  subroutine bicg_solver_algorithm(lhs, rhs, space)
    implicit none
    type(field_type), intent(inout)    :: lhs
    type(field_type), intent(in)       :: rhs
    integer, intent(in)                :: space

    character(len=str_def)             :: cmessage
    ! The temporary fields
    type(field_type)                   :: res, cr, p, v, s, t, cs

    ! the scalars 
    ! the greeks - standard BiCGStab
    real(kind=r_def)                   :: rho,alpha,omega,beta, norm
    real(kind=r_def)                   :: ts,tt
    ! others
    real(kind=r_def)                   :: err,sc_err, init_err
    integer                            :: iter
    integer                            :: rhs_fs, rhs_gq
    type(function_space_type)          :: fs
    type( gaussian_quadrature_type )   :: gq 

    ! compute the residual this is a global sum to the PSy ---
    !PSY call invoke ( inner_prod(rhs,rhs,sc_err))
    call invoke_inner_prod(rhs,rhs,sc_err)
    sc_err = max(sqrt(sc_err), 0.1_r_def)
    write(cmessage,'("solver_algorithm: starting ... ||b|| = ",E15.8)') sc_err
    call log_event(trim(cmessage), LOG_LEVEL_INFO)
    !PSY call invoke ( set_field_scalar(0.0_r_def, lhs))
    call invoke_set_field_scalar(0.0_r_def, lhs)

    rhs_fs = rhs%which_function_space()
    rhs_gq = rhs%which_gaussian_quadrature()
    v = field_type(vector_space = fs%get_instance(rhs_fs), &
                   gq = gq%get_instance(rhs_gq) )

    call mat_ax( v, lhs, space )

    !PSY call invoke ( inner_prod(v,v,err))
    call invoke_inner_prod(v,v,err)

    res = field_type(vector_space = fs%get_instance(rhs_fs), &
                     gq = gq%get_instance(rhs_gq) )
    !PSY call invoke ( minus_field_data(rhs,v,res))
    call invoke_minus_field_data(rhs,v,res)

    !PSY call invoke ( inner_prod(res,res,err))
    call invoke_inner_prod(res,res,err)
    err = sqrt(err)/sc_err
    init_err=err
    if (err < solver_tol) then 
      write(cmessage,'("solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') 0,init_err,err
      call log_event(trim(cmessage),LOG_LEVEL_INFO)
      return
   end if  

    alpha  = 1.0_r_def
    omega  = 1.0_r_def
    norm   = 1.0_r_def

    cr = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )
    !PSY call invoke ( copy_field_data(res,cr))
    call invoke_copy_field_data(res,cr)

    p = field_type(vector_space = fs%get_instance(rhs_fs), &
                   gq = gq%get_instance(rhs_gq) )
    !PSY call invoke ( set_field_scalar(0.0_r_def, p))
    call invoke_set_field_scalar(0.0_r_def, p)

    t = field_type(vector_space = fs%get_instance(rhs_fs), &
                   gq = gq%get_instance(rhs_gq) )
    s = field_type(vector_space = fs%get_instance(rhs_fs), &
                   gq = gq%get_instance(rhs_gq) )
    cs = field_type(vector_space = fs%get_instance(rhs_fs), &
                   gq = gq%get_instance(rhs_gq) )
                   
    !PSY call invoke ( set_field_scalar(0.0_r_def, v))
    call invoke_set_field_scalar(0.0_r_def, v)

    do iter = 1, max_iter

      !PSY call invoke ( inner_prod(cr,res,rho))
      call invoke_inner_prod(cr,res,rho)
      beta = (rho/norm) * (alpha/omega)
      !PSY call invoke ( axpy((-beta*omega),v,res,t))
      call invoke_axpy((-beta*omega),v,res,t)

      call preconditioner( s, t, no_pre_cond )
      !PSY call invoke ( axpy(beta,p,s,p))
      call invoke_axpy(beta,p,s,p)
      call mat_ax( v, p, space )
      !PSY call invoke ( inner_prod(cr,v,norm))
      call invoke_inner_prod(cr,v,norm)
      alpha = rho/norm
      !PSY call invoke ( axpy(-alpha,v,res,s))
      call invoke_axpy(-alpha,v,res,s)

      call preconditioner( cs, s, no_pre_cond )
      call mat_ax(t, cs, space )

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

      if (err < solver_tol) then 
        write(cmessage,'("solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') iter,init_err,err
        call log_event(trim(cmessage),LOG_LEVEL_INFO)
        exit
      end if
    end do

    if(iter >= max_iter) then
      write(cmessage,'("solver_algorithm: NOT converged in", I3," iters, Res=",E15.8)') &
            iter, err
      call log_event(trim(cmessage),LOG_LEVEL_ERROR)
      call log_event(" ... time to flee, bye.",LOG_LEVEL_ERROR)
      stop
    end if

  end subroutine bicg_solver_algorithm
!--------------------------------------------------

!> @brief CG solver for the system A.x = b with no preconditioning. 
!! @details solves A.x = b where the operation A.x is encoded in a kernel using
!! the conjugate gradient method. The choice of matrix is 
!! encoded in the matrix vector kernel that is called. 
!! @param[in] rhs_field The input b
!! @param[inout] lhs_field The answer, x
!! @param[in] space The function space lhs and rhs are defined on
  subroutine cg_solver_algorithm(lhs, rhs, space)
    implicit none
    type(field_type), intent(inout)    :: lhs
    type(field_type), intent(in)       :: rhs
    integer, intent(in)                :: space

    character(len=str_def)             :: cmessage
    ! The temporary fields
    type(field_type)                   :: res, p, Ap

    ! the scalars 
    real(kind=r_def)                   :: alpha, beta
    real(kind=r_def)                   :: rs_new, rs_old
    ! others
    real(kind=r_def)                   :: err,sc_err, init_err
    integer                            :: iter
    integer                            :: rhs_fs, rhs_gq
    type(function_space_type)          :: fs
    type( gaussian_quadrature_type )   :: gq 

    call invoke_inner_prod( rhs, rhs, rs_old )

    ! compute the residual this is a global sum to the PSy ---

    rhs_fs = rhs%which_function_space()
    rhs_gq = rhs%which_gaussian_quadrature()


    res = field_type(vector_space = fs%get_instance(rhs_fs), &
                     gq = gq%get_instance(rhs_gq) )
    p   = field_type(vector_space = fs%get_instance(rhs_fs), &
                     gq = gq%get_instance(rhs_gq) )
    Ap   = field_type(vector_space = fs%get_instance(rhs_fs), &
                     gq = gq%get_instance(rhs_gq) )  

!! First guess: lhs = rhs  
!    !PSY call invoke ( copy_field_data(rhs,lhs))
!    call invoke_copy_field_data(rhs,lhs)
! First guess: lhs = 0  
    !PSY call invoke ( set_field_scalar(0.0_r_def, lhs))
    call invoke_set_field_scalar(0.0_r_def, lhs)

    call mat_ax( Ap, lhs, space )

    !PSY call invoke ( minus_field_data(rhs,Ap,res))
    call invoke_minus_field_data( rhs, Ap, res )
    !PSY call invoke ( copy_field_data(res,p))
    call invoke_copy_field_data( res, p )
    !PSY call invoke ( inner_prod(res,res,rs_old))
    call invoke_inner_prod( res, res, rs_old )

    err = sqrt(rs_old)
    sc_err = max(err, 0.1_r_def)
    init_err=sc_err    
    write(cmessage,'("cg solver_algorithm: starting ... ||b|| = ",E15.8)') sc_err
    call log_event(trim(cmessage),LOG_LEVEL_INFO)
    if (err < solver_tol) then 
      write(cmessage,'("cg solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') 0,init_err,err
      call log_event(trim(cmessage),LOG_LEVEL_INFO)
      return
    end if 
    
    do iter = 1, max_iter
      call mat_ax( Ap, p, space )
      !PSY call invoke ( inner_prod(p,Ap,rs_new))
      call invoke_inner_prod( p, Ap, rs_new )
      alpha = rs_old/rs_new 

      !PSY call invoke ( axpy(alpha,p,lhs,lhs))
      call invoke_axpy( alpha, p, lhs, lhs )
      !PSY call invoke ( axpy(-alpha,Ap,res,res))
      call invoke_axpy( -alpha, Ap, res, res )

      ! check for convergence
      !PSY call invoke ( inner_prod(res,res,err))
      call invoke_inner_prod(res,res,rs_new)
      err = sqrt(rs_new)/sc_err

      write(cmessage,'("cg solver_algorithm[",I2,"]: res = ", E15.8)')        &
            iter, err
      call log_event(trim(cmessage), LOG_LEVEL_DEBUG)
      if (err < solver_tol) then 
        write(cmessage,'("cg solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') iter,init_err,err
        call log_event(trim(cmessage),LOG_LEVEL_INFO)
        exit
      end if

      beta = rs_new/rs_old
      rs_old = rs_new
      !PSY call invoke ( axpy(beta,p,res,p))
      call invoke_axpy(beta,p,res,p)

    end do
    
    if(iter.ge.max_iter) then
      write(cmessage,'("cg solver_algorithm: NOT converged in", I3," iters, Res=",E15.8)') &
      iter, err
      call log_event(trim(cmessage),LOG_LEVEL_ERROR)
      call log_event(" ... time to flee, bye.",LOG_LEVEL_ERROR)
      stop
    end if

  end subroutine cg_solver_algorithm

!--------------------------------------------------

!> @brief Jacobi solver for the system A.x = b. 
!! @details solves A.x = b where the operation A.x is encoded in a kernel using
!! a fixed (n_iter) number of iterations. The choice of matrix is 
!! encoded in the matrix vector kernel that is called. No measure of convergence
!! is used instead the algorithm is assumed to have converged sufficiently
!! after (n_iter) iterations
!! @param[in] rhs_field The input b
!! @param[inout] lhs_field The answser, x
!! @param[in] space The function space lhs and rhs are defined on
!! @param[in] n_iter The number of Jacobi iterations to perform
  subroutine jacobi_solver_algorithm(lhs, rhs, space, n_iter)


  implicit none

  integer,          intent(in)    :: space, n_iter
  type(field_type), intent(inout) :: lhs, rhs
  type(field_type)                :: Ax, lumped_weight, res

  real(kind=r_def), parameter :: mu = 0.9_r_def

  integer :: iter
  integer                            :: rhs_fs, rhs_gq
  type( function_space_type )        :: fs
  type( gaussian_quadrature_type )   :: gq

  rhs_fs = rhs%which_function_space()
  rhs_gq = rhs%which_gaussian_quadrature()

  Ax = field_type(vector_space = fs%get_instance(rhs_fs), &
                  gq = gq%get_instance(rhs_gq) )
  lumped_weight = field_type(vector_space = fs%get_instance(rhs_fs), &
                        gq = gq%get_instance(rhs_gq) )
  res = field_type(vector_space = fs%get_instance(rhs_fs), &
                   gq = gq%get_instance(rhs_gq) )

! Compute mass lump
  !PSY call invoke ( set_field_scalar(1.0_r_def, Ax))
  call invoke_set_field_scalar(1.0_r_def, Ax)
  call mat_ax( lumped_weight, Ax, space )

  !PSY call invoke ( divide_field(rhs, lumped_weight, lhs))
  call invoke_divide_field( rhs, lumped_weight, lhs )

! initial guess
  !PSY call invoke ( set_field_scalar(0.0_r_def, lhs))
  call invoke_set_field_scalar(0.0_r_def, lhs)
!  !PSY call invoke ( copy_field_data(lhs,Ax))
!  call invoke_copy_field_data( lhs, Ax )
  do iter = 1,n_iter
    call mat_ax( Ax, lhs, space )
    !PSY call invoke ( minus_field_data(rhs,Ax,res))
    call invoke_minus_field_data( rhs, Ax, res )
    !PSY call invoke ( divide_field(res, lumped_weight, res))
    call invoke_divide_field( res, lumped_weight, res )
    !PSY call invoke ( axpy(mu,res,lhs,lhs))
    call invoke_axpy( mu, res, lhs, lhs )
  
! Ready for next iteration  
  end do

  end subroutine jacobi_solver_algorithm

!--------------------------------------------------

!> @brief GMRes solver for the system A.x = b. 
!! @details solves A.x = b where the operation A.x is encoded in a kernel using
!! GMRes algorithm. The choice of matrix is 
!! encoded in the matrix vector kernel that is called. No measure of convergence
!! is used instead the algorithm is assumed to have converged sufficiently
!! after (n_iter) iterations
!! @param[in] rhs_field The input b
!! @param[inout] lhs_field The answser, x
!! @param[in] space The function space lhs and rhs are defined on
!! @param[in] n_iter The number of Jacobi iterations to perform
  subroutine gmres_solver_algorithm(lhs, rhs, space)

    use constants_mod, only: gcrk
   
    implicit none
    type(field_type), intent(inout)    :: lhs
    type(field_type), intent(in)       :: rhs
    integer, intent(in)                :: space

    character(len=str_def)             :: cmessage
    ! The temporary fields
    type(field_type)                   :: Ax, r, s, w, v(gcrk) 

    ! the scalars 
    real(kind=r_def)                   :: h(gcrk+1, gcrk), u(gcrk), g(gcrk+1)
    real(kind=r_def)                   :: beta, si, ci, nrm, h1, h2, p, q
    ! others
    real(kind=r_def)                   :: err, sc_err, init_err
    integer                            :: iter, i, j, k, m
    integer                            :: rhs_fs, rhs_gq
    type(function_space_type)          :: fs
    type( gaussian_quadrature_type )   :: gq

    rhs_fs = rhs%which_function_space()
    rhs_gq = rhs%which_gaussian_quadrature()

    Ax = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )
    r  = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )
    s  = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )
    w   = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )
    do iter = 1,gcrk
      v(iter) = field_type(vector_space = fs%get_instance(rhs_fs), &
                           gq = gq%get_instance(rhs_gq) )        
    end do

    !PSY call invoke ( inner_prod(rhs,rhs,err))
    call invoke_inner_prod( rhs, rhs, err )
    sc_err = max( sqrt(err), 0.01_r_def )
    init_err = sc_err

    if (err < solver_tol) then 
      write(cmessage,'("gmres solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') 0,init_err,err
      call log_event(trim(cmessage),LOG_LEVEL_INFO)
      return
    end if

    call mat_ax( Ax, lhs, space )

    !PSY call invoke ( minus_field_data(rhs,Ax,r))
    call invoke_minus_field_data( rhs, Ax, r )

!    call Precon(s,r,preit,precnd)
    !PSY call invoke ( copy_field_data(r,s))
    call invoke_copy_field_data(r, s)

    !PSY call invoke ( inner_prod(s,s,err))
    call invoke_inner_prod( s, s, err )
    beta = sqrt(err)    

    !PSY call invoke ( copy_field_data(s,v(1)))
    call invoke_copy_scaled_field_data( 1.0_r_def/beta, s, v(1) )


    h(:,:) = 0.0_r_def
    g(:)   = 0.0_r_def
    g(1)   = beta

    do iter = 1, max_iter

      do j = 1, GCRk

! This is the correct settings => call Precon(w,v(:,:,j),pstit,pstcnd)
        call preconditioner( w, v(j), no_pre_cond )
        call mat_ax( s, w, space)
! This is the correct settings => call Precon(w,s,preit,precnd)
        call preconditioner( w, s, no_pre_cond )

        do k = 1, j
          !PSY call invoke ( inner_prod(v(k),w,h(k,j)))
          call invoke_inner_prod( v(k), w, h(k,j)  )
          !PSY call invoke ( axpy(-h(k,j),v(k),w,w))
          call invoke_axpy( -h(k,j), v(k), w, w )
        end do
        !PSY call invoke ( inner_prod(w,w,err))
        call invoke_inner_prod( w, w, err  )
        h(j+1,j) = sqrt( err )
        if( j < GCRk ) then
          !PSY call invoke ( copy_scaled_field_data(1.0_r_def/h(j+1,j),w,v(j+1)))
          call invoke_copy_scaled_field_data(1.0_r_def/h(j+1,j), w, v(j+1))       
        end if
      end do
 
! Solve (7.23) of Wesseling (see Saad's book)
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
!  This is the correct settings => call Precon(s,v(:,:,i),pstit,pstcnd)
        call preconditioner( s, v(i), no_pre_cond )
        !PSY call invoke ( axpy(u(i), s, lhs, lhs))
        call invoke_axpy( u(i), s, lhs, lhs )
      end do

! Check for convergence
      call mat_ax( Ax, lhs, space ) 
     !PSY call invoke ( minus_field_data(rhs,Ax,r))
      call invoke_minus_field_data( rhs, Ax, r )    
     
     !PSY call invoke ( inner_prod(r,r,err))
      call invoke_inner_prod(r, r, err )
      beta = sqrt(err)

      err = beta/sc_err
      if( err <  solver_tol ) then
        write(cmessage,'("GMRES solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') iter,init_err,err
        call log_event(trim(cmessage),LOG_LEVEL_INFO)
        exit
      end if

!  This is the correct settings => call Precon(s,r,preit,precnd)
      call preconditioner( s, r, no_pre_cond )
      !PSY call invoke ( copy_scaled_field_data(1.0_r_def/beta,s,v(1)))
      call invoke_copy_scaled_field_data(1.0_r_def/beta, s, v(1))
     
      g(:) = 0.0_r_def
      g(1) = beta

    end do

    if( iter >= max_iter .and. err >  solver_tol ) then
      write(cmessage,'("GMRES solver_algorithm: NOT converged in", I3," iters, Res=",E15.8)') &
      iter, err
      call log_event(trim(cmessage),LOG_LEVEL_ERROR)
      call log_event(" ... time to flee, bye.",LOG_LEVEL_ERROR)
      stop
    end if

  end subroutine gmres_solver_algorithm

!--------------------------------------------------

!> @brief GCR solver for the system A.x = b. 
!! @details solves A.x = b where the operation A.x is encoded in a kernel using
!! the Preconditioned GCR(k) algorithm from Wesseling. The choice of matrix is 
!! encoded in the matrix vector kernel that is called. No measure of convergence
!! is used instead the algorithm is assumed to have converged sufficiently
!! after (n_iter) iterations
!! @param[in] rhs_field The input b
!! @param[inout] lhs_field The answser, x
!! @param[in] space The function space lhs and rhs are defined on
!! @param[in] n_iter The number of Jacobi iterations to perform
  subroutine gcr_solver_algorithm(lhs, rhs, space)

    use constants_mod, only: gcrk
   
    implicit none
    type(field_type), intent(inout)    :: lhs
    type(field_type), intent(in)       :: rhs
    integer, intent(in)                :: space

    character(len=str_def)             :: cmessage
    ! The temporary fields
    type(field_type)                   :: Ax, r, s(gcrk), v(gcrk) 

    ! the scalars 
    real(kind=r_def)                   :: alpha 
    ! others
    real(kind=r_def)                   :: err, sc_err, init_err
    integer                            :: iter, m, n
    integer                            :: rhs_fs, rhs_gq
    type(function_space_type)          :: fs
    type( gaussian_quadrature_type )   :: gq

    rhs_fs = rhs%which_function_space()
    rhs_gq = rhs%which_gaussian_quadrature()

    Ax = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )
    r  = field_type(vector_space = fs%get_instance(rhs_fs), &
                    gq = gq%get_instance(rhs_gq) )

    do iter = 1,gcrk
      s(iter)  = field_type(vector_space = fs%get_instance(rhs_fs), &
                            gq = gq%get_instance(rhs_gq) )
      v(iter)  = field_type(vector_space = fs%get_instance(rhs_fs), &
                            gq = gq%get_instance(rhs_gq) )   
    end do
    !PSY call invoke ( inner_prod(rhs,rhs,err))
    call invoke_inner_prod( rhs, rhs, err )
    sc_err = max( sqrt(err), 0.01_r_def )
    init_err = sc_err

    if (err < solver_tol) then 
      write(cmessage,'("gcr solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') 0,init_err,err
      call log_event(trim(cmessage),LOG_LEVEL_INFO)
      return
    end if

    call mat_ax(Ax, lhs, space)

    !PSY call invoke ( minus_field_data(rhs,Ax,r))
    call invoke_minus_field_data( rhs, Ax, r )

    do iter = 1, max_iter
      do m = 1, GCRk
! This is the correct settings -> call Precon(s(:,:,m),r,prit,prec)
        call preconditioner( s(m), r, no_pre_cond )
        call mat_ax( v(m), s(m), space )

        do n = 1, m-1
          !PSY call invoke ( inner_prod(v(m),v(n),alpha) )
          call invoke_inner_prod( v(m), v(n), alpha )
          !PSY call invoke ( axpy( -alpha, v(n), v(m), v(m))
          call invoke_axpy( -alpha, v(n), v(m), v(m))
          !PSY call invoke ( axpy( -alpha, s(n), s(m), s(m))
          call invoke_axpy( -alpha, s(n), s(m), s(m))
        end do
        !PSY call invoke ( inner_prod(v(m),v(m),err) )
        call invoke_inner_prod( v(m), v(m), err )
        alpha = sqrt(err)
        !PSY call invoke ( copy_scaled_field_data(1.0_r_def/alpha, v(m), v(m)) )
        call invoke_copy_scaled_field_data( 1.0_r_def/alpha, v(m), v(m) ) 
        !PSY call invoke ( copy_scaled_field_data(1.0_r_def/alpha, s(m), s(m)) )
        call invoke_copy_scaled_field_data( 1.0_r_def/alpha, s(m), s(m) )

        !PSY call invoke ( inner_prod(r,v(m),alpha) )
        call invoke_inner_prod( r, v(m), alpha )
        !PSY call invoke ( axpy( alpha, s(m), lhs, lhs)
        call invoke_axpy( alpha, s(m), lhs, lhs )
        !PSY call invoke ( axpy( -alpha, v(m), r, r)
        call invoke_axpy( -alpha, v(m), r, r )
      end do

      !PSY call invoke ( inner_prod(r,r,err) )
      call invoke_inner_prod( r, r, err )
      err = sqrt( err )/sc_err
      if( err <  solver_tol ) then
        write(cmessage,'("GCR solver_algorithm:converged in ", I2," iters, init=",E12.4," final=",E15.8)') iter,init_err,err
        call log_event(trim(cmessage),LOG_LEVEL_INFO)
        exit
      end if
    end do

    if( iter >= max_iter .and. err >  solver_tol ) then
      write(cmessage,'("GCR solver_algorithm: NOT converged in", I3," iters, Res=",E15.8)') &
      iter, err
      call log_event(trim(cmessage),LOG_LEVEL_ERROR)
      call log_event(" ... time to flee, bye.",LOG_LEVEL_ERROR)
      stop
    end if

end subroutine gcr_solver_algorithm

!--------------------------------------------------

!> @brief wrapper for computing  A.x  
!! @details Wrapper routine for calling the appropriate matrix
!! vector routine for computing A*x for matrix A and vector x
!! @param[in]    x The input field
!! @param[inout] Ax The output field
!! @param[in] space The function space that x is in, defines the matrix_vector 
!! routine to use
  subroutine mat_Ax(Ax, x, space) 
   
    implicit none
    type(field_type), intent(inout) :: Ax
    type(field_type), intent(in)    :: x
    integer,          intent(in)    :: space

    !PSY call invoke ( set_field_scalar(0.0_r_def, Ax))
    call invoke_set_field_scalar(0.0_r_def, Ax)

    select case ( space )
      case ( w0 )
        call invoke_matrix_vector_w0( Ax, x )
      case ( w1 )       
        call invoke_matrix_vector_w1( Ax, x )
      case ( w2 )
        call invoke_matrix_vector_w2( Ax, x )
    end select

    return
  end subroutine mat_Ax

!--------------------------------------------------

!> @brief Applies a selected prconditioner to a vector x  
!! @details Applies one of s number of preconditioners to a field x
!! and returns the preconditioned field y. Currently no preconditioner
!! is applied and y = x.
!! @param[in]    x The input field
!! @param[inout] y The output field
!! @param[in] pre_cond_type The type of preconditioner to be used
!! routine to use
  subroutine preconditioner(y, x, pre_cond_type) 
    use constants_mod, only: diagonal_pre_cond
    implicit none
    type(field_type), intent(inout) :: y
    type(field_type), intent(in)    :: x
    integer,          intent(in)    :: pre_cond_type
    character(len=str_def)          :: cmessage

    select case ( pre_cond_type )
      case ( diagonal_pre_cond )
! Diagonal preconditioner
        write(cmessage,'("Diagonal preconditioner not implemented yet")')
        call log_event(trim(cmessage),LOG_LEVEL_ERROR)
      case default
! Default - do nothing
        !PSY call invoke ( copy_field_data(x, y) )
        call invoke_copy_field_data( x, y )   
    end select

    return
  end subroutine preconditioner 
  
end module solver_mod
