!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @mainpage Dynamo
!> Illustration of the PSyKAl (Parallel-system/Kernel/Algorithm) architecture
!> for Gung Ho. Whilst the computational and optimisation infrastructure is
!> being developed, the science code is being developed using
!> a hand-rolled PSy layer, PSy-lite. A PSyKAl-lite needs a dynamo!
!> Eventually, PSyKAl-lite will be replaced with the real PSy and Dynamo
!> will be the implementation of the Gung Ho dynamical core.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Calls various init subroutines to create a mesh, function spaces
!> and then prognostic fields on those function spaces.
!> Controls the main timestepping loop, initialises and calls
!> timestepping algorithms and also handles output of diagnostics at
!> specified intervals as well as checkpoint restart files.

program gungho

  use cli_mod,                        only : get_initial_filename
  use constants_mod,                  only : i_def
  use gungho_mod,                     only : load_configuration
  use ESMF
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use formulation_config_mod,         only : transport_only, use_moisture, use_physics
  use init_mesh_mod,                  only : init_mesh
  use init_dynamo_mod,                only : init_dynamo
  use init_physics_mod,               only : init_physics
  use iter_timestep_alg_mod,          only : iter_alg_init, &
                                             iter_alg_step
  use runge_kutta_init_mod,           only : runge_kutta_init
  use rk_alg_timestep_mod,            only : rk_alg_init, &
                                             rk_alg_step
  use transport_config_mod,           only : scheme, &
                                             transport_scheme_method_of_lines, &
                                             transport_scheme_bip_cosmic, &
                                             transport_scheme_cusph_cosmic
  use rk_transport_mod,               only : rk_transport_init, &
                                             rk_transport_step, &
                                             rk_transport_final
  use cosmic_transport_alg_mod,       only : cosmic_transport_init, &
                                             cosmic_transport_step
  use cusph_cosmic_transport_alg_mod, only : cusph_cosmic_transport_init, &
                                             cusph_cosmic_transport_step
  use conservation_algorithm_mod,     only : conservation_algorithm
  use log_mod,                        only : log_event,         &
                                             log_set_level,     &
                                             log_scratch_space, &
                                             LOG_LEVEL_ERROR,   &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_DEBUG,   &
                                             LOG_LEVEL_TRACE
  use restart_config_mod,             only : restart_filename => filename
  use restart_control_mod,            only : restart_type
  use output_config_mod,              only : diagnostic_frequency, &
                                             subroutine_timers, &
                                             write_nodal_output, &
                                             write_xios_output
  use io_mod,                         only : output_nodal, &
                                             output_xios_nodal, &
                                             xios_domain_init
  use timestepping_config_mod,        only : method, &
                                           timestepping_method_semi_implicit, &
                                           timestepping_method_rk
  use derived_config_mod,             only : set_derived_config
  use checksum_alg_mod,               only : checksum_alg
  use diagnostic_alg_mod,             only : divergence_diagnostic_alg, &
                                             pressure_diagnostic_alg,   &
                                             density_diagnostic_alg
  use mr_indices_mod,                 only : imr_v, imr_c, imr_r, imr_nc, &
                                             imr_nr, nummr
  use runtime_constants_mod,          only : get_cell_orientation, get_detj_at_w2
  use global_mesh_collection_mod,     only : global_mesh_collection, &
                                             global_mesh_collection_type
  use timer_mod,                      only : timer, output_timer

  use xios
  use mpi
  use mod_wait

  implicit none

  character(:), allocatable :: filename

  type(ESMF_VM)      :: vm
  integer(i_def)     :: rc
  integer(i_def)     :: total_ranks, local_rank
  integer(i_def)     :: petCount, localPET, ierr

  integer(i_def)     :: comm = -999

  character(len=*), parameter   :: xios_id   = "lfric_client"
  character(len=*), parameter   :: xios_ctx  = "gungho_atm"

  type(restart_type) :: restart

  integer(i_def)     :: mesh_id

  ! prognostic fields
  type( field_type ) :: u, rho, theta, xi, mr(nummr), rho_in_wth

  ! prognostic fields in non-native spaces (e.g. for physics)
  type( field_type ) :: u1_in_w3, u2_in_w3, u3_in_w3, theta_in_w3, p_in_w3, p_in_wth

  ! Array to hold cell orientation values in W3 field
  type( field_type )               :: cell_orientation
  ! Array to hold detj values at W2 dof locations
  type( field_type )               :: detj_at_w2

  integer(i_def)                   :: timestep, ts_init, dtime

  integer(i_def) :: i

  character(5) :: name

  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise MPI
 
  call mpi_init(ierr)

  ! initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = comm)

  ! Initialise ESMF using mpi communicator initialised by XIOS
  ! and get the rank information from the virtual machine
  call ESMF_Initialize(vm=vm, &
                       defaultlogfilename="gungho.Log", &
                       logkindflag=ESMF_LOGKIND_MULTI, &
                       mpiCommunicator=comm, &
                       rc=rc)

  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                          LOG_LEVEL_ERROR )

  call ESMF_VMGet(vm, localPet=localPET, petCount=petCount, rc=rc)
  if (rc /= ESMF_SUCCESS) &
    call log_event( 'Failed to get the ESMF virtual machine.', LOG_LEVEL_ERROR )

  total_ranks = petCount
  local_rank  = localPET

  ! Currently log_event can only use ESMF so it cannot be used before ESMF is
  ! initialised.
  call log_event( 'gungho running...', LOG_LEVEL_INFO )

  call get_initial_filename( filename )
  call load_configuration( filename )
  call set_derived_config()
  deallocate( filename )

  restart = restart_type( restart_filename, local_rank, total_ranks )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------
  if ( subroutine_timers ) call timer('gungho')

  allocate( global_mesh_collection, &
            source = global_mesh_collection_type() )

  ! Create the mesh and function space collection
  call init_mesh(local_rank, total_ranks, mesh_id)

  ! Full global meshes no longer required, so reclaim
  ! the memory from global_mesh_collection
  write(log_scratch_space,'(A)') &
      "Purging global mesh collection."
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  deallocate(global_mesh_collection)

  ! Create and initialise prognostic fields
  timestep = 0
  call init_dynamo(mesh_id, u, rho, theta, rho_in_wth, mr, xi, restart)

  if (use_physics)then
    call init_physics(mesh_id, u, rho, theta, rho_in_wth,         &
                      u1_in_w3, u2_in_w3, u3_in_w3, theta_in_w3,  &
                      p_in_w3, p_in_wth)
  end if

  if ( transport_only .and. scheme == transport_scheme_cusph_cosmic) then
    cell_orientation = get_cell_orientation()
  end if

  ! Create detj values at W2 dof locations
  detj_at_w2 = get_detj_at_w2()

  !-----------------------------------------------------------------------------
  ! IO init
  !-----------------------------------------------------------------------------

  ! If xios output then set up XIOS domain and context
  if (write_xios_output) then

    dtime = 1

    call xios_domain_init(xios_ctx, comm, dtime, mesh_id, vm, local_rank, total_ranks)

  end if

  ! Initial output
  ! We only want these once at the beginning of a run 
  ts_init = max( (restart%ts_start() - 1), 0 ) ! 0 or t previous.

  if (ts_init == 0) then 

    ! Original nodal output
    if ( write_nodal_output)  then

      call output_nodal('theta', ts_init, theta, mesh_id)
      call output_nodal('xi',    ts_init, xi,    mesh_id)
      call output_nodal('u',     ts_init, u,     mesh_id)
      call output_nodal('rho',   ts_init, rho,   mesh_id)

      if (use_moisture)then
       call output_nodal('m_v',   ts_init, mr(imr_v),   mesh_id)
       call output_nodal('m_c',   ts_init, mr(imr_c),   mesh_id)
       call output_nodal('m_r',   ts_init, mr(imr_r),   mesh_id)
       call output_nodal('m_nc',   ts_init, mr(imr_nc),   mesh_id)
       call output_nodal('m_nr',   ts_init, mr(imr_nr),   mesh_id)
      end if

    end if

    ! XIOS output
    if (write_xios_output) then

      ! Make sure XIOS calendar is set to timestep 1 as it starts there
      ! not timestep 0.
      call xios_update_calendar(1)
 
      ! Output scalar fields
      call output_xios_nodal("init_theta", theta, mesh_id)
      call output_xios_nodal("init_rho", rho, mesh_id)

      ! Output vector fields
      call output_xios_nodal("init_u", u, mesh_id)
      call output_xios_nodal("init_xi", xi, mesh_id)

      if (use_moisture)then
       call output_xios_nodal('init_m_v',  mr(imr_v),   mesh_id)
       call output_xios_nodal('init_m_c',  mr(imr_c),   mesh_id)
       call output_xios_nodal('init_m_r',  mr(imr_r),   mesh_id)
       call output_xios_nodal('init_m_nc', mr(imr_nc),   mesh_id)
       call output_xios_nodal('init_m_nr', mr(imr_nr),   mesh_id)
      end if

    end if

    call pressure_diagnostic_alg(rho, theta, ts_init, mesh_id)
    call divergence_diagnostic_alg(u, ts_init, mesh_id)

  end if
  
  !-----------------------------------------------------------------------------
  ! model step 
  !-----------------------------------------------------------------------------
  do timestep = restart%ts_start(),restart%ts_end()

    ! Update XIOS calendar
    if (write_xios_output) then
      call log_event( "Gungho: Updating XIOS timestep", LOG_LEVEL_INFO )
      call xios_update_calendar(timestep)
    end if

    call log_event( &
    "/****************************************************************************\ ", &
     LOG_LEVEL_TRACE )
     write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
     call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if ( transport_only ) then

      select case( scheme )
        case ( transport_scheme_method_of_lines)
          if (timestep == restart%ts_start()) then 
            ! Initialise and output initial conditions on first timestep
            call runge_kutta_init()
            call rk_transport_init( mesh_id, u, rho, theta)
          end if
          call rk_transport_step( u, rho, theta)
        case ( transport_scheme_bip_cosmic)
          if (timestep == restart%ts_start()) then 
            ! Initialise and output initial conditions on first timestep
            call cosmic_transport_init(mesh_id, u)
          end if
          call cosmic_transport_step(rho,detj_at_w2)
        case ( transport_scheme_cusph_cosmic)
          if (timestep == restart%ts_start()) then 
            ! Initialise and output initial conditions on first timestep
            call cusph_cosmic_transport_init(mesh_id, u)
          end if
          call cusph_cosmic_transport_step(mesh_id, rho, cell_orientation, detj_at_w2)
          call density_diagnostic_alg(rho, timestep)
          call conservation_algorithm(timestep, rho, u, theta, xi)
        case default
         call log_event("Dynamo: Incorrect transport option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR)
         stop
      end select
    else
       select case( method )
         case( timestepping_method_semi_implicit )  ! Semi-Implicit 
           ! Initialise and output initial conditions on first timestep
           if (timestep == restart%ts_start()) then 
             call runge_kutta_init()
             call iter_alg_init(mesh_id, u, rho, theta)
             call log_event( "Dynamo: Outputting initial fields", LOG_LEVEL_INFO )
             call conservation_algorithm(timestep, rho, u, theta, xi)
           end if
           call iter_alg_step(u, rho, theta, mr, xi,           &
                              u1_in_w3, u2_in_w3, u3_in_w3, &
                              rho_in_wth, theta_in_w3,      &
                              p_in_w3, p_in_wth, timestep)

         case( timestepping_method_rk )             ! RK
           ! Initialise and output initial conditions on first timestep
           if (timestep == restart%ts_start()) then
             call runge_kutta_init()
             call rk_alg_init( mesh_id, u, rho, theta)
             call log_event( "Dynamo: Outputting initial fields", LOG_LEVEL_INFO )
             call conservation_algorithm(timestep, rho, u, theta, xi)
           end if
           call rk_alg_step( u, rho, theta, xi)
         case default
           call log_event("Dynamo: Incorrect time stepping option chosen, "// &
                          "stopping program! ",LOG_LEVEL_ERROR)
           stop
       end select

       call conservation_algorithm(timestep, rho, u, theta, xi)

    end if

    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( &
    '\****************************************************************************/ ', &
     LOG_LEVEL_INFO)

    ! Use diagnostic output frequency to determine whether to write diagnostics
    ! on this timestep

    if ( mod(timestep, diagnostic_frequency) == 0 ) then

      ! Original nodal output
      if ( write_nodal_output)  then

         call log_event("Gungho: writing nodal output", LOG_LEVEL_INFO)
         call output_nodal("theta", timestep, theta, mesh_id)
         call output_nodal("rho", timestep, rho, mesh_id)
         call output_nodal("u", timestep, u, mesh_id)
         call output_nodal("xi", timestep, xi, mesh_id)

         if (use_moisture)then
           call output_nodal('m_v',    timestep, mr(imr_v),   mesh_id)
           call output_nodal('m_c',    timestep, mr(imr_c),   mesh_id)
           call output_nodal('m_r',    timestep, mr(imr_r),   mesh_id)
           call output_nodal('m_nc',    timestep, mr(imr_nc),   mesh_id)
           call output_nodal('m_nr',    timestep, mr(imr_nr),   mesh_id)
         end if

      end if

      ! XIOS output
      if (write_xios_output) then
 
         ! Output scalar fields
         call log_event("Gungho: writing xios output", LOG_LEVEL_INFO)
         call output_xios_nodal("theta", theta, mesh_id)
         call output_xios_nodal("rho", rho, mesh_id)

         if (use_moisture)then
           call output_xios_nodal('m_v', mr(imr_v),   mesh_id)
           call output_xios_nodal('m_c', mr(imr_c),   mesh_id)
           call output_xios_nodal('m_r', mr(imr_r),   mesh_id)
           call output_xios_nodal('m_nc', mr(imr_nc),   mesh_id)
           call output_xios_nodal('m_nr', mr(imr_nr),   mesh_id)
         end if
 
         ! Output vector fields
         call output_xios_nodal("u", u, mesh_id)
         call output_xios_nodal("xi", xi, mesh_id)

      end if

      call pressure_diagnostic_alg(rho, theta, timestep, mesh_id)
      call divergence_diagnostic_alg(u, timestep, mesh_id)
    end if

  end do ! end ts loop

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------

  ! Log fields
  call rho%log_field(   LOG_LEVEL_DEBUG, 'rho' )
  call theta%log_field( LOG_LEVEL_DEBUG, 'theta' )
  call u%log_field(     LOG_LEVEL_DEBUG, 'u' )

  ! Write checksums to file
  call checksum_alg('dynamo', rho, 'rho', theta, 'theta', u, 'u')

  ! Write checkpoint/restart files if required
  if( restart%write_file() ) then 

    write(log_scratch_space,'(A,A)') "writing file:",  &
                            trim(restart%endfname("rho"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     call rho%write_checkpoint(trim(restart%endfname("rho")))

     write(log_scratch_space,'(A,A)') "writing file:",  &
                             trim(restart%endfname("u"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     call u%write_checkpoint(trim(restart%endfname("u")))

     write(log_scratch_space,'(A,A)') "writing file:",  &
                             trim(restart%endfname("theta"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     call theta%write_checkpoint(trim(restart%endfname("theta")))
  
     write(log_scratch_space,'(A,A)') "writing file:",  &
                             trim(restart%endfname("xi"))
     call log_event(log_scratch_space,LOG_LEVEL_INFO)
     call xi%write_checkpoint(trim(restart%endfname("xi")))

     if (use_moisture)then
       do i=1,nummr
         write(name, '(A,I2.2)') 'mr_', i
         write(log_scratch_space,'(A,A)') "writing file:",  &
            trim(restart%endfname(name))
         call log_event(log_scratch_space,LOG_LEVEL_INFO)
         call mr(i)%write_checkpoint(trim(restart%endfname(name)))
       end do
     end if
  end if

  ! Call timestep finalizers
  if ( transport_only .and. scheme == transport_scheme_method_of_lines) then
    call rk_transport_final( rho, theta)
  end if

  if ( subroutine_timers ) then
    call timer('gungho')
    call output_timer()
  end if
  
  call log_event( 'gungho completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Finalise XIOS context if we used it for IO
  if (write_xios_output) then
    call xios_context_finalize()
  end if

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise ESMF
  call ESMF_Finalize(endflag=ESMF_END_KEEPMPI,rc=rc)

  ! Finalise mpi
  call mpi_finalize(ierr)

end program gungho
