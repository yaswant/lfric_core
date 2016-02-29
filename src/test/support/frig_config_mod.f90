module frig_config_mod

  use constants_mod, only : i_def, i_native, l_def, r_def

  implicit none

  private
  public :: frig_base_mesh_config,           frig_initial_density_config, &
            frig_extrusion_config,           frig_finite_element_config,  &
            frig_initial_temperature_config, frig_initial_wind_config,    &
            frig_idealised_config,           frig_planet_config

  integer, parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_base_mesh_config( filename, geometry, partitioner, f_lat_deg )

    use base_mesh_config_mod, only : read_base_mesh_namelist, &
                                     key_from_geometry, key_from_partitioner

    implicit none

    character(*),      intent(in) :: filename
    integer(i_native), intent(in) :: geometry
    integer(i_native), intent(in) :: partitioner
    real(r_def),       intent(in) :: f_lat_deg

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("frig_base_mesh_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&base_mesh")')
    write( temporary_unit, '("filename = ''", A, "''")') filename
    write( temporary_unit, '("geometry = ''", A, "''")') key_from_geometry(geometry)
    write( temporary_unit, &
           '("partitioner = ''", A, "''")') key_from_partitioner(partitioner)
    write( temporary_unit, '("fplane = ", L)') .False.
    write( temporary_unit, '("f_lat_deg = ", E10.3)') f_lat_deg
    write( temporary_unit, '("/")')

    rewind(temporary_unit)
    call read_base_mesh_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_base_mesh_config: Unable to close temporary file'

  end subroutine frig_base_mesh_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_extrusion_config( method,     &
                                    domain_top, &
                                    number_of_layers )

    use extrusion_config_mod, only : read_extrusion_namelist, &
                                     key_from_method

    implicit none

    integer(i_native), intent(in) :: method
    real(r_def),       intent(in) :: domain_top
    integer(i_def),    intent(in) :: number_of_layers

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("frig_extrusion_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&extrusion")' )
    write( temporary_unit, &
           '("extrusion_method = ''", A, "''")' ) key_from_method(method)
    write( temporary_unit, '("domain_top = ", E10.3)' ) domain_top
    write( temporary_unit, '("number_of_layers = ", I0)' ) number_of_layers
    write( temporary_unit, '("/")')

    rewind(temporary_unit)
    call read_extrusion_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_base_mesh_config: Unable to close temporary file'

  end subroutine frig_extrusion_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_idealised_config( test )

    use idealised_config_mod, only : read_idealised_namelist, key_from_test

    implicit none

    integer, intent(in) :: test

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) stop 'frig_idealised_config: Unable to open temporary file'

    write( temporary_unit, '("&idealised")' )
    write( temporary_unit, '("test = ''", A, "''")' ) key_from_test(test)
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_idealised_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_idealised_config: Unable to close temporary file'
  end subroutine frig_idealised_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_finite_element_config( shape, element_order, rehabilitate )

    use finite_element_config_mod, only : key_from_shape, &
                                          read_finite_element_namelist

    implicit none

    integer(i_native), intent(in) :: shape
    integer(i_def),    intent(in) :: element_order
    logical(l_def),    intent(in) :: rehabilitate

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) stop 'frig_finite_element_config: Unable to open temporary file'

    write( temporary_unit, '("&finite_element")' )
    write( temporary_unit, '("shape = ", A)' ) key_from_shape( shape )
    write( temporary_unit, '("element_order = ", I0)' ) element_order
    write( temporary_unit, '("rehabilitate = ", L)' ) rehabilitate
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_finite_element_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_finite_element_config: Unable to close temporary file'

  end subroutine frig_finite_element_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_initial_density_config( tracer_max, tracer_background, &
                                          y1, y2 )

    use initial_density_config_mod, only : read_initial_density_namelist

    implicit none

    real(r_def), intent(in) :: tracer_max
    real(r_def), intent(in) :: tracer_background
    real(r_def), intent(in) :: y1
    real(r_def), intent(in) :: y2

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) stop 'frig_initial_density_config: Unable to open temporary file'

    write( temporary_unit, '("&initial_density")' )
    write( temporary_unit, '("tracer_max = ", E10.3)' ) tracer_max
    write( temporary_unit, '("tracer_background = ", E10.3)' ) tracer_background
    write( temporary_unit, '("y1 = ", E10.3)' ) y1
    write( temporary_unit, '("y2 = ", E10.3)' ) y2
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_initial_density_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_initial_density_config: Unable to close temporary file'

  end subroutine frig_initial_density_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_initial_temperature_config( bvf_square )

    use initial_temperature_config_mod, only : read_initial_temperature_namelist

    implicit none

    real(r_def), intent(in) :: bvf_square

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) stop 'frig_initial_temperature_config: Unable to open temporary file'

    write( temporary_unit, '("&initial_temperature")' )
    write( temporary_unit, '("bvf_square = ", E10.3)' ) bvf_square
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_initial_temperature_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_initial_temperature_config: Unable to close temporary file'

  end subroutine frig_initial_temperature_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_initial_wind_config( profile, u0, v0, rotation_angle )

    use initial_wind_config_mod, only : read_initial_wind_namelist, &
                                        key_from_profile

    implicit none

    integer,     intent(in) :: profile
    real(r_def), intent(in) :: u0
    real(r_def), intent(in) :: v0
    real(r_def), intent(in) :: rotation_angle

    integer :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) stop 'frig_initial_temperature_config: Unable to open temporary file'

    write( temporary_unit, '("&initial_wind")' )
    write( temporary_unit, '("profile = ''", A, "''")' ) key_from_profile(profile)
    write( temporary_unit, '("u0 = ", E10.3)' ) u0
    write( temporary_unit, '("v0 = ", E10.3)' ) v0
    write( temporary_unit, '("rotation_angle = ", E10.3)' ) rotation_angle
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_initial_wind_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_initial_temperature_config: Unable to close temporary file'

  end subroutine frig_initial_wind_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_planet_config( radius,  &
                                 gravity, &
                                 omega,   &
                                 rd,      &
                                 cp,      &
                                 p_zero,  &
                                 scaling )

    use planet_config_mod, only : read_planet_namelist

    implicit none

    real(r_def),           intent(in) :: radius
    real(r_def),           intent(in) :: gravity
    real(r_def),           intent(in) :: omega
    real(r_def),           intent(in) :: rd
    real(r_def),           intent(in) :: cp
    real(r_def),           intent(in) :: p_zero
    real(r_def), optional, intent(in) :: scaling

    integer     :: condition
    real(r_def) :: factor

    if (present(scaling)) then
      factor = scaling
    else
      factor = 1.0_r_def
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("frig_planet_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&planet")')
    write( temporary_unit, '("gravity = ", E14.7)') gravity
    write( temporary_unit, '("radius = ", E14.7)') radius
    write( temporary_unit, '("omega = ", E14.7)') omega
    write( temporary_unit, '("rd = ", E13.6)') rd
    write( temporary_unit, '("cp = ", E13.6)') cp
    write( temporary_unit, '("p_zero = ", E10.3)') p_zero
    write( temporary_unit, '("scaling_factor = ", E10.3)') factor
    write( temporary_unit, '("/")')

    rewind(temporary_unit)
    call read_planet_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_planet_config: Unable to close temporary file'

  end subroutine frig_planet_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine frig_subgrid_config( rho_approximation,        &
                                  transport_stencil_length, &
                                  rho_stencil_length )

    use subgrid_config_mod, only : read_subgrid_namelist,      &
                                   key_from_rho_approximation

    implicit none

    integer(i_native), intent(in) :: rho_approximation
    integer(i_def),    intent(in) :: transport_stencil_length
    integer(i_def),    intent(in) :: rho_stencil_length

    integer     :: condition

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("frig_subgrid_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&subgrid")')
    write( temporary_unit, '("rho_approximation = ", A)') &
                                key_from_rho_approximation( rho_approximation )
    write( temporary_unit, &
                '("transport_stencil_length = ", I0)') transport_stencil_length
    write( temporary_unit, '("rho_stencil_length = ", I0)') rho_stencil_length
    write( temporary_unit, '("/")')

    rewind(temporary_unit)
    call read_subgrid_namelist( temporary_unit )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop 'frig_subgrid_config: Unable to close temporary file'

  end subroutine frig_subgrid_config

end module frig_config_mod