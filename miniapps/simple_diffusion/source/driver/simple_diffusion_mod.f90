!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> simple_diffusion knows what configuration it needs.
!>
module simple_diffusion_mod

  implicit none

  private

  character(*), public, parameter ::                                &
      simple_diffusion_required_namelists(5) =  [ 'base_mesh     ', &
                                                  'extrusion     ', &
                                                  'finite_element', &
                                                  'partitioning  ', &
                                                  'planet        ']

end module simple_diffusion_mod
