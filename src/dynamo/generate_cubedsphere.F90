!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!  Program to generate a cubed-sphere mesh and write this in ugrid format
!  to the specified file.
!  Passes command-line arguments to gencube_type and uses ncdf_quad_mod to
!  write resulting mesh.
!  Invocation without arguments, or omission of any one or more
!  arguments, leads to the default output of: -o ugrid_quads_2d.nc -ndivs 4
!-------------------------------------------------------------------------------
program generate_cubedsphere
!-------------------------------------------------------------------------------
use gencube_mod,              only : gencube_ps_type
use generate_cubedsphere_mod, only : parse_args
use ugrid_2d_mod,             only : ugrid_2d_type
use ugrid_file_mod,           only : ugrid_file_type
use ncdf_quad_mod,            only : ncdf_quad_type
use constants_mod,            only : i_def, r_def, str_def
use iso_fortran_env,          only : stdout => output_unit

implicit none
!-------------------------------------------------------------------------------
  type(gencube_ps_type)                  :: csgen
  type(ugrid_2d_type)                    :: ugrid_2d
  class(ugrid_file_type), allocatable    :: ugrid_file
  character(len=str_def)                 :: filename, sztext
  integer(kind=i_def)                    :: ndivs
  logical                                :: nowrite
  integer                                :: fsize
  

  call parse_args(filename, ndivs, nowrite)

  allocate(ncdf_quad_type::ugrid_file)
  call ugrid_2d%set_file_handler(ugrid_file)

  csgen = gencube_ps_type(ndivs)

  write(stdout, "(A)") "Generating cubed-sphere mesh with..."
  write(stdout, "(A,I5)") "  ndivs: ", ndivs

  call ugrid_2d%set_by_generator(csgen)
  write(stdout, "(A)") "...generation complete."

  if(.not.nowrite) then
    write(stdout, "(A)", advance="NO") "Writing ugrid mesh to "//trim(adjustl(filename))//" ..."
    call ugrid_2d%write_to_file(trim(filename))
    inquire(file=filename, size=fsize)
    write(sztext, *) fsize
    write(stdout, "(A)") "... "//trim(adjustl(sztext))//" bytes written."
  else
    write(stdout, "(A)") "-nowrite selected, no output written."
  end if

  stop

end program generate_cubedsphere
