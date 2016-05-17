!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
module generate_cubedsphere_mod

  implicit none

  private
  public :: parse_args

contains
!-------------------------------------------------------------------------------
!>  @brief      Displays program's usage information on the specified output.
!!
!!  @details    Displays available range of parameters and syntax for generation
!!              and writing of cubed-sphere mesh.
!!
!!  @param[in]  dest  Destination stream to which usage should be written.
!-------------------------------------------------------------------------------
subroutine write_usage(dest)
  implicit none

  integer, intent(in)                    :: dest

  write(dest, "(A)") "Usage: generate_cubedsphere -h | -r <input_file> | "//&
                     "[[-o <output_file>] [-ndivs <divs_per_panel>] "//&
                     "[-nowrite]]"
  write(dest, "(A)") "   -h                      Print this help."
  write(dest, "(A)") "   -o  <output_file>       Write ugrid data to "//&
                                                 "<output_file>."
  write(dest, "(A)") "   -ndivs <divs_per_panel> Generate mesh with "//&
                     "<divs_per_panel> subdivisions per panel of "//&
                     "the cubed-sphere."
  write(dest, "(A)") "   -nowrite                Generate mesh in "//&
                     "memory only, do not write ugrid file. Useful "//&
                     "for benchmarking."
  write(dest, "(A)") "   -r  <input_file>        Read existing mesh"//&
                                                 " <input_file>."
  write(dest, *)
  write(dest, "(A)") "Defaults: -o ugrid_quads_2d.nc -ndivs 4"

end subroutine write_usage
!-------------------------------------------------------------------------------
!>  @brief      Parses any command-line arguments and assigns resulting
!!              values to subroutine arguments.
!!  @details    Handles erroneous input by printing usage to stderr and
!!              exiting.  Assigns default values to any argument that is not
!!              provided by the user.
!!
!!  @param[out]  filename  Filename to which ugrid output is to be written.
!!  @param[out]  ndivs  Number of subdivisions per panel of cubed-sphere mesh.
!!  @param[out]  nowrite  Whether mesh should be written to a file or discarded.
!-------------------------------------------------------------------------------
subroutine parse_args(filename, ndivs, nowrite)
  use constants_mod,       only : i_def, r_def, str_def
  use iso_fortran_env,     only : stdout => output_unit, &
                                  stderr => error_unit
  implicit none

  character(len=*), intent(out)          :: filename
  integer(kind=i_def), intent(out)       :: ndivs
  logical                                :: nowrite
  
  integer                                :: argc, arg
  character(len=str_def)                 :: argv(8)


  filename = "ugrid_quads_2d.nc"
  ndivs = 4_i_def
  argc = command_argument_count()
  nowrite = .FALSE.

  do arg = 1, argc
      call get_command_argument(arg, argv(arg))
  end do

  arg = 1
  do while(arg <= argc)
      select case(argv(arg))
        case("-h")
            call write_usage(stdout)
            stop
        case("-o")
            arg = arg + 1
            filename = trim(adjustl(argv(arg)))
        case("-r")
            arg = arg + 1
            filename = trim(adjustl(argv(arg)))
            call read_file(filename)
            stop
        case("-ndivs")
            arg = arg + 1
            read(argv(arg), *) ndivs
        case("-nowrite")
            nowrite = .TRUE.
        case default
            write(stderr, "(A)") "Unrecognised option: "//argv(arg)
            call write_usage(stderr)
            stop
      end select
      arg = arg + 1
  end do

end subroutine parse_args
!-------------------------------------------------------------------------------
!>  @brief      Test routine to read and verify existing ugrid mesh file.
!!
!!  @details    Accepts name of file containing ugrid mesh, reads and displays
!!              mesh dimensions to stdout.
!!
!!  @param[in]  filename  Name of file containing ugrid mesh to read.
!-------------------------------------------------------------------------------
subroutine read_file(filename)
  use ugrid_2d_mod,        only : ugrid_2d_type
  use ugrid_file_mod,      only : ugrid_file_type
  use ncdf_quad_mod,       only : ncdf_quad_type
  use constants_mod,       only : i_def
  use iso_fortran_env,     only : stdout => output_unit
  implicit none

  character(len=*), intent(in)           :: filename

  type(ugrid_2d_type)                    :: infile
  class(ugrid_file_type), allocatable    :: ugrid_file

  integer(kind=i_def)                    :: nodes, edges, faces
  integer(kind=i_def)                    :: nodes_per_face, edges_per_face
  integer(kind=i_def)                    :: nodes_per_edge, max_faces_per_node


  allocate(ncdf_quad_type::ugrid_file)

  call infile%set_file_handler(ugrid_file)
  call infile%read_from_file(trim(adjustl(filename)))

  call infile%get_dimensions(nodes, edges, faces, nodes_per_face, &
                             edges_per_face, nodes_per_edge, max_faces_per_node)

  write(stdout, "(A)") "File "//trim(adjustl(filename))//" contains a ugrid "//&
                       "mesh with dimensions: "

  write(stdout, "(A,19X,I7)") " Nodes: ", nodes
  write(stdout, "(A,19X,I7)") " Edges: ", edges
  write(stdout, "(A,19X,I7)") " Faces: ", faces
  write(stdout, "(A,10X,I7)") " Nodes per face: ", nodes_per_face
  write(stdout, "(A,10X,I7)") " Edges per face: ", edges_per_face
  write(stdout, "(A,10X,I7)") " Nodes per edge: ", nodes_per_edge
  write(stdout, "(A,2X,I7)") " Maximum faces per node: ", max_faces_per_node

end subroutine read_file
!-------------------------------------------------------------------------------
end module generate_cubedsphere_mod
