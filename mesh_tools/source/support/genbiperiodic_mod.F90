!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief   Module to define the genbiperiodic_type, a subclass of the
!>          ugrid_generator_type which generates a biperiodic mesh in a
!>          format suitable for storage as a ugrid file.
!> @details Type implements the ugrid_generator_type interface to
!>          construct a biperiodic mesh.  All required connectivity is
!>          calculated and made availabel to the ugrid writer.
!>
!-------------------------------------------------------------------------------
module genbiperiodic_mod
!-------------------------------------------------------------------------------

  use calc_global_cell_map_mod,       only: calc_global_cell_map
  use constants_mod,                  only: r_def, i_def, str_def, str_long
  use global_mesh_map_collection_mod, only: global_mesh_map_collection_type
  use global_mesh_map_mod,            only: generate_global_mesh_map_id
  use log_mod,                        only: log_event, log_scratch_space, &
                                            LOG_LEVEL_ERROR
  use reference_element_mod,          only: W, S, E, N, SWB, SEB, NWB, NEB
  use ugrid_generator_mod,            only: ugrid_generator_type

  implicit none

  private

  ! Mesh Vertex directions: local aliases for reference_element_mod values
  integer(i_def), parameter :: NW = NWB
  integer(i_def), parameter :: NE = NEB
  integer(i_def), parameter :: SE = SEB
  integer(i_def), parameter :: SW = SWB

  integer(i_def), parameter :: NPANELS = 1

  ! Prefix for error messages
  character(len=*),   parameter :: PREFIX = "[Biperiodic Mesh] "
  character(str_def), parameter :: MESH_CLASS = "plane"

  type, extends(ugrid_generator_type), public :: genbiperiodic_type

    private

    character(str_def)          :: mesh_name
    character(str_def)          :: mesh_class
    character(str_def)          :: coord_units_x
    character(str_def)          :: coord_units_y
    character(str_long)         :: constructor_inputs
    real(r_def)                 :: dx, dy
    integer(i_def)              :: edge_cells_x, edge_cells_y
    integer(i_def)              :: npanels
    integer(i_def)              :: nmaps

    character(str_def), allocatable :: target_mesh_names(:)
    integer(i_def),     allocatable :: target_edge_cells_x(:)
    integer(i_def),     allocatable :: target_edge_cells_y(:)
    type(global_mesh_map_collection_type), allocatable :: global_mesh_maps

    integer(i_def), allocatable :: cell_next(:,:)     ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: verts_on_cell(:,:) ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: edges_on_cell(:,:) ! (4, edge_cells_x*edge_cells_y)
    integer(i_def), allocatable :: verts_on_edge(:,:) ! (2, edge_cells_x*edge_cells_y)
    real(r_def),    allocatable :: vert_coords(:,:)   ! (2, edge_cells_x*edge_cells_y)
    real(r_def),    allocatable :: cell_coords(:,:)   ! (2, edge_cells_x*edge_cells_y)

  contains
    procedure :: calc_adjacency
    procedure :: calc_face_to_vert
    procedure :: calc_edges
    procedure :: calc_coords
    procedure :: calc_cell_centres
    procedure :: generate
    procedure :: get_metadata
    procedure :: get_dimensions
    procedure :: get_coordinates
    procedure :: get_connectivity
    procedure :: get_global_mesh_maps


    procedure :: write_mesh
  end type genbiperiodic_type

!-------------------------------------------------------------------------------
  interface genbiperiodic_type
    module procedure genbiperiodic_constructor
  end interface genbiperiodic_type

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> @brief   Constructor for genbiperiodic_type
!> @details Accepts mesh dimension and optional domain size arguments
!>          for initialisation.
!>
!> @param[in] mesh_name       Name of this mesh topology
!> @param[in] edge_cells_x    Number of cells in biperiodic mesh x axis
!> @param[in] edge_cells_y    Number of cells in biperiodic mesh y axis
!> @param[in] domain_x        Domain size in x-axis
!> @param[in] domain_y        Domain size in y-axis
!> @param[in, optional] target_mesh_names    
!>                            Names of mesh(es) to map to
!> @param[in, optional] target_edge_cells_x
!>                            Number of cells in x axis of
!>                            target mesh(es) to map to
!> @param[in, optional] target_edge_cells_y
!>                            Number of cells in y axis of
!>                            target mesh(es) to map to
!> @return    self  Instance of genbiperiodic_type
!-------------------------------------------------------------------------------
function genbiperiodic_constructor( mesh_name,                  &
                                    edge_cells_x, edge_cells_y, &
                                    domain_x, domain_y,         &
                                    target_mesh_names,          &
                                    target_edge_cells_x,        &
                                    target_edge_cells_y )       &
                            result( self )

  implicit none

  character(len=*), intent(in) :: mesh_name
  integer(i_def),   intent(in) :: edge_cells_x, edge_cells_y
  real(r_def),      intent(in) :: domain_x, domain_y
  character(str_def), optional, intent(in) :: target_mesh_names(:)
  integer(i_def),     optional, intent(in) :: target_edge_cells_x(:)
  integer(i_def),     optional, intent(in) :: target_edge_cells_y(:)



  type( genbiperiodic_type ) :: self
  character(str_def)         :: rchar1
  character(str_def)         :: rchar2

  character(str_long) :: target_mesh_names_str
  character(str_long) :: target_edge_cells_x_str
  character(str_long) :: target_edge_cells_y_str

  integer(i_def) ::i

  integer(i_def) :: remainder = 0_i_def

  if ( edge_cells_x < 2 .or. &
       edge_cells_y < 2 ) then
    call log_event( PREFIX//"Invalid dimension argument.", LOG_LEVEL_ERROR )
  end if

  self%mesh_name    = trim(mesh_name)
  self%mesh_class   = trim(MESH_CLASS)
  self%edge_cells_x = edge_cells_x
  self%edge_cells_y = edge_cells_y
  self%npanels      = NPANELS
  self%nmaps        = 0_i_def

  if (domain_x < 0.0_r_def)                                               &
      call log_event( PREFIX//" x-domain argument must be non-negative.", &
                      LOG_LEVEL_ERROR )

  if (domain_y < 0.0_r_def)                                               &
      call log_event( PREFIX//" y-domain argument must be non-negative.", &
                      LOG_LEVEL_ERROR )

  self%dx = domain_x / self%edge_cells_x
  self%dy = domain_y / self%edge_cells_y


  write(rchar1,'(F10.2)') domain_x
  write(rchar2,'(F10.2)') domain_y
  write(self%constructor_inputs,'(2(A,I0),(A))')  &
      'edge_cells_x=', self%edge_cells_x,  ';' // &
      'edge_cells_y=', self%edge_cells_y,  ';' // &
      'domain_x='//trim(adjustl(rchar1))// ';' // &
      'domain_y='//trim(adjustl(rchar2))


  if ( present(target_edge_cells_x) .and. &
       present(target_edge_cells_y) .and. &
       present(target_mesh_names) ) then

    if ( size(target_edge_cells_x) == size(target_mesh_names) .and. &
         size(target_edge_cells_y) == size(target_mesh_names) ) then

      self%nmaps = size(target_mesh_names)
      allocate(self%target_mesh_names(self%nmaps))
      allocate(self%target_edge_cells_x(self%nmaps))
      allocate(self%target_edge_cells_y(self%nmaps))


      self%target_mesh_names(:)   = target_mesh_names(:)
      self%target_edge_cells_x(:) = target_edge_cells_x(:)
      self%target_edge_cells_y(:) = target_edge_cells_y(:)

      do i=1, self%nmaps

        ! Check that mesh is not being mapped to itself
        if (self%edge_cells_x == target_edge_cells_x(i) .and. &
            self%edge_cells_y == target_edge_cells_y(i)) then
          write(log_scratch_space, '(A)') &
               'Invalid target while attempting to map mesh to itself'
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

        ! for x-axis
        if (target_edge_cells_x(i) < self%edge_cells_x) then
          remainder = mod(self%edge_cells_x, target_edge_cells_x(i))
          write(log_scratch_space,'(2(A,I0),A)')           &
              'Target edge_cells[',target_edge_cells_x(i), &
              '] must be a factor of source edge_cells[',  &
               self%edge_cells_x, ']'
        else if (target_edge_cells_x(i) > self%edge_cells_x) then
          remainder = mod(target_edge_cells_x(i), self%edge_cells_x)
          write(log_scratch_space,'(2(A,I0),A)')           &
              'Source edge_cells[',target_edge_cells_x(i), &
              '] must be a factor of target edge_cells[',  &
               self%edge_cells_x, ']'
        end if

        if (remainder == 0_i_def) then
          self%target_edge_cells_x(i) = target_edge_cells_x(i)
        else
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

        ! for y-axis
        if (target_edge_cells_y(i) < self%edge_cells_y) then
          remainder = mod(self%edge_cells_y, target_edge_cells_y(i))
          write(log_scratch_space,'(2(A,I0),A)')             &
               'Target edge_cells[',target_edge_cells_y(i),  &
               '] must be a factor of source edge_cells[',   &
               self%edge_cells_y, ']'
        else if (target_edge_cells_y(i) > self%edge_cells_y) then
          remainder = mod(target_edge_cells_y(i), self%edge_cells_y)
          write(log_scratch_space,'(2(A,I0),A)')             &
               'Source edge_cells[',target_edge_cells_y(i),  &
               '] must be a factor of target edge_cells[',   &
               self%edge_cells_y, ']'
        end if

        if (remainder == 0_i_def) then
          self%target_edge_cells_y(i) = target_edge_cells_y(i)
        else
          call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
        end if

      end do

    end if

    write(target_mesh_names_str, '(A)')   "'"//trim(adjustl(target_mesh_names(1)))//"'"
    write(target_edge_cells_x_str,'(I0)') target_edge_cells_x(1)
    write(target_edge_cells_y_str,'(I0)') target_edge_cells_y(1)
    if (size(target_mesh_names) > 1) then
      do i=2, self%nmaps
        write(target_mesh_names_str,'(A)')                  & 
            trim(adjustl(target_mesh_names_str)) // ",'" // &
            trim(adjustl(target_mesh_names(i)))  // "'"
        write(target_edge_cells_x_str,'(A,I0)')             &
            trim(adjustl(target_edge_cells_x_str)) // ',',  &
            target_edge_cells_x(i)
        write(target_edge_cells_y_str,'(A,I0)')             &
            trim(adjustl(target_edge_cells_y_str)) // ',',  &
            target_edge_cells_y(i)
      end do
    end if

    write(target_mesh_names_str,'(A)')    &
        'target_mesh_names=['//trim(adjustl(target_mesh_names_str))//']'
    write(target_edge_cells_x_str,'(A)')  &
        'target_edge_cells_x=['//trim(adjustl(target_edge_cells_x_str))//']'
    write(target_edge_cells_y_str,'(A)')  &
        'target_edge_cells_y=['//trim(adjustl(target_edge_cells_y_str))//']'


    write(self%constructor_inputs,'(A)')        &
        trim(self%constructor_inputs) // ';' // &
        trim(target_mesh_names_str)   // ';' // &
        trim(target_edge_cells_x_str) // ';' // &
        trim(target_edge_cells_y_str)

  end if

  return
end function genbiperiodic_constructor
!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the set of cells to which it is
!>          adjacent.
!> @details Allocates and populates the instance's cell_next(:,:) array
!>          with the id of each cell to which the index cell is adjacent.
!>
!> @param[in]   self       The genbiperiodic_type instance reference.
!> @param[out]  cell_next  A rank 2 (4,ncells)-sized array containing the
!>                         adjacency map.
!-------------------------------------------------------------------------------
subroutine calc_adjacency(self, cell_next)

  implicit none

  class(genbiperiodic_type),   intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: cell_next(:,:)

  integer(i_def) :: edge_cells_x, edge_cells_y, ncells
  integer(i_def) :: cell, astat

  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells = self%edge_cells_x * self%edge_cells_y

  allocate(cell_next(4, ncells), stat=astat)

  if (astat /= 0_i_def)                                         &
      call log_event( PREFIX//"Failure to allocate cell_next.", &
                      LOG_LEVEL_ERROR )

  do cell=1, ncells
    ! Cell default values
    cell_next(W, cell) = cell-1
    cell_next(S, cell) = cell+edge_cells_x
    cell_next(E, cell) = cell+1
    cell_next(N, cell) = cell-edge_cells_x

    ! Top row
    if (cell <= edge_cells_x) then
      cell_next(N, cell) = (edge_cells_y-1)*edge_cells_x + cell
    end if
    ! Bottom row
    if (cell > ncells-edge_cells_x) then
      cell_next(S, cell) = cell - (edge_cells_y-1)*edge_cells_x
    end if
    ! Left edge
    if (mod(cell, edge_cells_x) == 1) then
      cell_next(W, cell) = cell + edge_cells_x-1
    end if
    ! Right edge
    if (mod(cell, edge_cells_x) == 0) then
      cell_next(E, cell) = cell - (edge_cells_x-1)
    end if
  end do

end subroutine calc_adjacency

!-------------------------------------------------------------------------------
!> @brief   For each cell, calculates the four vertices whih comprise it.
!> @details Allocates and populates the instance's mesh(:,:) array with
!>          the vertices which form each cell.
!>
!> @param[in]   self           The genbiperiodic_type instance reference.
!> @param[out]  verts_on_cell  A rank 2 (4,ncells)-sized integer array of 
!>                             vertices which constitute each cell.
!-------------------------------------------------------------------------------
subroutine calc_face_to_vert(self, verts_on_cell)

  implicit none

  class(genbiperiodic_type),   intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: verts_on_cell(:,:)

  integer(i_def) :: edge_cells_x, edge_cells_y, ncells
  integer(i_def) :: y, vert, cell, nxf, astat



  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells = self%edge_cells_x * self%edge_cells_y

  cell = 1
  nxf  = 1

  allocate(verts_on_cell(4, ncells), stat=astat)

  if (astat /= 0)                                          &
      call log_event( PREFIX//"Failure to allocate mesh.", &
                      LOG_LEVEL_ERROR )

  do vert = 1, 4
    verts_on_cell(vert, cell) = nxf
    nxf = nxf+1
  end do

  ! East neighbour
  verts_on_cell(NW , self%cell_next(E, cell)) = verts_on_cell(NE, cell)
  verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)

  ! South neighbour
  verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
  verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)

  ! First row
  do cell = 2, edge_cells_x-1
    verts_on_cell(SE, cell) = nxf
    verts_on_cell(NE, cell) = nxf+1
    nxf = nxf + 2

    ! East neighbour
    verts_on_cell(NW , self%cell_next(E, cell)) = verts_on_cell(NE, cell)
    verts_on_cell(SW , self%cell_next(E, cell)) = verts_on_cell(SE, cell)

    ! South neighbour
    verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
    verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)
  end do

  ! Inner rows
  do y = 1, edge_cells_y-2
    ! First cell in row
    cell = y*edge_cells_x+1
    verts_on_cell(SW, cell) = nxf
    verts_on_cell(SE, cell) = nxf+1
    nxf = nxf+2

    ! South neighbour
    verts_on_cell(NW , self%cell_next(S, cell)) = verts_on_cell(SW, cell)
    verts_on_cell(NE , self%cell_next(S, cell)) = verts_on_cell(SE, cell)

    ! Remainder of row, every other cell
    do cell = y*edge_cells_x+3, (y+1)*edge_cells_x-1, 2
      verts_on_cell(SW, cell) = nxf
      verts_on_cell(SE, cell) = nxf+1
      nxf = nxf+2

      ! South neighbour
      verts_on_cell(NW, self%cell_next(S, cell)) = verts_on_cell(SW, cell)
      verts_on_cell(NE, self%cell_next(S, cell)) = verts_on_cell(SE, cell)
    end do

    ! Special case at end of row for odd edge_cells_x
    if (mod(edge_cells_x, 2) == 1) then
      cell = (y+1)*edge_cells_x
      verts_on_cell(SW, cell) = nxf
      nxf = nxf+1

      ! South neighbour
      verts_on_cell(NW, self%cell_next(S, cell)) = verts_on_cell(SW, cell)

      ! West neighbour
      verts_on_cell(SE, self%cell_next(W, cell)) = verts_on_cell(SW, cell)
    end if
  end do

  ! Right edge
  do cell = edge_cells_x, ncells, edge_cells_x
    ! Copy from left edge
    verts_on_cell(NE, cell) = verts_on_cell(NW, cell-edge_cells_x+1)
    verts_on_cell(SE, cell) = verts_on_cell(SW, cell-edge_cells_x+1)
  end do

  ! Inner step 2 cells
  do y = 1, edge_cells_y-2
    do cell = y*edge_cells_x+2, (y+1)*edge_cells_x, 2
      verts_on_cell(NW, cell) = verts_on_cell(SW, self%cell_next(N, cell))
      verts_on_cell(NE, cell) = verts_on_cell(SE, self%cell_next(N, cell))
      verts_on_cell(SE, cell) = verts_on_cell(SW, self%cell_next(E, cell))
      verts_on_cell(SW, cell) = verts_on_cell(SE, self%cell_next(W, cell))
    end do
  end do

  ! Special case at end of row for odd edge_cells_x
  if (mod(edge_cells_x, 2) == 1) then
    do cell = 2*edge_cells_x, ncells, edge_cells_x
      verts_on_cell(NW, cell) = verts_on_cell(SW, self%cell_next(N, cell))
    end do
  end if

  ! Last row
  do cell = ncells-edge_cells_x+1, ncells
    ! Copy from first row
    verts_on_cell(SE, cell) = verts_on_cell(NE, cell-(edge_cells_y-1)*edge_cells_x)
    verts_on_cell(SW, cell) = verts_on_cell(NW, cell-(edge_cells_y-1)*edge_cells_x)

    ! Copy from N
    verts_on_cell(NW, cell) = verts_on_cell(SW, self%cell_next(N, cell))
    verts_on_cell(NE, cell) = verts_on_cell(SE, self%cell_next(N, cell))
  end do

  return
end subroutine calc_face_to_vert

!-------------------------------------------------------------------------------
!> @brief   Calculates the edges which are found on each cell and the
!>          pair of vertices which are found on each edge.
!> @details Allocates and populates both the edges_on_cell and
!>          verts_on_edge arrays for the instance.
!>
!> @param[in]   self           The genbiperiodic_type instance reference.
!> @param[out]  edges_on_cell  A rank-2 (4,ncells)-sized integer array of
!>                             the edges found on each cell.
!> @param[out]  verts_on_edge  A rank-2 (2,2*ncells)-sized integer array
!>                             of the vertices found on each edge.
!>                             Vertex IDs listed with nearest vertex to
!>                             NW corner.
!-------------------------------------------------------------------------------
subroutine calc_edges(self, edges_on_cell, verts_on_edge)

  implicit none

  class(genbiperiodic_type),   intent(in)  :: self
  integer(i_def), allocatable, intent(out) :: edges_on_cell(:,:)
  integer(i_def), allocatable, intent(out) :: verts_on_edge(:,:)

  integer(i_def) :: edge_cells_x, edge_cells_y, ncells
  integer(i_def) :: cell, nxf, astat

  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells = self%edge_cells_x * self%edge_cells_y

  cell = 1
  nxf  = 1

  allocate(edges_on_cell(4, edge_cells_x*edge_cells_y), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate edges_on_cell.", &
                      LOG_LEVEL_ERROR )

  allocate(verts_on_edge(2, 2*edge_cells_x*edge_cells_y), stat=astat)

  if (astat /= 0)                                                   &
      call log_event( PREFIX//"Failure to allocate verts_on_edge.", &
                      LOG_LEVEL_ERROR )


  ! The assumed working layout for the folowing code with respect to 
  ! cell/panel is
  !
  ! ID Origin     Top (North)
  !    (NW) O--------------------O (NE)
  !         |                    +
  !         |   edge/vertex id   +
  !         |  + anti-clockwise  +
  !   Left  |   \   Numbering    + Right 
  !  (West) |    \               + (East)
  !         |     +------->      +
  !         |                    +
  !         |                    +
  !    (SE) O+-+-+-+-+-+-+-+-+-+-O (SE)
  !             Bottom (South)
  !
  ! This is a 2D grid, and compass points are used to orientate
  ! around the panel/cell. (top, left, right, bottom) could
  ! have been used, though top and bottom could have been
  ! confused with 3D cells. However, for readability,
  ! the terms (left, bottom, top, right) are used interchangably
  ! with the cardinal compass directions in the above diagram.
  ! It should also be noted that the cardinal compass directions
  ! are only used as an alternative way to reference the elements
  ! on the the panel/cells, i.e. the Top (North) edge of a cell
  ! may not actually be aligned in the same direction as geographic
  ! north.
  !
  ! When considered as a panel, the Eastern and Southern edges
  ! are taken as the `ghost` edges/vertices in the case of any
  ! periodicity.
  !  
  ! So that vertices on edges are consistent, listed vertices
  ! connected to edges will go from N-S, or W-E.

  !--------------------------------------------
  ! Top panel row, cell@left panel edge (ID:1)
  !--------------------------------------------
  ! Cell western edge
  edges_on_cell(W, cell)  = nxf
  verts_on_edge(1, nxf)   = self%verts_on_cell(NW, cell)
  verts_on_edge(2, nxf)   = self%verts_on_cell(SW, cell)


  ! Cell southern edge
  edges_on_cell(S, cell)  = nxf+1
  verts_on_edge(1, nxf+1) = self%verts_on_cell(SW, cell)
  verts_on_edge(2, nxf+1) = self%verts_on_cell(SE, cell)


  ! Cell eastern edge
  edges_on_cell(E, cell)  = nxf+2
  verts_on_edge(1, nxf+2) = self%verts_on_cell(NE, cell)
  verts_on_edge(2, nxf+2) = self%verts_on_cell(SE, cell)


  ! Cell northern edge
  edges_on_cell(N, cell)  = nxf+3
  verts_on_edge(1, nxf+3) = self%verts_on_cell(NW, cell)
  verts_on_edge(2, nxf+3) = self%verts_on_cell(NE, cell)

  nxf = nxf+4


  !-----------------------------------------------------------
  ! Top panel row, remaining cells, i.e. IDs = 2:edge_cells_x
  !-----------------------------------------------------------
  do cell = 2, edge_cells_x

    ! Cell western edge
    edges_on_cell(W, cell) = edges_on_cell(E,self%cell_next(W,cell))

    ! Cell southern edge
    edges_on_cell(S, cell) = nxf
    verts_on_edge(1, nxf)  = self%verts_on_cell(SW, cell)
    verts_on_edge(2, nxf)  = self%verts_on_cell(SE, cell)


    if (cell ==  edge_cells_x) then
      ! This cell is on the right-hand edge of the panel.

      ! For a biperiodic panel, if the cell is on the right-hand edge
      ! of the panel then the cell's eastern neighbour is actually
      ! the left-most cell on this row (because of the periodicity).
      ! At this point, left-most cell on this row has already
      ! had its edge ids assigned.
 
      ! In other words, the eastern edge of this cell exists on
      ! another cell where it has already been assigned an id. 

      ! Cell eastern edge
      edges_on_cell(E, cell)  = edges_on_cell(W,self%cell_next(E,cell))

      ! Cell northern edge
      edges_on_cell(N, cell)  = nxf+1
      verts_on_edge(1, nxf+1) = self%verts_on_cell(NW, cell)
      verts_on_edge(2, nxf+1) = self%verts_on_cell(NE, cell)

      nxf = nxf+2

    else 

      ! Cell eastern edge
      edges_on_cell(E, cell)  = nxf+1
      verts_on_edge(1, nxf+1) = self%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf+1) = self%verts_on_cell(SE, cell)


      ! Cell northern edge
      edges_on_cell(N, cell)  = nxf+2
      verts_on_edge(1, nxf+2) = self%verts_on_cell(NW, cell)
      verts_on_edge(2, nxf+2) = self%verts_on_cell(NE, cell)

      nxf = nxf+3
    end if
  
  end do

  !-----------------------------------------
  ! Internal panel rows
  !-----------------------------------------
  do cell = edge_cells_x+1, ncells-edge_cells_x

    if (mod(cell,edge_cells_x) == 1) then
      ! This cell is on the left-hand panel edge
      !-------------------------------------------

      ! Cell western edge
      edges_on_cell(W, cell)  = nxf
      verts_on_edge(1, nxf)   = self%verts_on_cell(NW, cell)
      verts_on_edge(2, nxf)   = self%verts_on_cell(SW, cell)


      ! Cell southern edge
      edges_on_cell(S, cell)  = nxf+1
      verts_on_edge(1, nxf+1) = self%verts_on_cell(SW, cell)
      verts_on_edge(2, nxf+1) = self%verts_on_cell(SE, cell)


      ! Cell eastern edge
      edges_on_cell(E, cell)  = nxf+2
      verts_on_edge(1, nxf+2) = self%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf+2) = self%verts_on_cell(SE, cell)

      nxf = nxf+3

    else if (mod(cell,edge_cells_x) == 0) then
      ! This cell is on the right-hand panel edge
      !--------------------------------------------

      ! Cell western edge
      edges_on_cell(W, cell) = edges_on_cell(E,self%cell_next(W,cell))

      ! Cell southern edge
      edges_on_cell(S, cell) = nxf
      verts_on_edge(1, nxf)  = self%verts_on_cell(SW, cell)
      verts_on_edge(2, nxf)  = self%verts_on_cell(SE, cell)

      nxf=nxf+1

      ! Cell eastern edge
      edges_on_cell(E, cell) = edges_on_cell(W,self%cell_next(E,cell))

    else
      ! This cell is internal to the panel
      !-------------------------------------------

      ! Cell western edge
      edges_on_cell(W, cell) = edges_on_cell(E,self%cell_next(W,cell))

      ! Cell southern edge
      edges_on_cell(S, cell) = nxf
      verts_on_edge(1, nxf)  = self%verts_on_cell(SW, cell)
      verts_on_edge(2, nxf)  = self%verts_on_cell(SE, cell)


      ! Cell eastern edge
      edges_on_cell(E, cell)  = nxf+1
      verts_on_edge(1, nxf+1) = self%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf+1) = self%verts_on_cell(SE, cell)

      nxf=nxf+2

    end if

    ! Cell northern edge
    ! The northern edges on these cells exist on the cells in
    ! the row above/previous to this one. Those edges have
    ! already been assigned ids.
    edges_on_cell(N, cell) = edges_on_cell(S,self%cell_next(N,cell))

  end do ! Panel inner rows


  ! Panel bottom row
  do cell = ncells-edge_cells_x+1, ncells

    if (mod(cell,edge_cells_x) == 1) then
      ! This cell is on the left-hand panel edge
      !-------------------------------------------

      ! Cell western edge
      edges_on_cell(W, cell) = nxf
      verts_on_edge(1, nxf)  = self%verts_on_cell(NW, cell)
      verts_on_edge(2, nxf)  = self%verts_on_cell(SW, cell)

      ! Cell southern edge
      ! For a biperiodic panel, if the cell is on the bottom edge
      ! of the panel then the cell's southern neighbour is actually
      ! the on the top-edge of the panel (because of the periodicity).
      ! At this point, all cells on the top-edge of the panel have
      ! had their edge ids assigned.
 
      ! In other words, the southern edge of this cell exists on
      ! another cell where it has already been assigned an id. 
      edges_on_cell(S, cell) = edges_on_cell(N,self%cell_next(S,cell))

      ! Cell eastern edge
      edges_on_cell(E, cell)  = nxf+1
      verts_on_edge(1, nxf+1) = self%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf+1) = self%verts_on_cell(SE, cell)

      nxf = nxf+2

    else if (mod(cell,edge_cells_x) == 0) then
      ! This cell is on the right-hand panel edge
      !--------------------------------------------

      ! Cell western edge
      edges_on_cell(W, cell) = edges_on_cell(E,self%cell_next(W,cell))
      
      ! Cell southern edge
      ! For a biperiodic panel, if the cell is on the bottom edge
      ! of the panel then the cell's southern neighbour is actually
      ! the on the top-edge of the panel (because of the periodicity).
      ! At this point, all cells on the top-edge of the panel have
      ! had their edge ids assigned.
 
      ! In other words, the southern edge of this cell exists on
      ! another cell where it has already been assigned an id. 
       edges_on_cell(S, cell) = edges_on_cell(N,self%cell_next(S,cell))
      
      ! Cell eastern edge
      ! For a biperiodic panel, if the cell is on the right-hand edge
      ! of the panel then the cell's eastern neighbour is actually
      ! the left-most cell on this row (because of the periodicity).
      ! At this point, left-most cell on this row has already
      ! had its edge ids assigned.
 
      ! In other words, the eastern edge of this cell exists on
      ! another cell where it has already been assigned an id.
      edges_on_cell(E, cell) = edges_on_cell(W,self%cell_next(E,cell))

    else
      ! Cell western edge
      edges_on_cell(W, cell) = edges_on_cell(E,self%cell_next(W,cell))

      ! Cell southern edge
      ! For a biperiodic panel, if the cell is on the bottom edge
      ! of the panel then the cell's southern neighbour is actually
      ! the on the top-edge of the panel (because of the periodicity).
      ! At this point, all cells on the top-edge of the panel have
      ! had their edge ids assigned.
      edges_on_cell(S, cell) = edges_on_cell(N,self%cell_next(S,cell))

      ! Cell eastern edge
      verts_on_edge(1, nxf) = self%verts_on_cell(NE, cell)
      verts_on_edge(2, nxf) = self%verts_on_cell(SE, cell)
      edges_on_cell(E, cell) = nxf
      nxf = nxf+1
    end if

    ! Cell north edge
    ! The northern edges on these cells exist on the cells in
    ! the row above/previous to this one. Those edges have
    ! already been assigned ids.
    edges_on_cell(N, cell) = edges_on_cell(S,self%cell_next(N,cell))

  end do

  return
end subroutine calc_edges

!-------------------------------------------------------------------------------
!> @brief   Calculates the coordinates of vertices in the mesh.
!> @details Assigns an (x,y) coordinate in units of dx and dy to each mesh
!>          vertex according to its Cartesian position in the mesh.
!>
!> @param[in]   self         The genbiperiodic_type instance reference.
!> @param[out]  vert_coords  A rank 2 (2,ncells)-sized real array of x and y
!>                           coordinates for each vertex.
!-------------------------------------------------------------------------------
subroutine calc_coords(self, vert_coords, coord_units_x, coord_units_y)

  implicit none

  class(genbiperiodic_type), intent(in)  :: self
  real(r_def), allocatable,  intent(out) :: vert_coords(:,:)
  character(str_def), intent(out) :: coord_units_x
  character(str_def), intent(out) :: coord_units_y

  integer(i_def) :: ncells, edge_cells_x, edge_cells_y
  integer(i_def) :: cell, x, y, astat, vert

  edge_cells_x = self%edge_cells_x
  edge_cells_y = self%edge_cells_y
  ncells = edge_cells_x*edge_cells_y

  allocate(vert_coords(2, ncells), stat=astat)

  if (astat /= 0)                                                 &
      call log_event( PREFIX//"Failure to allocate vert_coords.", &
                      LOG_LEVEL_ERROR )

  ! Cells begin numbering in rows from NW corner of panel
  do cell = 1, ncells
   y = 1+(cell-1)/edge_cells_x
   x = cell-(y-1)*edge_cells_x
   vert = self%verts_on_cell(NW, cell)
   vert_coords(1, vert) = real(x-1 - edge_cells_x/2,   r_def) * self%dx
   vert_coords(2, vert) = real(edge_cells_y/2 - (y-1), r_def) * self%dy
  end do

  coord_units_x = 'm'
  coord_units_y = 'm'

  return
end subroutine calc_coords

!-------------------------------------------------------------------------------
!> @brief   Calculates the mesh cell centres.
!> @details The face centres for the mesh are calculated based on the current
!>          node coordinates of the node on the NW corner of the cell.
!>          The node_cordinates are assumed to be in [m] in the
!>          x,y plane. Resulting face centre coordinates are in [m].
!>
!> @param[in,out]  self  The genbiperiodic_type instance reference.
!-------------------------------------------------------------------------------
subroutine calc_cell_centres(self)

  implicit none

  class(genbiperiodic_type), intent(inout) :: self

  integer(i_def) :: ncells
  real(r_def)    :: ratio

  ! Counters
  integer(i_def) :: cell, base_vert


  ncells = self%npanels*self%edge_cells_x*self%edge_cells_y

  ! 1.0 Initialise the face centres
  if (.not. allocated(self%cell_coords)) allocate( self%cell_coords(2,ncells) )
  self%cell_coords(:,:) = 0.0_r_def

  ! 2.0 Open cells have `ghost` nodes/edges and are located
  ! on the Eastern/Southern edges of the panel. Identify the open cells
  ! assuming that the numbering is along rows beginning from the NW corner
  ! of the panel. All cells, including the open cells have a unique NW
  ! vertex. Use this NW and self%dx, self%dy to calculate the face centre,
  ! assuming the cells are parallel in the x and y directions.
  do cell=1, ncells
    base_vert = self%verts_on_cell(NW, cell)
    self%cell_coords(1,cell) = self%vert_coords(1, base_vert) + self%dx/2.0_r_def
    self%cell_coords(2,cell) = self%vert_coords(2, base_vert) - self%dy/2.0_r_def
  end do

end subroutine calc_cell_centres

!-------------------------------------------------------------------------------
!> @brief Populates the arguments with the dimensions defining
!>        the biperiodic mesh.
!>
!> @param[in]   self                The genbiperiodic_type instance reference.
!> @param[out]  num_nodes           The number of nodes on the mesh.
!> @param[out]  num_edges           The number of edges on the mesh.
!> @param[out]  num_faces           The number of faces on the mesh.
!> @param[out]  num_nodes_per_face  The number of nodes around each face.
!> @param[out]  num_edges_per_face  The number of edges around each face.
!> @param[out]  num_nodes_per_face  The number of nodes around each edge.
!-------------------------------------------------------------------------------
subroutine get_dimensions( self,                                   &
                           num_nodes, num_edges, num_faces,        &
                           num_nodes_per_face, num_edges_per_face, &
                           num_nodes_per_edge )
  implicit none

  class(genbiperiodic_type), intent(in) :: self

  integer(i_def), intent(out) :: num_nodes
  integer(i_def), intent(out) :: num_edges
  integer(i_def), intent(out) :: num_faces
  integer(i_def), intent(out) :: num_nodes_per_face
  integer(i_def), intent(out) :: num_edges_per_face
  integer(i_def), intent(out) :: num_nodes_per_edge

  num_nodes = self%edge_cells_x*self%edge_cells_y
  num_edges = 2*self%edge_cells_x*self%edge_cells_y
  num_faces = self%edge_cells_x*self%edge_cells_y

  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!> @brief   Populates the argument array with the coordinates of the
!>          mesh's vertices.
!> @details Exposes the instance's vert_coords array to the caller.
!>
!> @param[in]   self              The genbiperiodic_type instance reference.
!> @param[out]  node_coordinates  The argument to receive the vert_coords data.
!> @param[out]  cell_coordinates  The argument to receive cell centre coordinates.
!> @param[out]  coord_units_x     Units for x-coordinates
!> @param[out]  coord_units_y     Units for y-coordinates
!-------------------------------------------------------------------------------
subroutine get_coordinates(self, node_coordinates, &
                                 cell_coordinates, &
                                 coord_units_x,    &
                                 coord_units_y)

  implicit none

  class(genbiperiodic_type), intent(in)  :: self
  real(r_def),               intent(out) :: node_coordinates(:,:)
  real(r_def),               intent(out) :: cell_coordinates(:,:)
  character(str_def),        intent(out) :: coord_units_x
  character(str_def),        intent(out) :: coord_units_y


  node_coordinates = self%vert_coords
  cell_coordinates = self%cell_coords
  coord_units_x    = self%coord_units_x
  coord_units_y    = self%coord_units_y

  return
end subroutine get_coordinates

!-------------------------------------------------------------------------------
!> @brief   Populates the argument arrays with the corresponding mesh
!>          connectivity information.
!> @details Implements the connectivity-providing interface required
!>          by the ugrid writer.
!>
!> @param[in]   self
!> @param[out]  face_node_connectivity  Face-node connectivity.
!> @param[out]  edge_node_connectivity  Edge-node connectivity.
!> @param[out]  face_edge_connectivity  Face-edge connectivity.
!> @param[out]  face_face_connectivity  Face-face connectivity.
!-------------------------------------------------------------------------------
subroutine get_connectivity( self,                   &
                             face_node_connectivity, &
                             edge_node_connectivity, &
                             face_edge_connectivity, &
                             face_face_connectivity )

  implicit none

  class(genbiperiodic_type), intent(in) :: self
  integer(i_def), intent(out) :: face_node_connectivity(:,:)
  integer(i_def), intent(out) :: edge_node_connectivity(:,:)
  integer(i_def), intent(out) :: face_edge_connectivity(:,:)
  integer(i_def), intent(out) :: face_face_connectivity(:,:)

  face_node_connectivity = self%verts_on_cell
  edge_node_connectivity = self%verts_on_edge
  face_edge_connectivity = self%edges_on_cell
  face_face_connectivity = self%cell_next

  return
end subroutine get_connectivity

!-------------------------------------------------------------------------------
!> @brief   Generates the mesh and connectivity.
!> @details Calls each of the instance methods which calculate the
!>          specified mesh and populate the arrays.
!>
!> @param[in,out]  self  The genbiperiodic_type instance reference.
!-------------------------------------------------------------------------------
subroutine generate(self)

  implicit none

  class(genbiperiodic_type), intent(inout) :: self

  call calc_adjacency(self, self%cell_next)
  call calc_face_to_vert(self, self%verts_on_cell)
  call calc_edges(self, self%edges_on_cell, self%verts_on_edge)
  if (self%nmaps > 0) call calc_global_mesh_maps(self)
  call calc_coords(self,               &
                   self%vert_coords,   &
                   self%coord_units_x, &
                   self%coord_units_y)
  call calc_cell_centres(self)

  return
end subroutine generate


!-------------------------------------------------------------------------------
!> @brief   PRIVATE subroutine to generate the reqeusted global mesh maps.
!> @details A map is generated for each requested target mesh based on this
!>          mesh objects mesh details, and those of the requested target
!>          meshes using calc_global_cell_map.
!>
!> @param[in,out]  self  The genbiperiodic_type instance reference.
!-------------------------------------------------------------------------------
subroutine calc_global_mesh_maps(self)

  implicit none

  class(genbiperiodic_type), intent(inout) :: self

  integer(i_def) :: source_id, source_cpp, source_ncells, &
                    target_edge_cells_x, target_edge_cells_y, target_cpp, &
                    target_ncells, target_cells_per_source_cell,i
  integer(i_def), allocatable :: cell_map(:,:)

  allocate( self%global_mesh_maps, source=global_mesh_map_collection_type())

  source_id  = 1
  source_cpp = self%edge_cells_x*self%edge_cells_y
  source_ncells = source_cpp*self%npanels

  do i=1, size(self%target_mesh_names)

    target_edge_cells_x = self%target_edge_cells_x(i)
    target_edge_cells_y = self%target_edge_cells_y(i)
    target_cpp          = target_edge_cells_x*target_edge_cells_y
    target_ncells       = target_cpp*self%npanels
    target_cells_per_source_cell = max(1,target_ncells/source_ncells)
    allocate(cell_map(target_cells_per_source_cell,source_ncells))

    call calc_global_cell_map(self, target_edge_cells_x, target_edge_cells_y, cell_map )
    call self%global_mesh_maps%add_global_mesh_map( source_id, i+1, cell_map )
    deallocate(cell_map)

  end do
  
  return
end subroutine calc_global_mesh_maps


!-------------------------------------------------------------------------------
!> @brief Writes out the mesh and connectivity for debugging purposes.
!>
!> @param[in,out]  self  The genbiperiodic_type instance reference.
!-------------------------------------------------------------------------------
subroutine write_mesh(self)

  use iso_fortran_env,     only : stdout => output_unit

  implicit none

  class(genbiperiodic_type), intent(in) :: self

  integer(i_def) :: i, cell, ncells

  ncells = self%edge_cells_x * self%edge_cells_y

  write(stdout,*) "cell_next"
  do cell=1, ncells
      write(stdout,"(I3,T8,4I4)") cell, self%cell_next(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "verts_on_cell"
  do cell=1, ncells
    write(stdout,"(I3,T8,4I4)") cell, self%verts_on_cell(:,cell)
  end do

  write(stdout,*)
  write(stdout,*) "verts_on_edge"
  do i=1, size(self%verts_on_edge, 2)
    write(stdout,"(I3,T8,2I4)") i, self%verts_on_edge(:,i)
  end do

  write(stdout,*)
  write(stdout,*) "edges_on_cell"
  do cell=1, ncells
    write(stdout,"(I3,T8,4I4)") cell, self%edges_on_cell(:,cell)
  end do

  write(*,*)
  write(*,*) "vert_coords"
  do cell=1, ncells
    write(*,"(I3,T8,2F11.4)") cell, self%vert_coords(:,cell)
  end do

  return
end subroutine write_mesh

!-----------------------------------------------------------------------------
!> @brief Returns mesh metadata information.
!>
!> @param[in]             self               The generator strategy object.
!> @param[out, optional]  mesh_name          Name of mesh instance to generate
!> @param[out, optional]  mesh_class         Primitive shape, i.e. sphere, plane
!> @param[out, optional]  npanels            Number of panels use to describe mesh
!> @param[out, optional]  edge_cells_x       Number of panel edge cells (x-axis).
!> @param[out, optional]  edge_cells_y       Number of panel edge cells (y-axis).
!> @param[out, optional]  constructor_inputs Inputs used to create this mesh from
!>                                           the this ugrid_generator_type
!> @param[out, optional]  nmaps              Number of maps to create with this mesh
!>                                           as source mesh
!> @param[out, optional]  maps_mesh_names    Mesh names of the target meshes that
!>                                           this mesh has maps for.
!> @param[out, optional]  maps_edge_cells_x  Number of panel edge cells (x-axis) of
!>                                           target mesh(es) to create map(s) for.
!> @param[out, optional]  maps_edge_cells_y  Number of panel edge cells (y-axis) of
!>                                           target mesh(es) to create map(s) for.
!-----------------------------------------------------------------------------
subroutine get_metadata( self,               &
                         mesh_name,          &
                         mesh_class,         &
                         npanels,            &
                         edge_cells_x,       &
                         edge_cells_y,       &
                         constructor_inputs, &
                         nmaps,              &
                         maps_mesh_names,    &
                         maps_edge_cells_x,  &
                         maps_edge_cells_y )
  implicit none

  class(genbiperiodic_type),     intent(in)  :: self
  character(str_def),  optional, intent(out) :: mesh_name
  character(str_def),  optional, intent(out) :: mesh_class
  integer(i_def),      optional, intent(out) :: npanels
  integer(i_def),      optional, intent(out) :: edge_cells_x
  integer(i_def),      optional, intent(out) :: edge_cells_y
  integer(i_def),      optional, intent(out) :: nmaps
  character(str_long), optional, intent(out) :: constructor_inputs

  character(str_def), optional, allocatable, intent(out) :: maps_mesh_names(:)
  integer(i_def),     optional, allocatable, intent(out) :: maps_edge_cells_x(:)
  integer(i_def),     optional, allocatable, intent(out) :: maps_edge_cells_y(:)

  if (present(mesh_name))          mesh_name          = self%mesh_name
  if (present(mesh_class))         mesh_class         = self%mesh_class
  if (present(npanels))            npanels            = self%npanels
  if (present(edge_cells_x))       edge_cells_x       = self%edge_cells_x
  if (present(edge_cells_y))       edge_cells_y       = self%edge_cells_y
  if (present(constructor_inputs)) constructor_inputs = trim(self%constructor_inputs)
  if (present(nmaps))              nmaps              = self%nmaps

  if (self%nmaps > 0) then
    if (present(maps_mesh_names))   maps_mesh_names   = self%target_mesh_names
    if (present(maps_edge_cells_x)) maps_edge_cells_x = self%target_edge_cells_x
    if (present(maps_edge_cells_y)) maps_edge_cells_y = self%target_edge_cells_y
  end if

  return
end subroutine get_metadata


!-------------------------------------------------------------------------------
!> @brief  Gets the global mesh map collection which uses
!>         this mesh as the source mesh
!>
!> @return global_mesh_maps global_mesh_map_collection_type
!-------------------------------------------------------------------------------
function get_global_mesh_maps(self) result(global_mesh_maps)

  implicit none

  class(genbiperiodic_type), target, intent(in) :: self

  type(global_mesh_map_collection_type), pointer :: global_mesh_maps

  nullify(global_mesh_maps)
  global_mesh_maps => self%global_mesh_maps

  return
end function get_global_mesh_maps

!-------------------------------------------------------------------------------
end module genbiperiodic_mod
