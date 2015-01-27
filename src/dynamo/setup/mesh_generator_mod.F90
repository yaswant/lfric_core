!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Generates a few simple meshes and the associated conntectivity relations
! This would be expected to be replaced with a preprocessor step/read from file
!-------------------------------------------------------------------------------
module mesh_generator_mod

use reference_element_mod, only : nfaces,   nedges,   nverts, &
                                 nfaces_h, nedges_h, nverts_h
use constants_mod,         only : r_def

implicit none

private get_lid_from_gid, binary_search

! Local parameters to describe the faces of the cube
!> Describes the face on the west side of the cell
integer, parameter :: W=1     
!> Describes the face on the south side of the cell
integer, parameter :: S=2
!> Describes the face on the east side of the cell
integer, parameter :: E=3
!> Describes the face on the north side of the cell
integer, parameter :: N=4
!> Describes the face on the bottom of the cell
integer, parameter :: B=5
!> Describes the face on the top of the cell
integer, parameter :: T=6

!> Describes the vertex at the south west bottom corner of the cell
integer, parameter :: SWB=1
!> Describes the vertex at the south east bottom corner of the cell
integer, parameter :: SEB=2
!> Describes the vertex at the north east bottom corner of the cell
integer, parameter :: NEB=3
!> Describes the vertex at the north west bottom corner of the cell
integer, parameter :: NWB=4
!> Describes the vertex at the south west top corner of the cell
integer, parameter :: SWT=5
!> Describes the vertex at the south east top corner of the cell
integer, parameter :: SET=6
!> Describes the vertex at the north east top corner of the cell
integer, parameter :: NET=7
!> Describes the vertex at the north west top corner of the cell
integer, parameter :: NWT=8

!> Describes the west bottom edge of the cell 
integer, parameter :: WB=1
!> Describes the south bottom edge of the cell 
integer, parameter :: SB=2
!> Describes the east bottom edge of the cell 
integer, parameter :: EB=3
!> Describes the north bottom edge of the cell 
integer, parameter :: NB=4
!> Describes the south west edge of the cell 
integer, parameter :: SW=5
!> Describes the south east edge of the cell 
integer, parameter :: SE=6
!> Describes the north east edge of the cell 
integer, parameter :: NE=7
!> Describes the north west edge of the cell 
integer, parameter :: NW=8
!> Describes the west top edge of the cell 
integer, parameter :: WT=9
!> Describes the south top edge of the cell 
integer, parameter :: ST=10
!> Describes the east top edge of the cell 
integer, parameter :: ET=11
!> Describes the north top edge of the cell 
integer, parameter :: NT=12

!> Global number of edges in a single 2D layer
integer :: nedge_h_g
!> Global number of vertices in a single 2D layer
integer :: nvert_h_g

!> Global number of faces in the full 3D domain
integer :: nface_g
!> Global number of edges in the full 3D domain
integer :: nedge_g
!> Global number of vertices in the full 3D domain
integer :: nvert_g

! In terminology of Logg 08 these are:
! cell_next    => MeshConnectivity(3,3) ( cells incident to cells )
! vert_on_cell => MeshConnectivity(3,0) ( vertices incident to cells )
! mesh_vertex  => MeshGeometry

! This is the minimal set of information, from which all other connectivity can be computed

!> Lookup table of the local indices of all the cells next to a given cell 
integer,          allocatable :: cell_next(:,:)
!> Lookup table of the local indices of all the vertices on a given cell 
integer,          allocatable :: vert_on_cell(:,:)
!> The x-, y- and z-coordinates of all the vertices on the local partition
real(kind=r_def), allocatable :: mesh_vertex(:,:)

! Extra connectivity for easy dofmap computation
! In terminology of Logg 08 these are:
! face_on_cell => MeshConnectivity(3,2) ( faces incident to cells )
! edge_on_cell => MeshConnectivity(3,1) ( edges incident to cells )

!> Lookup table of the local indices of all the faces on a given cell 
integer, allocatable ::  face_on_cell(:,:)
!> Lookup table of the local indices of all the edges on a given cell 
integer, allocatable ::  edge_on_cell(:,:)

! together this gives all the cell -> d connectivity (3,d), d=0,1,2

! types definitions for computing/storing the domain size
  type coordinate
    real(kind=r_def) :: x,y,z
  end type coordinate

  type domain_limits
    type(coordinate) :: minimum, maximum
  end type

!> The domain limits (x,y,z) for Cartesian domains
  type(domain_limits) :: domain_size

!> The z-coordinate of the top of the domain 
  real(kind=r_def) :: domain_top

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> Subroutine to allocate mesh arrays for connectivity
!> @param[in] ncells the number of cells in a horizontal layer on this partition
!> @param[in] nlayers the number of vertical layers
subroutine mesh_generator_init(ncells,nlayers)
!-----------------------------------------------------------------------------
! Subroutine to allocate connectivity
!-----------------------------------------------------------------------------

! number of cells in a horizontal layer on this partition
  integer, intent(in) :: ncells 
! number of vertical layers
  integer, intent(in) :: nlayers   


  allocate ( cell_next(nfaces,ncells*nlayers) )
  allocate ( vert_on_cell(nverts,ncells*nlayers) )

end subroutine mesh_generator_init

!> Generate a biperiodic domain (local to the partition).
!> @param nx, ny Number of cells in X and Y.
!> @param nlayers Number of vertical layers
!> @param dx, dy, dz Cell width in X, Y and Z direction.
!>
subroutine mesh_generator_biperiodic( nx, ny, nlayers, dx, dy, dz )

  use constants_mod, only: earth_radius
  use log_mod,       only: log_event, log_scratch_space, &
                           LOG_LEVEL_DEBUG, LOG_LEVEL_ERROR
  use mesh_mod, only : partitioned_cells, num_cells

  integer,              intent( in ) :: nx, ny
  integer,              intent( in ) :: nlayers
  real( kind = r_def ), intent( in ) :: dx, dy, dz

  integer         :: i, j, k, jd  ! Loop indices
  integer         :: lid, gid     ! local and global ids
  integer         :: gidx, gidy   ! global id converted to x,y coords

! List of edges horizontally around base of cell
  integer, allocatable :: edge_on_cell_h(:,:)

  allocate ( edge_on_cell_h(nedges_h,num_cells) )

! reset earth radius to 1 to avoid problems with routines 
! multiplying position by earth_radius
  earth_radius = 1.0_r_def

  do lid=1,num_cells

    gid=partitioned_cells(lid)

    !convert from a 1d id to x,y coord of cell
    gidy=1+(gid-1)/nx
    gidx=gid-nx*(gidy-1)

! Calculate the global ids of the cells next to this one
! j-1 cell (South face)
    cell_next(S,lid) = gid - nx
! i+1 cell (East face)
    cell_next(E,lid) = gid + 1
! j+1 cell (North face)
    cell_next(N,lid) = gid + nx      
! i-1 cell (West face)
    cell_next(W,lid) = gid - 1
! k-1 cell (bottom face)
    cell_next(B,lid) = gid - nx*ny
! k+1 cell (top face)
    cell_next(T,lid) = gid + nx*ny

! Now do periodicity/connectivity along edges 
! South
    if (gidy == 1)then
      cell_next(S,lid) = cell_next(S,lid)+nx*ny
    end if
! North  
    if (gidy == ny)then
      cell_next(N,lid) = cell_next(N,lid)-nx*ny
    endif
! West
    if (gidx == 1)then
      cell_next(W,lid) = cell_next(W,lid)+nx
    end if
! East  
    if (gidx == nx)then
      cell_next(E,lid) = cell_next(E,lid)-nx
    endif

  end do
  
  ! perform vertical extrusion for connectivity 
  do k=1,nlayers-1
    do i=1,num_cells
      lid = i + k*(num_cells)
      jd = lid - (num_cells)        
      do j=1,nfaces
        cell_next(j,lid) = cell_next(j,jd) + nx*ny
      end do
    end do
  end do
  
  ! cell_next(lid) still contains global ids so convert to local ids
  do k=0,nlayers-1
    do i=1,num_cells
      lid = i + k*(num_cells)
      do j=1,nfaces
        cell_next(j,lid)=get_lid_from_gid(cell_next(j,lid))
      end do
    end do
  end do

  ! Set connectivity at lower/upper boundary to some dummy cell
  do i=1,num_cells
    cell_next(B,i) = 0
    cell_next(T,i+(nlayers-1)*(num_cells)) = 0
  end do

! compute vertices on cells

! vertices around bottom face of cube
  nvert_h_g=0
  vert_on_cell(:,:) = 0
  do lid=1,num_cells
! 1. south west corner of cell
    if(vert_on_cell(SWB,lid)==0)then 
      nvert_h_g=nvert_h_g+1
      vert_on_cell(SWB,lid) = nvert_h_g
      if(cell_next(W,lid)>0)then                     ! and south east corner of cell to west 
        vert_on_cell(SEB,cell_next(W,lid)) = nvert_h_g
        if(cell_next(S,cell_next(W,lid))>0)then      ! and north east corner of cell to south west
          vert_on_cell(NEB,cell_next(S,cell_next(W,lid))) = nvert_h_g
        end if
      end if
      if(cell_next(S,lid)>0)then                     ! and north west corner of cell to south
        vert_on_cell(NWB, cell_next(S,lid)) = nvert_h_g
        if(cell_next(W,cell_next(S,lid))>0)then      ! and again north east corner of cell to south west (in case other route to southwest goes through a missing cell)
          vert_on_cell(NEB,cell_next(W,cell_next(S,lid))) = nvert_h_g
        end if
      end if
    end if
! 2. south east corner of cell
    if(vert_on_cell(SEB,lid)==0)then 
      nvert_h_g=nvert_h_g+1
      vert_on_cell(SEB,lid) = nvert_h_g
      if(cell_next(E,lid)>0)then                     ! and south west corner of cell to east 
        vert_on_cell(SWB,cell_next(E,lid)) = nvert_h_g
        if(cell_next(S,cell_next(E,lid))>0)then      ! and north west corner of cell to south east
          vert_on_cell(NWB,cell_next(S,cell_next(E,lid))) = nvert_h_g
        end if
      end if
      if(cell_next(S,lid)>0)then                     ! and north east corner of cell to south
        vert_on_cell(NEB,cell_next(S,lid)) = nvert_h_g
        if(cell_next(E,cell_next(S,lid))>0)then      ! and again north west corner of cell to south east (in case other route to southeast goes through a missing cell)
          vert_on_cell(NWB,cell_next(E,cell_next(S,lid))) = nvert_h_g
        end if
      end if
    end if
! 3. north east corner of cell
    if(vert_on_cell(NEB,lid)==0)then 
      nvert_h_g=nvert_h_g+1
      vert_on_cell(NEB,lid) = nvert_h_g
      if(cell_next(E,lid)>0)then                     ! and north west corner of cell to east 
        vert_on_cell(NWB,cell_next(E,lid)) = nvert_h_g
        if(cell_next(N,cell_next(E,lid))>0)then      ! and south west corner of cell to north east
          vert_on_cell(SWB,cell_next(N,cell_next(E,lid))) = nvert_h_g
        end if
      end if
      if(cell_next(N,lid)>0)then                     ! and south east corner of cell to north
        vert_on_cell(SEB,cell_next(N,lid)) = nvert_h_g
        if(cell_next(E,cell_next(N,lid))>0)then      ! and again south west corner of cell to north east (in case other route to northeast goes through a missing cell)
          vert_on_cell(SWB,cell_next(E,cell_next(N,lid))) = nvert_h_g
        end if
      end if
    end if
! 4. north west corner of cell
    if(vert_on_cell(NWB,lid)==0)then 
      nvert_h_g=nvert_h_g+1
      vert_on_cell(NWB,lid) = nvert_h_g
      if(cell_next(W,lid)>0)then                     ! and north east corner of cell to west 
        vert_on_cell(NEB,cell_next(W,lid)) = nvert_h_g
        if(cell_next(N,cell_next(W,lid))>0)then      ! and south east corner of cell to north west
          vert_on_cell(SEB,cell_next(N,cell_next(W,lid))) = nvert_h_g
        end if
      end if
      if(cell_next(N,lid)>0)then                     ! and south west corner of cell to north
        vert_on_cell(SWB,cell_next(N,lid)) = nvert_h_g
        if(cell_next(W,cell_next(N,lid))>0)then      ! and again south east corner of cell to north west (in case other route to northwest goes through a missing cell)
          vert_on_cell(SEB,cell_next(W,cell_next(N,lid))) = nvert_h_g
        end if
      end if
    end if
  end do
! vertices around top face of cube
  do lid=1,num_cells
    vert_on_cell(SWT,lid)=vert_on_cell(SWB,lid)+nvert_h_g
    vert_on_cell(SET,lid)=vert_on_cell(SEB,lid)+nvert_h_g
    vert_on_cell(NET,lid)=vert_on_cell(NEB,lid)+nvert_h_g
    vert_on_cell(NWT,lid)=vert_on_cell(NWB,lid)+nvert_h_g
  end do

! perform vertical extrusion for connectivity 
  do k=1,nlayers-1
    do i=1,num_cells
      lid = i + k*(num_cells)
      jd = lid - (num_cells)
      vert_on_cell(SWB,lid)=vert_on_cell(SWT,jd)
      vert_on_cell(SEB,lid)=vert_on_cell(SET,jd)
      vert_on_cell(NEB,lid)=vert_on_cell(NET,jd)
      vert_on_cell(NWB,lid)=vert_on_cell(NWT,jd)
      vert_on_cell(SWT,lid)=vert_on_cell(SWB,lid)+nvert_h_g
      vert_on_cell(SET,lid)=vert_on_cell(SEB,lid)+nvert_h_g
      vert_on_cell(NET,lid)=vert_on_cell(NEB,lid)+nvert_h_g
      vert_on_cell(NWT,lid)=vert_on_cell(NWB,lid)+nvert_h_g
    end do
  end do

! compute edges horizontally on cell
  nedge_h_g=0
  edge_on_cell_h(:,:)=0
  do lid=1,num_cells
! 1. south edge of cell
    if(edge_on_cell_h(SB,lid)==0)then
      nedge_h_g=nedge_h_g+1
      edge_on_cell_h(SB,lid) = nedge_h_g
      if(cell_next(S,lid)>0)then                     ! and north edge of cell to south
        edge_on_cell_h(NB,cell_next(S,lid)) = nedge_h_g
      end if
    end if
! 2. east edge of cell
    if(edge_on_cell_h(EB,lid)==0)then
      nedge_h_g=nedge_h_g+1
      edge_on_cell_h(EB,lid) = nedge_h_g
      if(cell_next(E,lid)>0)then                     ! and west edge of cell to east
        edge_on_cell_h(WB,cell_next(E,lid)) = nedge_h_g
      end if
    end if
! 3. north edge of cell
    if(edge_on_cell_h(NB,lid)==0)then
      nedge_h_g=nedge_h_g+1
      edge_on_cell_h(NB,lid) = nedge_h_g
      if(cell_next(N,lid)>0)then                     ! and south edge of cell to north
        edge_on_cell_h(SB,cell_next(N,lid)) = nedge_h_g
      end if
    end if
! 4. west edge of cell
    if(edge_on_cell_h(WB,lid)==0)then
      nedge_h_g=nedge_h_g+1
      edge_on_cell_h(WB,lid) = nedge_h_g
      if(cell_next(W,lid)>0)then                     ! and east edge of cell to west
        edge_on_cell_h(EB,cell_next(W,lid)) = nedge_h_g
      end if
    end if
  end do

  deallocate ( edge_on_cell_h )

  nface_g = nedge_h_g*nlayers + (num_cells)*(nlayers + 1)
  nedge_g = nedge_h_g*(nlayers + 1) + nlayers*nvert_h_g
  nvert_g = nvert_h_g*(nlayers + 1)
  
! allocate coordinate array
  allocate ( mesh_vertex(3,nvert_g) )

! Compute vertices
  do lid=1,num_cells

    gid=partitioned_cells(lid)

    !convert from a 1d id to x,y coord of cell
    gidy=1+(gid-1)/nx
    gidx=gid-nx*(gidy-1)

    do k=1,nlayers+1
      mesh_vertex(1,vert_on_cell(SWB,lid)+(k-1)*(num_cells)) = &
                                                     real(gidx-1 - nx/2)*dx 
      mesh_vertex(2,vert_on_cell(SWB,lid)+(k-1)*(num_cells)) = &
                                                     real(gidy-1 - ny/2)*dy 
      mesh_vertex(3,vert_on_cell(SWB,lid)+(k-1)*(num_cells)) = &
                                                     real(k-1)*dz
    end do

  end do
  domain_top = dz * real(nlayers)

  call set_domain_size()

! Diagnostic information  
  call log_event( 'grid connectivity', LOG_LEVEL_DEBUG )
  do i = 1, (num_cells) * nlayers
    write( log_scratch_space,'(7i6)' ) i, &
                            cell_next(S,i), cell_next(E,i), cell_next(N,i), &
                            cell_next(W,i), cell_next(B,i), cell_next(T,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do
  call log_event( 'verts on cells', LOG_LEVEL_DEBUG )
  do i = 1, (num_cells) * nlayers
    write( log_scratch_space, '(9i6)' ) i, &
                             vert_on_cell(SWB,i), vert_on_cell(SEB,i), &
                             vert_on_cell(NEB,i), vert_on_cell(NWB,i), &
                             vert_on_cell(SWT,i), vert_on_cell(SET,i), &
                             vert_on_cell(NET,i), vert_on_cell(NWT,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call log_event( 'vert coords', LOG_LEVEL_DEBUG )
  do i = 1, nvert_g
    write( log_scratch_space, '(i6,4f12.4)' ) &
         i, mesh_vertex(1,i), mesh_vertex(2,i), mesh_vertex(3,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

end subroutine mesh_generator_biperiodic

!> Generate a biperiodic domain of size nx * ny where nx * ny = ncells.
!>
!> <pre>
!> Location of panels:
!>             .....
!>            :  3  |
!>            :     |
!>             -----
!>      -----  .....  .....  -----
!>     |  5  :|  1  :|  2  :|  4  :
!>     |     :|     :|     :|     :
!>      .....  -----  -----  .....
!>             .....
!>            |  6  :
!>            |     :
!>             -----
!>
!>     Solid lines: left and bottom edges of panel
!>     Dotted lines: right and top edges of panel
!>     Currently reads in data created from John Thuburns gengrid_cube program
!> </pre>
!>
!> @param filename Data file for cubed sphere data.
!> @param ncells Total number of cells in horizontal layer.
!> @param nlayers Number of vertical layers.
!> @param dz Vertical grid spacing.
!>

subroutine mesh_generator_cubedsphere( filename, ncells, nlayers, dz )

  use log_mod,              only: log_event, log_scratch_space, &
                                  LOG_LEVEL_DEBUG, LOG_LEVEL_ERROR
  use ugrid_2d_mod,         only: ugrid_2d_type
  use ugrid_file_mod,       only: ugrid_file_type 
  use ncdf_quad_mod,        only: ncdf_quad_type 
  use coord_algorithms_mod, only: llr2xyz
  
  implicit none

  character(*),     intent(in) :: filename
  integer,          intent(in) :: ncells
  integer,          intent(in) :: nlayers
  real(kind=r_def), intent(in) :: dz

! Loop indices
  integer :: i, j, k, vert, id, jd
! data from file
  integer :: nvert_in, nface_in, nedge_in
  integer :: num_nodes_per_face, num_nodes_per_edge, num_edges_per_face
! lat/long coordinates
  !2D ugrid and generator strategy
  type(ugrid_2d_type)                 :: ugrid_2d
  class(ugrid_file_type), allocatable :: file_handler
  real(kind=r_def) :: long, lat, r

! topologically a cube
  nedge_h_g = 2*ncells
  nvert_h_g = (ncells + 2)
  
  nface_g = nedge_h_g*nlayers + ncells*(nlayers + 1)
  nedge_g = nedge_h_g*(nlayers + 1) + nvert_h_g*nlayers
  nvert_g = nvert_h_g*(nlayers + 1)

! allocate coordinate array
  allocate ( mesh_vertex(3,nvert_g) )

  allocate(ncdf_quad_type :: file_handler)
  call ugrid_2d%set_file_handler(file_handler)
  call ugrid_2d%read_from_file(trim(filename))

  call ugrid_2d%get_dimensions(                &
     num_nodes          = nvert_in,            &
     num_edges          = nedge_in,            &
     num_faces          = nface_in,            &
     num_nodes_per_face = num_nodes_per_face,  &
     num_edges_per_face = num_edges_per_face,  &
     num_nodes_per_edge = num_nodes_per_edge)

  if ( nface_in /= ncells ) then
    write( log_scratch_space, '(A, I0, A, I0)' ) &
         'Number of cells in ugrid file does not match: ', &
         nface_in, ' vs. ', ncells
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    stop
  end if
  if ( nvert_in /= nvert_h_g ) then
    write( log_scratch_space, '(A, I0, A, I0)' ) &
         'Number of vertices in ugrid file does not match: ', &
         nvert_in, ' vs. ', nvert_h_g
    stop
  end if  

  !Get coordinates of vertices
  call ugrid_2d%get_node_coords( mesh_vertex(1:3,1:nvert_h_g) )
  call ugrid_2d%get_face_node_connectivity(vert_on_cell(1:num_nodes_per_face,1:nface_in))
  call ugrid_2d%get_face_face_connectivity(cell_next(1:num_edges_per_face,1:nface_in))
  
   !mesh_vertex(3,1:nvert_h_g) = earth_radius

! add connectivity for up/down
  do j=1,ncells
    cell_next(5,j) = j - ncells
    cell_next(6,j) = j + ncells
  end do
       
! perform vertical extrusion for connectivity 
  do k=1,nlayers-1
    do i=1,ncells
      id = i + k*ncells      
      jd = id - ncells        
      do j=1,nfaces
        cell_next(j,id) = cell_next(j,jd) + ncells
      end do
    end do
  end do
    
! Set connectivity at lower/upper boundary to some dummy cell
  do i=1,ncells
    cell_next(5,i) = 0
    j =  i+(nlayers-1)*ncells
    cell_next(6,j) = 0
  end do
    
! perform vertical extrusion for vertices
  do j=1,nvert_in
    do k=1,nlayers                       
      mesh_vertex(1,j+k*nvert_in) = mesh_vertex(1,j)
      mesh_vertex(2,j+k*nvert_in) = mesh_vertex(2,j) 
      mesh_vertex(3,j+k*nvert_in) = mesh_vertex(3,j) + dz*real(k)
    end do
  end do
! Depth of atmosphere
  domain_top = dz * real(nlayers)

! Convert (long,lat,r) -> (x,y,z)
  do j=1,nvert_in
    do k=0,nlayers
      long = mesh_vertex(1,j+k*nvert_in)
      lat  = mesh_vertex(2,j+k*nvert_in)
      r    = mesh_vertex(3,j+k*nvert_in)
      call llr2xyz(long,lat,r,mesh_vertex(1,j+k*nvert_in),                  &
                              mesh_vertex(2,j+k*nvert_in),                  &
                              mesh_vertex(3,j+k*nvert_in))
    end do
  end do
  
! assign vertices to cells
  do j=1,ncells
    vert_on_cell(5,j) = vert_on_cell(1,j) + nvert_in
    vert_on_cell(6,j) = vert_on_cell(2,j) + nvert_in
    vert_on_cell(7,j) = vert_on_cell(3,j) + nvert_in
    vert_on_cell(8,j) = vert_on_cell(4,j) + nvert_in
! do vertical extrusion
    do k=1,nlayers-1
      do vert=1,8
        vert_on_cell(vert,j+k*ncells) = vert_on_cell(vert,j) + k*nvert_in
      end do
    end do
  end do 

  call set_domain_size()
 
! Diagnostic information  
  call log_event( 'grid connectivity', LOG_LEVEL_DEBUG )
  do i = 1, ncells * nlayers
    write( log_scratch_space, '(7i6)' ) i, &
                              cell_next(1,i), cell_next(2,i), cell_next(3,i), &
                              cell_next(4,i), cell_next(5,i), cell_next(6,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do
  call log_event( 'verts on cells', LOG_LEVEL_DEBUG )
  do i = 1, ncells * nlayers
    write( log_scratch_space, '(9i6)' ) i, &
                              vert_on_cell(1,i), vert_on_cell(2,i), &
                              vert_on_cell(3,i), vert_on_cell(4,i), &
                              vert_on_cell(5,i), vert_on_cell(6,i), &
                              vert_on_cell(7,i), vert_on_cell(8,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call log_event( 'vert coords', LOG_LEVEL_DEBUG )
  do i = 1, nvert_g
    write( log_scratch_space, '(i6,4e16.8)' ) &
         i, mesh_vertex(1,i), mesh_vertex(2,i), mesh_vertex(3,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do 

end subroutine mesh_generator_cubedsphere

!> Compute mesh connectivity.
!>
!> <pre>
!> cells->faces (3,2)
!> cells->edges (3,1)
!> </pre>
!>
!> @param ncells Total number of cells in horizontal layer.
!>
subroutine mesh_connectivity( ncells )

  use log_mod, only : log_event, log_scratch_space, &
                      LOG_LEVEL_DEBUG, LOG_LEVEL_ERROR

  integer, intent( in ) :: ncells
 
  integer :: i            ! loop counter used to write out connectivities
  integer :: lid          ! loop couonter over local ids
  integer :: nedge_layer  ! Number of edges for a single layer  
  integer :: nface_layer  ! Number of faces for a single layer  

  allocate( face_on_cell(nfaces,ncells) )
  allocate( edge_on_cell(nedges,ncells) )
  face_on_cell(:,:) = 0
  edge_on_cell(:,:) = 0

! compute faces on cells
  nface_layer=0
  do lid=1,ncells
! 1. south facing face of cell
    if(face_on_cell(S,lid)==0)then
      nface_layer=nface_layer+1
      face_on_cell(S,lid) = nface_layer
      if(cell_next(S,lid)>0)then                     ! and north facing face of cell to south
        face_on_cell(N,cell_next(S,lid)) = nface_layer
      end if
    end if
! 2. east facing face of cell
    if(face_on_cell(E,lid)==0)then
      nface_layer=nface_layer+1
      face_on_cell(E,lid) = nface_layer
      if(cell_next(E,lid)>0)then                     ! and west facing face of cell to east
        face_on_cell(W,cell_next(E,lid)) = nface_layer
      end if
    end if
! 3. north facing face of cell
    if(face_on_cell(N,lid)==0)then
      nface_layer=nface_layer+1
      face_on_cell(N,lid) = nface_layer
      if(cell_next(N,lid)>0)then                     ! and south facing face of cell to north
        face_on_cell(S,cell_next(N,lid)) = nface_layer
      end if
    end if
! 4. west facing face of cell
    if(face_on_cell(W,lid)==0)then
      nface_layer=nface_layer+1
      face_on_cell(W,lid) = nface_layer
      if(cell_next(W,lid)>0)then                     ! and east facing face of cell to west
        face_on_cell(E,cell_next(W,lid)) = nface_layer
      end if
    end if
! 5. bottom facing face of cell
     nface_layer=nface_layer+1
     face_on_cell(B,lid) = nface_layer
! 6. top facing face of cell
     nface_layer=nface_layer+1
     face_on_cell(T,lid) = nface_layer
  end do
  
! compute edges on cells
  nedge_layer=0
  do lid=1,ncells
! 1. top and bottom edge on south facing face of cell
    if(edge_on_cell(SB,lid)==0)then
      edge_on_cell(SB,lid) = nedge_layer+1
      edge_on_cell(ST,lid) = nedge_layer+2
      if(cell_next(S,lid)>0)then                     ! and north facing face of cell to south
        edge_on_cell(NB,cell_next(S,lid)) = nedge_layer+1
        edge_on_cell(NT,cell_next(S,lid))= nedge_layer+2
      end if
      nedge_layer=nedge_layer+2
    end if
! 2. top and bottom egde on east facing face of cell
    if(edge_on_cell(EB,lid)==0)then
      edge_on_cell(EB,lid) = nedge_layer+1
      edge_on_cell(ET,lid)= nedge_layer+2
      if(cell_next(E,lid)>0)then                     ! and west facing face of cell to east
        edge_on_cell(WB,cell_next(E,lid)) = nedge_layer+1
        edge_on_cell(WT,cell_next(E,lid))= nedge_layer+2
      end if
      nedge_layer=nedge_layer+2
    end if
! 3. top and bottom edge on north facing face of cell
    if(edge_on_cell(NB,lid)==0)then
      edge_on_cell(NB,lid) = nedge_layer+1
      edge_on_cell(NT,lid)= nedge_layer+2
      if(cell_next(N,lid)>0)then                     ! and south facing face of cell to north
        edge_on_cell(SB,cell_next(N,lid)) = nedge_layer+1
        edge_on_cell(ST,cell_next(N,lid)) = nedge_layer+2
      end if
      nedge_layer=nedge_layer+2
    end if
! 4. top and bottom edge on west facing face of cell
    if(edge_on_cell(WB,lid)==0)then
      edge_on_cell(WB,lid) = nedge_layer+1
      edge_on_cell(WT,lid)= nedge_layer+2
      if(cell_next(W,lid)>0)then                     ! and east facing face of cell to west
        edge_on_cell(EB,cell_next(W,lid)) = nedge_layer+1
        edge_on_cell(ET,cell_next(W,lid))= nedge_layer+2
      end if
      nedge_layer=nedge_layer+2
    end if
! 5. vertical edge at south west corner of cell
    if(edge_on_cell(SW,lid)==0)then 
      nedge_layer=nedge_layer+1
      edge_on_cell(SW,lid) = nedge_layer
      if(cell_next(W,lid)>0)then                     ! and at south east corner of cell to west 
        edge_on_cell(SE,cell_next(W,lid)) = nedge_layer
        if(cell_next(S,cell_next(W,lid))>0)then      ! and at north east corner of cell to south west
          edge_on_cell(NE,cell_next(S,cell_next(W,lid))) = nedge_layer
        end if
      end if
      if(cell_next(S,lid)>0)then                     ! and at north west corner of cell to south
        edge_on_cell(NW,cell_next(S,lid)) = nedge_layer
        if(cell_next(W,cell_next(S,lid))>0)then      ! and again at north east corner of cell to south west (in case other route to southwest goes through a missing cell)
          edge_on_cell(NE,cell_next(W,cell_next(S,lid))) = nedge_layer
        end if
      end if
    end if
! 6. vertical edge at south east corner of cell
    if(edge_on_cell(SE,lid)==0)then 
      nedge_layer=nedge_layer+1
      edge_on_cell(SE,lid) = nedge_layer
      if(cell_next(E,lid)>0)then                     ! and at south west corner of cell to east 
        edge_on_cell(SW,cell_next(E,lid)) = nedge_layer
        if(cell_next(S,cell_next(E,lid))>0)then      ! and north west corner of cell to south east
          edge_on_cell(NW,cell_next(S,cell_next(E,lid))) = nedge_layer
        end if
      end if
      if(cell_next(S,lid)>0)then                     ! and at north east corner of cell to south
        edge_on_cell(NE,cell_next(S,lid)) = nedge_layer
        if(cell_next(E,cell_next(S,lid))>0)then      ! and again at north west corner of cell to south east (in case other route to southeast goes through a missing cell)
          edge_on_cell(NW,cell_next(E,cell_next(S,lid))) = nedge_layer
        end if
      end if
    end if
! 7. vertical edge at north east corner of cell
    if(edge_on_cell(NE,lid)==0)then 
      nedge_layer=nedge_layer+1
      edge_on_cell(NE,lid) = nedge_layer
      if(cell_next(E,lid)>0)then                     ! and at north west corner of cell to east 
        edge_on_cell(NW,cell_next(E,lid)) = nedge_layer
        if(cell_next(N,cell_next(E,lid))>0)then      ! and at south west corner of cell to north east
          edge_on_cell(SW,cell_next(N,cell_next(E,lid))) = nedge_layer
        end if
      end if
      if(cell_next(N,lid)>0)then                     ! and at south east corner of cell to north
        edge_on_cell(SE,cell_next(N,lid)) = nedge_layer
        if(cell_next(E,cell_next(N,lid))>0)then      ! and again at south west corner of cell to north east (in case other route to northeast goes through a missing cell)
          edge_on_cell(SW,cell_next(E,cell_next(N,lid))) = nedge_layer
        end if
      end if
    end if
! 8. vertical edge at north west corner of cell
    if(edge_on_cell(NW,lid)==0)then 
      nedge_layer=nedge_layer+1
      edge_on_cell(NW,lid) = nedge_layer
      if(cell_next(W,lid)>0)then                     ! and at north east corner of cell to west 
        edge_on_cell(NE,cell_next(W,lid)) = nedge_layer
        if(cell_next(N,cell_next(W,lid))>0)then      ! and at south east corner of cell to north west
          edge_on_cell(SE,cell_next(N,cell_next(W,lid))) = nedge_layer
        end if
      end if
      if(cell_next(N,lid)>0)then                     ! and at south west corner of cell to nortth
        edge_on_cell(SW,cell_next(N,lid)) = nedge_layer
        if(cell_next(W,cell_next(N,lid))>0)then      ! and again at south east corner of cell to north west (in case other route to northwest goes through a missing cell)
          edge_on_cell(SE,cell_next(W,cell_next(N,lid))) = nedge_layer
        end if
      end if
    end if
  end do

  call log_event( 'faces on cells', LOG_LEVEL_DEBUG )
  do i = 1, ncells
    write( log_scratch_space, '(7i6)' ) i, &
                              face_on_cell(S,i), face_on_cell(E,i), &
                              face_on_cell(N,i), face_on_cell(W,i), &
                              face_on_cell(B,i), face_on_cell(T,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call log_event( 'edges on cells', LOG_LEVEL_DEBUG )
  do i = 1, ncells
    write( log_scratch_space, '(13i6)' ) i, &
                              edge_on_cell(SB,i), edge_on_cell(EB,i), &
                              edge_on_cell(NB,i), edge_on_cell(WB,i), &
                              edge_on_cell(SW,i), edge_on_cell(SE,i), &
                              edge_on_cell(NE,i), edge_on_cell(NW,i), &
                              edge_on_cell(ST,i), edge_on_cell(ET,i),&
                              edge_on_cell(NT,i),edge_on_cell(WT,i)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

end subroutine mesh_connectivity

!-------------------------------------------------------------------------------
      
subroutine xyz2llr(x,y,z,long,lat,r)
  
!-------------------------------------------------------------------------------
!  Subroutine to convert cartesian coordinates to longitude and latitude
!------------------------------------------------------------------------------- 

  use constants_mod, only: PI


  real(kind=r_def), intent(in)  :: x, y, z
  real(kind=r_def), intent(out) :: long, lat, r
  real(kind=r_def)              :: tln, tlt
  real(kind=r_def)              :: tol = 10.0e-8_r_def

  if ( abs(x) < tol ) then
    if ( y >= 0.0_r_def ) then
      long = 0.5_r_def*PI
    else
      long = 1.5_r_def*PI
    end if
  else
    tln = y/x
    long = atan(tln)
    if ( x < 0.0_r_def ) then
      long = long + PI
    end if
    if ( long < 0.0_r_def ) then
      long = long + 2.0_r_def*PI
    end if
  end if

  r = sqrt(x*x + y*y)
  if ( abs(r) < tol ) then
    if (z > 0.0_r_def ) then
      lat =  0.5_r_def*PI
    else
      lat = -0.5_r_def*PI
    end if
  else
    tlt = z/r
    lat = atan(tlt)
  end if
  
  r = sqrt(x*x + y*y + z*z)  

end subroutine xyz2llr

!-------------------------------------------------------------------------------
!> @brief Subroutine to return the vertex coordinates for a column of cells
!> @param[in] cell the horizontal cell index
!> @param[in] ncells the number of horizontal cells
!> @param[in] nlayers the number of vertical layers
!> @param[out] vert_coords the coordinates of the vertices
subroutine get_cell_coords(cell, ncells, nlayers, vert_coords)

!-------------------------------------------------------------------------------
!  Subroutine to return coordinates of vertices on cell 
!------------------------------------------------------------------------------- 

  integer, intent(in) :: cell, ncells, nlayers
  real(kind=r_def), intent(out) :: vert_coords(3,nverts,nlayers)
  
  integer :: i, j, k, cell_idx, v
  
  do k = 1,nlayers
    cell_idx = cell+(k-1)*ncells
    do j = 1,nverts
      v = vert_on_cell(j,cell_idx)
      do i = 1,3      
        vert_coords(i,j,k) = mesh_vertex(i,v)
      end do
    end do
  end do

end subroutine get_cell_coords

!-------------------------------------------------------------------------------
!> @brief Subroutine to compute the domain limits (x,y,z) for Cartesian domains,
!>        and (lambda,phi,r) for spherical domains
subroutine set_domain_size()

  use constants_mod, only: PI
  use mesh_mod,      only: l_spherical
  
  if ( l_spherical ) then
    domain_size%minimum%x =  0.0_r_def 
    domain_size%maximum%x =  2.0_r_def*PI
    domain_size%minimum%y = -0.5_r_def*PI
    domain_size%maximum%y =  0.5_r_def*PI
    domain_size%minimum%z =  0.0_r_def
    domain_size%maximum%z =  domain_top
  else
    domain_size%minimum%x =  minval(mesh_vertex(1,:))
    domain_size%maximum%x =  maxval(mesh_vertex(1,:))
    domain_size%minimum%y =  minval(mesh_vertex(2,:))
    domain_size%maximum%y =  maxval(mesh_vertex(2,:))
    domain_size%minimum%z =  minval(mesh_vertex(3,:))
    domain_size%maximum%z =  maxval(mesh_vertex(3,:))  
  end if

!> @todo Need to do a global reduction of maxs and mins when the code is parallel

end subroutine set_domain_size

!-------------------------------------------------------------------------------
!> @brief Function to convert a flux vector in cartesian coordinates (x,y,z) 
!> to one in spherical coodinates (lambda,phi,r)
!> @detail convert 3d cartesian velocity u = (u1,u2,u3) in cartesian coordinates at (x,y,z) 
!> to spherical velocity v = (v1,v2,v3) (in m/s) in spherical coordinates at (lambda,phi,r)
!> @param[in]  x_vec vector (x,y,z) location in cartesian coodinates
!> @param[in]  cartesian_vec components of a flux vector (u,v,w) in cartesian coordinates
!> @param[out] spherical_vec components of a flux vector (u,v,w) on spherical coordinates
pure function cart2sphere_vector(x_vec, cartesian_vec) result ( spherical_vec )
     use constants_mod, only: r_def, PI

     implicit none

     real(kind=r_def), intent(in)  :: x_vec(3)
     real(kind=r_def), intent(in)  :: cartesian_vec(3)
     real(kind=r_def)              :: spherical_vec(3)

     real(kind=r_def) :: r, t, phi

     t = x_vec(1)**2 + x_vec(2)**2
     r = sqrt(t + x_vec(3)**2)

     spherical_vec(1) = (- x_vec(2)*cartesian_vec(1) &
                         + x_vec(1)*cartesian_vec(2) ) / t
     spherical_vec(2) = (-x_vec(1)*x_vec(3)*cartesian_vec(1) &
                         -x_vec(2)*x_vec(3)*cartesian_vec(2) &
                                        + t*cartesian_vec(3))/(r*r*sqrt(t))
     spherical_vec(3) = (x_vec(1)*cartesian_vec(1) &
                       + x_vec(2)*cartesian_vec(2) &
                       + x_vec(3)*cartesian_vec(3)) / r

! convert from (dlambda/dt,dphi/dt,dr/dt) to (u,v,w) in m/s
     phi = 0.5_r_def*PI - acos(x_vec(3)/r)
     spherical_vec(1) = spherical_vec(1)*r*cos(phi)
     spherical_vec(2) = spherical_vec(2)*r

end function cart2sphere_vector

!-------------------------------------------------------------------------------
!> @brief Function to convert a vector in sperical coordinates to one in
!> cartesian coodinates
!> @param[in]  llr location in spherical coodinates
!> @param[in]  dlambda input vector in spherical coordinates
!> @param[return] dx output vector in cartesian coordinates
function sphere2cart_vector( dlambda, llr ) result ( dx )
  use constants_mod,     only: r_def
  use matrix_invert_mod, only: matrix_invert_3x3
  implicit none

  real(kind=r_def), intent(in)  :: dlambda(3)
  real(kind=r_def), intent(in)  :: llr(3)
  real(kind=r_def)              :: dx(3)

  real(kind=r_def)              :: A(3,3), A_inv(3,3)

! Form transformation matrix
  A(1,:) = (/ -sin(llr(1)),              cos(llr(1)),              0.0_r_def   /)
  A(2,:) = (/ -sin(llr(2))*cos(llr(1)), -sin(llr(2))*sin(llr(1)), -cos(llr(2)) /)
  A(3,:) = (/  cos(llr(2))*cos(llr(1)),  cos(llr(2))*sin(llr(1)), -sin(llr(2)) /)
! form inverse
  A_inv = matrix_invert_3x3(A)
  dx(:) = matmul(A_inv, dlambda)
  return

end function sphere2cart_vector

!> @brief Returns the local index of the cell on the local
!> partition that corresponds to the given global index.
!> @param[in] gid The global index to search for on the local partition
!> @param[return] lid The local index that corresponds to the given global index
!> or -1 if the cell with the given global index is not present of the local partition
pure function get_lid_from_gid(gid) result (lid)
!
! Perform a search through the global cell lookup table looking for the
! required global index.
!
! The partitioned_cells array holds global indices in three groups: core cells,
! followed by owned cells and finally, the indices of the halo cells. The
! cells are numerically ordered within the different groups so a binary search
! can be used, but not between groups, so need to do separate binary searches
! through the core, owned and halo cells and exit if a match is found
!
  use mesh_mod, only : partitioned_cells, num_core, num_owned, num_halo, &
                       num_cells_x, num_cells_y
  implicit none

  integer, intent(in) :: gid           ! global index
  integer             :: lid           ! local index
  integer             :: nlayer        ! layer of supplied gid
  integer             :: gid_in_layer  ! supplied gid projected to bottom layer

  ! Set the default return code
  lid=-1
  ! If the supplied gid is not valid just return
  if(gid < 1) return

  ! The global index lookup table (partitioned_cells) only has the indices for
  ! a single layer, so convert the full 3d global index into the global index
  ! within the layer and a layer number
  gid_in_layer=modulo(gid-1,(num_cells_x*num_cells_y))+1
  nlayer=(gid-1)/(num_cells_x*num_cells_y)

  ! Search though core cells - partitioned_cells(1:num_core) - looking for the gid
  lid=binary_search(partitioned_cells( 1:num_core ), gid)
  if(lid /= -1)then
    lid=lid+nlayer*(num_core+num_owned+num_halo)  !convert back to 3d lid
    return
  end if

  ! Search though owned cells - partitioned_cells(num_core+1:num_core+num_owned) - looking for the gid
  lid=binary_search(partitioned_cells( num_core+1:num_core+num_owned ), gid)
  if(lid /= -1)then
    lid=lid+num_core+nlayer*(num_core+num_owned+num_halo)  !convert back to 3d lid
    return
  end if

  ! Search though halo cells - partitioned_cells(num_core+num_owned+1:num_core+num_owned+num_halo) - looking for the gid
  lid=binary_search(partitioned_cells( num_core+num_owned+1:num_core+num_owned+num_halo ), gid)
  if(lid /= -1)then
    lid=lid+num_core+num_owned+nlayer*(num_core+num_owned+num_halo)  !convert back to 3d lid
    return
  end if

  ! No lid has been found in either the core, owned or halo cells on this partition, so return with lid=-1
  return
  
end function get_lid_from_gid

!> @brief Performs a binary search through the given array looking for
!> a particular entry and returns the index of the entry found
!> or -1 if no matching entry can be found. The values held in "array_to_be_searched"
!> must be in numerically increasing order
!> @param[in] array_to_be_searched The array that will be searched for the given entry
!> @param[in] value_to_find The entry that is to be searched for
pure function binary_search(array_to_be_searched, value_to_find) result (id)

  implicit none

  integer, intent(in) :: array_to_be_searched( : )
  integer, intent(in) :: value_to_find
  integer :: bot, top  ! Lower and upper index limits between which to search for the value
  integer :: id        ! Next index for refining the search. If an entry is found this will
                       ! contain the index of the matching entry

  ! Set bot and top to be the whole array to begin with
  bot=1
  top=size(array_to_be_searched)

  search: do
    ! If top is lower than bot then there is no more array to be searched
    if(top < bot) exit search
    ! Refine the search
    id=(bot+top)/2
    if(array_to_be_searched(id) == value_to_find)then  ! found matching entry
      return
    else if(array_to_be_searched(id) < value_to_find)then ! entry has to be between id and top
      bot = id + 1
    else ! entry has to be between bot and id
      top = id - 1
    endif
  end do search

  ! Didn't find a match - return failure code
  id=-1

end function binary_search

end module mesh_generator_mod
