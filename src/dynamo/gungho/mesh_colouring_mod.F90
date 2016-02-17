!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!------------------------------------------------------------------------------
!> @brief Computes mesh colouring for vector spaces
!>
!> @details Contains algorithms for colouring of meshes according to different
!>          traversal orders and different colouring policies.
!>          At present, greedy and balanced colouring are available for each
!>          mesh type.
!------------------------------------------------------------------------------
module mesh_colouring_mod
!------------------------------------------------------------------------------
  use reference_element_mod,   only : W, S, E, N, nfaces_h
  use log_mod,                 only : log_event, LOG_LEVEL_ERROR, &
                                                 LOG_LEVEL_INFO
  use constants_mod,           only : i_def
!------------------------------------------------------------------------------
  implicit none
!------------------------------------------------------------------------------
  integer, parameter               :: MAXCOLS = 50 ! Temporary hardcode until
                                                   ! dynamic palette
  integer(kind=i_def), target      :: cells_per_colour(MAXCOLS)
!------------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------
!>  @brief  Set colouring pattern for vector spaces
!>
!>  @detail Acquires details of cells in partition and applies colouring
!>          to this space according to the specified colouring policy.
!>
!>  @param[in]  num_cells   Number of cells defined on the mesh.
!>  @param[in]  cell_next   Adjacency array for the mesh.
!>  @param[out] num_colours   Number of colours used.
!>  @param[out] num_cell_per_colour   Count of cells in each colour
!>  @param[out] cells_in_colour   List of cell indices in each colour
!------------------------------------------------------------------------------
subroutine set_colours(num_cells,              &
                       cell_next,              &
                       num_colours,            &
                       num_cell_per_colour,    &
                       cells_in_colour)
  implicit none
  integer(i_def), intent(in)                    :: num_cells
  integer(i_def), allocatable, intent(in)       :: cell_next(:,:)
  integer, intent(out)                          :: num_colours
  integer(kind=i_def), allocatable, intent(out) :: num_cell_per_colour(:)
  integer(kind=i_def), allocatable, intent(out) :: cells_in_colour(:,:)


  ! Array holding the colour of each cell. Index 0 cell holds 0 value for
  ! elements having fewer than maximum count of neighbours.
  integer(kind=i_def), allocatable        :: colour_map(:)

  ! Array for marking used (unavailable) colours for a cell
  integer(kind=i_def)                     :: used_colours(0:MAXCOLS)
  ! The next available colour
  integer(kind=i_def)                     :: free_colour
  ! Loop and status variables
  integer                                 :: i, astat
  integer(kind=i_def)                     :: colour, test_cell, cell, maxcell
  ! Prefix for allocation error messages
  character(len=*), parameter :: prefix="[Compute Colours] Failure to allocate "

  ! Determine largest neighbour element that will be a subscript to cell_next
  maxcell = maxval(cell_next(:nfaces_h, :num_cells))

  allocate(colour_map(0:maxcell), stat=astat)
  if(astat/=0) then
    call log_event(prefix//"colour_map.", LOG_LEVEL_ERROR)
    ! Although the logger will call "stop" for errors GFortran is unable to
    ! perform inter-file analysis so believes it is possible to use colour_map
    ! without initialising it.
    ! The following "stop" will never be reached but it lets GFortran know
    ! that colour_map will always be initialised.
    stop
    !> @todo Keep an eye on GFortran development and remove uncalled "stop"
    !>       when it is no longer necessary to satisfy "Uninitialised variable
    !>       may be used" warning. This warning is only seen for the
    !>       "production" target.
  end if

  colour_map = 0_i_def
  cells_per_colour = 0_i_def

  do cell = 1, num_cells
    do i=0, MAXCOLS
      used_colours(i) = 0_i_def
    end do

    test_cell = cell_next(N, cell)
    if(test_cell > 0) then
      used_colours(colour_map(test_cell)) = 1

      ! Cardinal point orientation may vary when crossing a panel boundary
      ! for cubed-sphere mesh types.
      if(cell == cell_next(S, test_cell)) then
        used_colours(colour_map(cell_next(E, test_cell))) = 1
        used_colours(colour_map(cell_next(W, test_cell))) = 1
      else
        used_colours(colour_map(cell_next(N, test_cell))) = 1
        used_colours(colour_map(cell_next(S, test_cell))) = 1
      end if
    end if

    test_cell = cell_next(E, cell)
    if(test_cell > 0) used_colours(colour_map(test_cell)) = 1

    test_cell = cell_next(W, cell)
    if(test_cell > 0) used_colours(colour_map(test_cell)) = 1

    test_cell = cell_next(S, cell)
    if(test_cell > 0) then
      used_colours(colour_map(test_cell)) = 1

      ! Cardinal point orientation may vary when crossing a panel boundary
      ! for cubed-sphere mesh types.
      if(cell == cell_next(N, test_cell)) then
        used_colours(colour_map(cell_next(E, test_cell))) = 1
        used_colours(colour_map(cell_next(W, test_cell))) = 1
      else
        used_colours(colour_map(cell_next(N, test_cell))) = 1
        used_colours(colour_map(cell_next(S, test_cell))) = 1
      end if
    end if

    free_colour = choose_colour_greedy(used_colours)
    ! Alternate colour choice procedures may be applied here, e.g...
    ! free_colour = choose_colour_balanced(used_colours)
    colour_map(cell) = free_colour
    cells_per_colour(free_colour) = cells_per_colour(free_colour) + 1
  end do

  num_colours = MAXCOLS
  ! Allocate return data and populate
  do i = 1, MAXCOLS
    if(cells_per_colour(i) == 0) then
      num_colours = i-1
      exit
    end if
  end do

  allocate(num_cell_per_colour(num_colours), stat=astat)

  if(astat/=0) call log_event(prefix//"num_cell_per_colour.", &
                                LOG_LEVEL_ERROR)

  do colour = 1, num_colours
    num_cell_per_colour(colour) = cells_per_colour(colour)
  end do

  allocate(cells_in_colour(num_colours, maxval(num_cell_per_colour)), &
           stat=astat)
  if(astat/=0) call log_event(prefix//"cells_in_colour.", &
                                LOG_LEVEL_ERROR)

  cells_in_colour = 0

  do colour = 1, num_colours
    i = 0
    do cell = 1, num_cells
      if (colour_map(cell) == colour) then 
        i = i+1
        cells_in_colour(colour, i) = cell
      end if
    end do
  end do


end subroutine set_colours
!-----------------------------------------------------------------------------
!> @brief  Display colour map on console
!>
!> @detail Subroutine to print colour map to stdout; useful for debugging.
!> 
!> @param[in] colour_map     array of colours assigned by set_colours()
!> @param[in] display_width  Optional integer number of cells per line
!------------------------------------------------------------------------------
subroutine write_colours(colour_map, display_width)

  integer, intent(in)             :: colour_map(0:)
  integer, intent(in), optional   :: display_width

  integer                         :: i, num_cells, width
  character(len=32)               :: varfmt
  character(len=512)              :: log_buf

  if(present(display_width)) then
    width = display_width
  else
    width = 10
  end if

  write(varfmt,*) width
  num_cells = size(colour_map)-1

  write(log_buf, "(A)") "Colours(0:num_cells):"
  call log_event(log_buf, LOG_LEVEL_INFO)
  do i=0, num_cells, width
    write(log_buf, "(I5,T16,"//trim(adjustl(varfmt))//"I5)") i/width*width, &
                    colour_map(i:min(num_cells, i+width-1))
    call log_event(log_buf, LOG_LEVEL_INFO)
  end do

end subroutine write_colours
!-----------------------------------------------------------------------------
!> @brief  Choose next available colour using greedy algorithm
!> @detail Function to select the next available colour; creates an imbalanced
!>         palette using a minimal number of colours.
!>
!> @param[in] used_colours array of colours already used
!> @return    colour       integer signifying the chosen colour
!------------------------------------------------------------------------------
pure function choose_colour_greedy(used_colours)  result(colour)

  integer, intent(in)     :: used_colours(0:MAXCOLS)
  integer                 :: colour
  integer                 :: idx

  colour = 0
  do idx=1, MAXCOLS
    if(used_colours(idx) == 0) then
      colour = idx
      exit
    endif
  end do

end function choose_colour_greedy
!-----------------------------------------------------------------------------
!> @brief   Choose next available colour using balanced algorithm
!-----------------------------------------------------------------------------
!> @detail  Function to select the next available colour; creates a balanced
!>          palette using all available colours. 
!>
!> @param[in] used_colours array of colours already used
!> @return    colour       integer signifying the chosen colour
!------------------------------------------------------------------------------
pure function choose_colour_balanced(used_colours)  result(colour)

  integer, intent(in)     :: used_colours(0:MAXCOLS)
  integer                 :: colour

  colour = minloc(cells_per_colour, 1, mask=(used_colours(1:MAXCOLS)==0))

end function choose_colour_balanced
!------------------------------------------------------------------------------
end module mesh_colouring_mod
!------------------------------------------------------------------------------
