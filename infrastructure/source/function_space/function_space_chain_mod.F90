!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> Module for function_space_chain_type which defines an ordered list of pointers
!> to function_space_type objects.
!>
!> Typical usage is to instantiate a function_space_chain_type,
!> then add pointers to function_space_type objects in the chain order
!> required.
!>
!> e.g. MultiGrid, Physics
!>
!> Entry function space (Head of Chain)<br>
!> Primal mesh <-> Coarsen(Primal)x1 <-> Coarsen(Primal)x2
!> <-> Coarsen(Primal)x3<br>
!> Primal mesh <-> Physics mesh.
!>

module function_space_chain_mod

use constants_mod,              only: i_def, l_def
use linked_list_mod,            only: linked_list_type, &
                                      linked_list_item_type
use linked_list_int_mod,        only: linked_list_int_type
use linked_list_data_mod,       only: linked_list_data_type
use log_mod,                    only: log_event, log_scratch_space, &
                                      LOG_LEVEL_ERROR, LOG_LEVEL_DEBUG
use function_space_mod,         only: function_space_type
use mesh_mod,                   only: mesh_type
use mesh_map_collection_mod,    only: mesh_map_collection_type
use function_space_pointer_mod, only: function_space_pointer_type

implicit none

private

!===============================================================================
! Type which is an ordered collection of function space pointers held in
! a linked list

type, public :: function_space_chain_type
  private

  !> Linked list of function_space_pointer_type objects
  type(linked_list_type) :: function_space_chain_list

  ! Despite not having a user finaliser this still seems to be necessary.
  integer(i_def), allocatable :: dummy_for_gnu

contains

  !> @brief Adds a given function space to the next location in this chain.
  !> @param[in] function_space Pointer to function space to add to this chain.
  procedure, public :: add
  !>
  !> @brief   Returns a pointer to the function space at the start of this
  !>          chain.
  !> @details This function will also update the current function space
  !>          in this chain to be the function space at the head of the chain.
  !> @return  start_function_space Pointer to function space at
  !>                               Entry/Exit point of this chain.
  procedure, public :: get_start
  !>
  !> @brief   Returns pointer to the next function space in this chain.
  !> @details This function will increment the current position in this
  !>          chain tailwards (forward) by 1.
  !> @return  next_function_space Pointer to the next function space in
  !>                              this chain.
  procedure, public :: get_next
  !>
  !> @brief   Returns pointer to the previous function space before the current
  !>          function space in this chain.
  !> @details This function will increment the current position in this chain
  !>          headwards (backward) by 1.
  !> @return  previous_function_space Pointer to the previous function
  !>                                  space in this chain
  procedure, public :: get_previous

  !> @brief Check if a chain item exists
  !> @param[in] fs_id ID of function space to check
  !> @return answer Flag showing if item exists
  procedure, public :: exists

  !> @brief Set the current item in the list to point to a given ID
  !> @param[in] fs_id ID of function space to set to current
  procedure, public :: set_current
  !> @brief   Manually release memory used by this object.
  procedure, public :: clear

  final :: function_space_chain_destructor

end type function_space_chain_type

interface function_space_chain_type
  module procedure function_space_chain_constructor
end interface

type(function_space_chain_type), public, allocatable :: &
                                           single_layer_function_space_chain
type(function_space_chain_type), public, allocatable :: &
                                           multigrid_function_space_chain
type(function_space_chain_type), public, allocatable :: &
                                           W2_multigrid_function_space_chain
type(function_space_chain_type), public, allocatable :: &
                                           Wtheta_multigrid_function_space_chain

!===============================================================================
contains ! Module procedures


!===============================================================================
!> @brief Constructor for function space chain objects
!>        (function_space_chain_type).
!> @return instance A function space chain object
!===============================================================================
function function_space_chain_constructor() result(instance)

implicit none

type(function_space_chain_type) :: instance

instance%function_space_chain_list = linked_list_type()

return
end function function_space_chain_constructor


!===============================================================================
subroutine add(self, function_space)

implicit none

class(function_space_chain_type) :: self
type(function_space_type), intent(inout), pointer :: function_space

type(function_space_pointer_type)          :: this_function_space_pointer
type(function_space_pointer_type), pointer :: previous_function_space_pointer => null()

type(linked_list_item_type), pointer :: previous_list_item      => null()
type(function_space_type),   pointer :: previous_function_space => null()


previous_list_item => self%function_space_chain_list%get_tail()

if ( associated(previous_list_item) ) then

  ! A mesh map "may" need to be created.
  ! 'cast' to the function_space_pointer_type
  select type(v => previous_list_item%payload)
  type is (function_space_pointer_type)
    previous_function_space_pointer => v
  end select

  previous_function_space => previous_function_space_pointer%get_target()

  call create_function_space_chain_mesh_maps( previous_function_space, &
                                              function_space )
end if

this_function_space_pointer = function_space_pointer_type(function_space)
call self%function_space_chain_list%insert_item(this_function_space_pointer)

nullify(previous_function_space_pointer)
nullify(previous_list_item)
nullify(previous_function_space)

return
end subroutine add


!===============================================================================
function get_start(self) result (start_function_space)

implicit none

class(function_space_chain_type) :: self
type(function_space_type), pointer :: start_function_space

class(linked_list_item_type), pointer :: head => null()

type(function_space_pointer_type), pointer :: function_space_pointer => null()

start_function_space => null()
head => self%function_space_chain_list%get_head()

if (associated(head)) then
  call self%function_space_chain_list%set_current(head)
else
  ! The head of this linked list is not associated, so there are no
  ! items in this list. Given that the chain must have at least two
  ! different function spaces something must have gone wrong.
  write(log_scratch_space,'(A)') "No items in this function space chain."
  call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  return
end if

! 'cast' to the function_space_type
select type(v => head%payload)
type is (function_space_pointer_type)
  function_space_pointer => v
end select

start_function_space => function_space_pointer%get_target()

nullify( head )
nullify( function_space_pointer )

return
end function get_start


!===============================================================================
function get_next(self) result (next_function_space)

implicit none

class(function_space_chain_type) :: self
type(function_space_type), pointer :: next_function_space

class(linked_list_item_type), pointer :: next => null()
type(linked_list_item_type),  pointer :: tmp  => null()

type(function_space_pointer_type), pointer :: function_space_pointer => null()

next_function_space => null()
tmp  => self%function_space_chain_list%get_current()
next => tmp%next

if (associated(next)) then
  call self%function_space_chain_list%set_current(next)
else
  ! There is no "next" function space so this is most likely at
  ! the end of the function space chain
  write(log_scratch_space,'(A)') &
      "End of function_space_chain; No more function spaces."
  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
  return
end if


! 'cast' to the function_space_type
select type(v => next%payload)
type is (function_space_pointer_type)
  function_space_pointer => v
end select

next_function_space => function_space_pointer%get_target()

nullify( next )
nullify( tmp )
nullify( function_space_pointer )
return
end function get_next


!==============================================================================
function get_previous(self) result (previous_function_space)

implicit none

class(function_space_chain_type) :: self
type(function_space_type), pointer :: previous_function_space

class(linked_list_item_type), pointer :: previous => null()
type(linked_list_item_type),  pointer :: head => null()
type(linked_list_item_type),  pointer :: tmp  => null()

type(function_space_pointer_type), pointer :: function_space_pointer => null()

previous_function_space => null()
head => self%function_space_chain_list%get_head()
tmp  => self%function_space_chain_list%get_current()
previous => tmp%prev


if (associated(previous)) then
  call self%function_space_chain_list%set_current(previous)
else
  ! There is no "previous" function space so this is most likely at
  ! the start of the function space chain
  if ( associated( tmp, head) ) then
    write(log_scratch_space,'(A)') &
        "At head of function_space_chain; No previous function spaces."
  else
    write(log_scratch_space,'(A)') "No previous function spaces."
  end if

  call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

  return
end if


! 'cast' to the function_space_type
select type(v => previous%payload)
type is (function_space_pointer_type)
  function_space_pointer => v
end select

previous_function_space => function_space_pointer%get_target()

nullify( previous )
nullify( head )
nullify( tmp  )
nullify( function_space_pointer )

return
end function get_previous

function exists(self, fs_id) result (answer)

  implicit none

  class(function_space_chain_type), intent(in) :: self
  integer(i_def), intent(in) :: fs_id
  logical(l_def) :: answer

  answer = self%function_space_chain_list%item_exists(fs_id)

  return

end function exists
!==============================================================================


subroutine set_current(self, fs_id)

  implicit none

  class(function_space_chain_type), intent(inout) :: self
  integer(i_def), intent(in) :: fs_id
  class(linked_list_item_type), pointer   :: item => null()

  item => self%function_space_chain_list%get_item(fs_id)
  if(associated( item) ) then
     call self%function_space_chain_list%set_current(item)
  else
     write(log_scratch_space,'(A,I0,A)') &
     "function_space_chain_mod:set_current: fs_id=",fs_id, " not found in chain"
     call log_event(log_scratch_space,LOG_LEVEL_ERROR)
  end if

  return
end subroutine set_current


!==============================================================================
!  @brief Subroutine to clear up objects - called by destructor
!         Explcitly deallocates any allocatable arrays in the object
!         to avoid memory leaks
!> @param[inout] self, The calling function_space_chain instance
subroutine clear(self)

implicit none

class(function_space_chain_type), intent(inout) :: self

call self%function_space_chain_list%clear()

if (allocated(self%dummy_for_gnu)) deallocate(self%dummy_for_gnu)

return
end subroutine clear



!-------------------------------------------------------------------------------
!> @brief Finalizer routine which should automatically call clear
!>        when object is out of scope.
!> @param[inout] self, The calling function_space_chain instance
subroutine function_space_chain_destructor(self)

  implicit none

  type(function_space_chain_type), intent(inout) :: self

  call self%clear()

  return
end subroutine function_space_chain_destructor




!==============================================================================
!> @brief   Generates required mesh maps for the given function space source
!>          and target.
!> @details This routine is provide to generate mesh maps for partitioned meshes
!>          as they are added to the chain. In future, these maps will be read in
!>          from the NetCDF UGRID file.
!> @param[in] source_function_space Pointer to the source function space.
!> @param[in] target_function_space Pointer to the target function space.
!==============================================================================
subroutine create_function_space_chain_mesh_maps( source_function_space, &
                                                  target_function_space )

implicit none

type(function_space_type), intent(in), pointer :: source_function_space
type(function_space_type), intent(in), pointer :: target_function_space

type(mesh_type), pointer :: source_mesh => null()
type(mesh_type), pointer :: target_mesh => null()

integer(i_def) :: source_mesh_id
integer(i_def) :: target_mesh_id


if ( associated(source_function_space) .and. &
     associated(target_function_space) ) then

  ! Get the local meshes associated with consecutive function spaces
  source_mesh => source_function_space%get_mesh()
  target_mesh => target_function_space%get_mesh()

  ! The local and global mesh ids associated with the source and
  ! target meshes
  source_mesh_id = source_mesh%get_id()
  target_mesh_id = target_mesh%get_id()

  ! Mapping is only necessary if consecutive function spaces are on
  ! different meshes.
  if (source_mesh_id /= target_mesh_id) then

    call source_mesh%add_mesh_map(target_mesh)
    call target_mesh%add_mesh_map(source_mesh)

  end if

else
  write(log_scratch_space,'(A)') &
      'Unassociated pointers to Source or Target function spaces.'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

nullify(source_mesh)
nullify(target_mesh)

return
end subroutine create_function_space_chain_mesh_maps

!===============================================================================
end module function_space_chain_mod
