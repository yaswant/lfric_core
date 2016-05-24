!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Linked list data integer type

!> @details A type to hold ints for putting in a linked list
!>         inherits from linked_list_data_type

module linked_list_int_mod

  use linked_list_data_mod, only : linked_list_data_type

  implicit none

  private

  type, extends(linked_list_data_type), public        :: linked_list_int_type
    private
    contains
    ! Nothing in here - it's all in the base class
  end type linked_list_int_type

  interface linked_list_int_type
       module procedure int_constructor
  end interface linked_list_int_type

contains

! constructor
type(linked_list_int_type) function int_constructor(id)

   integer, intent(in)           :: id
   call int_constructor%set_id(id)

end function int_constructor


end module linked_list_int_mod
