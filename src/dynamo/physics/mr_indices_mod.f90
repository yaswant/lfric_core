!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!> @brief Define indices for the mixing ratio field vectors
!>
!> @details Define and set indices for the mixing ratio field vectors
module mr_indices_mod

  implicit none

  integer, parameter :: imr_v  = 1  ! vapour
  integer, parameter :: imr_c  = 2  ! cloud mass
  integer, parameter :: imr_r  = 3  ! rain mass
  integer, parameter :: imr_nc = 4  ! cloud number concentration
  integer, parameter :: imr_nr = 5  ! rain number concentration
  integer, parameter :: nummr  = 5  ! Total number of mixing ratio variables

end module mr_indices_mod
