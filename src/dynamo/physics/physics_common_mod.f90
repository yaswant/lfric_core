!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Collection of routines that are needed for physics. 

!> @detail Collection of routines that are needed for physics. These may be 
!>         replaced/overloaded by specific schemes as they are brought in from
!>         the UM

module physics_common_mod

  use constants_mod,                 only: r_def
  implicit none
  private

  public qsaturation

contains

  ! Function to return the saturation mr over water
  ! Based on tetans formula
  ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  function qsaturation (T, p)    
    real(r_def), intent(IN) :: T, p
    real(r_def) :: Qsaturation
    ! Temperature in Kelvin
    ! Pressure in mb
    real(r_def), parameter ::tk0c = 273.15, qsa1 = 3.8, &
       qsa2 = - 17.2693882, qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !
    if (T > qsa3 .and. p * exp (qsa2 * (t - tk0c) / (T - qsa3)) > qsa4) then
      qsaturation=qsa1/(p*exp(qsa2*(t-tk0c)/(T-qsa3))-qsa4)
    else
      qsaturation=999.0
    end if
  end function qsaturation

end module physics_common_mod
