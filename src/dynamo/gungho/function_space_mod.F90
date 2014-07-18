!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A type which holds information about the function space.

!> @details A container which holds type definition of the function space and 
!> has holds a number of static copies of the function spaces require by the
!> model. It provides accessor functions (getters) to various information weld
!> in the type

module function_space_mod

use constants_mod, only:dp

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
  
type, public :: function_space_type
  private
  integer              :: ndf, ncell, undf, ngp_h, ngp_v, fs
  integer              :: dim_space, dim_space_diff
  !> A two dimensional, allocatable array which holds the indirection map 
  !! or dofmap for the whole function space over the bottom level of the domain.
  integer, allocatable :: dofmap(:,:)
  !> 4-dim allocatable array of reals which hold the values of the basis function
  real(kind=dp), allocatable :: basis(:,:,:,:)
  !> 4-dim allocatable array of reals which hold the values of the basis function
  !! for the differential  functions space
  real(kind=dp), allocatable :: diff_basis(:,:,:,:)

  !> A two dimensional, allocatable array of reals which holds the coordinates
  !! of the function_space degrees of freedom
  real(kind=dp), allocatable :: nodal_coords(:,:)

contains
  !final :: destructor

  !> Function returns a pointer to a function space. If the required function
  !> space had not yet been created, it creates one before returning the pointer
  !> to it
  procedure, nopass :: get_instance
!> Function to get the total unique degrees of freedom for this space
!! returns an integer
!! @param[in] self the calling function space
  procedure :: get_undf

!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of cells
  procedure :: get_ncell

!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling functions_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the dofmap
  procedure :: get_cell_dofmap

!> Function which obtains the number of dofs per cell
!! @param[in] self The calling functions space
!! return an integer, the number of dofs per cell
  procedure :: get_ndf

!>  Accessor procedure for the basis function
!! @param[in] self the calling function space
!! @return A pointer to the array to hold the values of the basis function 
  procedure :: get_basis

!> Accessor procedure for the basis function for the
!! differential function space
!! @param[in] self the calling function space
!! @return A pointer to the real array to hold the values of the
!!  basis function 
  procedure :: get_diff_basis

!> Accessor function to get the coordinates of the function space
!! @return A pointer to the two dimensional array, (xyz,ndf)
  procedure :: get_nodes
  
  !> function returns the enumerated integer for the functions_space which
  !! is this function_space
  procedure :: which

end type function_space_type

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------
!> integer that defines the type of function space required
integer, public, parameter      :: V0 = 100
integer, public, parameter      :: V1 = 101
integer, public, parameter      :: V2 = 102
integer, public, parameter      :: V3 = 103

!> These are static copies of all the function spaces that will be required 
type(function_space_type), target, allocatable, save :: v0_function_space
type(function_space_type), target, allocatable, save :: v1_function_space
type(function_space_type), target, allocatable, save :: v2_function_space
type(function_space_type), target, allocatable, save :: v3_function_space

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public get_ncell, get_cell_dofmap, which
contains

function get_instance(function_space) result(instance)
  use basis_function_mod,         only : &
              v0_basis, v1_basis, v2_basis, v3_basis, &
              v0_diff_basis, v1_diff_basis, v2_diff_basis, v3_diff_basis, &
              v0_nodal_coords, v1_nodal_coords, v2_nodal_coords, v3_nodal_coords

  use dofmap_mod,                 only : &
              v0_dofmap, v1_dofmap, v2_dofmap, v3_dofmap
  use gaussian_quadrature_mod,    only : ngp_h, ngp_v
  use mesh_mod, only : num_cells, v_unique_dofs

  implicit none

  integer :: function_space
  type(function_space_type), pointer :: instance

  select case (function_space)
  case (V0)
    if(.not.allocated(v0_function_space)) then
      allocate(v0_function_space)   
      call init_function_space(self=v0_function_space, &
         num_cells = num_cells ,num_dofs = v_unique_dofs(1,2), &
         num_unique_dofs = v_unique_dofs(1,1) ,  &
         dim_space = 1, dim_space_diff = 3,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=v0_dofmap, &
         basis=v0_basis, diff_basis=v0_diff_basis, &
         nodal_coords=v0_nodal_coords, fs=V0) 
    end if
    instance => v0_function_space
  case (V1)
    if(.not.allocated(v1_function_space)) then
      allocate(v1_function_space)  
      call init_function_space(self=v1_function_space, &
         num_cells = num_cells ,num_dofs = v_unique_dofs(2,2), &
         num_unique_dofs = v_unique_dofs(2,1) ,  &
         dim_space = 3, dim_space_diff = 3,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=v1_dofmap, &
         basis=v1_basis, diff_basis=v1_diff_basis, &
         nodal_coords=v1_nodal_coords,fs=V1 )
    end if
    instance => v1_function_space
  case (V2)
    if(.not.allocated(v2_function_space)) then 
      allocate(v2_function_space)
      call init_function_space(self=v2_function_space, &
         num_cells = num_cells ,num_dofs = v_unique_dofs(3,2), &
         num_unique_dofs = v_unique_dofs(3,1) ,  &
         dim_space = 3, dim_space_diff = 1,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=v2_dofmap, &
         basis=v2_basis, diff_basis=v2_diff_basis, &
         nodal_coords=v2_nodal_coords,fs=V2 )
    end if
    instance => v2_function_space
  case (V3)
    if(.not.allocated(v3_function_space)) then
      allocate(v3_function_space)
      call init_function_space(self=v3_function_space, &
         num_cells = num_cells ,num_dofs = v_unique_dofs(4,2), &
         num_unique_dofs = v_unique_dofs(4,1) ,  &
         dim_space = 1, dim_space_diff = 1,  &
         ngp_h = ngp_h, ngp_v = ngp_v, &
         dofmap=v3_dofmap, &
         basis=v3_basis, diff_basis=v3_diff_basis, &
         nodal_coords=v3_nodal_coords, fs=V3 )
    end if
    instance => v3_function_space
  case default
    !not a recognised function space - return a null pointer
    instance => null()
  end select

  return
end function get_instance


!> Subroutine initialises a function space.
!! @param[in] num_cells
!! @param[in] num_dofs
!! @param[in] num_unique_dofs
!! @param[in] dim_space The dimension of this function space
!! @param[in] dim_space_diff The dimension of the differentiated function space
!! @param[in] ngp_h The number of guassian quadrature points in the horizonal
!! @param[in] ngp_v The number of guassian quadrature points in the vertical
subroutine init_function_space(self, &
                               num_cells,num_dofs, &
                               num_unique_dofs,  &
                               dim_space, dim_space_diff,  &
                               ngp_h,ngp_v, &
                               dofmap, &
                               basis, diff_basis, &
                               nodal_coords, fs)
  implicit none

  class(function_space_type) :: self
  integer, intent(in) :: num_cells, num_dofs, num_unique_dofs
  integer, intent(in) :: dim_space, dim_space_diff
  integer, intent(in) :: ngp_h,ngp_v
! The following four arrays have intent inout because the move_allocs in the
! code need access to the arrays to free them in their original locations
  integer, intent(inout), allocatable  :: dofmap(:,:)
  real(kind=dp), intent(inout), allocatable  :: basis(:,:,:,:)
  real(kind=dp), intent(inout), allocatable  :: diff_basis(:,:,:,:)
  real(kind=dp), intent(inout), allocatable  :: nodal_coords(:,:)
  integer, intent(in) :: fs

  self%ncell           =  num_cells
  self%ndf             =  num_dofs
  self%undf            =  num_unique_dofs
  self%dim_space       =  dim_space
  self%dim_space_diff  =  dim_space_diff
  self%ngp_h           =  ngp_h
  self%ngp_v           =  ngp_v  
  call move_alloc(dofmap, self%dofmap)
  call move_alloc(basis , self%basis)
  call move_alloc(diff_basis , self%diff_basis)
  call move_alloc(nodal_coords , self%nodal_coords) 
  self%fs              = fs
  return
end subroutine init_function_space

!-----------------------------------------------------------------------------
! Get total unique dofs for this space
!-----------------------------------------------------------------------------

!> Function to get the total unique degrees of freedom for this space
!! returns an integer
!! @param[in] self the calling function space
integer function get_undf(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_undf=self%undf

  return
end function get_undf

!-----------------------------------------------------------------------------
! Get the number of cells for this function space
!-----------------------------------------------------------------------------
!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of cells
integer function get_ncell(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_ncell=self%ncell

  return
end function get_ncell


!-----------------------------------------------------------------------------
! Get the number of dofs for a single cell 
!-----------------------------------------------------------------------------
integer function get_ndf(self)
  implicit none
  class(function_space_type), intent(in) :: self

  get_ndf=self%ndf

  return
end function get_ndf

!-----------------------------------------------------------------------------
! Get the dofmap for a single cell
!-----------------------------------------------------------------------------
!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling functions_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the dofmap
function get_cell_dofmap(self,cell) result(map)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer,                            intent(in) :: cell
  integer, pointer                               :: map(:)

  map => self%dofmap(cell,:)
  return
end function get_cell_dofmap

!-----------------------------------------------------------------------------
! Get the basis function
!-----------------------------------------------------------------------------
function get_basis(self)  result(basis)
  implicit none
  class(function_space_type), target, intent(in)  :: self  
  real(kind=dp),              pointer             :: basis(:,:,:,:)

  basis => self%basis

  return
end function get_basis

!-----------------------------------------------------------------------------
! Get the differential of the basis function
!-----------------------------------------------------------------------------
function get_diff_basis(self) result(diff_basis)
  implicit none
  class(function_space_type), target, intent(in)  :: self  
  real(kind=dp),              pointer             :: diff_basis(:,:,:,:)

  diff_basis => self%diff_basis

  return
end function get_diff_basis

! ----------------------------------------------------------------
! Get the nodal coordinates of the function_space
! ----------------------------------------------------------------
function get_nodes(self) result(nodal_coords)
  implicit none
  class(function_space_type), target, intent(in)  :: self
  real(kind=dp),              pointer             :: nodal_coords(:,:)
  
  nodal_coords => self%nodal_coords
  
  return
end function get_nodes

function which(self) result(fs)
  implicit none
  class(function_space_type),  intent(in) :: self
  integer :: fs
  
  fs = self%fs
  return
end function which

end module function_space_mod
