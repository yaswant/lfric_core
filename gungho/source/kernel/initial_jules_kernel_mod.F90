!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2019.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
! @brief Initialise Jules fields

module initial_jules_kernel_mod

  use argument_mod, only: arg_type, &
    GH_FIELD, GH_INTEGER, GH_WRITE, GH_READ, &
    ANY_SPACE_1, ANY_SPACE_2, ANY_SPACE_3, ANY_SPACE_4, &
    CELLS
  use constants_mod, only: r_def, i_def
  use kernel_mod, only: kernel_type

  implicit none

  !> Kernel metadata for Psyclone 
  type, public, extends(kernel_type) :: initial_jules_kernel_type
      private
      type(arg_type) :: meta_args(19) = (/              &
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_fraction
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_2), & ! leaf_area_index
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_2), & ! canopy_height
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! sd_orog
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! soil_albedo
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! soil_roughness
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! albedo_obs_sw
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! albedo_obs_vis
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! albedo_obs_nir
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_temperature
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_snow_mass
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_1), & ! tile_snow_rgrain
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! snow_soot
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_3), & ! chloro_sea
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! sea_ice_thickness
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! sea_ice_pond_frac
          arg_type(GH_FIELD,   GH_WRITE,  ANY_SPACE_4), & ! sea_ice_pond_depth
          arg_type(GH_INTEGER, GH_READ               ), & ! n_surf_tile
          arg_type(GH_INTEGER, GH_READ               )  & ! n_sea_ice_tile
          /)
      integer :: iterates_over = CELLS

  contains
      procedure, nopass :: initial_jules_code
  end type

  ! overload the default structure constructor for function space
  interface initial_jules_kernel_type
      module procedure initial_jules_kernel_constructor
  end interface

  public initial_jules_code
contains

  function initial_jules_kernel_constructor() result(self)
    implicit none
    type(initial_jules_kernel_type) :: self
    return
  end function initial_jules_kernel_constructor

  ! @param[in]  nlayers            The number of layers
  ! @param[out] tile_fraction      Surface tile fractions
  ! @param[out] leaf_area_index    Leaf Area Index
  ! @param[out] canopy_height      Canopy height
  ! @param[out] sd_orog            Standard deviation of orography
  ! @param[out] soil_albedo        Snow-free soil albedo
  ! @param[out] soil_roughness     Bare soil surface roughness length
  ! @param[out] albedo_obs_sw      Observed snow-free shortwave albedo
  ! @param[out] albedo_obs_vis     Observed snow-free visible albedo
  ! @param[out] albedo_obs_nir     Observed snow-free near-IR albedo
  ! @param[out] tile_temperature   Surface tile temperatures
  ! @param[out] tile_snow_mass     Snow mass on tiles (kg/m2)
  ! @param[out] tile_snow_rgrain   Snow grain size on tiles (microns)
  ! @param[out] snow_soot          Snow soot content (kg/kg)
  ! @param[out] chloro_sea         Chlorophyll content of the sea
  ! @param[out] sea_ice_thickness  Sea ice thickness (m)
  ! @param[out] sea_ice_pond_frac  Meltpond fraction on sea ice
  ! @param[out] sea_ice_pond_depth Meltpond depth on sea ice (m)
  ! @param[in]  n_surf_tile        Number of surface tiles
  ! @param[in]  n_sea_ice_tile     Number of sea ice tiles
  ! @param[in]  ndf_tile           Number of DOFs per cell for tiles
  ! @param[in]  undf_tile          Number of total DOFs for tiles
  ! @param[in]  map_tile           Dofmap for cell at the base of the column
  ! @param[in]  ndf_pft            Number of DOFs per cell for PFTs
  ! @param[in]  undf_pft           Number of total DOFs for PFTs
  ! @param[in]  map_pft            Dofmap for cell at the base of the column
  ! @param[in]  ndf_2d             Number of DOFs per cell for 2D fields
  ! @param[in]  undf_2d            Number of total DOFs for 2D fields
  ! @param[in]  map_2d             Dofmap for cell at the base of the column
  ! @param[in]  ndf_sice           Number of DOFs per cell for sea ice tiles
  ! @param[in]  undf_sice          Number of total DOFs for sea ice tiles
  ! @param[in]  map_sice           Dofmap for cell at the base of the column
  subroutine initial_jules_code(nlayers,                       &
                                tile_fraction,                 &
                                leaf_area_index,               &
                                canopy_height,                 &
                                sd_orog,                       &
                                soil_albedo,                   &
                                soil_roughness,                &
                                albedo_obs_sw,                 &
                                albedo_obs_vis,                &
                                albedo_obs_nir,                &
                                tile_temperature,              &
                                tile_snow_mass,                &
                                tile_snow_rgrain,              &
                                snow_soot,                     &
                                chloro_sea,                    &
                                sea_ice_thickness,             &
                                sea_ice_pond_frac,             &
                                sea_ice_pond_depth,            &
                                n_surf_tile,                   &
                                n_sea_ice_tile,                &
                                ndf_tile, undf_tile, map_tile, &
                                ndf_pft, undf_pft, map_pft,    &
                                ndf_2d, undf_2d, map_2d,       &
                                ndf_sice, undf_sice, map_sice)

      implicit none

      !Arguments
      integer(kind=i_def), intent(in) :: nlayers, n_surf_tile, n_sea_ice_tile
      integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
      integer(kind=i_def), intent(in) :: map_tile(ndf_tile)
      integer(kind=i_def), intent(in) :: ndf_pft, undf_pft
      integer(kind=i_def), intent(in) :: map_pft(ndf_pft)
      integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
      integer(kind=i_def), intent(in) :: map_2d(ndf_2d)
      integer(kind=i_def), intent(in) :: ndf_sice, undf_sice
      integer(kind=i_def), intent(in) :: map_sice(ndf_sice)

      real(kind=r_def), intent(out) :: tile_fraction(undf_tile)
      real(kind=r_def), intent(out) :: tile_temperature(undf_tile)
      real(kind=r_def), intent(out) :: tile_snow_mass(undf_tile)
      real(kind=r_def), intent(out) :: tile_snow_rgrain(undf_tile)

      real(kind=r_def), intent(out) :: leaf_area_index(undf_pft)
      real(kind=r_def), intent(out) :: canopy_height(undf_pft)

      real(kind=r_def), intent(out) :: sd_orog(undf_2d)
      real(kind=r_def), intent(out) :: soil_albedo(undf_2d)
      real(kind=r_def), intent(out) :: soil_roughness(undf_2d)
      real(kind=r_def), intent(out) :: albedo_obs_sw(undf_2d)
      real(kind=r_def), intent(out) :: albedo_obs_vis(undf_2d)
      real(kind=r_def), intent(out) :: albedo_obs_nir(undf_2d)
      real(kind=r_def), intent(out) :: snow_soot(undf_2d)
      real(kind=r_def), intent(out) :: chloro_sea(undf_2d)

      real(kind=r_def), intent(out) :: sea_ice_thickness(undf_sice)
      real(kind=r_def), intent(out) :: sea_ice_pond_frac(undf_sice)
      real(kind=r_def), intent(out) :: sea_ice_pond_depth(undf_sice)

      !Internal variables
      integer(kind=i_def) :: i


      ! Tile fraction (100% bare soil)
      do i = 1, n_surf_tile
        if (i == 8) then
          tile_fraction(map_tile(i)) = 1.0_r_def
        else
          tile_fraction(map_tile(i)) = 0.0_r_def
        end if
      end do

      ! Leaf area index (UKV values: lai_pft)
      leaf_area_index(map_pft(1)) = 5.0_r_def
      leaf_area_index(map_pft(2)) = 4.0_r_def
      leaf_area_index(map_pft(3)) = 2.0_r_def
      leaf_area_index(map_pft(4)) = 4.0_r_def
      leaf_area_index(map_pft(5)) = 1.0_r_def

      ! Canopy height (UKV values: canht_pft)
      canopy_height(map_pft(1)) = 19.01_r_def
      canopy_height(map_pft(2)) = 16.38_r_def
      canopy_height(map_pft(3)) =  1.46_r_def
      canopy_height(map_pft(4)) =  1.26_r_def
      canopy_height(map_pft(5)) =  1.59_r_def

      ! Standard deviation of orography
      sd_orog(map_2d(1)) = 0.0_r_def

      ! Snow-free soil albedo (Cardington UKV value: albsoil_soilt)
      soil_albedo(map_2d(1)) = 0.1065_r_def

      ! Bare soil surface roughness length
      soil_roughness(map_2d(1)) = 0.0_r_def
      
      ! Observed snow-free shortwave albedo
      albedo_obs_sw(map_2d(1)) = 0.0_r_def

      ! Observed snow-free visible albedo
      albedo_obs_vis(map_2d(1)) = 0.0_r_def

      ! Observed snow-free near infra-red albedo
      albedo_obs_nir(map_2d(1)) = 0.0_r_def


      ! Tile temperature
      do i = 1, n_surf_tile
        tile_temperature(map_tile(i)) = 295.0_r_def
      end do

      ! Snow mass on tiles (kg/m2)
      do i = 1, n_surf_tile
        tile_snow_mass(map_tile(i)) = 0.0_r_def
      end do

      ! Snow grain size on tiles (microns)
      do i = 1, n_surf_tile
        tile_snow_rgrain(map_tile(i)) = 300.0_r_def
      end do

      ! Snow soot content (kg/kg)
      snow_soot(map_2d(1)) = 0.0_r_def

      ! Chlorophyll content of the sea near-surface (kg/m3)
      chloro_sea(map_2d(1)) = 5.0e-7_r_def

      ! Variables on sea ice tiles
      do i = 1, n_sea_ice_tile
        ! Sea ice thickness (m)
        sea_ice_thickness(map_sice(i)) = 0.0_r_def
        ! Meltpond fraction
        sea_ice_pond_frac(map_sice(i)) = 0.0_r_def
        ! Meltpond depth (m)
        sea_ice_pond_depth(map_sice(i)) = 0.0_r_def
      end do

  end subroutine initial_jules_code

end module initial_jules_kernel_mod
