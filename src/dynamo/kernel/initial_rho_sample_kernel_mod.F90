!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes initial rho field

module initial_rho_sample_kernel_mod

    use argument_mod,               only : arg_type, func_type,      &
                                           GH_FIELD, GH_REAL,        &
                                           GH_READ, GH_WRITE,        &
                                           ANY_SPACE_9, W3, CELLS
    use constants_mod,              only : r_def
    use idealised_config_mod,       only : test
    use kernel_mod,                 only : kernel_type

    implicit none

    !-------------------------------------------------------------------------------
    ! Public types
    !-------------------------------------------------------------------------------
    !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
    type, public, extends(kernel_type) :: initial_rho_sample_kernel_type
        private
        type(arg_type) :: meta_args(3) = (/                                 &
            arg_type(GH_FIELD,   GH_WRITE,  W3),                            &
            arg_type(GH_FIELD*3, GH_READ, ANY_SPACE_9),                     &
            arg_type(GH_REAL,    GH_READ)                                   &
            /)
        integer :: iterates_over = CELLS
    contains
        procedure, nopass ::initial_rho_sample_code
    end type

    !-------------------------------------------------------------------------------
    ! Constructors
    !-------------------------------------------------------------------------------

    ! overload the default structure constructor for function space
    interface initial_rho_sample_kernel_type
        module procedure initial_rho_sample_kernel_constructor
    end interface

    !-------------------------------------------------------------------------------
    ! Contained functions/subroutines
    !-------------------------------------------------------------------------------
    public initial_rho_sample_code
contains

    type(initial_rho_sample_kernel_type) function initial_rho_sample_kernel_constructor() result(self)
        return
    end function initial_rho_sample_kernel_constructor

    !> @brief Computes the initial rho field
    !! @param[in] nlayers Number of layers
    !! @param[in] ndf_w3 Number of degrees of freedom per cell
    !! @param[in] undf_w3 Total number of degrees of freedom
    !! @param[in] map_w3 Dofmap for the cell at the base of the column
    !! @param[inout] rho Density
    !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
    !! @param[in] undf_chi Number of degrees of freedom for chi
    !! @param[in] map_chi Dofmap for the cell at the base of the column for chi
    !! @param[in] chi_basis Basis functions evaluated at gaussian quadrature points
    !! @param[inout] chi_1 X component of the chi coordinate field
    !! @param[inout] chi_2 Y component of the chi coordinate field
    !! @param[inout] chi_3 Z component of the chi coordinate field
    !! @param[in] time Current time of the model run
    subroutine initial_rho_sample_code(nlayers, ndf_w3, undf_w3, map_w3, rho, &
        ndf_chi, undf_chi, map_chi, chi_basis, chi_1, chi_2, chi_3, time)

        use analytic_density_profiles_mod, only : analytic_density

        implicit none

        !Arguments
        integer,                                        intent(in)  :: nlayers, ndf_w3, ndf_chi, undf_w3, undf_chi
        integer, dimension(ndf_w3),                     intent(in)  :: map_w3
        integer, dimension(ndf_chi),                    intent(in)  :: map_chi
        real(kind=r_def), dimension(undf_w3),           intent(out) :: rho
        real(kind=r_def), dimension(undf_chi),          intent(in)  :: chi_1, chi_2, chi_3
        real(kind=r_def), dimension(1,ndf_chi, ndf_w3), intent(in)  :: chi_basis
        real(kind=r_def),                               intent(in)  :: time
        !Internal variables
        integer               :: df1, df, k

        real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
        real(kind=r_def)                             :: x(3)

        ! compute the RHS & LHS integrated over one cell and solve
        do k = 0, nlayers-1
            do df1 = 1, ndf_chi
                chi_1_e(df1) = chi_1( map_chi(df1) + k)
                chi_2_e(df1) = chi_2( map_chi(df1) + k)
                chi_3_e(df1) = chi_3( map_chi(df1) + k)
            end do


            do df = 1, ndf_w3
              x(:) = 0.0_r_def
              do df1 = 1, ndf_chi
                x(1) = x(1) + chi_1_e(df1)*chi_basis(1,df1,df)
                x(2) = x(2) + chi_2_e(df1)*chi_basis(1,df1,df)
                x(3) = x(3) + chi_3_e(df1)*chi_basis(1,df1,df)
              end do

              rho(map_w3(df) + k) = analytic_density(x, test, time)

            end do
        end do


    end subroutine initial_rho_sample_code

end module initial_rho_sample_kernel_mod
