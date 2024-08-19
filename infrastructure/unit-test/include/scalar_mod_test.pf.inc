! Preprocessor is used to template test code.

#define TMPL_PASTE(a) a
#define TMPL_MODULE_NAME TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_mod_test
#define TMPL_TEST_CLASS TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_test_type
#define TMPL_SCALAR_MOD TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_mod
#define TMPL_SCALAR_TYPE TMPL_PASTE(TMPL_PASTE(scalar_)TMPL_KIND)_type
#define TMPL_LITERAL(v) TMPL_PASTE(TMPL_PASTE(v)_)TMPL_KIND

module TMPL_MODULE_NAME

  use, intrinsic :: iso_fortran_env, only : TMPL_KIND
  use constants_mod, only: i_def
  use mpi_mod, only : global_mpi
  use pFUnit_Mod

  implicit none

  private
  public :: test_all

  @TestCase
  type, extends(MPITestCase), public :: TMPL_TEST_CLASS
    private
  contains
    procedure setUp
    procedure tearDown
    procedure test_all
  end type TMPL_TEST_CLASS

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setUp( this )

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: this

    !Store the MPI communicator for later use
    call global_mpi%initialise(this%getMpiCommunicator())

  end subroutine setUp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tearDown( this )

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: this

    ! Clear the stored MPI communicator
    call global_mpi%finalise()

  end subroutine tearDown

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  @Test( npes=[2] )
  subroutine test_all( this )

    use TMPL_SCALAR_MOD, only: TMPL_SCALAR_TYPE

    implicit none

    class(TMPL_TEST_CLASS), intent(inout) :: this

    type(TMPL_SCALAR_TYPE) :: s1, s2

    TMPL_TYPE(TMPL_KIND) :: s1_val, s2_val
    TMPL_TYPE(TMPL_KIND) :: sum_val, min_val, max_val

    integer(i_def) :: local_rank

    local_rank  = global_mpi%get_comm_rank()

    s1_val = TMPL_LITERAL(10.0) + TMPL_CAST(local_rank, kind=TMPL_KIND)

    s1%value = s1_val

    ! Check the scalar type has been correctly constructed and that
    ! we have access to the data
    @assertEqual( s1_val, s1%value, TMPL_LITERAL(1.0e-2) )

    ! Test the global reductions
    ! - a bit trivial until we have parallel tests
    sum_val = s1%get_sum()
    @assertEqual( TMPL_LITERAL(21.0), sum_val, TMPL_LITERAL(1.0e-2) )

    min_val = s1%get_min()
    @assertEqual( TMPL_LITERAL(10.0), min_val, TMPL_LITERAL(1.0e-2) )

    max_val = s1%get_max()
    @assertEqual( TMPL_LITERAL(11.0), max_val, TMPL_LITERAL(1.0e-2) )

    ! Test alternative way to create a scalar - by calling the initialiser
    s2_val = TMPL_LITERAL(20.0) + TMPL_CAST(local_rank, kind=TMPL_KIND)
    call s2%initialise(s2_val, global_mpi)
    ! Make sure we can call one of the reductions on this scalar
    sum_val = s2%get_sum()
    @assertEqual( TMPL_LITERAL(41.0), sum_val, TMPL_LITERAL(1.0e-2) )

  end subroutine test_all

end module TMPL_MODULE_NAME
