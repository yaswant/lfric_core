! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_setup_io_mod

USE clock_mod,                     ONLY: clock_type
USE constants_mod,                 ONLY: i_def, str_max_filename
USE lfric_xios_file_mod,           ONLY: xios_file_type
USE linked_list_mod,               ONLY: linked_list_type,                     &
                                         linked_list_item_type

! Configuration modules
USE files_config_mod,              ONLY: ancil_directory,                      &
                                         checkpoint_stem_name,                 &
                                         land_area_ancil_path,                 &
                                         start_dump_filename,                  &
                                         start_dump_directory
USE initialization_config_mod,     ONLY: init_option,                          &
                                         init_option_fd_start_dump,            &
                                         ancil_option,                         &
                                         ancil_option_aquaplanet
USE io_config_mod,                 ONLY: diagnostic_frequency,                 &
                                         checkpoint_write,                     &
                                         checkpoint_read,                      &
                                         write_dump

IMPLICIT NONE

PRIVATE
PUBLIC :: init_lfricinp_files

CONTAINS

SUBROUTINE init_lfricinp_files(files_list, clock)

! lfricinp modules
USE lfricinp_ancils_mod, ONLY: l_land_area_fraction

IMPLICIT NONE

TYPE(linked_list_type), INTENT(OUT) :: files_list
TYPE(clock_type),       INTENT(IN)  :: clock

TYPE(xios_file_type)            :: tmp_file
INTEGER(i_def)                  :: tmp_freq
CHARACTER(LEN=str_max_filename) :: checkpoint_write_fname,                     &
                                   checkpoint_read_fname,                      &
                                   dump_fname,                                 &
                                   ancil_fname

! Setup diagnostic output file
CALL tmp_file%init_xios_file("lfric_diag", freq=diagnostic_frequency)
CALL files_list%insert_item(tmp_file)

! Setup dump-writing context information
IF ( write_dump ) THEN
  ! Create dump filename from base name and end timestep
  WRITE(dump_fname,'(A,A,I6.6)')                                               &
     TRIM(start_dump_directory)//'/'//TRIM(start_dump_filename),"_",           &
     clock%get_last_step()

  ! Setup dump file for end timestep
  CALL tmp_file%init_xios_file( "lfric_fd_dump", dump_fname,                   &
                                clock%get_last_step())
  CALL files_list%insert_item(tmp_file)
END IF

! Setup dump-reading context information
IF ( init_option == init_option_fd_start_dump .OR.                             &
    ancil_option == ancil_option_aquaplanet ) THEN
  ! Create dump filename from stem
  WRITE(dump_fname,'(A)') TRIM(start_dump_directory)//'/'//                    &
                          TRIM(start_dump_filename)

  ! Setup dump file
  CALL tmp_file%init_xios_file( "read_lfric_fd_dump", path=dump_fname)
  CALL files_list%insert_item(tmp_file)
END IF

IF (l_land_area_fraction) THEN
  ! Set land area ancil filename from namelist
  WRITE(ancil_fname,'(A)') TRIM(ancil_directory)//'/'//                        &
                           TRIM(land_area_ancil_path)
  CALL tmp_file%init_xios_file("land_area_ancil", path=ancil_fname)
  CALL files_list%insert_item(tmp_file)
END IF

! Setup checkpoint writing context information
IF ( checkpoint_write ) THEN
  ! Create checkpoint filename from stem and end timestep
  WRITE(checkpoint_write_fname,'(A,A,I6.6)')                                   &
                       TRIM(checkpoint_stem_name),"_", clock%get_last_step()

  CALL tmp_file%init_xios_file( "lfric_checkpoint_write", &
                                checkpoint_write_fname, clock%get_last_step(), &
                                field_group_id="checkpoint_fields" )
  CALL files_list%insert_item(tmp_file)
END IF

! Setup checkpoint reading context information
IF ( checkpoint_read ) THEN
  ! Create checkpoint filename from stem and (start - 1) timestep
  WRITE(checkpoint_read_fname,'(A,A,I6.6)')                                    &
               TRIM(checkpoint_stem_name),"_", (clock%get_first_step() - 1)

  CALL tmp_file%init_xios_file( "lfric_checkpoint_read", checkpoint_read_fname,&
                                clock%get_first_step() - 1,                    &
                                field_group_id="checkpoint_fields" )
  CALL files_list%insert_item(tmp_file)
END IF

END SUBROUTINE init_lfricinp_files

END MODULE lfricinp_setup_io_mod
