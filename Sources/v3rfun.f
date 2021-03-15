!-------------------------------------------------------------------------------
!  The @header2, @begin_table, @item3 and @end_table commands are custom defined
!  commands in Doxygen.in. They are defined under ALIASES. For the page created
!  here, the 80 column limit is exceeded. Arguments of aliases are separated by
!  ','. If you intended ',' to be a string you must use an escaped comma '\,'.
!
!>  @page v3rfun_cl_parsing_sec Command Line Arguments
!>
!>  @tableofcontents
!>
!>  @section v3rfun_cl_parsing_intro Introduction
!>  This contains a description of the command line arguments. All arguments
!>  take the form of
!>
!>  @fixed_width{xv3rfun [-arg] filename}
!>
!>  @section v3rfun_cl_parsing_arg_sec Command Line Arguments
!>  @header2{Argument, Takes Value, Discription}
!>  @begin_table
!>     @item3{@fixed_width{-h},     N, Displays the help text and exits the program.}
!>  @end_table
!>
!>  @section v3rfun_cl_pasring_prog_ref_sec Programmers Reference
!>  Reference material for the coding to implement command line parsing is found
!>  in the @ref v3rfun.f
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file v3rfun.f
!>  @brief Contains the main routines for V3RFUN.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  The V3RFUN code is part of a suite of codes, called V3FIT, designed for
!>  equilibrium reconstruction of stellarator plasmas. It takes as input
!>  descriptions of the field coils and diagnostic coils, and generates files
!>  containing the coil and plasma response functions and the plasma response
!>  functions. The output files from the V3RFUN code are designed to be used as
!>  input to the V3FIT code. Physical quantities in the V3FIT project are
!>  measured in SI units, unless otherwise specified.
!>
!>  @author James D. Hanson
!>  @author Mark Cianciosa
!>
!>  Below is a brief discription of the major top level objects of the code. For
!>  discriptions of lower level objects consult the referenced top level
!>  objects.
!*******************************************************************************
!*******************************************************************************
!  MAIN PROGRAM
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief V3RFUN main program.
!>
!>  Highest level V3RFUN routine.
!-------------------------------------------------------------------------------
      PROGRAM v3rfun
      USE v3rfun_context
      USE profiler
      USE v3_utilities
!$    USE omp_lib

      IMPLICIT NONE

!  local variables
!  Forward declare the v3rfun_parse_command_line function.
      INTERFACE
         FUNCTION v3rfun_parse_command_line()
         USE v3rfun_input
         IMPLICIT NONE
         CHARACTER (len=v3rfun_file_length) :: v3rfun_parse_command_line
         END FUNCTION
      END INTERFACE

      TYPE (v3rfun_context_class), POINTER :: context => null()
      TYPE (diagnostic_dot_coil), POINTER  :: head
      INTEGER                              :: id_num
#if defined(MPI_OPT)
      INTEGER                              :: error

!  Start of executable code
!  The Intel version of the MPI libraries does not provide a correct value for
!  time until after MPI_INIT is called. Make sure this is the first think called
!  so that correct timing information can be used.
      CALL MPI_INIT(error)
#endif

      CALL profiler_construct

!$    CALL OMP_SET_NUM_THREADS(OMP_GET_NUM_PROCS())
!$OMP PARALLEL
!$    error = OMP_GET_MAX_THREADS()
!$OMP END PARALLEL
!$    WRITE (*,*) error

      context =>                                                               &
     &   v3rfun_context_construct(TRIM(v3rfun_parse_command_line()))

!  Compute the coil response functions.
      head => context%coils
      id_num = 1
      DO WHILE(ASSOCIATED(head))
         SELECT CASE (head%d_type)

            CASE ('flux_loop', 'flux_loop_circular', 'magnetic_probe',         &
     &            'magnetic_probe_tokamak', 'b_rogowski', 'i_rogowski',        &
     &            'f_rogowski','s_rogowski')
               CALL v3rfun_context_write_mrf(context, head, id_num)

            CASE ('b_point_probe')
               CALL v3rfun_context_write_point(context, head, id_num)

            CASE DEFAULT
               CALL err_fatal('Unknown magnetic diagnostic type.')

         END SELECT

         head => head%next
         id_num = id_num + 1
      END DO

      CALL v3rfun_context_destruct(context)

      CALL profiler_destruct

#if defined(MPI_OPT)
      CALL MPI_FINALIZE(error)
#endif

      END PROGRAM v3rfun

!-------------------------------------------------------------------------------
!>  @brief Parse command line arguments.
!>
!>  The command line arguments are simple. Expect either the -h flag or read the
!>  namelist input filename. If any other argument is encountered display the
!>  help and terminate the program.
!>
!>  @returns The namelist input filename.
!-------------------------------------------------------------------------------
      FUNCTION v3rfun_parse_command_line()
      USE v3rfun_input

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=v3rfun_file_length) :: v3rfun_parse_command_line

!  local variables
      INTEGER                            :: num_args, i
      REAL (rprec)                       :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Read the zeroith arg to get the number of arguments. This should also be the
!  command name.
      CALL getcarg(0, v3rfun_parse_command_line, num_args)

      IF (num_args .le. 0 .or. num_args .gt. 2) THEN
         CALL v3rfun_print_help
      END IF

      DO i = 1, num_args
         CALL getcarg(i, v3rfun_parse_command_line, num_args)
         IF (TRIM(v3rfun_parse_command_line) .eq. '-h') THEN
            CALL v3rfun_print_help
         END IF
      END DO

      CALL profiler_set_stop_time('v3rfun_parse_command_line',                 &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Print out help text.
!>
!>  Command line help text should take the form of '-flag'. A Y or N indicating
!>  the flag takes a value or not followed by a space. A short message
!>  describing the flag. A reference for all the command line arguments can be
!>  found at @ref FIXME.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_print_help()
      USE v3rfun_context
      USE magnetic_response

      IMPLICIT NONE

!  Start of executable code
!  All command line messages need to fit within this width.
!                '           s              c                          '
      WRITE(*,*) '                                                     '
      WRITE(*,1000) code_name
      WRITE(*,1100) magnetic_response_current
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Usage: xv3rfun [-arg] filename ...                   '
      WRITE(*,*) '                                                     '
      WRITE(*,*) 'Options:                                             '
      WRITE(*,*) 'All options are displayed as [arg][takesoption][text]'
      WRITE(*,*) '  -h       N Display this information                '
      WRITE(*,*) '                                                     '

      STOP

1000  FORMAT('Code: ',a)
1100  FORMAT('Version: ',a)

      END SUBROUTINE
