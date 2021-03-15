!-------------------------------------------------------------------------------
!  The @fixed_width, @begin_table, @item2 and @end_table commands are custom
!  defined commands in Doxygen.in. They are defined under ALIASES. For the page
!  created here, the 80 column limit is exceeded. Arguments of aliases are
!  separated by ','. If you intended ',' to be a string you must use an escaped
!  comma '\,'.
!
!>  @page diagnostic_dot_sec Magnetic Diagnostic Dot File
!>
!>  @tableofcontents
!>  @section diagnostic_dot_intro_sec Introduction
!>  This document contains information about the magnetic diagnostic
!>  specifications. Diagnostics are specified in a structured text file. This
!>  file is read in by the @ref v3rfun_input::name_diagnostic_dot variable of
!>  the name list input file.
!>
!>  @section  diagnostic_dot_coil_spec_sec Diagnostic Coil Specification
!>  Every diagnostic coil specification looks like:\n
!>  @fixed_width{key_word}\n
!>  @fixed_width{name}\n
!>  The remaining specification is keyword dependent.
!>
!>  @section diagnostic_dot_key_sec Keywords
!>  @begin_table
!>  @item2{@fixed_width{flux_loop},              An arbitrary shaped flux loop.\n\n
!>                                               Lines after ID:\n
!>                                               Arbitrary number of nodes\, each consists of
!>                                               three reals.\n\n
!>                                               Description: A @fixed_width{flux_loop} is a
!>                                               closed filamentary polygon. Each line after
!>                                               the ID line contains the x\, y\, and z
!>                                               coordinates of one node. The last node need
!>                                               not be the same as the first node: the parser
!>                                               will ensure that the flux loop is closed.}
!>  @item2{@fixed_width{flux_loop_circular},     A circular flux loop.\n\n
!>                                               Lines after ID:\n
!>                                               First line after ID: one real\, the radius of
!>                                               the circle.\n Second line after ID: three
!>                                               reals\, the x\, y\, z coordinates of the
!>                                               center of the circle.\n Third line after ID:
!>                                               three reals\, the x\, y\, z components of the
!>                                               normal to the circle.\n\n
!>                                               Description: A flux_loop_circular is a
!>                                               filamentary circle. The parser takes care of
!>                                               the normalization of the normal vector.}
!>  @item2{@fixed_width{magnetic_probe},         Point measurement of the magnetic field.\n\n
!>                                               Lines after ID:\n
!>                                               First line after ID: one real\, the radius
!>                                               of the magnetic probe.\n Second line after ID:
!>                                               three reals\, the x\, y\, z coordinates of
!>                                               probe position.\n Third line after ID: three
!>                                               reals\, the x\, y\, z components of the
!>                                               normal to the probe.\n\n
!>                                               Description: A magnetic_probe is very similar
!>                                               to a small flux_loop_circular diagnostic.
!>                                               However\, the signal for a magnetic_probe is
!>                                               a B field component (units in Tesla) rather
!>                                               than a magnetic flux (units Webers). The
!>                                               parser takes care of the normalization of the
!>                                               normal vector.}
!>  @item2{@fixed_width{magnetic_probe_tokamak}, Point measurement of the magnetic field.\n\n
!>                                               Lines after ID:\n
!>                                               First line after ID: five reals: probe
!>                                               radius\, major radius\, toroidal angle (in
!>                                               degrees)\, z position\, and orientation angle
!>                                               (in degrees).\n\n
!>                                               Description: A magnetic_probe_tokamak is very
!>                                               similar to a magnetic_probe. However\, the
!>                                               parameterization is different. The normal to
!>                                               a magnetic_probe_tokamak lies in the
!>                                               Rhat-zhat plane. It has no component in the
!>                                               phi-hat direction. The orientation angle
!>                                               specifies the angle from the Rhat vector.}
!>  @item2{@fixed_width{rogowski},               Rogowski coil.\n\n
!>                                               Lines after ID:\n
!>                                               First line after ID: two reals\, the number
!>                                               of turns\, and the cross-sectional area (in
!>                                               meter^2) of the Rogowski coil. Second through
!>                                               last lines after ID: arbitrary in number each
!>                                               containing three reals.\n\n
!>                                               Description: The guide curve for a rogowski
!>                                               is specified by giving node positions. Each
!>                                               of the second through last line after the ID
!>                                               line specifies a node position. The guide
!>                                               curve is made of straight line segments
!>                                               between the node positions. Note that the
!>                                               parser makes no assumption about whether or
!>                                               not a Rogowski coil guide curve is closed or
!>                                               not. Thus\, a rogowski can model either a
!>                                               partial Rogowski coil\, or a full Rogowski
!>                                               coil. The signal computed for a Rogowski is
!>                                               the line integrated B along the guide curve\,
!>                                               divided by the length of the guide curve. The
!>                                               units are Tesla.}
!>  @item2{@fixed_width{b_rogowski},             Same as a @fixed_width{rogowski}.}
!>  @item2{@fixed_width{i_rogowski},             Same as a @fixed_width{rogowski} but with
!>                                               measuing current in Amperes. If this is a
!>                                               full rogowski\, the signal is the total
!>                                               current enclosed.}
!>  @item2{@fixed_width{f_rogowski},             Same as a @fixed_width{rogowski} but with
!>                                               measuing flux in Webers.}
!>  @item2{@fixed_width{s_rogowski},             A rogowski defined by a series of segments
!>                                               connected in series to represent a single
!>                                               measurement.}
!>  @item2{@fixed_width{b_point_probe},          Point measurement of the magnetic field in
!>                                               tesla.\n\n
!>                                               Lines after ID:\n
!>                                               First line after ID: three reals\, the x\,
!>                                               y\, z coordinates of probe position.\n
!>                                               Second line after ID: three reals\, the xhat\,
!>                                               yhat\, zhat components of the direction the probe
!>                                               measures in.\n\n
!>                                               Description: These probes differ from the
!>                                               magnetic probes by modling the probe as an
!>                                               exact point measurement insead of modeling as
!>                                               a loop of wire. The second vector direction
!>                                               defines the components that the probe reads.}
!>  @item2{@fixed_width{b_point_probe_cyl},      Same as @fixed_width{b_point_probe} but defined
!>                                               in cylindrical coordinates.}
!>  @item2{@fixed_width{end_of_file},            End of the file.}
!>  @end_table
!>
!>  @section diagnostic_dot_prog_ref_sec Programmers Reference
!>  Reference material for the coding to parse these files can be found in the
!>  @ref diagnostic_dot module.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file diagnostic_dot.f
!>  @brief Contains module @ref diagnostic_dot
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Module for opening and reading a diagnostic dot file. The file format for
!>  these files are documented in @ref diagnostic_dot_sec.
!*******************************************************************************
      MODULE diagnostic_dot
      USE stel_kinds
      USE stel_constants
      USE bsc_T
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  diagnostic_dot module parameters
!*******************************************************************************
!>  Maximum line length.
      INTEGER, PARAMETER :: diagnostic_dot_line_len = 200
!>  Maximum id name.
      INTEGER, PARAMETER :: diagnostic_dot_id_name = 30

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) diagnostic_dot_coil
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  A single coil. A coil set is structured as a singly linked list.
!-------------------------------------------------------------------------------
      TYPE diagnostic_dot_coil
!>  Flux factor.
         REAL (rprec)                           :: factor
!>  Measurement units.
         CHARACTER (len=2)                      :: units
!>  Diagnostic Type
         CHARACTER (len=diagnostic_dot_id_name) :: d_type
!>  Diagnostic id name
         CHARACTER (len=diagnostic_dot_id_name) :: id_name

!>  Probe position.
         REAL (rprec), DIMENSION(3)             :: position
!>  Probe direction.
         REAL (rprec), DIMENSION(3)             :: direction

!>  Biotsavart Coil object.
         TYPE (bsc_coil), POINTER               :: coil => null()
!>  Reference to the next vertex.
         TYPE (diagnostic_dot_coil), POINTER    :: next => null()
      END TYPE

      CONTAINS

!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct the coil diagnostics.
!>
!>  Allocates memory and initializes all the diagnostic coils by parseing the
!>  diagnostic dot file.
!>
!>  @param[in] filename File name of the diagnostic dot file.
!>  @returns A pointer to the constructed coil.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_construct(filename)
      USE safe_open_mod

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_dot_coil), POINTER :: diagnostic_dot_construct
      CHARACTER (len=*), INTENT(in)           :: filename

!  local variables
      TYPE (diagnostic_dot_coil), POINTER     :: current_coil
      TYPE (diagnostic_dot_coil), POINTER     :: new_coil
      INTEGER                                 :: status
      INTEGER                                 :: iou
      CHARACTER (len=diagnostic_dot_line_len) :: line
      INTEGER                                 :: line_number
      REAL (rprec)                            :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      diagnostic_dot_construct => null()
      new_coil => null()
      current_coil => null()

!  Open the diagnostic dot file.
      CALL safe_open(iou, status, filename, 'old', 'formatted')

!  Read lines until a keywood is reached. Read the first line out side the loop.
!  Once a keyword is reached, hand off to the appropriate parsing function. The
!  parsing functions will update the next line when keyword is reached.
      READ (iou, '(a)', iostat=status) line
      line_number = 1

      DO
         SELECT CASE (TRIM(line))

            CASE DEFAULT
               READ (iou, '(a)', iostat=status) line
               line_number = line_number + 1

            CASE ('end_of_file')
               EXIT

            CASE ('flux_loop')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_filaments(iou, line, line_number,       &
     &                                           .false., 'Wb')
               new_coil%d_type = 'flux_loop'

            CASE ('flux_loop_circular')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_probe(iou, line, line_number,           &
     &                                       .false., 'Wb')
               new_coil%d_type = 'flux_loop_circular'

            CASE ('magnetic_probe')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_probe(iou, line, line_number,           &
     &                                       .false., 'T')
               new_coil%d_type = 'magnetic_probe'

            CASE ('magnetic_probe_tokamak')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_probe(iou, line, line_number,           &
     &                                       .true., 'T')
               new_coil%d_type = 'magnetic_probe_tokamak'

            CASE ('rogowski', 'b_rogowski')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_filaments(iou, line, line_number,       &
     &                                           .true., 'T')
               new_coil%d_type = 'b_rogowski'

            CASE ('i_rogowski')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_filaments(iou, line, line_number,       &
     &                                           .true., 'A')
               new_coil%d_type = 'i_rogowski'

            CASE ('f_rogowski')
               new_coil =>                                                     &
     &            diagnostic_dot_parse_filaments(iou, line, line_number,       &
     &                                           .true., 'Wb')
               new_coil%d_type = 'f_rogowski'

            CASE ('s_rogowski')
               new_coil => diagnostic_dot_parse_s_rogowski(iou, line,          &
     &                                                     line_number)
               new_coil%d_type = 's_rogowski'

            CASE ('b_point_probe')
                new_coil =>                                                    &
     &             diagnostic_dot_parse_point(iou, line, line_number,          &
     &                                        'T', .true.)
                new_coil%d_type = 'b_point_probe'

            CASE ('b_point_probe_cyl')
                new_coil =>                                                    &
     &             diagnostic_dot_parse_point(iou, line, line_number,          &
     &                                        'T', .false.)
                new_coil%d_type = 'b_point_probe'

         END SELECT

!  Check that the current coil was associtaed before assiging the pointer. When
!  a diagnostic dot file line is read but doens't contain a keyword, a new
!  coil does not get allocated. If a new coil was created check to see if it is
!  the first coils. Otherwise assign it to the next coil and move the current
!  coil to the new coil.
         IF (ASSOCIATED(new_coil)) THEN
            IF (.not.ASSOCIATED(diagnostic_dot_construct)) THEN
               current_coil => new_coil
               diagnostic_dot_construct => current_coil
            ELSE
               current_coil%next => new_coil
               current_coil => current_coil%next
            END IF
         END IF

         new_coil => null()
      END DO

      CLOSE(iou)

      CALL profiler_set_stop_time('diagnostic_dot_construct',                  &
     &                            start_time)

      END FUNCTION

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref diagnostic_dot_coil list.
!>
!>  Deallocates memory and uninitializes a @ref diagnostic_dot_coil list.
!>
!>  @param[inout] this A @ref diagnostic_dot_coil list.
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_dot_destruct(list)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_dot_coil), POINTER :: list

!  local variables
      TYPE (diagnostic_dot_coil), POINTER :: current_coil

!  Start of executable code
      DO WHILE(ASSOCIATED(list))

!  Check if the end of the list has been reached yet. If this is not the last
!  node, point current_coil to the start of the linked list. Move the start of
!  the list to next node. current_coil still points to the last node. With the
!  list pointing to the next node, erase the current node. Otherwise erase the
!  list.
         IF (ASSOCIATED(list%next)) THEN
            current_coil => list
            list => list%next

            IF (ASSOCIATED(current_coil%coil)) THEN
               DEALLOCATE(current_coil%coil)
               current_coil%coil => null()
            END IF

            DEALLOCATE(current_coil)
         ELSE
            DEALLOCATE(list)
         END IF

      END DO

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Report error.
!>
!>  Errors are formatted the the subject, line, line_number and calling
!>  function. This subroutine does not return so there is no need for profiling
!>  information.
!>
!>  @param[in] subject     Description of what caused the error.
!>  @param[in] line        The current line of the file.
!>  @param[in] line_number The current line number of the file.
!>  @param[in] caller      Name of the calling function.
!-------------------------------------------------------------------------------
      SUBROUTINE diagnostic_dot_error(subject, line, line_number,              &
     &                                caller)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: subject
      CHARACTER (len=*), INTENT(in) :: line
      INTEGER, INTENT(in)           :: line_number
      CHARACTER (len=*), INTENT(in) :: caller

!  Start of executable code
      IF (LEN(subject) .eq. 4) THEN
         WRITE (*,1000) subject, line_number
      ELSE
         WRITE (*,1100) subject, line_number
      END IF
      WRITE (*,1200) line
      WRITE (*,1300) caller
      CALL EXIT(1)

1000  FORMAT('Failed to read ', a4,' on line: ',i6)
1100  FORMAT('Failed to read ', a6,' on line: ',i6)
1200  FORMAT(a)
1300  FORMAT('Called from: ',a)

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Parse a name.
!>
!>  Reads a name and from the diagnostic dot file.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[in]    caller      Name of the calling function.
!>  @returns A name.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_name(iou, line, line_number, caller)

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=diagnostic_dot_id_name) ::                                &
     &   diagnostic_dot_parse_name
      INTEGER                                             :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                              :: line_number
      CHARACTER (len=*), INTENT(in)                       :: caller

!  local variables
      INTEGER                                             :: status
      REAL (rprec)                                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      status = 0
      READ (iou,1000,iostat=status) line
      line_number = line_number + 1

      IF (status .ne. 0) THEN
         CALL diagnostic_dot_error('name', line, line_number, caller)
      END IF

      diagnostic_dot_parse_name = TRIM(line)

      CALL profiler_set_stop_time('diagnostic_dot_parse_name',                 &
     &                            start_time)

1000  FORMAT(a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse a real.
!>
!>  Reads a single real and from the diagnostic dot file.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[in]    caller      Name of the calling function.
!>  @returns A real.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_real(iou, line, line_number, caller)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec) :: diagnostic_dot_parse_real
      INTEGER                                             :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                              :: line_number
      CHARACTER (len=*), INTENT(in)                       :: caller

!  local variables
      INTEGER                                             :: status
      REAL (rprec)                                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      status = 0
      READ (iou,1000,iostat=status) line
      line_number = line_number + 1
      READ (line,*,iostat=status) diagnostic_dot_parse_real

      IF (status .ne. 0) THEN
         CALL diagnostic_dot_error('real', line, line_number, caller)
      END IF

      CALL profiler_set_stop_time('diagnostic_dot_parse_real',                 &
     &                            start_time)

1000  FORMAT(a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse a node.
!>
!>  Reads a single node and from the diagnostic dot file. Node may be arbitrary
!>  in length. Leave error handling up to the caller.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[out]   status      The result of the read attempt.
!>  @returns A real.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_node(iou, line, line_number, status)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3) :: diagnostic_dot_parse_node
      INTEGER                                             :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                              :: line_number
      INTEGER, INTENT(out)                                :: status

!  local variables
      REAL (rprec)                                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      status = 0
      READ (iou,1000,iostat=status) line
      line_number = line_number + 1
      READ (line,*,iostat=status) diagnostic_dot_parse_node

      CALL profiler_set_stop_time('diagnostic_dot_parse_node',                 &
     &                            start_time)

1000  FORMAT(a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse a node with sensitivity.
!>
!>  Reads a single node and from the diagnostic dot file. Node may be arbitrary
!>  in length. Leave error handling up to the caller. The node contains a
!>  sensisitivity value for the calibration as well.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[out]   status      The result of the read attempt.
!>  @returns A real.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_node4(iou, line, line_number,              &
     &                                    status)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(4) :: diagnostic_dot_parse_node4
      INTEGER                                             :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                              :: line_number
      INTEGER, INTENT(out)                                :: status

!  local variables
      REAL (rprec)                                        :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      status = 0
      READ (iou,1000,iostat=status) line
      line_number = line_number + 1
      READ (line,*,iostat=status) diagnostic_dot_parse_node4

      CALL profiler_set_stop_time('diagnostic_dot_parse_node4',                &
     &                            start_time)

1000  FORMAT(a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse coils defined by filaments connected in series.
!>
!>  Parses a diagnostic rogowski defined as a collection of partial rogowskis
!>  connected in series.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @returns A constructed @ref diagnostic_dot_coil instance.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_s_rogowski(iou, line, line_number)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_dot_coil), POINTER ::                                   &
     &   diagnostic_dot_parse_s_rogowski
      INTEGER                                :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                 :: line_number

!  local variables
      REAL (rprec), DIMENSION(:,:), POINTER  :: coil
      REAL (rprec), DIMENSION(:,:), POINTER  :: coil_temp
      REAL (rprec), DIMENSION(4)             :: node
      INTEGER                                :: num_coil
      CHARACTER (len=diagnostic_dot_id_name) :: id_name
      INTEGER                                :: status
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

! Start with 8 coil nodes. The nodes will get reallocated as the coil fills up.
      ALLOCATE(coil(4,8))
      num_coil = 0
      coil_temp => null()

!  Read the id name.
      id_name = diagnostic_dot_parse_name(iou, line, line_number,              &
     &             'diagnostic_dot_parse_s_rogowski')

      status = 0
      DO WHILE(status .eq. 0)
         node = diagnostic_dot_parse_node4(iou, line, line_number,             &
     &                                     status)

         IF (status .eq. 0) THEN
!  If the coil array is full, double it's size in a new array and copy all the
!  data. Delete the coil array and point it to the new array.
            IF (num_coil + 1 .gt. SIZE(coil, 2)) THEN
               ALLOCATE(coil_temp(4,2*SIZE(coil, 2)))
               coil_temp(:,1:SIZE(coil, 2)) = coil
               DEALLOCATE(coil)
               coil => coil_temp
               coil_temp => null()
            END IF

            num_coil = num_coil + 1
            coil(:,num_coil) = node
         END IF

      END DO

      ALLOCATE(diagnostic_dot_parse_s_rogowski)
!  Construct the coil.
      ALLOCATE(diagnostic_dot_parse_s_rogowski)
      ALLOCATE(diagnostic_dot_parse_s_rogowski%coil)
      diagnostic_dot_parse_s_rogowski%next => null()
      diagnostic_dot_parse_s_rogowski%units = TRIM('wb')
      diagnostic_dot_parse_s_rogowski%id_name = id_name

      CALL bsc_construct(diagnostic_dot_parse_s_rogowski%coil,                 &
     &                   'fil_rogo_s', TRIM(id_name), '', 1.0_dp,              &
     &                   coil(1:3,1:num_coil),                                 &
     &                   sen=coil(4,1:num_coil - 1))
      diagnostic_dot_parse_s_rogowski%factor = 1.0

      DEALLOCATE(coil)

      CALL profiler_set_stop_time('diagnostic_dot_parse_s_rogowski',           &
     &                            start_time)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse coils defined by filaments.
!>
!>  Parses a diagnostic flux loop. Flux loops are defined by an id and an
!>  arbitrary number of nodes. Rogowski coils are defined by an id. Next is the
!>  number of turn and area. The guild center is defined as an arbitrary
!>  number of nodes.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[in]    is_rogowski Parse the filaments as a rogowski coil.
!>  @param[in]    units       Units the coil measures in.
!>  @returns A constructed @ref diagnostic_dot_coil instance.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_filaments(iou, line, line_number,          &
     &                                        is_rogowski, units)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_dot_coil), POINTER ::                                   &
     &   diagnostic_dot_parse_filaments
      INTEGER                                :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                 :: line_number
      LOGICAL, INTENT(in)                    :: is_rogowski
      CHARACTER (len=*), INTENT(in)          :: units

!  local variables
      REAL (rprec), DIMENSION(:,:), POINTER  :: coil
      REAL (rprec), DIMENSION(:,:), POINTER  :: coil_temp
      REAL (rprec), DIMENSION(3)             :: node
      INTEGER                                :: num_coil
      INTEGER                                :: status
      CHARACTER (len=diagnostic_dot_id_name) :: id_name
      REAL (rprec)                           :: num_turns
      REAL (rprec)                           :: area
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

! Start with 8 coil nodes. The nodes will get reallocated as the coil fills up.
      ALLOCATE(coil(3,8))
      num_coil = 0
      coil_temp => null()

!  Read the id name.
      id_name = diagnostic_dot_parse_name(iou, line, line_number,              &
     &             'diagnostic_dot_parse_filaments')

      IF (is_rogowski) THEN
!  Read the radius and area
         line = diagnostic_dot_parse_name(iou, line, line_number,              &
     &             'diagnostic_dot_parse_filaments')
         READ(line,*,iostat=status) num_turns, area
         IF (status .ne. 0) THEN
            CALL diagnostic_dot_error('spec', line, line_number,               &
     &              'diagnostic_dot_parse_filaments')
         END IF
         IF (num_turns*area .eq. 0.0) THEN
            WRITE (*,1000) num_turns, area
         END IF
      END IF

      status = 0
      DO WHILE(status .eq. 0)

         node = diagnostic_dot_parse_node(iou, line, line_number,              &
     &                                    status)

         IF (status .eq. 0) THEN
!  If the coil array is full, double it's size in a new array and copy all the
!  data. Delete the coil array and point it to the new array.
            IF (num_coil + 1 .gt. SIZE(coil, 2)) THEN
               ALLOCATE(coil_temp(3,2*SIZE(coil, 2)))
               coil_temp(:,1:SIZE(coil, 2)) = coil
               DEALLOCATE(coil)
               coil => coil_temp
               coil_temp => null()
            END IF

            num_coil = num_coil + 1
            coil(:,num_coil) = node
         END IF

      END DO

!  Check to see if the coil is closed. Borrow the node variable to avoid
!  defining new ones since it is not used anymore. Rogowski coils do not need to
!  be closed.
      IF (.not.is_rogowski .and. num_coil .gt. 2) THEN
         node(1) = SUM((coil(1,2:num_coil) - coil(1,1:num_coil - 1))**2)
         node(2) = SUM((coil(2,2:num_coil) - coil(2,1:num_coil - 1))**2)
         node(3) = SUM((coil(3,2:num_coil) - coil(3,1:num_coil - 1))**2)

!  Store the average length^2 in node(1)
         node(1) = SUM(node)/num_coil

!  Store the length^2 between the last and first nodes in node(2)
         node(2) = SUM((coil(:,1) - coil(:,num_coil))**2)

         IF (node(2)/node(1) .lt. 1.0E-12) THEN
            num_coil = num_coil - 1
         END iF
      END IF

!  Construct the coil.
      ALLOCATE(diagnostic_dot_parse_filaments)
      ALLOCATE(diagnostic_dot_parse_filaments%coil)
      diagnostic_dot_parse_filaments%next => null()
      diagnostic_dot_parse_filaments%units = TRIM(units)
      diagnostic_dot_parse_filaments%id_name = id_name

      IF (is_rogowski) THEN
         CALL bsc_construct(diagnostic_dot_parse_filaments%coil,               &
     &                      'fil_rogo', TRIM(id_name), '', 1.0_dp,             &
     &                      coil(:,1:num_coil), anturns=num_turns,             &
     &                      xsarea=area)
      ELSE
         CALL bsc_construct(diagnostic_dot_parse_filaments%coil,               &
     &                      'fil_loop', TRIM(id_name), '', 1.0_dp,             &
     &                      coil(:,1:num_coil))
      END IF

      SELECT CASE (TRIM(units))

         CASE ('T')
            diagnostic_dot_parse_filaments%factor = 1.0/(num_turns*area)

         CASE ('A')
            diagnostic_dot_parse_filaments%factor =                            &
     &         1.0/(diagnostic_dot_parse_filaments%coil%ave_n_area*            &
     &              2.0E-7*twopi)

         CASE DEFAULT
            diagnostic_dot_parse_filaments%factor = 1.0

      END SELECT

      DEALLOCATE(coil)

      CALL profiler_set_stop_time('diagnostic_dot_parse_filaments',            &
     &                            start_time)

1000  FORMAT('Number of turns (',i6,') x cross section area (',i6,             &
     &       ') cannot be zero.')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse a circular flux loop or magnetic probe.
!>
!>  Parses a diagnostic flux loop or magnetic probe. Flux loops and probes are
!>  defined by an id, the radius, the center and the normal vector. Since this
!>  has a fixed number of reads, the next line must be read before the function
!>  returns.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[in]    is_tokamak  Parse the probe as a tokamak coil.
!>  @param[in]    units       Units the coil measures in.
!>  @returns A constructed @ref diagnostic_dot_coil instance.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_probe(iou, line, line_number,              &
     &                                    is_tokamak, units)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_dot_coil), POINTER ::                                   &
     &   diagnostic_dot_parse_probe
      INTEGER                                :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                 :: line_number
      LOGICAL, INTENT(in)                    :: is_tokamak
      CHARACTER (len=*), INTENT(in)          :: units

!  local variables
      CHARACTER (len=diagnostic_dot_id_name) :: id_name
      REAL (rprec)                           :: radius
      REAL (rprec)                           :: pitch
      REAL (rprec), DIMENSION(3)             :: center
      REAL (rprec), DIMENSION(3)             :: normal
      INTEGER                                :: status
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Parse the id.
      id_name = diagnostic_dot_parse_name(iou, line, line_number,              &
     &             'diagnostic_dot_parse_probe')

      IF (is_tokamak) THEN
!  Parse the magnetic probe spec. Borrow the normal to hold the r, phi, z
!  position.
         status = 0
         READ (iou,1000, iostat=status) line

         IF (status .ne. 0) THEN
            CALL diagnostic_dot_error('spec', line, line_number,               &
     &              'diagnostic_dot_parse_magnetic_probe_tokamak')
         END IF

         READ (line,*,iostat=status) radius, normal, pitch
         IF (status .ne. 0) THEN
            CALL diagnostic_dot_error('spec', line, line_number,               &
     &              'diagnostic_dot_parse_magnetic_probe_tokamak')
         END IF
!  Translate into a bsc_coil. Once all the elements are set over write the
!  values
         normal(2) = normal(2)*degree
         pitch = pitch*degree
         center(1) = normal(1)*COS(normal(2))
         center(2) = normal(1)*SIN(normal(2))
         center(3) = normal(3)
         normal(1) = COS(pitch)*COS(normal(2))
         normal(2) = COS(pitch)*SIN(normal(2))
         normal(3) = SIN(pitch)
      ELSE
!  Parse the radius.
         radius = diagnostic_dot_parse_real(iou, line, line_number,            &
     &               'diagnostic_dot_parse_probe')

!  Parse the center.
         center = diagnostic_dot_parse_node(iou, line, line_number,            &
     &                                      status)
         IF (status .ne. 0) THEN
            CALL diagnostic_dot_error('center', line, line_number,             &
     &              'diagnostic_dot_parse_probe')
         END IF

!  Parse the normal.
         normal = diagnostic_dot_parse_node(iou, line, line_number,            &
     &                                      status)
         IF (status .ne. 0) THEN
            CALL diagnostic_dot_error('normal', line, line_number,             &
     &              'diagnostic_dot_parse_probe')
         END IF
      END IF

!  Construct the coil.
      ALLOCATE(diagnostic_dot_parse_probe)
      ALLOCATE(diagnostic_dot_parse_probe%coil)
      diagnostic_dot_parse_probe%next => null()
      diagnostic_dot_parse_probe%units = TRIM(units)
      diagnostic_dot_parse_probe%id_name = id_name

      SELECT CASE (TRIM(units))

         CASE ('T')
            diagnostic_dot_parse_probe%factor = 1.0/(pi*radius*radius)

         CASE DEFAULT
            diagnostic_dot_parse_probe%factor = 1.0

      END SELECT

      CALL bsc_construct(diagnostic_dot_parse_probe%coil, 'fil_circ',          &
     &                   id_name, '', 1.0_dp, rcirc=radius,                    &
     &                   xcent=center(1:3), enhat=normal(1:3))

!  The function that called this expects a next line to already be read.
      READ (iou,1000,iostat=status) line
      line_number = line_number + 1

      CALL profiler_set_stop_time('diagnostic_dot_parse_probe',                &
     &                            start_time)

1000  FORMAT(a)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Parse a point field measurement.
!>
!>  Parses point measurement probe. Point measurements are defined as a point
!>  and direction. The direction should be normalized. The phi_offset model is
!>  easiest to use when dealing with cyclindical coordinates. Convert all points
!>  to cyclindical coordinates for internal purposes.
!>
!>  @param[in]    iou         File handle of diagnostic dot file.
!>  @param[out]   line        On exit contains the current line of the file.
!>  @param[inout] line_number The current line number of the file.
!>  @param[in]    units       Units the coil measures in.
!>  @param[in]    is_cart     True is point and vector defined in cylindrical
!>                            coordinates.
!>  @returns A constructed @ref diagnostic_dot_coil instance.
!-------------------------------------------------------------------------------
      FUNCTION diagnostic_dot_parse_point(iou, line, line_number,              &
     &                                    units, is_cart)
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (diagnostic_dot_coil), POINTER ::                                   &
     &   diagnostic_dot_parse_point
      INTEGER                                :: iou
      CHARACTER (len=diagnostic_dot_line_len), INTENT(out) :: line
      INTEGER, INTENT(inout)                 :: line_number
      CHARACTER (len=*), INTENT(in)          :: units
      LOGICAL, INTENT(in)                    :: is_cart

!  local variables
      CHARACTER (len=diagnostic_dot_id_name) :: id_name
      REAL (rprec), DIMENSION(3)             :: position
      REAL (rprec), DIMENSION(3)             :: direction
      INTEGER                                :: status
      REAL (rprec)                           :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      id_name = diagnostic_dot_parse_name(iou, line, line_number,              &
     &             'diagnostic_dot_parse_point')

      position = diagnostic_dot_parse_node(iou, line, line_number,             &
     &                                     status)
      IF (status .ne. 0) THEN
         CALL diagnostic_dot_error('position', line, line_number,              &
     &           'diagnostic_dot_parse_point')
      END IF

      direction = diagnostic_dot_parse_node(iou, line, line_number,            &
     &                                      status)
      IF (status .ne. 0) THEN
         CALL diagnostic_dot_error('direction', line, line_number,             &
     &           'diagnostic_dot_parse_point')
      END IF

!  Normalize to a unit vector.
      direction = direction/SQRT(DOT_PRODUCT(direction, direction))

!  Convert to internal cylindical coordinates. If already cylindrical convert
!  the phi direction to radians.
      IF (is_cart) THEN
         direction = cart_to_cyl_vec(position, direction)
         position = cart_to_cyl(position)
      ELSE
         position(2) = position(2)*degree
      END IF

!  Construct the coil.
      ALLOCATE(diagnostic_dot_parse_point)

      diagnostic_dot_parse_point%id_name = id_name
      diagnostic_dot_parse_point%next => null()
      diagnostic_dot_parse_point%position = position
      diagnostic_dot_parse_point%direction = direction
      diagnostic_dot_parse_point%units = TRIM(units)

!  The function that called this expects a next line to already be read.
      READ (iou,1000,iostat=status) line
      line_number = line_number + 1

      CALL profiler_set_stop_time('diagnostic_dot_parse_point',                &
     &                            start_time)

1000  FORMAT(a)

      END FUNCTION

      END MODULE
