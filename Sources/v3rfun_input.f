!-------------------------------------------------------------------------------
!  The @header, @table_section, @table_subsection, @item and @end_table commands
!  are custom defined commands in Doxygen.in. They are defined under ALIASES.
!  For the page created here, the 80 column limit is exceeded. Arguments of
!  aliases are separated by ','. If you intended ',' to be a string you must use
!  an escaped comma '\,'.
!
!>  @page v3rfun_namelist_sec Namelist v3r_coils definition
!>
!>  @tableofcontents
!>  @section v3rfun_namelist_intro_sec Introduction
!>  This page documents the contents of a namelist input file. V3RFUN namelist
!>  variables are defined in the @fixed_width{v3r_coils} common block.
!>
!>  @section v3rfun_namelist_var_sec Namelist Variables
!>  @header{Input variable, Description, Code Reference}
!>
!>  @table_section{v3rfun_filename_sec, Filename Variables}
!>    @item{name_coils_dot,      Filename for the field coil discription.,            v3rfun_input::name_coils_dot}
!>    @item{name_diagnostic_dot, Filename for the field diagnostic coil discription., v3rfun_input::name_diagnostic_dot}
!>  @end_table
!>
!>  @table_section{v3rfun_control_sec, Control Variables}
!>    @item{idrfun,            V3RFUN identification for the run.,                     v3rfun_input::idrfun}
!>    @item{lstell_sym,        Control for stellarator symmetry.,                      v3rfun_input::lstell_sym}
!>    @item{l_read_coils_dot,  Control to ignore the coils dot file.,                  v3rfun_input::l_read_coils_dot}
!>    @item{use_con_shell,     Computes the response function for a conducting shell., v3rfun_input::use_con_shell}
!>    @item{len_integrate_ddc, Integration length in meters.,                          v3rfun_input::len_integrate_ddc}
!>  @end_table
!>
!>  @table_section{v3rfun_super_sec, Super Conducting Coils}
!>    @item{is_super_con, Flags a coil as a super conductor. Super conducting
!>                        are steady state so integated magnetic won't pick
!>                        up an induced signal.,                              v3rfun_input::is_super_con}
!>  @end_table
!>
!>  @table_section{v3rfun_grid_sec, Response Grid Variables}
!>    @item{ir,                  Number of radial grid points.,         v3rfun_input::ir}
!>    @item{jz,                  Number of vertical grid points.,       v3rfun_input::jz}
!>    @item{kp,                  Number of toroidal grid points.,       v3rfun_input::kp}
!>    @item{kp_shell,            Number of shell toroidal grid points., v3rfun_input::kp_shell}
!>    @item{n_field_periods_nli, Number of field periods.,              v3rfun_input::n_field_periods_nli}
!>    @item{rmin,                Minimum R for plasma grid.,            v3rfun_input::rmin}
!>    @item{rmax,                Maximum R for plasma grid.,            v3rfun_input::rmax}
!>    @item{zmin,                Minimum Z for plasma grid.,            v3rfun_input::zmin}
!>    @item{zmax,                Maximum Z for plasma grid.,            v3rfun_input::zmax}
!>    @item{major_radius,        Shell major radius for shell grid.,    v3rfun_input::major_radius}
!>    @item{minor_radius,        Shell minor radius for shell grid.,    v3rfun_input::minor_radius}
!>  @end_table
!>
!>  @table_section{v3rfun_shift_sec, Tilt and Shift Variables}
!>    @item{cg_shift_1,        Vector to shift all the coils before rotation.,                           v3rfun_input::cg_shift_1}
!>    @item{cg_shift_2,        Vector to shift all the coils after rotation.,                            v3rfun_input::cg_shift_2}
!>    @item{cg_rot_xcent,      Position of center of rotation.,                                          v3rfun_input::cg_rot_xcent}
!>    @item{cg_rot_theta,      Spherical polar angle to specify axis of rotation.,                       v3rfun_input::cg_rot_theta}
!>    @item{cg_rot_phi,        Spherical azimuthal angle to specify axis of rotation.,                   v3rfun_input::cg_rot_phi}
!>    @item{cg_rot_angle,      Angle to rotate about axis of rotation. Left hand convention. Put left
!>                             thumb along axis of rotation\, fingers indicate direction of positive
!>                             rotation.,                                                                v3rfun_input::cg_rot_angle}
!>    @item{l_rot_coil_center, Controls the center of rotation.
!>                             * True  Use current-averaged center of coil-group for center of rotation.
!>                             * False Use position specified in cg_rot_xcent for center of rotation.,   v3rfun_input::l_rot_coil_center}
!>  @end_table
!>
!>  @section v3r_coils Example File
!>  @code
!>  ! Example V3RFUN input file.
!>  &v3fit_main_nli
!>    idrfun = 'example'
!>    name_coils_dot = 'coils.example'
!>    name_diagnostic_dot = 'diagnostic.example.dot'
!>    rmin=0.98
!>    rmax=2.02
!>    zmin=-0.52
!>    zmax=0.52
!>    ir=101
!>    jz=101
!>    kp=32
!>    lstell_sym =  F
!>
!>    use_con_shell = T
!>    major_radius = 1.5
!>    minor_radius = 0.525
!>    kp_shell = 128
!>  /
!>  @endcode
!>
!>  @section v3rfun_namelist_prog_ref_sec Programmers Reference
!>  Reference material for the coding to implement this namelist is found in the
!>  @ref v3rfun_input module.
!-------------------------------------------------------------------------------
!*******************************************************************************
!>  @file v3rfun_input.f
!>  @brief Contains module @ref v3rfun_input.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  This file contains all the variables and maximum sizes of the inputs for a
!>  V3RFUN namelist input file. The module contained within does not represent
!>  an object instance. Instead all variables are contained in a global context.
!>  This is required due to limitations of FORTRAN 95 and namelist inputs. Any
!>  information needed by a V3RFUN task should be copied to a
!>  @ref v3rfun_context object. All non-parameters are inputs to the namelist
!>  input file.
!>
!>  @ref v3rfun_namelist_sec "Namelist v3r_coils definition"
!>
!>  @note Some of the references are missing here. This is due to a bug in
!>  Doxygen when variable decalarations span multiple lines.
!*******************************************************************************
      MODULE v3rfun_input
      USE stel_kinds
      USE vsvd0, ONLY: nigroup
      USE profiler

      IMPLICIT NONE

!*******************************************************************************
!  v3rfun input module parameters
!*******************************************************************************
!>  Filename length
      INTEGER, PARAMETER :: v3rfun_file_length = 120
!>  Name length
      INTEGER, PARAMETER :: v3rfun_name_length = 25

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3rfun_namelist_class
!
!*******************************************************************************
!>  Filename for the field coils.
      CHARACTER (len=v3rfun_file_length) :: name_coils_dot = ''
!>  Filename for the diagnostic coils.
      CHARACTER (len=v3rfun_file_length) :: name_diagnostic_dot = ''

!>  v3rfun identification for the run.
      CHARACTER (len=v3rfun_name_length) :: idrfun = ''

!>  Control for stellarator symmetry.
      LOGICAL :: lstell_sym = .false.
!>  Control to ignore the coils dot file.
      LOGICAL :: l_read_coils_dot = .true.
!>  Computes the response function for a conducting shell.
      LOGICAL :: use_con_shell = .false.

!>  Number of radial grid points.
      INTEGER :: ir = 101
!>  Number of vertical grid points.
      INTEGER :: jz = 101
!>  Number of toroidal grid points.
      INTEGER :: kp = 1
!>  Number of shell toroidal grid points.
      INTEGER :: kp_shell = 1
!>  Number of field periods.
      INTEGER :: n_field_periods_nli = 0

!>  Minimum R for plasma grid.
      REAL (rprec) :: rmin = 0.0
!>  Maximum R for plasma grid.
      REAL (rprec) :: rmax = 0.0
!>  Minimum Z for plasma grid.
      REAL (rprec) :: zmin = 0.0
!>  Maximum Z for plasma grid.
      REAL (rprec) :: zmax = 0.0

!>  Shell major radius for shell grid.
      REAL (rprec) :: major_radius = 0.0
!>  Shell minor radius for shell grid.
      REAL (rprec) :: minor_radius = 0.0

!>  Integration length in meters.
      REAL (rprec) :: len_integrate_ddc = 0.001

!>  Tag super conducting coils.
      LOGICAL, DIMENSION(nigroup) :: is_super_con = .false.

!>  Vector to shift all the coils before rotation.
      REAL (rprec), DIMENSION(nigroup,3) :: cg_shift_1 = 0.0
!>  Vector to shift all the coils after rotation.
      REAL (rprec), DIMENSION(nigroup,3) :: cg_shift_2 = 0.0
!>  Position of center of rotation.
      REAL (rprec), DIMENSION(nigroup,3) :: cg_rot_xcent = 0.0
!>  Spherical polar angle to specify axis of rotation.
      REAL (rprec), DIMENSION(nigroup)   :: cg_rot_theta = 0.0
!>  Spherical azimuthal angle to specify axis of rotation.
      REAL (rprec), DIMENSION(nigroup)   :: cg_rot_phi = 0.0
!>  Angle to rotate about axis of rotation.
!>  Left hand convention. Put left thumb along axis of rotation, fingers
!>  indicate direction of positive rotation.
      REAL (rprec), DIMENSION(nigroup)   :: cg_rot_angle = 0.0
!>  Controls the center of rotation.
!>  * True  Use current-averaged center of coil-group for center of rotation.
!>  * False Use position specified in cg_rot_xcent for center of rotation.
      LOGICAL, DIMENSION(nigroup)        :: l_rot_coil_center = .true.

!  Declare Namelist
      NAMELIST/v3r_coils/                                                      &
!  Files
     &   name_coils_dot, name_diagnostic_dot, idrfun,                          &
!  Control parameters
     &   lstell_sym, l_read_coils_dot, use_con_shell,                          &
!  Super conductors
     &   is_super_con,                                                         &
!  Grid sizes
     &   rmin, rmax, zmin, zmax, major_radius, minor_radius,                   &
     &   ir, jz, kp, kp_shell, n_field_periods_nli, len_integrate_ddc,         &
!  Coil tilt and shift
     &   cg_shift_1, cg_shift_2, cg_rot_xcent, cg_rot_theta, cg_rot_phi,       &
     &   cg_rot_angle, l_rot_coil_center

      CONTAINS

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Reads the namelist input file.
!>
!>  Reads the namelist input file. V3RFUN currently requires no processing of
!>  the namelist input.
!>
!>  @param[in] namelist_file The file name of the namelist input file.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_input_read_namelist(namelist_file)
      USE safe_open_mod
      USE v3_utilities

      IMPLICIT NONE

!  Declare Arguments
      CHARACTER (len=*), INTENT(in) :: namelist_file

!  local variables
      INTEGER                       :: iou_mnli
      INTEGER                       :: status
      REAL (rprec)                  :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Initalize a default value of the I\O unit. V3FIT increments from there.
      iou_mnli = 0
      CALL safe_open(iou_mnli, status, TRIM(namelist_file),                    &
     &               'old', 'formatted')
      CALL assert_eq(0, status, 'v3rfun_input_read_namelist' //                &
     &   ': Safe_open of ' // TRIM(namelist_file) // ' failed')

!  Read the namelist input file.
      READ (iou_mnli,nml=v3r_coils)
      CLOSE (iou_mnli,iostat=status)
      CALL assert_eq(0, status, 'v3rfun_input_read_namelist' //                &
     &   ': Error closing ' // TRIM(namelist_file) // ' failed')

      IF (l_read_coils_dot .and. name_coils_dot .eq. '') THEN
         WRITE(*,*) 'l_read_coils_dot is TRUE, but'
         WRITE(*,*) 'name_coils_dot = blank. '
         WRITE(*,*) 'Setting l_read_coils_dot to FALSE'
         l_read_coils_dot = .false.
      END IF

      WRITE (*,*) ' *** V3RFUN namelist input read from ' //                   &
     &            TRIM(namelist_file)

      CALL profiler_set_stop_time('v3rfun_input_read_namelist',                &
     &                            start_time)

      END SUBROUTINE

      END MODULE
