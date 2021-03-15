!*******************************************************************************
!>  @file v3rfun_context.f
!>  @brief Contains module @ref v3rfun_context.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Defines a @ref v3rfun_context_class object to contain all the memory for
!>  running v3rfun.
!*******************************************************************************
      MODULE v3rfun_context
      USE v3rfun_input
      USE diagnostic_dot
      USE biotsavart

      IMPLICIT NONE

!*******************************************************************************
!  v3rfun context module parameters
!*******************************************************************************
!>  Name of this code.
      CHARACTER (len=*), PARAMETER :: code_name = 'V3RFUN'

!*******************************************************************************
!  DERIVED-TYPE DECLARATIONS
!  1) v3rfun context
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Base class representing a v3rfun context. This contains all memory needed to
!>  operate v3rfun.
!-------------------------------------------------------------------------------
      TYPE v3rfun_context_class
!>  The date that the V3RFUN code was run.
         CHARACTER (len=v3rfun_name_length)    :: date_run

!>  The number of phi planes to store the plasma response function.
         INTEGER                               :: kp_store
!>  The number of phi planes to store the conducting shell responce function.
         INTEGER                               :: kp_shell_store

!>  mdsig list input/output file.
         INTEGER                               :: mdsig_list_iou
!>  Mutual inductance list input/output file.
         INTEGER                               :: mi_iou

!>  Linked list containing all the diagnostic coils.
         TYPE (diagnostic_dot_coil), POINTER   :: coils => null()

!>  Array of external field coil currents
         REAL (rprec), DIMENSION(:), POINTER   :: extcur => null()
!>  Array of external field coil responces.
         REAL (rprec), DIMENSION(:), POINTER   :: rdiag => null()

!>  Array of cartesian coordinates for the plasma response function grid.
         REAL (rprec), DIMENSION(:,:,:,:,:), POINTER ::                        &
     &      xcart_grid_e => null()
!>  Array of phi values for the plasma response function grid.
         REAL(rprec), DIMENSION(:,:), POINTER :: phi_grid_e => null()
!>  Array of plasma response functions.
         REAL(rprec), DIMENSION(:,:,:,:), POINTER    ::                        &
     &      pl_response => null()
!>  Array of stellarator symmetric plasma response functions.
         REAL(rprec), DIMENSION(:,:,:,:), POINTER    ::                        &
     &      pl_response_ss => null()

!>  Array of cartesian coordinates for the conducting shell response function
!>  grid.
         REAL (rprec), DIMENSION(:,:,:,:), POINTER ::                          &
     &      xcart_s_grid_e => null()
!>  Array of phi values for the conducting shell response function grid.
         REAL(rprec), DIMENSION(:,:), POINTER      ::                          &
     &      phi_grid_shell => null()
!>  Array of conducting shell response functions.
         REAL(rprec), DIMENSION(:,:,:), POINTER    ::                          &
     &      s_response => null()
!>  Array of stellarator symmetric conducting shell response functions.
         REAL(rprec), DIMENSION(:,:,:), POINTER    ::                          &
     &      s_response_ss => null()
      END TYPE

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
      PRIVATE :: cart2cyl_v

      CONTAINS
!*******************************************************************************
!  CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Construct a @ref v3rfun_context_class object.
!>
!>  Allocates memory and initializes a @ref v3rfun_context_class object.
!>
!>  @returns A pointer to a constructed @ref v3rfun_context_class object.
!-------------------------------------------------------------------------------
      FUNCTION v3rfun_context_construct(filename)
      USE v3_utilities
      USE safe_open_mod
      USE magnetic_response

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3rfun_context_class), POINTER :: v3rfun_context_construct
      CHARACTER (len=*), INTENT(in)        :: filename

!  local variables
      CHARACTER (len=v3rfun_file_length)   :: temp_filename
      CHARACTER (len=v3rfun_file_length)   :: base_name
      CHARACTER (len=8)                    :: date
      CHARACTER (len=10)                   :: time
      CHARACTER (len=5)                    :: zone
      INTEGER                              :: i
      INTEGER                              :: status
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      ALLOCATE(v3rfun_context_construct)

      CALL date_and_time(date, time, zone)
      v3rfun_context_construct%date_run = date // ' ' // time //               &
     &                                    ' ' // zone

!  Read namelist input file.
      CALL v3rfun_input_read_namelist(filename)

!  Open the output file.
      base_name = name_diagnostic_dot(INDEX(name_diagnostic_dot, '/',          &
     &                                      .true.) + 1:)
      temp_filename = TRIM(base_name) // '_mdsig.LIST'
      v3rfun_context_construct%mdsig_list_iou = 0
      WRITE(*,*) 'list_filename is ', temp_filename
      CALL safe_open(v3rfun_context_construct%mdsig_list_iou, status,          &
     &               TRIM(temp_filename), 'replace', 'formatted',              &
     &               delim_in='none')
      CALL assert_eq(0, status, 'v3rfun_context_construct' //                  &
     &   ': Safe_open of ' // TRIM(temp_filename) // ' failed')

      IF (l_read_coils_dot) THEN
         temp_filename = TRIM(base_name) // '_MI'
         v3rfun_context_construct%mi_iou = 0
         WRITE(*,*) 'mi_filename is ', TRIM(temp_filename)
         CALL safe_open(v3rfun_context_construct%mi_iou, status,               &
     &                  TRIM(temp_filename), 'replace', 'formatted',           &
     &                  delim_in='none')
         CALL assert_eq(0, status, 'v3rfun_context_construct' //               &
     &      ': Safe_open of ' // TRIM(temp_filename) // ' failed')

!  Write mutual inductance file header.
         WRITE(v3rfun_context_construct%mi_iou,1000) code_name
         WRITE(v3rfun_context_construct%mi_iou,1100)                           &
     &      magnetic_response_current
         WRITE(v3rfun_context_construct%mi_iou,1200)                           &
     &      TRIM(v3rfun_context_construct%date_run)
         WRITE(v3rfun_context_construct%mi_iou,1300)                           &
     &      TRIM(name_coils_dot)
         WRITE(v3rfun_context_construct%mi_iou,1400)                           &
     &      TRIM(name_diagnostic_dot)
         WRITE(v3rfun_context_construct%mi_iou,*)
      ELSE
         nfp_bs = n_field_periods_nli
      END IF

!  Check stellarator symmetric grid sizes.
      IF (lstell_sym) THEN
         zmax = MAX(ABS(zmax), ABS(zmin))
         zmin = -zmax
         IF (MOD(kp, 2) .ne. 0) THEN
            WRITE (*,1500) 'kp', kp
            CALL EXIT(1)
         END IF
         v3rfun_context_construct%kp_store = kp/2 + 1

         IF (use_con_shell .and. MOD(kp_shell, 2) .ne. 0) THEN
            WRITE (*,1500) 'kp_shell', kp_shell
            CALL EXIT(1)
         END IF
         v3rfun_context_construct%kp_shell_store = kp_shell/2 + 1
      ELSE
         v3rfun_context_construct%kp_store = kp
         v3rfun_context_construct%kp_shell_store = kp_shell
      END IF

      v3rfun_context_construct%coils => null()
      v3rfun_context_construct%coils =>                                        &
     &   diagnostic_dot_construct(name_diagnostic_dot)

      CALL profiler_set_stop_time('v3rfun_context_construct',                  &
     &                            start_time)

      IF (l_read_coils_dot) THEN
         CALL v3rfun_context_construct_field_coils(                            &
     &           v3rfun_context_construct)

         ALLOCATE(v3rfun_context_construct%rdiag(SIZE(coil_group)))
         ALLOCATE(v3rfun_context_construct%extcur(SIZE(coil_group)))

         DO i = 1, SIZE(coil_group)
!  Find extcur array element for computing scaled inductance matrix lines taken
!  out of makegrid code.
            status = coil_group(i)%ncoil
            status = MAXLOC(ABS(coil_group(i)%coils(1:status)%current),        &
     &                      1)
            v3rfun_context_construct%extcur(i) =                               &
     &         coil_group(i)%coils(status)%current
         END DO

         WHERE (v3rfun_context_construct%extcur .eq. 0.0)
            v3rfun_context_construct%extcur = 1.0
         END WHERE


      END IF

      CALL v3rfun_context_construct_responce_grids(                            &
     &        v3rfun_context_construct)

      CALL profiler_set_stop_time('v3rfun_context_construct',                  &
     &                            start_time)

1000  FORMAT('  MUTUAL INDUCTANCES computed by ',a)
1100  FORMAT('  code version is ',a)
1200  FORMAT('  Date run: ',a)
1300  FORMAT('  Field coil information from ',a)
1400  FORMAT('  Diagnostic information from ',a)
1500  FORMAT(a,' (',i4,                                                        &
     &       ') must be even when using stellarator symmetry.')

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Construct and initialize the field coil objects.
!>
!>  Allocates memory and the field coil objects then apply the rotations and
!>  shits as defined in the namelist inout file.
!>
!>  @params[in] this An instance of a @ref v3rfun_context_class.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_context_construct_field_coils(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3rfun_context_class), POINTER :: this

!  local variables
      INTEGER                              :: i, j
      REAL (rprec)                         :: current
      REAL (rprec), DIMENSION(3)           :: center
      REAL (rprec), DIMENSION(3)           :: mean_r
      TYPE (bsc_rs)                        :: rotation_shift
      REAL (rprec)                         :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      IF (name_coils_dot .eq. '') THEN
         WRITE (*,1000) l_read_coils_dot
         CALL EXIT(1)
      END IF

      CALL parse_coils_file(name_coils_dot)
      WRITE (*,1100) TRIM(name_coils_dot), SIZE(coil_group)

!  Rotate and shift the constructed coils.
      IF (SIZE(coil_group) .gt. nigroup) THEN
         WRITE (*,1200) nigroup
         CALL EXIT(1)
      END IF

      WRITE (*,*)
      WRITE (*,1300)

      DO i = 1, SIZE(coil_group)
         WRITE (*,*)
         WRITE (*,1400) i, TRIM(coil_group(i)%s_name)

!  Generate rotation shift for first shift and apply.
         center = 0.0
         CALL bsc_construct_rs(rotation_shift, 0.0_dp, 0.0_dp, 0.0_dp,         &
     &                         center, cg_shift_1(i,:))
         CALL bsc_rot_shift(coil_group(i), rotation_shift)
         WRITE (*,1500) cg_shift_1(i,:)

         IF (l_rot_coil_center(i)) THEN
!  Compute current-averaged center of coil group.
            current = 0.0
            DO j = 1, coil_group(i)%ncoil
               current = current + coil_group(i)%coils(j)%current
               CALL bsc_mean_r(coil_group(i)%coils(j), mean_r)
               center = center                                                 &
     &                + mean_r*coil_group(i)%coils(j)%current
            END DO

            IF (current .ne. 0) THEN
               center = center/current
            END IF
            WRITE (*,1600) center
         ELSE
            center = cg_rot_xcent(i,:)
         END IF
         WRITE (*,1700) center

!  Generate rotation shift for second shift and apply.
         CALL bsc_construct_rs(rotation_shift, cg_rot_theta(i),                &
     &                         cg_rot_phi(i), cg_rot_angle(i),                 &
     &                         center, cg_shift_2(i,:))
         CALL bsc_rot_shift(coil_group(i), rotation_shift)
         WRITE (*,1800) cg_rot_theta(i),  cg_rot_phi(i),                       &
     &                  cg_rot_angle(i)
         WRITE (*,1900) cg_shift_2(i,:)
      END DO

      CALL profiler_set_stop_time(                                             &
     &        'v3rfun_context_construct_field_coils', start_time)

1000  FORMAT('Expected field coil when l_read_coils_dot = ',l)
1100  FORMAT('Coils file ',a,' read, number of coil groups is ',i4)
1200  FORMAT('Number of coil groups exceeds max size (',i3,').',               &
     &       'Increase nigroup in LIBSTELL/Sources/Modules/vsvd0.f')
1300  FORMAT('Rotate and Shift of the Coil Groups')
1400  FORMAT('Coil Group ',i4,' with s_name ',a)
1500  FORMAT('                     First Shift = ',3(2x,es12.5))
1600  FORMAT('   Current-Averaged center of cg = ',3(2x,es12.5))
1700  FORMAT('         Center of Rotation Used = ',3(2x,es12.5))
1800  FORMAT('      Rotation theta, phi, angle = ',3(2x,es12.5))
1900  FORMAT('                    Second Shift = ',3(2x,es12.5))

      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Allocate and initialize the response function arrays
!>
!>  Allocates memory and for the plasma and conducting shell response function
!>  arrays. This does not compute the actual response yet.
!>
!>  @params[in] this An instance of a @ref v3rfun_context_class.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_context_construct_responce_grids(this)
      USE coordinate_utilities

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3rfun_context_class), POINTER :: this

!  local variables
      REAL (rprec)                          :: fperiod
      REAL (rprec)                          :: delr
      REAL (rprec)                          :: delz
      REAL (rprec)                          :: delp
      REAL(rprec), DIMENSION(3)             :: xcyl
      REAL(rprec)                           :: phi0
      INTEGER                               :: k, l, i, j
      REAL (rprec)                          :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

!  Set the grid differental elements
      fperiod = twopi/nfp_bs
      delr = (rmax - rmin)/(ir - 1)
      delz = (zmax - zmin)/(jz - 1)
      delp = fperiod/kp

      ALLOCATE(this%xcart_grid_e(ir,jz,kp,nfp_bs,3))
      ALLOCATE(this%phi_grid_e(kp,nfp_bs))
      ALLOCATE(this%pl_response(ir,jz,kp,3))
      IF (lstell_sym) THEN
         ALLOCATE(this%pl_response_ss(ir,jz,this%kp_store,3))
      ENDIF

!  Compute the coordinates, store them in xcart_grid_e. Also compute and store
!  (in phi_grid_e) the phi values
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,k,l,xcyl,phi0)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO k = 1, kp
         phi0 = (k - 1)*delp
         DO l = 1, nfp_bs
            xcyl(2) = phi0 + (l - 1)*fperiod
            this%phi_grid_e(k,l) = xcyl(2)
            DO i = 1, ir
               xcyl(1) = rmin + (i - 1)*delr
               DO j = 1, jz
                  xcyl(3) = zmin + (j - 1)*delz
                  this%xcart_grid_e(i,j,k,l,1:3) = cyl_to_cart(xcyl)
               END DO ! over j
            END DO ! over i
         END DO ! over l
      END DO ! over k
!$OMP END DO
!$OMP END PARALLEL

      IF (use_con_shell) THEN
!  Set the conducting shell grid differental elements. Reuse the delz as the
!  grid spaceing in the u direction, Reuse delr as value of u. Reuse i as the
!  number of u grid points.
         i = MAX(INT(nfp_bs*kp_shell*minor_radius/major_radius), 10)
         delz = twopi/i
         delp = fperiod/kp_shell

         ALLOCATE(this%xcart_s_grid_e(i,kp_shell,nfp_bs,3))
         ALLOCATE(this%phi_grid_shell(kp_shell,nfp_bs))
         ALLOCATE(this%s_response(i,kp_shell,3))
         IF (lstell_sym) THEN
            ALLOCATE(this%s_response_ss(i,this%kp_shell_store,3))
         END IF

!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,k,l,xcyl,delr,phi0)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO k = 1, kp_shell
            phi0 = (k - 1)*delp
            DO l = 1, nfp_bs
               xcyl(2) = phi0 + (l - 1)*fperiod
               this%phi_grid_shell(k,l) = xcyl(2)
               DO j = 1, i
                  delr = (j - 1)*delz
                  xcyl(1) = major_radius + minor_radius*cos(delr)
                  xcyl(3) = minor_radius*sin(delr)
                  this%xcart_s_grid_e(j,k,l,1:3) = cyl_to_cart(xcyl)
               END DO ! over j
            END DO ! over l
         END DO ! over k
!$OMP END DO
!$OMP END PARALLEL
      END IF

      CALL profiler_set_stop_time(                                             &
     &        'v3rfun_context_construct_responce_grids', start_time)

      END SUBROUTINE

!*******************************************************************************
!  DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Deconstruct a @ref v3rfun_context_class object.
!>
!>  Deallocates memory and uninitializes a @ref v3rfun_context_class object.
!>  This all the v3rfun constructed objects are destroyed here.
!>
!>  @param[inout] this A @ref v3rfun_context_class instance.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_context_destruct(this)

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3rfun_context_class), POINTER :: this

!  Start of executable code
      CLOSE(this%mdsig_list_iou)

      IF (l_read_coils_dot) THEN
         CLOSE(this%mi_iou)
      END IF

      CALL cleanup_biotsavart

      IF (ASSOCIATED(this%coils)) THEN
         CALL diagnostic_dot_destruct(this%coils)
         this%coils => null()
      END IF

      IF (ASSOCIATED(this%xcart_grid_e)) THEN
         DEALLOCATE(this%xcart_grid_e)
         this%xcart_grid_e => null()
      END IF

      IF (ASSOCIATED(this%phi_grid_e)) THEN
         DEALLOCATE(this%phi_grid_e)
         this%phi_grid_e => null()
      END IF

      IF (ASSOCIATED(this%pl_response)) THEN
         DEALLOCATE(this%pl_response)
         this%pl_response => null()
      END IF

      IF (ASSOCIATED(this%pl_response_ss)) THEN
         DEALLOCATE(this%pl_response_ss)
         this%pl_response_ss => null()
      END IF

      IF (ASSOCIATED(this%xcart_s_grid_e)) THEN
         DEALLOCATE(this%xcart_s_grid_e)
         this%xcart_s_grid_e => null()
      END IF

      IF (ASSOCIATED(this%phi_grid_shell)) THEN
        DEALLOCATE(this%phi_grid_shell)
        this%phi_grid_shell => null()
      END IF

      IF (ASSOCIATED(this%s_response)) THEN
         DEALLOCATE(this%s_response)
         this%s_response => null()
      END IF

      IF (ASSOCIATED(this%s_response_ss)) THEN
         DEALLOCATE(this%s_response_ss)
         this%s_response_ss => null()
      END IF

      DEALLOCATE(this)

      END SUBROUTINE

!*******************************************************************************
!  UTILITY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Write out a magnetic coil response function.
!>
!>  Computes the magnetic response function for each a coil and writes the
!>  results to an mdsig file.
!>
!>  @param[inout] this   A @ref v3rfun_context_class instance.
!>  @param[in]    d_coil A @ref diagnostic_dot::diagnostic_dot_coil instance.
!>  @param[in]    id_num Diagnostic coil number.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_context_write_mrf(this, d_coil, id_num)
      USE magnetic_response
      USE v3_utilities
      USE bsc_cdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3rfun_context_class), INTENT(inout) :: this
      TYPE (diagnostic_dot_coil), INTENT(in)     :: d_coil
      INTEGER, INTENT(in)                        :: id_num

!  local variables
      INTEGER                                    :: i, j, k, l
      INTEGER                                    :: iss, jss, kss
      REAL(rprec), DIMENSION(3)                  :: acart
      TYPE (magnetic_response_class), POINTER    :: response
      CHARACTER (len=v3rfun_file_length)         :: temp_filename
      REAL(rprec), DIMENSION(:,:,:), POINTER     :: a_r
      REAL(rprec), DIMENSION(:,:,:), POINTER     :: a_f
      REAL(rprec), DIMENSION(:,:,:), POINTER     :: a_z
      REAL(rprec), DIMENSION(:,:), POINTER       :: a_s_r
      REAL(rprec), DIMENSION(:,:), POINTER       :: a_s_f
      REAL(rprec), DIMENSION(:,:), POINTER       :: a_s_z
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (*,1000) id_num

!  Coil Response ---------------------------------------------------------------
      IF (l_read_coils_dot) THEN
         DO i = 1, SIZE(coil_group)
!  Super conducting coils are steady state and will not induce a signal in
!  integrated magnetic diagnostics. For those coils set the mutual inductance to
!  zero.
            IF (.not.is_super_con(i)) THEN
               CALL bsc_fluxba(coil_group(i), d_coil%coil,                     &
     &                         len_integrate_ddc, this%rdiag(i))
            ELSE
               this%rdiag(i) = 0.0
            END IF
         END DO

         this%rdiag = this%rdiag*d_coil%factor
      END IF

!  Plasma Response -------------------------------------------------------------
!  Parallelize the outer most loops. Any local variables inside the loop must
!  be private.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(i,j,k,l,acart)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
      DO i = 1, ir
         DO j = 1, jz
            DO k = 1, kp
               this%pl_response(i,j,k,1:3) = 0.0
               DO l = 1, nfp_bs
                  CALL bsc_a(d_coil%coil,                                      &
     &                       this%xcart_grid_e(i,j,k,l,1:3), acart)
                  this%pl_response(i,j,k,1:3) =                                &
     &               this%pl_response(i,j,k,1:3) +                             &
     &               cart2cyl_v(acart, this%phi_grid_e(k,l))
               END DO ! l, field periods
            END DO ! k
         END DO ! j
      END DO ! i
!$OMP END DO
!$OMP END PARALLEL
      this%pl_response = this%pl_response*d_coil%factor

!  When using stellator symmetry, only half of the field period needs to be
!  stored due to the symmetry. A vertical point will have the same current as a
!  the symmetric point below the midplane. The counter point to point j is found
!  from by counting the grid in reverse.
      IF (lstell_sym) THEN
         DO i = 1, ir
            iss = i
            DO j = 1, jz
               jss = jz + 1 - j
               DO k = 1, SIZE(this%pl_response_ss, 3)
                  IF (k .eq. 1) THEN
                     kss = 1
                  ELSE
                     kss = kp + 2 - k
                  END IF
                  this%pl_response_ss(i,j,k,1) =                               &
     &               this%pl_response(i,j,k,1) -                               &
     &               this%pl_response(iss,jss,kss,1)
                  this%pl_response_ss(i,j,k,2) =                               &
     &               this%pl_response(i,j,k,2) +                               &
     &               this%pl_response(iss,jss,kss,2)
                  this%pl_response_ss(i,j,k,3) =                               &
     &               this%pl_response(i,j,k,3) +                               &
     &               this%pl_response(iss,jss,kss,3)
               END DO ! k
            END DO ! j
         END DO ! i
      END IF ! lstell_sym

!  Conducting Shell Response ---------------------------------------------------
      IF (use_con_shell) THEN
!  Parallelize the outer most loops. Any local variables inside the loop must
!  be private.
!$OMP PARALLEL
!$OMP& DEFAULT(SHARED)
!$OMP& PRIVATE(j,k,l,acart)
!$OMP DO
!$OMP& SCHEDULE(STATIC)
         DO j = 1, SIZE(this%s_response, 1)
            DO k = 1, kp_shell
               this%s_response(j,k,1:3) = 0.0
               DO l = 1, nfp_bs
                  CALL bsc_a(d_coil%coil,                                      &
     &                       this%xcart_s_grid_e(j,k,l,1:3), acart)
                  this%s_response(j,k,1:3) = this%s_response(j,k,1:3) +        &
     &               cart2cyl_v(acart, this%phi_grid_shell(k,l))
               END DO
            END DO
         END DO
!$OMP END DO
!$OMP END PARALLEL
         this%s_response = this%s_response*d_coil%factor

!  When using stellator symmetry, only half of the field period needs to be
!  stored due to the symmetry. A poloidal point will have the same current as a
!  the symmetric point below the midplane. The j=1 point corresponds to the
!  same u position on the opposite side of the field perior. For j!=0, the
!  symmetry point is found by counting indices backwards.
         IF (lstell_sym) THEN
            DO j = 1, SIZE(this%s_response, 1)
               IF (j .eq. 1) THEN
                  jss = 1
               ELSE
                  jss = SIZE(this%s_response, 1) + 2 - j
               END IF
               DO k = 1, SIZE(this%s_response_ss, 2)
               IF (k .eq. 1) THEN
                  kss = 1
               ELSE
                  kss = kp_shell + 2 - k
               END IF

               this%s_response_ss(j,k,1) = this%s_response(j,k,1)              &
     &                                   - this%s_response(jss,kss,1)
               this%s_response_ss(j,k,2) = this%s_response(j,k,2)              &
     &                                   + this%s_response(jss,kss,2)
               this%s_response_ss(j,k,3) = this%s_response(j,k,3)              &
     &                                   + this%s_response(jss,kss,3)
               END DO
            END DO
         END IF
      END IF

!  Assign the response function grids to pointers.
      a_r => null()
      a_f => null()
      a_z => null()
      a_s_r => null()
      a_s_f => null()
      a_s_z => null()
      IF (lstell_sym) THEN
         a_r => this%pl_response_ss(:,:,:,1)
         a_f => this%pl_response_ss(:,:,:,2)
         a_z => this%pl_response_ss(:,:,:,3)
         IF (use_con_shell) THEN
            a_s_r => this%s_response_ss(:,:,1)
            a_s_f => this%s_response_ss(:,:,2)
            a_s_z => this%s_response_ss(:,:,3)
         END IF
      ELSE
         a_r => this%pl_response(:,:,:,1)
         a_f => this%pl_response(:,:,:,2)
         a_z => this%pl_response(:,:,:,3)
         IF (use_con_shell) THEN
            a_s_r => this%s_response(:,:,1)
            a_s_f => this%s_response(:,:,2)
            a_s_z => this%s_response(:,:,3)
         END IF
      ENDIF

!  Construct a responce object.
      response => magnetic_response_construct(code_name, this%date_run,        &
     &                                        d_coil%id_name,                  &
     &                                        this%rdiag, this%extcur,         &
     &                                        kp, kp_shell,                    &
     &                                        rmin, rmax, zmin, zmax,          &
     &                                        nfp_bs, lstell_sym,              &
     &                                        a_r, a_f, a_z,                   &
     &                                        a_s_r, a_s_f, a_s_z,             &
     &                                        0.0_dp)

!  Write Coil Response ---------------------------------------------------------
!  Open up an mdsig file, for netcdf writing. Borrow the i variable to store the
!  status. Borrow the j variable as the netcdf input/output unit.
      temp_filename = TRIM(d_coil%id_name) // '_mdsig.nc'
      temp_filename = TRIM(ADJUSTL(temp_filename))
      i = 0
      CALL cdf_open(j, TRIM(temp_filename), 'w', i)
      CALL assert_eq(0, i, 'mdsig file ' // TRIM(temp_filename) //             &
     &               ' failed to open.')

!  Define NetCDF variables.
      CALL cdf_define(j, 'diagnostic_desc_d_type', d_coil%d_type)
      CALL cdf_define(j, 'diagnostic_desc_s_name', d_coil%id_name)
      CALL cdf_define(j, 'diagnostic_desc_l_name', d_coil%id_name)
      CALL cdf_define(j, 'diagnostic_desc_units', d_coil%units)
      CALL cdf_define(j, 'diagnostic_desc_sigma_default', 0.0)
      CALL magnetic_response_define(response, j)
      CALL bsc_cdf_define_coil(d_coil%coil, j, '')

!  Write NetCDF variables. The cdf_write subroutine switches from define mode to
!  write mode.
      CALL cdf_write(j, 'diagnostic_desc_d_type', d_coil%d_type)
      CALL cdf_write(j, 'diagnostic_desc_s_name', d_coil%id_name)
      CALL cdf_write(j, 'diagnostic_desc_l_name', d_coil%id_name)
      CALL cdf_write(j, 'diagnostic_desc_units', d_coil%units)
      CALL cdf_write(j, 'diagnostic_desc_sigma_default', 0.0)
      CALL magnetic_response_write(response, j)
      CALL bsc_cdf_write_coil(d_coil%coil, j, '')

!  Close NetCDF file.
      CALL cdf_close(j)

!  Write entry in list file
      WRITE(this%mdsig_list_iou,1100) id_num, TRIM(temp_filename)

!  Write to Mutual Inductance file
      IF (l_read_coils_dot) THEN
         WRITE (this%mi_iou,1200) id_num
         WRITE (this%mi_iou,1300) TRIM(d_coil%id_name)
         WRITE (this%mi_iou,1400) TRIM(d_coil%d_type)
         WRITE (this%mi_iou,1500) TRIM(d_coil%units)
         WRITE (this%mi_iou,1600)
         DO i = 1, SIZE(coil_group)
            WRITE (this%mi_iou,1700) i, TRIM(coil_group(i)%s_name),            &
     &                               this%rdiag(i)
         END DO
         WRITE (this%mi_iou,*)
      END IF

!  Cleanup
      CALL magnetic_response_destruct(response)

      CALL profiler_set_stop_time('v3rfun_context_write_mrf',                  &
     &                            start_time)

1000  FORMAT(5x,'Diagnostic #',i4)
1100  FORMAT(i4.4,x,a)
1200  FORMAT('Diagnostic Coil:',3x,i4)
1300  FORMAT('Short Name:',8x,a)
1400  FORMAT('MDDC Type:',9x,a)
1500  FORMAT('Signal units:',6x,a)
1600  FORMAT(3x,'i',1x,'ID',12x,'Inductance')
1700  FORMAT(i4,1x,a11,2x,es14.6)
      END SUBROUTINE

!-------------------------------------------------------------------------------
!>  @brief Write out a point probe response.
!>
!>  Computes the response function for the vacuum field at the point
!>  measurement. Point probes are assumed to be outside of the plasma.
!>
!>  @param[inout] this   A @ref v3rfun_context_class instance.
!>  @param[in]    d_coil A @ref diagnostic_dot::diagnostic_dot_coil instance.
!>  @param[in]    id_num Diagnostic coil number.
!-------------------------------------------------------------------------------
      SUBROUTINE v3rfun_context_write_point(this, d_coil, id_num)
      USE magnetic_response
      USE v3_utilities
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      TYPE (v3rfun_context_class), INTENT(inout) :: this
      TYPE (diagnostic_dot_coil), INTENT(in)     :: d_coil
      INTEGER, INTENT(in)                        :: id_num

!  local variables
      INTEGER                                    :: i
      INTEGER                                    :: mdsig_iou
      REAL (rprec), DIMENSION(3)                 :: b_vec
      TYPE (magnetic_response_class), POINTER    :: response
      CHARACTER (len=v3rfun_file_length)         :: temp_filename
      REAL (rprec)                               :: start_time

!  Start of executable code
      start_time = profiler_get_start_time()

      WRITE (*,1000) id_num

      IF (l_read_coils_dot) THEN
         DO i = 1, SIZE(coil_group)
           CALL bsc_b(coil_group(i), d_coil%position, b_vec)
           this%rdiag(i) = DOT_PRODUCT(b_vec, d_coil%direction)
         END DO
      END IF

      response => magnetic_response_construct(code_name, this%date_run,        &
     &                                        d_coil%id_name,                  &
     &                                        d_coil%position,                 &
     &                                        d_coil%direction,                &
     &                                        this%rdiag, this%extcur)

!  Write Coil Response ---------------------------------------------------------
!  Open up an mdsig file, for netcdf writing. Borrow the i variable to store the
!  status.
      temp_filename = TRIM(d_coil%id_name) // '_mdsig.nc'
      temp_filename = TRIM(ADJUSTL(temp_filename))
      mdsig_iou = 0
      CALL cdf_open(mdsig_iou, TRIM(temp_filename), 'w', i)
      CALL assert_eq(0, i, 'mdsig file ' // TRIM(temp_filename) //             &
     &               ' failed to open.')

!  Define NetCDF variables.
      CALL cdf_define(mdsig_iou, 'diagnostic_desc_d_type',                     &
     &                d_coil%d_type)
      CALL cdf_define(mdsig_iou, 'diagnostic_desc_s_name',                     &
     &                d_coil%id_name)
      CALL cdf_define(mdsig_iou, 'diagnostic_desc_l_name',                     &
     &                d_coil%id_name)
      CALL cdf_define(mdsig_iou, 'diagnostic_desc_units', d_coil%units)
      CALL cdf_define(mdsig_iou, 'diagnostic_desc_sigma_default', 0.0)
      CALL magnetic_response_define(response, mdsig_iou)

!  Write NetCDF variables. The cdf_write subroutine switches from define mode to
!  write mode.
      CALL cdf_write(mdsig_iou, 'diagnostic_desc_d_type', d_coil%d_type)
      CALL cdf_write(mdsig_iou, 'diagnostic_desc_s_name',                      &
     &               d_coil%id_name)
      CALL cdf_write(mdsig_iou, 'diagnostic_desc_l_name',                      &
     &               d_coil%id_name)
      CALL cdf_write(mdsig_iou, 'diagnostic_desc_units', d_coil%units)
      CALL cdf_write(mdsig_iou, 'diagnostic_desc_sigma_default', 0.0)
      CALL magnetic_response_write(response, mdsig_iou)

      CALL cdf_close(mdsig_iou)

!  Write entry in list file
      WRITE(this%mdsig_list_iou,1100) id_num, TRIM(temp_filename)

!  Cleanup
      CALL magnetic_response_destruct(response)

      CALL profiler_set_stop_time('rfun_context_write_point',                  &
     &                            start_time)

1000  FORMAT(5x,'Diagnostic #',i4)
1100  FORMAT(i4.4,x,a)
      END SUBROUTINE

!*******************************************************************************
!  PRIVATE
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Converts a vector from Cartesian coordinates to cylindrical coordinates
!>
!>    [vr  ]   [ Cos(phi) Sin(phi) 0]   [vx]
!>    [vPhi] = [-Sin(phi) Cos(phi) 0] * [vy]
!>    [vZ  ]   [        0        0 1]   [vz]
!>
!>  @param[in] vcart The vector to convert.
!>  @param[in] phi   The toroidal angle to perform the conversion at.
!>  @returns The converted vector.
!-------------------------------------------------------------------------------
      FUNCTION cart2cyl_v(vcart, phi)

      IMPLICIT NONE

!  Declare Arguments
      REAL (rprec), DIMENSION(3)             :: cart2cyl_v
      REAL (rprec), DIMENSION(3), INTENT(in) :: vcart
      REAL (rprec), INTENT(in)               :: phi

!  local variables
      REAL (rprec)                           :: cphi
      REAL (rprec)                           :: sphi

!  Start of executable code
      cphi = cos(phi)
      sphi = sin(phi)
      cart2cyl_v(1) = vcart(1)*cphi + vcart(2)*sphi
      cart2cyl_v(2) = -vcart(1)*sphi + vcart(2)*cphi
      cart2cyl_v(3) = vcart(3)

      END FUNCTION

      END MODULE
