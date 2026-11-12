#  Copy input file from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/coils.circle
                ${CMAKE_CURRENT_BINARY_DIR}/coils.circle
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/magnetic.dot
                ${CMAKE_CURRENT_BINARY_DIR}/magnetic.dot
                COPYONLY)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/v3rfun.in
                ${CMAKE_CURRENT_BINARY_DIR}/v3rfun.in
                COPYONLY)

#  Run V3RFUN.
add_test (NAME    v3rfun_single_phi_plane_v3rfun_test
          COMMAND $<TARGET_PROPERTY:xv3rfun,BINARY_DIR>/xv3rfun v3rfun.in)
set_tests_properties (v3rfun_single_phi_plane_v3rfun_test
                      PROPERTIES
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)

#  Check mdsig files to ensure that values match.
foreach (i RANGE 90 270 90)
    add_test (NAME    v3rfun_single_phi_plane_mddc_mrf_a_f_${i}_test
              COMMAND $<TARGET_PROPERTY:xnetcdf_diff,BINARY_DIR>/xnetcdf_diff -netcdf_file1=probe_0_mdsig.nc -netcdf_file2=probe_${i}_mdsig.nc -quantity=mddc_mrf_a_f -tol=6.0E-18 -sum)
    set_tests_properties (v3rfun_single_phi_plane_mddc_mrf_a_f_${i}_test
                          PROPERTIES
                          DEPENDS "v3rfun_single_phi_plane_v3rfun_test")
endforeach ()
