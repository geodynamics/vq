# Set the example directory path
SET(VC_EXAMPLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/)

# Copy the fault traces to the build directory for use in the tests
file(COPY ${VC_EXAMPLE_DIR}/fault_traces DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)

SET(RESOLUTIONS 12000 6000 4000 3000 2000)

SET(SETUP_SCRIPT ${VC_EXAMPLE_DIR}/single_fault/setup_mesh.sh)
SET(CHECK_SCRIPT ${VC_EXAMPLE_DIR}/check_results.py)

# Test the mesher program to ensure it creates faults at different resolutions
FOREACH(RES ${RESOLUTIONS})
    SET(EVENT_FILE ${CMAKE_CURRENT_BINARY_DIR}/none/events_${RES}.txt)
    SET(SWEEP_FILE ${CMAKE_CURRENT_BINARY_DIR}/none/sweeps_${RES}.txt)

    # Test that the mesher works
    ADD_TEST(NAME mesh_single_${RES} COMMAND ${SETUP_SCRIPT} ${RES} none 0.0)
    # Test that the simulation runs with the resulting mesh
    ADD_TEST(NAME run_single_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/none/ COMMAND vc params_${RES}.d)
    SET_TESTS_PROPERTIES (run_single_${RES} PROPERTIES DEPENDS mesh_single_${RES})
    # Check that the results are internally self consistent
    ADD_TEST(NAME test_consistent_single_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/none/ COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --check_consistent)
    SET_TESTS_PROPERTIES (test_consistent_single_${RES} PROPERTIES DEPENDS run_single_${RES})
    # Check that the mean slip is near our expectations
    ADD_TEST(NAME test_slip_single_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/none/ COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_slip 100.4)
    SET_TESTS_PROPERTIES (test_slip_single_${RES} PROPERTIES DEPENDS run_single_${RES})
    # Check that the mean interevent time is near our expectations
    ADD_TEST(NAME test_interevent_single_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/none/ COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_interevent 49.2)
    SET_TESTS_PROPERTIES (test_interevent_single_${RES} PROPERTIES DEPENDS run_single_${RES})
ENDFOREACH()

# Test meshing with fault tapering
FOREACH(RES ${RESOLUTIONS})
    SET(EVENT_FILE ${CMAKE_CURRENT_BINARY_DIR}/taper/events_${RES}.txt)
    SET(SWEEP_FILE ${CMAKE_CURRENT_BINARY_DIR}/taper/sweeps_${RES}.txt)

    # Test that the mesher works
    ADD_TEST(NAME mesh_taper_${RES} COMMAND ${SETUP_SCRIPT} ${RES} taper 0.0)
    # Test that the simulation runs with the resulting mesh
    ADD_TEST(NAME run_taper_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/taper/ COMMAND vc params_${RES}.d)
    SET_TESTS_PROPERTIES (run_taper_${RES} PROPERTIES DEPENDS mesh_taper_${RES})
    # Check that the results are internally self consistent
    ADD_TEST(NAME test_consistent_taper_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/taper/ COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --check_consistent)
    SET_TESTS_PROPERTIES (test_consistent_taper_${RES} PROPERTIES DEPENDS run_taper_${RES})
ENDFOREACH()

# Test meshing with fault tapering and velocity renormalization
FOREACH(RES ${RESOLUTIONS})
    SET(EVENT_FILE ${CMAKE_CURRENT_BINARY_DIR}/taper_renorm/events_${RES}.txt)
    SET(SWEEP_FILE ${CMAKE_CURRENT_BINARY_DIR}/taper_renorm/sweeps_${RES}.txt)

    # Test that the mesher works
    ADD_TEST(NAME mesh_taper_renorm_${RES} COMMAND ${SETUP_SCRIPT} ${RES} taper_renorm 0.0)
    # Test that the simulation runs with the resulting mesh
    ADD_TEST(NAME run_taper_renorm_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/taper_renorm/ COMMAND vc params_${RES}.d)
    SET_TESTS_PROPERTIES (run_taper_renorm_${RES} PROPERTIES DEPENDS mesh_taper_renorm_${RES})
    # Check that the results are internally self consistent
    ADD_TEST(NAME test_consistent_taper_renorm_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/taper_renorm/ COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --check_consistent)
    SET_TESTS_PROPERTIES (test_consistent_taper_renorm_${RES} PROPERTIES DEPENDS run_taper_renorm_${RES})
    # Check that the mean slip is near our expectations
    ADD_TEST(NAME test_slip_taper_renorm_${RES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/taper_renorm/ COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_slip 100.4)
    SET_TESTS_PROPERTIES (test_slip_taper_renorm_${RES} PROPERTIES DEPENDS run_taper_renorm_${RES})
ENDFOREACH()

