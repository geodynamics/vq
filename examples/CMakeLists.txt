# Set the example directory path
SET(VC_EXAMPLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/)

# Copy the fault traces to the build directory for use in the tests
FILE(COPY ${VC_EXAMPLE_DIR}/fault_traces DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)

SET(TAPER_METHODS none taper taper_renorm)
SET(METHOD_MEANSLIP 100.4 0 0)
SET(METHOD_INTEREVENT 49.2 0 0)
SET(RESOLUTIONS 12000 6000 4000 3000 2000)

SET(SETUP_PARAMS_SCRIPT ${VC_EXAMPLE_DIR}/single_fault/setup_params.sh)
SET(CHECK_SCRIPT ${VC_EXAMPLE_DIR}/check_results.py)
SET(ENV{PATH} $ENV{PATH}:${CMAKE_CURRENT_BINARY_DIR}/../src/)

LIST(LENGTH TAPER_METHODS N)
MATH(EXPR NUM_TAPER "${N}-1")

FOREACH(TAPER_IND RANGE ${NUM_TAPER})
    LIST(GET TAPER_METHODS ${TAPER_IND} TAPER)
    LIST(GET METHOD_MEANSLIP ${TAPER_IND} MEANSLIP)
    LIST(GET METHOD_INTEREVENT ${TAPER_IND} INTEREVENT)
    FOREACH(RES ${RESOLUTIONS})
        SET(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TAPER}/)
        FILE(MAKE_DIRECTORY ${TEST_DIR})
        SET(EVENT_FILE ${TEST_DIR}/events_${RES}.txt)
        SET(SWEEP_FILE ${TEST_DIR}/sweeps_${RES}.txt)
        FILE(COPY ${VC_EXAMPLE_DIR}/single_fault/params.d DESTINATION ${TEST_DIR})

        # Test that the mesher works and generate a parameter file
        ADD_TEST(
            NAME mesh_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND mesher
            --import_file=../../fault_traces/single_fault_trace.txt
            --import_file_type=trace --import_trace_element_size=${RES}
            --taper_fault_method=${TAPER}
            --export_file=single_fault_${RES}.txt
            --export_file_type=text
            --export_file=single_fault_${RES}.kml
            --export_file_type=kml
            --print_statistics=statistics_${RES}.txt
            )
        ADD_TEST(
            NAME param_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND ${SETUP_PARAMS_SCRIPT} ${RES} 0.0
            )
        SET_TESTS_PROPERTIES (
            param_${TAPER}_${RES} PROPERTIES DEPENDS mesh_${TAPER}_${RES}
            )

        # Test that the simulation runs with the resulting mesh
        ADD_TEST(
            NAME run_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND vc params_${RES}.d
            )
        SET_TESTS_PROPERTIES (
            run_${TAPER}_${RES} PROPERTIES DEPENDS param_${TAPER}_${RES}
            )

        # Check that the results are internally self consistent
        ADD_TEST(
            NAME test_consistent_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --check_consistent
            )
        SET_TESTS_PROPERTIES (test_consistent_${TAPER}_${RES} PROPERTIES DEPENDS run_${TAPER}_${RES})

        # Check that the mean slip is near our expectations
        IF (NOT ${MEANSLIP} EQUAL 0)
            ADD_TEST(
                NAME test_slip_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_slip ${MEANSLIP}
                )
            SET_TESTS_PROPERTIES (
                test_slip_${TAPER}_${RES} PROPERTIES DEPENDS run_${TAPER}_${RES}
                )
        ENDIF (NOT ${MEANSLIP} EQUAL 0)

        # Check that the mean interevent time is near our expectations
        IF (NOT ${INTEREVENT} EQUAL 0)
            ADD_TEST(
                NAME test_interevent_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_interevent ${INTEREVENT}
                )
            SET_TESTS_PROPERTIES (
                test_interevent_${TAPER}_${RES} PROPERTIES DEPENDS run_${TAPER}_${RES}
                )
        ENDIF(NOT ${INTEREVENT} EQUAL 0)
    ENDFOREACH()
ENDFOREACH()

