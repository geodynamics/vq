# Set the example directory path
SET(VC_EXAMPLE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/)
SET(VC_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/../src/)

# Copy the fault traces to the build directory for use in the tests
FILE(COPY ${VC_EXAMPLE_DIR}/fault_traces DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)
FILE(COPY ${VC_EXAMPLE_DIR}/ca_model DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)

SET(TAPER_METHODS none taper taper_renorm)
SET(METHOD_MEANSLIP 100.4 0 0)
SET(METHOD_INTEREVENT 49.2 0 0)

# If MPI is available, set up MPI based tests as well
IF(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)
    SET(NUM_PROCS 1 2 4)
ELSE(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)
    SET(NUM_PROCS 1)
ENDIF(DEFINED MPI_C_FOUND AND MPI_CXX_FOUND)

SET(SETUP_PARAMS_SCRIPT ${VC_EXAMPLE_DIR}/setup_params.sh)
SET(CHECK_SCRIPT ${VC_EXAMPLE_DIR}/check_results.py)
SET(ENV{PATH} $ENV{PATH}:${CMAKE_CURRENT_BINARY_DIR}/../src/)

LIST(LENGTH TAPER_METHODS N)
MATH(EXPR NUM_TAPER "${N}-1")
 
# --------------  Single Fault Section Tests ----------------------------
FOREACH(NPROC ${NUM_PROCS})
    # Don't use 12k meshing (i.e. single element) if more than 1 proc (otherwise VC will fail)
    IF (NPROC EQUAL 1)
        SET(RESOLUTIONS     12000 6000 4000 3000 2000)
    ELSE (NPROC EQUAL 1)
        SET(RESOLUTIONS     6000 4000 3000 2000)
    ENDIF (NPROC EQUAL 1)

    FOREACH(TAPER_IND RANGE ${NUM_TAPER})
        LIST(GET TAPER_METHODS ${TAPER_IND} TAPER)
        LIST(GET METHOD_MEANSLIP ${TAPER_IND} MEANSLIP)
        LIST(GET METHOD_INTEREVENT ${TAPER_IND} INTEREVENT)
        FOREACH(RES ${RESOLUTIONS})
            SET(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/PROCS${NPROC}/${TAPER}/)
            FILE(MAKE_DIRECTORY ${TEST_DIR})
            SET(EVENT_FILE     ${TEST_DIR}events_${RES}.txt)
            SET(SWEEP_FILE     ${TEST_DIR}sweeps_${RES}.txt)
            

            # Test that the mesher works and generate a parameter file
            ADD_TEST(
                NAME mesh_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                COMMAND mesher
                --import_file=../../../fault_traces/single_fault_trace.txt
                --import_file_type=trace --import_trace_element_size=${RES}
                --taper_fault_method=${TAPER}
                --export_file=single_fault_${RES}.txt
                --export_file_type=text
                --export_file=single_fault_${RES}.kml
                --export_file_type=kml
                --print_statistics=statistics_${RES}.txt
                )
            ADD_TEST(NAME param_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                COMMAND ${SETUP_PARAMS_SCRIPT} ${RES} 0.7 single_fault ${VC_EXAMPLE_DIR}/example_params.d)
            SET_TESTS_PROPERTIES (param_P${NPROC}_${TAPER}_${RES} PROPERTIES DEPENDS mesh_P${NPROC}_${TAPER}_${RES})

            # Test that the simulation runs with the resulting mesh
            IF (NPROC EQUAL 1)
                ADD_TEST(NAME run_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                    COMMAND ${VC_BINARY_DIR}/vc params_${RES}.d)
            ELSE (NPROC EQUAL 1)
                ADD_TEST(NAME run_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${VC_BINARY_DIR}/vc params_${RES}.d)
            ENDIF (NPROC EQUAL 1)
            SET_TESTS_PROPERTIES (run_P${NPROC}_${TAPER}_${RES} PROPERTIES DEPENDS param_P${NPROC}_${TAPER}_${RES})

            # Check that the results are internally self consistent
            ADD_TEST(NAME test_consistent_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --check_consistent)
            SET_TESTS_PROPERTIES (test_consistent_P${NPROC}_${TAPER}_${RES} PROPERTIES DEPENDS run_P${NPROC}_${TAPER}_${RES})

            # Check that the mean slip is near our expectations
            IF (NOT ${MEANSLIP} EQUAL 0)
                ADD_TEST(NAME test_slip_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                    COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_slip ${MEANSLIP})
                SET_TESTS_PROPERTIES (test_slip_P${NPROC}_${TAPER}_${RES} PROPERTIES DEPENDS run_P${NPROC}_${TAPER}_${RES})
            ENDIF (NOT ${MEANSLIP} EQUAL 0)

            # Check that the mean interevent time is near our expectations
            IF (NOT ${INTEREVENT} EQUAL 0)
                ADD_TEST(NAME test_interevent_P${NPROC}_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                    COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE} --sweep_file ${SWEEP_FILE} --mean_interevent ${INTEREVENT})
                SET_TESTS_PROPERTIES (test_interevent_P${NPROC}_${TAPER}_${RES} PROPERTIES DEPENDS run_P${NPROC}_${TAPER}_${RES})
            ENDIF(NOT ${INTEREVENT} EQUAL 0) 
        ENDFOREACH(RES ${RESOLUTIONS})
    ENDFOREACH(TAPER_IND RANGE ${NUM_TAPER})
ENDFOREACH(NPROC ${NUM_PROCS})




# --------------  Two Separate Faults  -------------------

SET(RESOLUTIONS_TWO 6000 3000)
SET(METHOD_INTEREVENT_TWO 50 0 0)

FOREACH(TAPER_IND RANGE ${NUM_TAPER})
    LIST(GET TAPER_METHODS ${TAPER_IND} TAPER)
    LIST(GET METHOD_MEANSLIP ${TAPER_IND} MEANSLIP)
    LIST(GET METHOD_INTEREVENT_TWO ${TAPER_IND} INTEREVENT)
    FOREACH(RES ${RESOLUTIONS_TWO})
        SET(TEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/two_fault/${TAPER}/)
        FILE(MAKE_DIRECTORY ${TEST_DIR})
        SET(EVENT_FILE_TWO ${TEST_DIR}events_${RES}.txt)
        SET(SWEEP_FILE_TWO ${TEST_DIR}sweeps_${RES}.txt)


        # Test that the mesher works and generate a parameter file
        ADD_TEST(
            NAME mesh_two_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND mesher
            --import_file=../../../fault_traces/two_fault_1.txt
            --import_file_type=trace --import_trace_element_size=${RES}
            --taper_fault_method=${TAPER}
            --import_file=../../../fault_traces/two_fault_2.txt
            --import_file_type=trace --import_trace_element_size=${RES}
            --taper_fault_method=${TAPER}
            --export_file=two_fault_${RES}.txt
            --export_file_type=text
            --export_file=two_fault_${RES}.kml
            --export_file_type=kml
            --print_statistics=statistics_${RES}.txt
            )
 
        ADD_TEST(NAME param_two_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND ${SETUP_PARAMS_SCRIPT} ${RES} 0.7 two_fault ${VC_EXAMPLE_DIR}/example_params.d)
        SET_TESTS_PROPERTIES (param_two_${TAPER}_${RES} PROPERTIES DEPENDS mesh_two_${TAPER}_${RES})

        # Test that the simulation runs with the resulting mesh
        ADD_TEST(NAME run_two_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND vc params_${RES}.d)
        SET_TESTS_PROPERTIES (run_two_${TAPER}_${RES} PROPERTIES DEPENDS param_two_${TAPER}_${RES})

        # Check that the results are internally self consistent
        ADD_TEST(NAME test_two_consistent_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
            COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE_TWO} --sweep_file ${SWEEP_FILE_TWO} --check_consistent)
        SET_TESTS_PROPERTIES (test_two_consistent_${TAPER}_${RES} PROPERTIES DEPENDS run_two_${TAPER}_${RES})

        # Check that the mean slip is near our expectations
        #    These two faults have same slip rate as single fault, so using
        #       same mean slip values
        IF (NOT ${MEANSLIP} EQUAL 0)
            ADD_TEST(NAME test_two_slip_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
                COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE_TWO} --sweep_file ${SWEEP_FILE_TWO} --mean_slip ${MEANSLIP})
            SET_TESTS_PROPERTIES (test_two_slip_${TAPER}_${RES} PROPERTIES DEPENDS run_two_${TAPER}_${RES})
        ENDIF (NOT ${MEANSLIP} EQUAL 0)

        # Check that the mean interevent time is near our expectations
        # TODO: These tests fail, Kasey needs to figure out why
        #IF (NOT ${INTEREVENT} EQUAL 0)
        #    ADD_TEST(NAME test_two_interevent_${TAPER}_${RES} WORKING_DIRECTORY ${TEST_DIR}
        #        COMMAND ${CHECK_SCRIPT} --event_file ${EVENT_FILE_TWO} --sweep_file ${SWEEP_FILE_TWO} --mean_interevent ${INTEREVENT})
        #    SET_TESTS_PROPERTIES (test_two_interevent_${TAPER}_${RES} PROPERTIES DEPENDS run_two_${TAPER}_${RES})
        #ENDIF(NOT ${INTEREVENT} EQUAL 0)             
    ENDFOREACH(RES ${RESOLUTIONS_TWO})
ENDFOREACH(TAPER_IND RANGE ${NUM_TAPER})







