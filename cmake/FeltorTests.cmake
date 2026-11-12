# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md
include(CTest)
include(Catch)

# ============================================
# Test environment
#
# Every TEST_CASE is registered on its own through catch_discover_tests. In an
# MPI build this runs the test on a single rank without mpirun. The variables
# below add whole-binary runs under mpirun and fix the number of OpenMP threads.

set(FELTOR_TEST_MPI_RANKS "" CACHE STRING
    "Semicolon-separated list of rank counts, e.g. \"2;4;6;8\". Each MPI test executable is additionally registered as one test per rank count, launched as a whole through MPIEXEC_EXECUTABLE. Empty: no additional tests.")
set(FELTOR_TEST_OMP_THREADS "" CACHE STRING
    "Value of OMP_NUM_THREADS for all tests (every rank in MPI tests) when CCCL_THRUST_DEVICE_SYSTEM is OMP. Empty: inherit from the environment of ctest. Ignored for other devices.")

foreach(np IN LISTS FELTOR_TEST_MPI_RANKS)
    if( NOT np MATCHES "^[1-9][0-9]*$")
        message(FATAL_ERROR "FELTOR_TEST_MPI_RANKS must contain positive integers, got \"${np}\"")
    endif()
endforeach()
if( NOT FELTOR_TEST_OMP_THREADS STREQUAL "" AND NOT FELTOR_TEST_OMP_THREADS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "FELTOR_TEST_OMP_THREADS must be empty or a positive integer, got \"${FELTOR_TEST_OMP_THREADS}\"")
endif()
if( FELTOR_WITH_MPI AND NOT FELTOR_TEST_MPI_RANKS STREQUAL "")
    message( STATUS "MPI tests run additionally with ${FELTOR_TEST_MPI_RANKS} ranks: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_PREFLAGS}")
endif()

# Number of OpenMP threads per process that tests run with ("" = not set)
if(CCCL_THRUST_DEVICE_SYSTEM STREQUAL "OMP")
    set(FELTOR_TEST_THREADS_ "${FELTOR_TEST_OMP_THREADS}")
else()
    set(FELTOR_TEST_THREADS_ "")
endif()

### @brief Given a src file named test_path (read-only) and a runtime folder create a
### target with name target_name (write-only)
###
### @param test_path SRC file name to be compiled. Executable will be named after its STEM.
###     Target will be named "dg_testfolder_stem" or "dg_stem".
### @param test_folder Folder name of built executable (i.e. exe will live in
###     build/inc/${test_folder}/${STEM}). If empty string "", then exe will live
###     in build/inc/dg/${STEM}, if "topology" then location is build/inc/dg/topology/${STEM},
###     and analogous "backend".
### @param with_MPI If on, then executable is compiled with MPI
### @param target_name (write-only) Contains target name (dg_testfolder_stem or dg_stem) on output
###     This allows to link libraries on target_name in calling scope
### @param RESOURCE_LOCK name (optional keyword argument) ctest never runs two
###     tests holding the same lock at the same time. Use it for tests that write
###     files with fixed names into the shared working directory, so that
###     ctest -j cannot make them overwrite each other's files.
function(add_dg_test test_path test_folder with_MPI target_name)
    cmake_parse_arguments(PARSE_ARGV 4 arg "" "RESOURCE_LOCK" "")
    if( DEFINED arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "add_dg_test: unknown arguments ${arg_UNPARSED_ARGUMENTS}")
    endif()
    if(CCCL_THRUST_DEVICE_SYSTEM STREQUAL "CUDA" OR CCCL_THRUST_DEVICE_SYSTEM STREQUAL "")
        set_source_files_properties("${test_path}" PROPERTIES LANGUAGE CUDA)
    endif()
    # Get a unique name for the test
    cmake_path(GET test_path STEM test_stem)
    if( test_folder STREQUAL "")
        set(test_alias "feltor::dg::${test_stem}")
        set(test_name "dg_${test_stem}")
    else()
        set(test_alias "feltor::dg::${test_folder}::${test_stem}")
        set(test_name "dg_${test_folder}_${test_stem}")
    endif()
    # Get custom main function
    if( with_MPI)
        set(test_main "${CMAKE_SOURCE_DIR}/tests/catch-tests-mpi.cpp")
    else()
        set(test_main "${CMAKE_SOURCE_DIR}/tests/catch-tests.cpp")
    endif()
    # Create the test executable and link dependencies
    add_executable(${test_name} "${test_main}" "${test_path}")
    add_executable(${test_alias} ALIAS ${test_name})
    # Change name of executable
    set_target_properties( ${test_name} PROPERTIES OUTPUT_NAME ${test_stem})
    # Create Alias name feltor::dg::...
    target_link_libraries(${test_name} PRIVATE dg_dg)
    target_link_libraries(${test_name} PRIVATE Catch2::Catch2)
    if( with_MPI)
        target_link_libraries(${test_name} PRIVATE MPI::MPI_CXX)
        target_compile_definitions(${test_name} PRIVATE WITH_MPI)
    endif()
    # Add as dependency to the dg_tests target
    add_dependencies(dg_tests ${test_name})
    # Make sure tests mirror folder structure in ./build/
    if( test_folder STREQUAL "")
        set_target_properties(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/inc/dg"
        )
    elseif( test_folder STREQUAL "topology")
        set_target_properties(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/inc/dg/topology"
        )
    elseif( test_folder STREQUAL "backend")
        set_target_properties(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/inc/dg/backend"
        )
    else() # matrix, file, geometries
        set_target_properties(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/inc/${test_folder}"
        )
    endif()
    # Register the tests with CTest: each TEST_CASE on its own (one rank if MPI)
    set(test_properties "")
    if( NOT FELTOR_TEST_THREADS_ STREQUAL "")
        list(APPEND test_properties
            ENVIRONMENT "OMP_NUM_THREADS=${FELTOR_TEST_THREADS_}"
            PROCESSORS ${FELTOR_TEST_THREADS_}
        )
    endif()
    if( DEFINED arg_RESOURCE_LOCK)
        list(APPEND test_properties RESOURCE_LOCK ${arg_RESOURCE_LOCK})
    endif()
    if( test_properties STREQUAL "")
        catch_discover_tests(${test_name})
    else()
        catch_discover_tests(${test_name} PROPERTIES ${test_properties})
    endif()
    # Register whole-binary runs for each rank count in FELTOR_TEST_MPI_RANKS
    if( with_MPI)
        if( FELTOR_TEST_THREADS_ STREQUAL "")
            set(num_threads 1)
        else()
            set(num_threads ${FELTOR_TEST_THREADS_})
        endif()
        foreach(np IN LISTS FELTOR_TEST_MPI_RANKS)
            add_test(NAME ${test_name}_np${np}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
                    ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
            )
            math(EXPR num_procs "${np}*${num_threads}")
            # PROCESSORS lets ctest -j avoid oversubscribing the machine
            set_tests_properties(${test_name}_np${np} PROPERTIES
                PROCESSORS ${num_procs}
                LABELS "mpi;np${np}"
            )
            if( NOT FELTOR_TEST_THREADS_ STREQUAL "")
                set_tests_properties(${test_name}_np${np} PROPERTIES
                    ENVIRONMENT "OMP_NUM_THREADS=${FELTOR_TEST_THREADS_}"
                )
            endif()
            set(locks "")
            if(CCCL_THRUST_DEVICE_SYSTEM STREQUAL "CUDA" OR CCCL_THRUST_DEVICE_SYSTEM STREQUAL "")
                # All ranks of a node share the GPU(s); run one multi-rank GPU test at a time
                list(APPEND locks gpu)
            endif()
            if( DEFINED arg_RESOURCE_LOCK)
                list(APPEND locks ${arg_RESOURCE_LOCK})
            endif()
            if( NOT locks STREQUAL "")
                set_tests_properties(${test_name}_np${np} PROPERTIES RESOURCE_LOCK "${locks}")
            endif()
        endforeach()
    endif()
    # Set function output
    set( target_name ${test_name} PARENT_SCOPE)
endfunction()
