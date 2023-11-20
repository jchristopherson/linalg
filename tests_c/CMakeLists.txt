include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# C API tests
if (${BUILD_C_API})
    message(STATUS "Building LINALG C API tests.")

    include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/include
        ${ferror_INCLUDE_DIR}
    )

    set(linalg_c_test_sources
        c_linalg_test.c
        c_test_core.c
        c_linalg_test_misc.c
        c_linalg_test_factor.c
        c_linalg_test_eigen.c
    )

    # Build the C API tests
    add_executable(linalg_c_test ${linalg_c_test_sources})
    link_library(linalg_c_test linalg ${PROJECT_INCLUDE_DIR})
    add_test(
        NAME linalg_c_test
        WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
        COMMAND $<TARGET_FILE:linalg_c_test>
    )
endif()