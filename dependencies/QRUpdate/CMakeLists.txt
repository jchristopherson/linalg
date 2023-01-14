# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# Use an installed version, if one is found; else, import the library
if (${qrupdate_FOUND})
    # Inform the user of what's going on
    message(STATUS "FERROR (${qrupdate_VERSION}) library found.")
else()
    # Inform the user of what's going on
    message(STATUS "QRUpdate library not found.  Downloading appropriate repository.")
    FetchContent_Declare(
        qrupdate
        GIT_REPOSITORY "https://github.com/mpimd-csc/qrupdate-ng"
        OVERRIDE_FIND_PACKAGE
    )

    FetchContent_MakeAvailable(qrupdate)

    if (WIN32)
        if (BUILD_SHARED_LIBS)
            add_custom_command(
                TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:qrupdate>
                $<TARGET_FILE_DIR:${PROJECT_NAME}
            )
        endif()
    endif()
endif()


# Make a parent-scope variable for the library
set(qrupdate_LIBRARY qrupdate)
set(qrupdate_LIBRARY ${qrupdate_LIBRARY} PARENT_SCOPE)