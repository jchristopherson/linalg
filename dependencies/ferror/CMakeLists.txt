# Get the macros and functions we'll need
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")
include(FetchContent)

# If found, use the installed version; else, import the library
if (${ferror_FOUND})
    # Inform the user of what's going on
    message(STATUS "FERROR (${ferror_VERSION}) library found.")

    # Get the mod file location
    get_target_property(ferror_INCLUDE_DIR ferror INTERFACE_INCLUDE_DIRECTORIES)
else()
    # Inform the user of what's going on
    message(STATUS "FERROR library not found.  Downloading appropriate repository.")
    
    # Fetch the proper content
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )

    FetchContent_MakeAvailable(ferror)

    if (WIN32)
        if (BUILD_SHARED_LIBS)
            add_custom_command(
                TARGET ${PROJECT_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:ferror>
                $<TARGET_FILE_DIR:${PROJECT_NAME}
            )
        endif()
    endif()
    
    set(ferror_INCLUDE_DIR ${ferror_BINARY_DIR}/include)
endif()

# Make a parent-scope variable for the library
set(ferror_LIBRARY ferror)
set(ferror_LIBRARY ${ferror_LIBRARY} PARENT_SCOPE)

# Make a parent-scope variable locating the include directory for ferror
set(ferror_INCLUDE_DIR ${ferror_INCLUDE_DIR} PARENT_SCOPE)