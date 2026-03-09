# --------------------------------------------------------------------------- #
#    CMake file for StochasticBlock/test_discrete                             #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Since we are using the block factory, objects from linked libraries
# may wrongly appear unused, and by default the linker does not include them,
# so we have to force the linking.
if (BUILD_SHARED_LIBS)
    if (UNIX)
        if (APPLE)
            add_link_options("-Wl,-dead_strip_dylibs")
        else ()
            add_link_options("-Wl,--no-as-needed")
        endif ()
    endif ()
else ()
    if (MSVC)
        add_link_options("/WHOLEARCHIVE" "/FORCE:MULTIPLE")
    else () # Unix
        if (APPLE)
            add_link_options("-Wl,-all_load")
        else ()
            add_link_options("-Wl,--whole-archive,--allow-multiple-definition")
        endif ()
    endif ()
endif ()

# ----- StochasticBlock_test_discrete --------------------------------------- #
add_executable(test_discrete test_discrete.cpp)
target_link_libraries(test_discrete PRIVATE SMS++::StochasticBlock)

if (TARGET SMS++::MILPSolver)
    target_link_libraries(test_discrete PRIVATE SMS++::MILPSolver)
endif ()

if (WIN32 AND BUILD_SHARED_LIBS)
    # copy runtime DLLs next to the exe after build
    add_custom_command(TARGET test_discrete POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:test_discrete>
            $<TARGET_FILE_DIR:test_discrete>
            COMMAND_EXPAND_LISTS
            COMMENT "Copying runtime DLLs next to test_discrete.exe"
    )
endif()

add_test(NAME test_discrete
         COMMAND test_discrete)

# --------------------------------------------------------------------------- #
