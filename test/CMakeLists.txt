# --------------------------------------------------------------------------- #
#    CMake file for StochasticBlock/tests                                     #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

add_executable(StochasticBlock_test test.cpp)
target_link_libraries(StochasticBlock_test PRIVATE SMS++::StochasticBlock)

if (WIN32 AND BUILD_SHARED_LIBS AND NOT DEFINED ENV{CONDA_BUILD})
    # copy runtime DLLs next to the exe after build
    add_custom_command(TARGET StochasticBlock_test POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy_if_different
                       $<TARGET_RUNTIME_DLLS:StochasticBlock_test>
                       $<TARGET_FILE_DIR:StochasticBlock_test>
                       COMMAND_EXPAND_LISTS
                       COMMENT "Copying runtime DLLs next to StochasticBlock_test.exe")
endif()

add_test(NAME StochasticBlock_test
         COMMAND StochasticBlock_test)

# --------------------------------------------------------------------------- #
