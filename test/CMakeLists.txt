# --------------------------------------------------------------------------- #
#    CMake file for unit tests                                                #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# FIXME: This test cannot be compiled by GCC
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_executable(StochasticBlock_test test.cpp)
    target_link_libraries(StochasticBlock_test PRIVATE SMS++::StochasticBlock)
endif ()
