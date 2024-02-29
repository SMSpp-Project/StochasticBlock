# --------------------------------------------------------------------------- #
#    CMake file for StochasticBlock/tests                                     #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#    Note: this file tries to respect the principles of Modern CMake,         #
#          keep it in mind when editing it.                                   #
#                                                                             #
#                              Niccolo' Iardella                              #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

add_executable(StochasticBlock_test test.cpp)
target_link_libraries(StochasticBlock_test PRIVATE SMS++::StochasticBlock)

add_test(NAME StochasticBlock_test
         COMMAND StochasticBlock_test)

# --------------------------------------------------------------------------- #
