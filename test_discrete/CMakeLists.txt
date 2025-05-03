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
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# Comprehensive test suite for DiscreteScenarioSet
add_executable(test_discretescenarioset test_discretescenarioset.cpp)
target_link_libraries(test_discretescenarioset PRIVATE SMS++::StochasticBlock)

# Register the test with CTest
add_test(NAME test_discretescenarioset
        COMMAND test_discretescenarioset)

# --------------------------------------------------------------------------- #
