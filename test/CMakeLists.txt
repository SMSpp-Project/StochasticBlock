# --------------------------------------------------------------------------- #
#    CMake file for StochasticBlock/test                                      #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# ----- StochasticBlock_test ------------------------------------------------ #

add_executable(${modName}_test test.cpp)
target_link_libraries(${modName}_test PRIVATE SMS++::StochasticBlock)

smspp_copy_runtime(${modName}_test)

add_test(NAME ${modName}_test
         COMMAND ${modName}_test)

# ----- test_discrete ------------------------------------------------------- #

add_executable(test_discrete test_discrete.cpp)
target_link_libraries(test_discrete PRIVATE SMS++::StochasticBlock)

smspp_copy_runtime(test_discrete)

add_test(NAME test_discrete
         COMMAND test_discrete)

# ----- MultiStageDiscreteScenarioSet_unit_test ----------------------------- #

# The view API of the scenario tree, clone() and the serialize() /
# deserialize() round trip of MultiStageDiscreteScenarioSet.
set(_t MultiStageDiscreteScenarioSet_unit_test)
add_executable(${_t} test_mssds.cpp)
target_link_libraries(${_t} PRIVATE SMS++::StochasticBlock)

smspp_copy_runtime(${_t})

# run in the build tree, since it writes (and removes) a netCDF file
add_test(NAME ${_t}
         COMMAND ${_t}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ----- IndependentMultiStageScenarioGenerator_unit_test -------------------- #

# The view API of IndependentMultiStageScenarioGenerator, whose stages do not
# depend on the history.
set(_t IndependentMultiStageScenarioGenerator_unit_test)
add_executable(${_t} test_imssg.cpp)
target_link_libraries(${_t} PRIVATE SMS++::StochasticBlock)

smspp_copy_runtime(${_t})

# run in the build tree, since it writes (and removes) a netCDF file
add_test(NAME ${_t}
         COMMAND ${_t}
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ----- the label of this module's tests ------------------------------------ #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}")
endif ()

# --------------------------------------------------------------------------- #
