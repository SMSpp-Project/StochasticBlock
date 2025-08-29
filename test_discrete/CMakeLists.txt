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
    if (MSVC)
        add_link_options("/OPT:NOREF")
    elseif (UNIX AND (NOT APPLE))
        add_link_options("-Wl,--no-as-needed")
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

# Comprehensive test suite for DiscreteScenarioSet
add_executable(test_discretescenarioset test_discretescenarioset.cpp)
target_link_libraries(test_discretescenarioset PRIVATE SMS++::StochasticBlock)

# Add CapacitatedFacilityLocationBlock dependency if available
if (TARGET ${modNamespace}::CapacitatedFacilityLocationBlock)
    target_link_libraries(test_discretescenarioset PRIVATE ${modNamespace}::CapacitatedFacilityLocationBlock)
    target_compile_definitions(test_discretescenarioset PRIVATE WITH_CAPACITATED_FACILITY_LOCATION)
endif()

# Register the test with CTest
add_test(NAME test_discretescenarioset
        COMMAND test_discretescenarioset)

# --------------------------------------------------------------------------- #