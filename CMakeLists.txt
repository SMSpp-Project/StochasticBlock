# --------------------------------------------------------------------------- #
#    CMake file for SMS++/StochasticBlock                                     #
# --------------------------------------------------------------------------- #

cmake_minimum_required( VERSION 3.20 )

project( StochasticBlock VERSION 1.0 LANGUAGES CXX )

set( modName StochasticBlock )
set( modNamespace SMS++ )

# --------------------------------------------------------------------------- #
#  Library target                                                              #
# --------------------------------------------------------------------------- #

add_library( ${modName} )
add_library( ${modNamespace}::${modName} ALIAS ${modName} )

set_target_properties( ${modName} PROPERTIES
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_compile_features( ${modName} PUBLIC cxx_std_17 )

if( MSVC AND BUILD_SHARED_LIBS )
    target_link_options( ${modName} INTERFACE
        "/INCLUDE:SMSpp_force_load_${modName}" )
endif()

target_sources( ${modName} PRIVATE
    src/StochasticBlock.cpp
    src/DiscreteScenarioSet.cpp
    src/ScenarioReductionBlock.cpp
)

target_include_directories( ${modName} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include/${modNamespace}>
)

# --------------------------------------------------------------------------- #
#  Dependencies                                                                #
# --------------------------------------------------------------------------- #

# Core SMS++ framework (Block, Solver, Configuration, etc.)
target_link_libraries( ${modName} PUBLIC SMS++::SMS++ )

# Eigen is used by DiscreteScenarioSet for distance calculations
find_package( Eigen3 REQUIRED )
target_link_libraries( ${modName} PUBLIC Eigen3::Eigen )

# netCDF for serialization/deserialization
find_package( netCDF REQUIRED )
target_link_libraries( ${modName} PUBLIC netCDF::netcdf )

# TwoStageStochasticBlock: DiscreteScenarioSet stores a Block* that callers
# provide as TwoStageStochasticBlock; we need the header for the forward
# declaration in set_Block(). Linked as PRIVATE since it is not part of
# our public API.
target_link_libraries( ${modName} PRIVATE SMS++::TwoStageStochasticBlock )

# --------------------------------------------------------------------------- #
#  Install                                                                     #
# --------------------------------------------------------------------------- #

install( TARGETS ${modName}
    EXPORT  ${modNamespace}Targets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install( DIRECTORY include/
    DESTINATION include/${modNamespace}
    FILES_MATCHING PATTERN "*.h"
)