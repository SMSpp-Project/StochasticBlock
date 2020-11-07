# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ StochasticBlock                                #
#                                                                             #
#                              Niccolo' Iardella                              #
#                          Operations Research Group                          #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.11)
project(StochasticBlock VERSION 0.3.0 LANGUAGES CXX)

# ----- Requirements -------------------------------------------------------- #
# This is to use the find modules we provide
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(SMSpp REQUIRED)

list(REMOVE_AT CMAKE_MODULE_PATH -1)

# ----- Library ------------------------------------------------------------- #
add_library(StochasticBlock)
target_compile_features(StochasticBlock PUBLIC cxx_std_17)

target_sources(StochasticBlock PRIVATE
               src/StochasticBlock.cpp)

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake), that must be different
# for local build and installation.
target_include_directories(
        StochasticBlock PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/SMS++>)

target_link_libraries(StochasticBlock PUBLIC
                      SMS++::SMSpp)

# Alias for testing
add_library(SMS++::StochasticBlock ALIAS StochasticBlock)

# ----- Install instructions ------------------------------------------------ #
# For installing the library and all its CMake modules on the system.
# Not required for local builds (see below).

include(GNUInstallDirs)

# Install library
install(TARGETS StochasticBlock
        EXPORT StochasticBlock-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SMS++)

# Install target import file
install(EXPORT StochasticBlock-targets
        FILE StochasticBlockTargets.cmake
        NAMESPACE SMS++::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/StochasticBlock)

# Generate package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/StochasticBlockConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Generate package configuration file
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/StochasticBlockConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/StochasticBlockConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/StochasticBlock)

# Install package version and configuration files
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/StochasticBlockConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/StochasticBlockConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/StochasticBlock)

# ----- Export from the build tree ------------------------------------------ #
# This is needed so other projects can find the library in a local build.
# It does not affect installation.

if (CMAKE_EXPORT_PACKAGE_REGISTRY)

    export(EXPORT StochasticBlock-targets
           FILE ${CMAKE_CURRENT_BINARY_DIR}/StochasticBlockTargets.cmake
           NAMESPACE SMS++::)

    # Register package in user's package registry (~/.cmake/registry)
    export(PACKAGE StochasticBlock)
endif ()

if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()
