cmake_minimum_required(VERSION 3.24)

project(singler_classic_markers
    VERSION 0.1.0
    DESCRIPTION "Classic marker detection in SingleR"
    LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Library
add_library(singler_classic_markers INTERFACE)
add_library(singler::singler_classic_markers ALIAS singler_classic_markers)

target_include_directories(singler_classic_markers INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/singler_classic_markers>)
target_compile_features(singler_classic_markers INTERFACE cxx_std_17)

# Dependencies
option(SINGLEPP_FETCH_EXTERN "Automatically fetch singler_classic_markers's external dependencies." ON)
if(SINGLEPP_FETCH_EXTERN)
    add_subdirectory(extern)
else()
    find_package(tatami_tatami 4.0.0 CONFIG REQUIRED)
    find_package(tatami_tatami_stats 3.0.0 CONFIG REQUIRED)
    find_package(ltla_sanisizer 0.1.5 CONFIG REQUIRED)
    find_package(libscran_topicks 0.2.2 CONFIG REQUIRED)
endif()

target_link_libraries(singler_classic_markers INTERFACE tatami::tatami tatami::tatami_stats ltla::sanisizer libscran::topicks)

# Tests
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(SINGLEPP_TESTS "Build singler_classic_markers's test suite." ON)
else()
    option(SINGLEPP_TESTS "Build singler_classic_markers's test suite." OFF)
endif()

if(SINGLEPP_TESTS)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif() 
endif()

# Install
install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/singler_classic_markers)

install(TARGETS singler_classic_markers
    EXPORT singler_classic_markersTargets)

install(EXPORT singler_classic_markersTargets
    FILE singler_singler_classic_markersTargets.cmake
    NAMESPACE singler::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singler_singler_classic_markers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/singler_singler_classic_markersConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singler_singler_classic_markers)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/singler_singler_classic_markersConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/singler_singler_classic_markersConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/singler_singler_classic_markersConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/singler_singler_classic_markers)
