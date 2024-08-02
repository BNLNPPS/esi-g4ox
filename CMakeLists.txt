cmake_minimum_required(VERSION 3.20)

project(g4ox VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set a debug postfix
set(CMAKE_DEBUG_POSTFIX "-dbg")

# Set installation destinations typical for *nix systems
include(GNUInstallDirs)

set(G4OX_INCLUDE_INSTALL_DIR        "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")
set(G4OX_RUNTIME_INSTALL_DIR        "${CMAKE_INSTALL_BINDIR}")
set(G4OX_LIBRARY_INSTALL_DIR        "${CMAKE_INSTALL_LIBDIR}")
set(G4OX_ARCHIVE_INSTALL_DIR        "${CMAKE_INSTALL_LIBDIR}")
set(G4OX_FRAMEWORK_INSTALL_DIR      "${CMAKE_INSTALL_LIBDIR}")
set(G4OX_DATA_INSTALL_DIR           "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")
set(G4OX_DATA_INSTALL_DIR_FULL      "${CMAKE_INSTALL_FULL_DATADIR}/${PROJECT_NAME}")
set(G4OX_CMAKE_CONFIG_INSTALL_DIR   "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
set(G4OX_CMAKE_CONFIG_TEMPLATE      "cmake/Config.cmake.in")
set(G4OX_CMAKE_CONFIG_DIR           "${CMAKE_CURRENT_BINARY_DIR}")
set(G4OX_CMAKE_VERSION_CONFIG_FILE  "${G4OX_CMAKE_CONFIG_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(G4OX_CMAKE_PROJECT_CONFIG_FILE  "${G4OX_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Config.cmake")

list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")

# Set up external dependencies
find_package(Geant4 REQUIRED)
find_package(g4cx REQUIRED)

add_subdirectory(src)

# Install files and export configs
include(CMakePackageConfigHelpers)

# Generate the config file that is includes the exports
configure_package_config_file(${G4OX_CMAKE_CONFIG_TEMPLATE}
    ${G4OX_CMAKE_PROJECT_CONFIG_FILE}
    INSTALL_DESTINATION ${G4OX_CMAKE_CONFIG_INSTALL_DIR}
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

# Create and install version file
write_basic_package_version_file(${G4OX_CMAKE_VERSION_CONFIG_FILE}
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

install(
    FILES ${G4OX_CMAKE_PROJECT_CONFIG_FILE} ${G4OX_CMAKE_VERSION_CONFIG_FILE}
    DESTINATION ${G4OX_CMAKE_CONFIG_INSTALL_DIR}
)

install(
    EXPORT g4oxTargets
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${G4OX_CMAKE_CONFIG_INSTALL_DIR}
)

# Set up tests
#include(CTest)
#
#if(BUILD_TESTING)
#    enable_testing()
#    add_subdirectory(test)
#endif()
