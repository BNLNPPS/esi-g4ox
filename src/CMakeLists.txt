set(source_files
    csg_tree.cpp
)

add_library(g4ox ${source_files})
target_link_libraries(g4ox PUBLIC Geant4::G4persistency Opticks::G4CX)

target_include_directories(g4ox
    INTERFACE
    $<INSTALL_INTERFACE:include>
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties(g4ox PROPERTIES
    VERSION ${PROJECT_VERSION}
    POSITION_INDEPENDENT_CODE ON)

install(TARGETS g4ox EXPORT g4oxTargets
    RUNTIME   DESTINATION ${G4OX_RUNTIME_INSTALL_DIR}
    LIBRARY   DESTINATION ${G4OX_LIBRARY_INSTALL_DIR}
    ARCHIVE   DESTINATION ${G4OX_ARCHIVE_INSTALL_DIR}
    FRAMEWORK DESTINATION ${G4OX_FRAMEWORK_INSTALL_DIR})

add_executable(make_csg_tree make_csg_tree.cpp)
target_link_libraries(make_csg_tree g4ox)
