set(G4OX_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

# ----------------------------------------------------------------------
# consgeo
# ----------------------------------------------------------------------
add_executable(consgeo
    consgeo.cpp
    ../include/argparse/argparse.hpp
)

target_link_libraries(consgeo
    PRIVATE
        Geant4::G4persistency      
        eic-opticks::G4CX
        eic-opticks::SysRap
)

target_include_directories(consgeo
    PRIVATE
        ${G4OX_INCLUDE_DIR}
)

# ----------------------------------------------------------------------
# simg4ox
# ----------------------------------------------------------------------
add_executable(simg4ox
    simg4ox.cpp
    g4app.h
    ../include/argparse/argparse.hpp
)

target_link_libraries(simg4ox
    PRIVATE
        Geant4::G4persistency
        eic-opticks::G4CX
        eic-opticks::SysRap
)

target_include_directories(simg4ox
    PRIVATE
        ${G4OX_INCLUDE_DIR}
)

# ----------------------------------------------------------------------
# simtox
# ----------------------------------------------------------------------
add_executable(simtox
    simtox.cpp
)

# This one only needs SysRap
target_link_libraries(simtox
    PRIVATE
        eic-opticks::SysRap
)
