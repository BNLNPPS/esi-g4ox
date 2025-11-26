cmake_minimum_required(VERSION 3.20)

project(g4ox LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(Geant4 REQUIRED ui_all vis_all)
find_package(Geant4 REQUIRED)
find_package(eic-opticks REQUIRED)

add_subdirectory(src)

# Install just the executables, like simphox does
install(TARGETS
    consgeo
    simg4ox
    simtox
)
