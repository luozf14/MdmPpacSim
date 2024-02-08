# ----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(MdmPpacSim)

enable_language(Fortran)
SET(CMAKE_Fortran_FLAGS "-fPIC -finit-local-zero -falign-commons -fno-automatic -ffixed-line-length-none -w")

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -Wall -g2 -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -Wall")

# ----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()

# Find ROOT
find_package(ROOT CONFIG REQUIRED RIO)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
include_directories(${ROOT_INCLUDE_DIRS})
# message(${ROOT_INCLUDE_DIRS})
# message(${ROOT_LIBRARIES})

# ----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${Geant4_INCLUDE_DIR})

# ----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc src/MDMTrace.cpp src/RAYTKIN1.F)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh include/MDMTrace.h)

# add_library(MDMTrace SHARED src/MDMTrace.cpp src/RAYTKIN1.F include/MDMTrace.h)
# target_link_libraries(MDMTrace c gfortran)

# ----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(MdmPpacSim MdmPpacSim.cc ${sources} ${add_sources} ${headers} ${add_headers})
target_link_libraries(MdmPpacSim ${Geant4_LIBRARIES} ${HBOOK_LIBRARIES} ${ROOT_LIBRARIES})

# ----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(MDMPPACSIM_SCRIPTS
    init_vis.mac
    run.mac
    vis.mac
    plotHisto.C
    rayin.dat
    waveFunction.txt
)

foreach(_script ${MDMPPACSIM_SCRIPTS})
    configure_file(
        ${PROJECT_SOURCE_DIR}/${_script}
        ${PROJECT_BINARY_DIR}/${_script}
        COPYONLY
    )
endforeach()
