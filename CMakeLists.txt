#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
set(name exampleB1)
project(${name} VERSION 0.1.0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#



option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()



### Opticks
message(STATUS "WITH_OPTICKS is set ON")

list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
set(OPTICKS_PREFIX "$ENV{OPTICKS_PREFIX}")
### we need to Make sure this OPTICKS_PREFIX matches with the one on OpticksBuildOptions
message(STATUS "---- OPTICKS ----")
include(OpticksBuildOptions NO_POLICY_SCOPE)
message(STATUS  "$ENV{OPTICKS_HOME}/cmake/Modules")
message(STATUS "CMAKE_MODULE_PATH : ${CMAKE_MODULE_PATH}")
message(STATUS "OPTICKS_HOME : $ENV{OPTICKS_HOME}")
message(STATUS "OPTICKS_PREFIX : $ENV{OPTICKS_PREFIX} ")

message(STATUS "---- OPTICKS END ----")
## These are needed for opticks
find_package(G4CX REQUIRED CONFIG)
find_package(U4       REQUIRED CONFIG)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${name} exampleB1.cc ${sources} ${headers} )
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories( ${CMAKE_MODULE_PATH})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(Opticks::G4CX Opticks::U4 )
#target_link_libraries(exampleB1 ${Geant4_LIBRARIES}  Opticks::G4CX Opticks::U4 Opticks::ExtG4 Opticks::CSG_GGeo )
target_link_libraries(${name} ${Geant4_LIBRARIES}  Opticks::G4CX Opticks::U4 )
#target_include_directories( ${name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

#target_compile_definitions( ${name} PUBLIC OPTICKS_G4CX )

#target_link_libraries(exampleB1 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  exampleB1.in
  exampleB1.out
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  tsg_offscreen.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()





#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(B1 DEPENDS ${name})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${name} DESTINATION bin)
