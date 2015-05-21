cmake_minimum_required( VERSION 2.8.3 )
project( kinematics_toolbox )

## Find catkin macros and libraries
find_package( catkin REQUIRED )

## System dependencies are found with CMake's conventions
find_package( Eigen3 REQUIRED )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS    include
  LIBRARIES       ${PROJECT_NAME}
  DEPENDS         EIGEN3
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories( include )
include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${catkin_INCLUDE_DIRS}
)

add_definitions( "-Wall -Wextra -std=c++11 -O2" )

#######################
## Build - Libraries ##
#######################

add_library( ${PROJECT_NAME}
  src/kinematics.cpp
)
target_link_libraries( ${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

#########################
## Build - Executables ##
#########################


#############
## Install ##
#############


#############
## Testing ##
#############
