cmake_minimum_required( VERSION 2.8.3 )
project( kinematics_toolbox )

## Find catkin macros and libraries
find_package( catkin REQUIRED 
  geometry_msgs 
)

## System dependencies are found with CMake's conventions
find_package( Eigen3 REQUIRED )

###################################
## catkin specific configuration ##
###################################

catkin_package(
  INCLUDE_DIRS    include
                  ${EIGEN3_INCLUDE_DIR} # TODO: why is this needed here with EIGEN3 as a depend later?
  LIBRARIES       ${PROJECT_NAME}
                  ${PROJECT_NAME}_ros_types
  CATKIN_DEPENDS  geometry_msgs
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

add_library( ${PROJECT_NAME}_ros_types
  src/kinematics_ros_types.cpp
)
target_link_libraries( ${PROJECT_NAME}_ros_types
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
