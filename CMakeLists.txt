cmake_minimum_required (VERSION 2.6)

SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0" CACHE STRING "debug" )
#SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=1" CACHE STRING "optimized release" )

IF(NOT CMAKE_BUILD_TYPE)
 SET(CMAKE_BUILD_TYPE Release CACHE STRING
           "Choose the type of build, options are: None Debug Release
           RelWithDebInfo MinSizeRel."
                 FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF(DEFINED CMAKE_BUILD_TYPE)
project (CHROMATIN_SIS_COARSE)

# The version number.
set (CHROMATIN_SIS_COARSE_VERSION_MAJOR 1)
set (CHROMATIN_SIS_COARSE_VERSION_MINOR 0)


#---------------------------------------------
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
"/home/gamze/Fall2013/igh/random_model/boost_1_54_0/")
find_package(Boost COMPONENTS system filesystem REQUIRED )
MESSAGE("Boost information:") 
MESSAGE("  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}") 
MESSAGE("
Boost_LIBRARIES:/home/gamze/Fall2013/igh/random_model/boost_1_54_0/stage/lib") 
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

#---------------------------------------------
if (EIGEN3_INCLUDE_DIR)
  # in cache already
  set(EIGEN3_FOUND TRUE)
else (EIGEN3_INCLUDE_DIR)
  find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
       PATH_SUFFIXES eigen3
       HINTS
       "/home/gamze/Fall2013/igh/random_model/eigen3"
       ${INCLUDE_INSTALL_DIR}
       ${KDE4_INCLUDE_DIR}
     )
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR )
  mark_as_advanced(EIGEN3_INCLUDE_DIR)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_INCLUDE_DIR)
#---------------------------------------------


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
add_subdirectory (include)
add_subdirectory (src bin)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_SOURCE_DIR}/test
    ${CMAKE_BINARY_DIR}/test)


