# ####################################################################
#  CMakeLists.txt
#  graphClustering
#
#  Created by 曹 巍 on 12/06/12.
#  Copyright (c) 2012年 生物情報工学研究室・東京大学. All rights reserved.
# ####################################################################

cmake_minimum_required(VERSION 2.6.4)

#since VERSION 2.8.0, CMP0014 -"Input directories must have CMakeLists.txt."
#cmake_policy(SET CMP0014 OLD) Not work

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.6)
   cmake_policy(VERSION 2.8)
endif("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 2.6)

#################################
# Compilers: 
# 
#################################
#set(CMAKE_C_COMPILER "gcc-4.2")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++-4.2")

project(gCluster)
####################
# System name
####################
#message( STATUS "System name:" ${CMAKE_SYSTEM_NAME})
#message( STATUS "Host system name:" ${CMAKE_HOST_SYSTEM_NAME})
#####################
# project version
######################
SET( ${PROJECT_NAME}_MAJOR_VERSION 0 ) 
SET( ${PROJECT_NAME}_MINOR_VERSION 1 ) 
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )

################################################
# Options:
################################################
option(USE_STATIC_BOOST "Build with static BOOST libraries instead of dynamic" OFF)

################################################
# Default build type: Release
# ----------------------------------------------
# Set a default build type for single-configuration
# CMake generators if no build type is set.
################################################
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
##############################################
# C++11 support
##############################################
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
   set(CPLUSPLUS_11 "-std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   set(CPLUSPLUS_11 "-std=c++0x")
else()
   set(CPLUSPLUS_11 "")
endif()
###############################################
# Compilation schemes: -DCMAKE_BUILD_TYPE=
# 
#  1  None          :(CMAKE_C_FLAGS or CMAKE_CXX_FLAGS used)
#  2. Debug         :(CMAKE_C_FLAGS_DEBUG or CMAKE_CXX_FLAGS_DEBUG)
#  3. Release       :(CMAKE_C_FLAGS_RELEASE or CMAKE_CXX_FLAGS_RELEASE)
#  4. RelWithDebInfo:(CMAKE_C_FLAGS_RELWITHDEBINFO or CMAKE_CXX_FLAGS_RELWITHDEBINFO)
#  5. Profile
#  6. MinSizeRel    :(CMAKE_C_FLAGS_MINSIZEREL or CMAKE_CXX_FLAGS_MINSIZEREL)  
###############################################
#----------------------------------------------
# Common
#----------------------------------------------
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} ${CPLUSPLUS_11} -Wall -Wno-unknown-pragmas -Wno-long-long -Wno-variadic-macros")

#----------------------------------------------
# 1. Debug scheme: C/C++
#----------------------------------------------
set(CMAKE_CXX_FLAGS_DEBUG     "-DNDEBUG -g3 -O0 -fno-inline")
set(CMAKE_CXX_FLAGS_DEBUGFULL "-g3 -O0-fno-inline")
set(CMAKE_C_FLAGS_DEBUG       "-DNDEBUG -g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
set(CMAKE_C_FLAGS_DEBUGFULL   "-g3 -fno-inline")

#-----------------------------------------------
# 2. Release scheme: C/C++
#-----------------------------------------------
#SET(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -DGOOGLE_STRIP_LOG=1 -Os -v -Xlinker -v")
#SET(CMAKE_CXX_FLAGS_RELEASE "-Os -v -Xlinker -v")
set(CMAKE_CXX_FLAGS_RELEASE "-Os")
set(CMAKE_C_FLAGS_RELEASE   "-Os")

#-----------------------------------------------
# 3. Release_With_Debug_Info scheme: C/C++
#-----------------------------------------------
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")

#-----------------------------------------------
# 4. Profile scheme: C/C++
#-----------------------------------------------

set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")

#-----------------------------------------------
# General
#-----------------------------------------------
set(CMAKE_CXX_FLAGS_ARRAY_CHECK    "-g3 -fno-inline -ftest-coverage -fprofile-arcs -fstack-protector")
set(CMAKE_C_FLAGS_ARRAY_CHECK      "-g3 -fno-inline -ftest-coverage -fprofile-arcs -fstack-protector")

######################################################
# set link option
######################################################
#set(BUILD_SHARED_LIBS FALSE)
if(USE_STATIC_BOOST)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIE -fPIC -static")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -pie -v -Xlinker -v")
else(USE_STATIC_BOOST)
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif(USE_STATIC_BOOST)

######################################################
# set default install prefix to project root directory
######################################################
IF( CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" ) 
	SET( CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" )
ENDIF()

#######################################################################
#  set boost and loki
#  Note: set the following variable before using find_package(Boost) 
#######################################################################
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    #SET(LOKI_INCLUDE_DIRS ${LOKI_INCLUDE_DIRS} "/Users/davecao/Downloads/loki/include")
    SET(BOOST_ROOT ${BOOST_ROOT} "/Users/davecao/Downloads/boost")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    #SET(LOKI_INCLUDE_DIRS ${LOKI_INCLUDE_DIRS} "/home/davecao/numerical_libs/loki/include")
    SET(BOOST_ROOT ${BOOST_ROOT} "/home/davecao/libs")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
message(STATUS "Boost:" ${BOOST_ROOT})
###################################
# Using local cmake modules: FindOpenCL
###################################
#set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)
######################
# find_package: 
######################
unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIRS CACHE)

if(USE_STATIC_BOOST)
   unset(Boost_LIBRARIES)
   message(WARNING "Linking against boost static libraries")
   set(Boost_USE_STATIC_LIBS ON)
   set(Boost_USE_MULTITHREADED ON)
	 find_package(Boost COMPONENTS program_options graph regex system REQUIRED)
else(USE_STATIC_BOOST)
   unset(Boost_LIBRARIES)
   message(WARNING "Linking against boost dynamic libraries")
   set(Boost_USE_STATIC_LIBS OFF)
endif(USE_STATIC_BOOST)
find_package(Boost COMPONENTS program_options graph regex system REQUIRED)
#find_package(GLUT REQUIRED)
#find_package(OpenGL REQUIRED)
#find_package(OpenCL REQUIRED)


#########################################################
# set binary output directory
##########################################################

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin" 
	CACHE PATH "Non-configurable executalbe directory")

##########################################################
# set libs output directory
##########################################################
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib" 
	CACHE PATH "Non-configurable executalbe directory") 
 
#################################################################
# Include directories:
##########################################################################

INCLUDE_DIRECTORIES(
    ${LOKI_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/bilab/
#    ${PROJECT_SOURCE_DIR}/Algorithm/
#    ${PROJECT_SOURCE_DIR}/BGL/
#    ${PROJECT_SOURCE_DIR}/cmdline/
#    ${PROJECT_SOURCE_DIR}/cmdline/gflags/
#    ${PROJECT_SOURCE_DIR}/console/
#    ${PROJECT_SOURCE_DIR}/FileIO/
#    ${PROJECT_SOURCE_DIR}/memory/
#    ${PROJECT_SOURCE_DIR}/MolecularManager/
#    ${PROJECT_SOURCE_DIR}/Ocl/
#    ${PROJECT_SOURCE_DIR}/StoragePolicy/
#    ${PROJECT_SOURCE_DIR}/utility/
#    ${EXECUTABLE_OUTPUT_PATH}
#    ${OPENCL_INCLUDE_DIRS}
#    ${OPENGL_INCLUDE_DIR}
#    ${GLUT_INCLUDE_DIR}
)
#if( OPENCL_HAS_CPP_BINDINGS )
#    message( STATUS "OpenCL has CPP bindings. Full include is:" ${OPENCL_INCLUDE_DIRS} )
#else( OPENCL_HAS_CPP_BINDINGS )
#    message( STATUS "NO OpenCL CPP bindings found." )
#endif( OPENCL_HAS_CPP_BINDINGS )

#if( GLUT_FOUND )
#    message( STATUS "OpenGL GLUT found. Full include is:" ${GLUT_INCLUDE_DIR} )
#else( GLUT_FOUND )
#    message( STATUS "NO GLUT LIB found." )
#endif( GLUT_FOUND )

#################################################################
# link directories:
##########################################################################
LINK_DIRECTORIES(
    ${Boost_LIBRARY_DIRS}
)

#########################################################
# Add source directory: build tasks
##########################################################
# subdirs is depreciated. Instead, using add_subdirectory
#subdirs(
#  ${PROJECT_SOURCE_DIR}/lib/loki_src
#  ${PROJECT_SOURCE_DIR}/src
#)

#add_subdirectory(${PROJECT_SOURCE_DIR}/cmdline) 
#add_subdirectory(${PROJECT_SOURCE_DIR}/Lokisrc)
#add_subdirectory(${PROJECT_SOURCE_DIR}/Ocl)

FILE(GLOB CC_SRCS TO_CMAKE_PATH ${PROJECT_SOURCE_DIR} [^.]*.cpp [^.]*.c)

################################################
# compile source files.
###############################################
add_executable( ${PROJECT_NAME}
	${CC_SRCS}
)

###############################################
#  Add dependencies
###############################################

#add_dependencies(${PROJECT_NAME} Loki gflags)

################################################
# library link:
################################################
target_link_libraries(${PROJECT_NAME}
#    m
#    Loki
#    gflags
#    OCL
    ${Boost_LIBRARIES}
#    boost_program_options
#		boost_graph
#		boost_regex
#    boost_system
#    ${OPENCL_LIBRARIES}
#    ${OPENGL_LIBRARIES}
#    ${GLUT_LIBRARIES}
)

