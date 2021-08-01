# Version $Id:$

cmake_minimum_required(VERSION 3.3)

option( ${PROJECT_NAME}_ENABLE_BOOST "Find Boost and, if successful, enable use in ${PROJECT_NAME}" true )
option( ${PROJECT_NAME}_ENABLE_COMPILER_INTRINSICS "Find compiler intrinsics headers and, if successful, enable use in ${PROJECT_NAME}" true )
option( ${PROJECT_NAME}_ENABLE_CUDA "Find Cuda and, if successful, enable use in ${PROJECT_NAME}" false )
