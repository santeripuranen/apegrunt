# Version $Id:$

cmake_minimum_required(VERSION 3.1)

#option( ${CMAKE_PROJECT_NAME}_ENABLE_CUDA "Find CUDA and, if successful, enable use in SuperDCA" true )

##############
## CUDA setup
###

set( ${CMAKE_PROJECT_NAME}_NO_CUDA true CACHE INTERNAL "Don't use CUDA, if true" ) # Initialize with default value 

if( ${CMAKE_PROJECT_NAME}_ENABLE_CUDA )
	setup_message( "check for CUDA" )
	find_package( CUDA QUIET )
	if( CUDA_FOUND )
		set( ${CMAKE_PROJECT_NAME}_NO_CUDA false CACHE INTERNAL "Don't use CUDA, if true" )

		## Find static CUDA runtime library
		#string(REGEX REPLACE "libcudart[._][a-z]*" "" CUDA_LIBRARIES_PATH ${CUDA_LIBRARIES} )
		#setup_message( "Look for libcudart_static.a in \"${CUDA_LIBRARIES_PATH}\"")
		#find_file( CUDA_LIBRARIES_STATIC libcudart_static.a ${CUDA_LIBRARIES_PATH} )
		#if( CUDA_LIBRARIES_STATIC )
		#	set( CUDA_LIBRARIES_STATIC_FOUND true CACHE INTERNAL "Static CUDA libraries found, if true")
		#else()
		#	set( CUDA_LIBRARIES_STATIC_FOUND false CACHE INTERNAL "Static CUDA libraries found, if true" )
		#endif()
		
		setup_message( "found CUDA v${CUDA_VERSION_STRING}" )
		setup_message( "include dir: ${CUDA_INCLUDE_DIRS}" INDENT )
		setup_message( "runtime library: ${CUDA_LIBRARIES}" INDENT )
		setup_message( "static runtime library: ${CUDA_cudart_static_LIBRARY}" INDENT )

		set( CUDA_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS} CACHE INTERNAL "CUDA include directories" )
		set( CUDA_LIBRARIES ${CUDA_LIBRARIES} CACHE INTERNAL "CUDA libraries" )
		set( CUDA_cudart_static_LIBRARY ${CUDA_cudart_static_LIBRARY} CACHE INTERNAL "CUDA runtime static libraries" )
		
		set( CUDA_COMPATIBLE_CC_BINARY_SET $ENV{CUDA_COMPATIBLE_CC_BINARY} )
		if( CUDA_COMPATIBLE_CC_BINARY_SET )
			list( APPEND CUDA_NVCC_FLAGS "--compiler-bindir $ENV{CUDA_COMPATIBLE_CC_BINARY}" )
			setup_message( "nvcc will use user-defined C-compiler for compiling CUDA code ('$ENV{CUDA_COMPATIBLE_CC_BINARY}')" INDENT )
		else()
			setup_message( "nvcc will use default C-compiler for CUDA code ('${CMAKE_C_COMPILER}')" INDENT )
		endif()
	else()
		setup_message( "WARNING: could not find CUDA" )
	endif()
endif()
if( ${CMAKE_PROJECT_NAME}_NO_CUDA )
	add_definitions( -D${CMAKE_PROJECT_NAME}_NO_CUDA )
	setup_message( "CUDA is DISABLED" )
else()
	setup_message( "CUDA is enabled" )
endif()	
