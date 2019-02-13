# Version $Id: $

# Aoegrunt src

############################
## Add sources and includes
###

include_directories(
	${CMAKE_CURRENT_BINARY_DIR}
	${APEGRUNT_INCLUDE_DIR}
	${Boost_INCLUDE_DIR}
	${TBB_INCLUDE_DIRS}
)

link_directories( ${Boost_LIBRARY_DIRS} )

set( APEGRUNT_LIBRARY_SOURCES
	Apegrunt_options.cpp
	aligned_allocator.cpp
) # *.cpp *.hpp *.cc

set( APEGRUNT_BINARY_SOURCES
	Apegrunt.cpp
)

#################################
## Add libraries and executables
###

set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib ) # dynamic libraries
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib ) # static libraries
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )


### General optimization flags	
set( APEGRUNT_OPTIMIZATION_FLAGS "${APEGRUNT_OPTIMIZATION_FLAGS} -O3 -mavx2 -static-libstdc++" )
set( APEGRUNT_OPTIMIZATION_FLAGS "${APEGRUNT_OPTIMIZATION_FLAGS} -ftree-vectorize -fwhole-program" )
set( APEGRUNT_OPTIMIZATION_FLAGS "${APEGRUNT_OPTIMIZATION_FLAGS} -flto -ffat-lto-objects" ) # -flto-report" )

# GCC-specific optimization flags
set( APEGRUNT_GCC_OPTIMIZATION_FLAGS "${APEGRUNT_GCC_OPTIMIZATION_FLAGS}" )

# Clang-specific optimization flags
set( APEGRUNT_CLANG_OPTIMIZATION_FLAGS "${APEGRUNT_CLANG_OPTIMIZATION_FLAGS}" )


### General release build flags
set( APEGRUNT_RELEASE_FLAGS "-w -DNDEBUG -fvisibility=hidden -fvisibility-inlines-hidden -Wl,--strip-all" )

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall ${APEGRUNT_OPTIMIZATION_FLAGS} ${APEGRUNT_GCC_OPTIMIZATION_FLAGS} ${APEGRUNT_RELEASE_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall ${APEGRUNT_OPTIMIZATION_FLAGS} ${APEGRUNT_CLANG_OPTIMIZATION_FLAGS} ${APEGRUNT_RELEASE_FLAGS}")
endif()

# debug build flags
set( APEGRUNT_DEBUG_FLAGS "-pg -g -ftree-vectorizer-verbose=2" )
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall ${APEGRUNT_OPTIMIZATION_FLAGS} ${APEGRUNT_GCC_OPTIMIZATION_FLAGS} ${APEGRUNT_DEBUG_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall ${APEGRUNT_OPTIMIZATION_FLAGS} ${APEGRUNT_CLANG_OPTIMIZATION_FLAGS} ${APEGRUNT_DEBUG_FLAGS}")
endif()

#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wignored-attributes" )

# set preferred linker
set( CMAKE_LINKER "ld.gold" )

# Add apegrunt library
if( APEGRUNT_LIBRARY_SOURCES )
	#add_library( libapegrunt SHARED ${APEGRUNT_LIBRARY_SOURCES} )
	#set_target_properties( libapegrunt PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" )

	add_library( libapegrunt STATIC ${APEGRUNT_LIBRARY_SOURCES} )
	#set_target_properties( libapegrunt PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" )

	set_target_properties( libapegrunt PROPERTIES ARCHIVE_OUTPUT_NAME "apegrunt" )
	set_target_properties( libapegrunt PROPERTIES COMPILE_FLAGS "--std=c++14" )
	set( APEGRUNT_LIBRARIES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/libapegrunt.a CACHE INTERNAL "The apegrunt libraries" )
endif()

# Add apegrunt standalone executable
if( APEGRUNT_BINARY_SOURCES )
	# need APEGRUNT_LIBRARY_SOURCES here for the APEGRUNT_STANDALONE_BUILD macro mechanism, rather than a "target_link_libraries( libapegrunt )"
	add_executable( apegrunt ${APEGRUNT_LIBRARY_SOURCES} ${APEGRUNT_BINARY_SOURCES} )
	set_target_properties( apegrunt PROPERTIES COMPILE_DEFINITIONS APEGRUNT_STANDALONE_BUILD )
	set_target_properties( apegrunt PROPERTIES COMPILE_FLAGS "--std=c++14" )
	#target_link_libraries( apegrunt libapegrunt )
endif()

# Add Boost libraries
if( NOT APEGRUNT_NO_BOOST )
	target_link_libraries( libapegrunt ${Boost_LIBRARIES} )
	target_link_libraries( apegrunt ${Boost_LIBRARIES} )
endif()

# Add TBB libraries
if( NOT APEGRUNT_NO_TBB )
	target_link_libraries( libapegrunt ${TBB_LIBRARIES} )
	target_link_libraries( apegrunt ${TBB_LIBRARIES} )
endif()

if( UNIX )
	target_link_libraries( libapegrunt pthread rt )
	target_link_libraries( apegrunt pthread rt )
endif( UNIX )

# Prevent linking against shared libraries on OS X;
# Apple gcc always links against a shared version of a library if present,
# regardless of -Bstatic or equivalent linker flags.
if(APPLE)
	set_target_properties( libapegrunt PROPERTIES LINK_SEARCH_END_STATIC TRUE )
	set_target_properties( apegrunt PROPERTIES LINK_SEARCH_END_STATIC TRUE )
endif(APPLE)
