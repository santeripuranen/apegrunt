# Version $Id: $

cmake_minimum_required(VERSION 3.3)

#project( ${PROJECT_NAME}_DEPENDENCIES )

# Setup message helper function
include( setup_message.cmake NO_POLICY_SCOPE )

# initialize CMAKE_MODULE_PATH
include( initialize_cmake_module_path.cmake NO_POLICY_SCOPE )

# Setup of external dependencies
include( boost.cmake NO_POLICY_SCOPE )
include( tbb.cmake NO_POLICY_SCOPE )
include( cuda.cmake NO_POLICY_SCOPE )
