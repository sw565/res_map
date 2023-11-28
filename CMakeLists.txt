cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

project(ACC)

set(INC_DIRS 
)


set (EXE_SPECS
  cmake.res_map
)

include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
