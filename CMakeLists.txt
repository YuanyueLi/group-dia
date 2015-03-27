cmake_minimum_required(VERSION 2.6)
project("group-dia") 

## Set default build type (if not set by user on command line)
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release)
endif()

## force build type into the cache (needs to be set beforehand)
set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)

## Set programme
set(group_dia
normalize_rt
analysis_identification
analysis_swath
analysis_validate
)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
aux_source_directory(${PROJECT_SOURCE_DIR}/library/ group_dia_lib)
aux_source_directory(${PROJECT_SOURCE_DIR}/alglib/ alglib)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++0x")
endif()
#add_definitions(/D_GLIBCXX_DEBUG)

find_package(OpenMS)
if (OpenMS_FOUND)
  include_directories(${OPENMS_INCLUDE_DIRS}  ./library ./alglib)
  add_definitions(${OPENMS_DEFINITIONS})  
  add_library(group_dia_lib SHARED ${group_dia_lib}) 
  add_library(alglib SHARED ${alglib})

  foreach(i ${group_dia})
      ADD_EXECUTABLE(${i} ./pro/${i}.cpp)
      target_link_libraries(${i} OpenMS group_dia_lib alglib)
  endforeach(i)

else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)

