cmake_minimum_required(VERSION 3.20)
project(Detatchment VERSION 0.1.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE)
IF(DEFINED CMAKE_BUILD_TYPE AND CMAKE_VERSION VERSION_GREATER "3.0")
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS  ${CMAKE_CONFIGURATION_TYPES})
ENDIF()

set(CMAKE_INSTALL_PREFIX .)

if(CMAKE_BUILD_TYPE MATCHES "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -O3 -pedantic -std=gnu++17 -fopenmp -DNDEBUG -lstdc++")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -Wall -O3 -fopenmp -lboost_date_time -std=gnu++17 -DNDEBUG -lstdc++")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -g3 -pedantic -fopenmp -std=gnu++17 -lstdc++")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -lboost_date_time -fopenmp -std=gnu++17 -g3 -lstdc++")
endif()

add_subdirectory(Graph)
add_subdirectory(Utils)
add_subdirectory(Tests)
add_subdirectory(Algorithms)
add_executable(run_experiments run_experiments.cpp)

install(TARGETS run_experiments DESTINATION ../${CMAKE_BUILD_TYPE})

target_link_libraries(run_experiments
  graph
  random 
  algorithms
  util
  stdc++fs
  stdc++ 
  tests
)