project(random_points_benchmark)

cmake_minimum_required(VERSION 2.4)

if (CMAKE_BUILD_TOOL MATCHES "make")
    add_definitions(-Wall -Werror -Wno-deprecated -pedantic-errors)
endif (CMAKE_BUILD_TOOL MATCHES "make")

set(CMAKE_BUILD_TYPE Release) # change to "Debug" for assert() checks etc.
MESSAGE(STATUS " CMAKE_BUILD_TYPE  = " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
endif (CMAKE_BUILD_TYPE MATCHES "Release")
set( CMAKE_CXX_FLAGS  ${CMAKE_CXX_FLAGS_RELEASE})

MESSAGE(STATUS " CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

find_package( Boost REQUIRED )
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "found Boost: " ${Boost_LIB_VERSION})
    MESSAGE(STATUS "boost-incude dirs are: " ${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "linking with boost: " ${Boost_LIBRARIES})
endif()

find_library(OVD_LIBRARY 
            NAMES openvoronoi
            PATHS /usr/local/lib/openvoronoi
            DOC "openvoronoi"
            REQUIRED
)
MESSAGE(STATUS "OVD_LIBRARY is now: " ${OVD_LIBRARY})
include_directories( /usr/local/include/openvoronoi )

set(OVD_TST_SRC
    main.cpp
)

add_executable(
    random_points_benchmark
    ${OVD_TST_SRC}
)
target_link_libraries(random_points_benchmark ${OVD_LIBRARY} ${Boost_LIBRARIES})


