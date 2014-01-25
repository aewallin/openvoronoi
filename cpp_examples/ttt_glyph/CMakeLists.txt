project(OVD_TST)

cmake_minimum_required(VERSION 2.4)

if (CMAKE_BUILD_TOOL MATCHES "make")
    add_definitions(-Wall -Werror -Wno-deprecated -pedantic-errors)
endif (CMAKE_BUILD_TOOL MATCHES "make")

# find BOOST 
find_package( Boost COMPONENTS program_options REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "found Boost: " ${Boost_LIB_VERSION})
    MESSAGE(STATUS "boost-incude dirs are: " ${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "linking with boost: " ${Boost_LIBRARIES})
endif()

# openvoronoi library
find_library(OVD_LIBRARY 
            NAMES libopenvoronoi openvoronoi
            PATHS /usr/local/lib/openvoronoi
            DOC "openvoronoi"
            REQUIRED
)
MESSAGE(STATUS "OVD_LIBRARY is now: " ${OVD_LIBRARY})
include_directories( /usr/local/include/openvoronoi )


# the treutype-tracer library - provides font input geometry
# https://github.com/aewallin/truetype-tracer
find_library(TTT_LIBRARY 
            NAMES libtruetypetracer truetypetracer
            PATHS /usr/local/lib/truetypetracer
            DOC "rpg"
            REQUIRED
)
include_directories( /usr/local/include/ttt )

find_package(Freetype REQUIRED)
if( FREETYPE_FOUND )
    MESSAGE(STATUS "FREETYPE_LIBRARY is: " ${FREETYPE_LIBRARY})
    MESSAGE(STATUS "FREETYPE_INCLUDE_DIRS is: " ${FREETYPE_INCLUDE_DIRS})
    include_directories(${FREETYPE_INCLUDE_DIRS})
endif()


set(OVD_TST_SRC
    main.cpp
)

add_executable(
    ovd_tst
    ${OVD_TST_SRC}
)
target_link_libraries(ovd_tst ${OVD_LIBRARY} ${TTT_LIBRARY} ${Boost_LIBRARIES})


