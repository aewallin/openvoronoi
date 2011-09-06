project(OpenVoronoi)

cmake_minimum_required(VERSION 2.4)



# some options to set what is built:
option(BUILD_CXX_LIB
  "Build/install the pure c++ library? " ON)

option(BUILD_PY_LIB 
  "Build/install the python ocl library? " ON)


if (NOT BUILD_CXX_LIB)
    MESSAGE(STATUS " Note: will NOT build pure c++ library")
endif(NOT BUILD_CXX_LIB)

if (NOT BUILD_PY_LIB)
    MESSAGE(STATUS " Note: will NOT build python library")
endif(NOT BUILD_PY_LIB)



#
# Turn compiler warnings up to 11, at least with gcc.  I dont know how to
# do this with other compilers we might support and I'm leaving it up to
# the relevant platform maintainers...
# #include'ing the boost graph-library creates deprecated warnings
# thus we use now use -Wno-deprecated here.
#
if (CMAKE_BUILD_TOOL MATCHES "make")
    MESSAGE(STATUS "setting gcc options: -Wall -Werror -Wno-deprecated -pedantic-errors")
    
    # optimization:
    add_definitions(-O3  -Wall  -Wno-deprecated -Werror -pedantic-errors)
    
    # NON-optimized build:
    #add_definitions(-O3 -Wall  -Wno-deprecated -Werror -pedantic-errors)
    
    add_definitions(-fPIC)
endif (CMAKE_BUILD_TOOL MATCHES "make")

#
# this figures out the Python include directories and adds them to the
# header file search path
#
execute_process(
    COMMAND python-config --includes
    COMMAND sed -r "s/-I//g; s/ +/;/g"
    COMMAND tr -d '\n'
    OUTPUT_VARIABLE Python_Includes
)
include_directories(${Python_Includes})


include_directories(${CMAKE_CURRENT_BINARY_DIR})


# find BOOST and boost-python
find_package( Boost COMPONENTS python REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "found Boost: " ${Boost_LIB_VERSION})
    MESSAGE(STATUS "boost-incude dirs are: " ${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "boost-python lib is: " ${Boost_PYTHON_LIBRARY})
    MESSAGE(STATUS "boost_LIBRARY_DIRS is: " ${Boost_LIBRARY_DIRS})
    MESSAGE(STATUS "Boost_LIBRARIES is: " ${Boost_LIBRARIES})    
endif()



# this defines the source-files

set(OVD_SRC
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram.cpp
    ${OpenVoronoi_SOURCE_DIR}/voronoivertex.cpp
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram_checker.cpp
    ${OpenVoronoi_SOURCE_DIR}/facegrid.cpp
    ${OpenVoronoi_SOURCE_DIR}/point.cpp

)


set( OVD_INCLUDE_FILES  
    ${OpenVoronoi_SOURCE_DIR}/point.hpp
    ${OpenVoronoi_SOURCE_DIR}/halfedgediagram.hpp
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram_graph.hpp
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram.hpp
    ${OpenVoronoi_SOURCE_DIR}/voronoivertex.hpp
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram_checker.hpp
    ${OpenVoronoi_SOURCE_DIR}/facegrid.hpp
    
)


# include directories
include_directories( ${OpenVoronoi_SOURCE_DIR} )

add_library(
    libopenvoronoi
    STATIC
    ${OVD_SRC}
)
set_target_properties(libopenvoronoi PROPERTIES PREFIX "") 
    

    # this makes the Python module
    add_library(
        openvoronoi
        MODULE
        open_voronoi_py.cpp
    )
    MESSAGE(STATUS "linking python binary with boost: " ${Boost_PYTHON_LIBRARY})
    target_link_libraries(openvoronoi libopenvoronoi ${Boost_LIBRARIES} ) #${Boost_PYTHON_LIBRARY})
    set_target_properties(openvoronoi PROPERTIES PREFIX "") 


#
# this figures out where to install the Python modules
#

execute_process(
    COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()"
    OUTPUT_VARIABLE Python_site_packages
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
MESSAGE(STATUS "Python libraries will be installed to: " ${Python_site_packages})

if (BUILD_PY_LIB)
    # this installs the python library
    install(
        TARGETS openvoronoi
        LIBRARY DESTINATION ${Python_site_packages}
    )
endif (BUILD_PY_LIB)

# this installs the c++ library
if (BUILD_CXX_LIB)
    install(
        TARGETS libopenvoronoi
        LIBRARY 
        DESTINATION lib/openvoronoi
        ARCHIVE DESTINATION lib/openvoronoi
        PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )

    # this installs the c++ include headers
    install(
        FILES ${OVD_INCLUDE_FILES}
        DESTINATION include/openvoronoi
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )
endif (BUILD_CXX_LIB)
