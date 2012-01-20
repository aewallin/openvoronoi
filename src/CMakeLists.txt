project(OpenVoronoi)

cmake_minimum_required(VERSION 2.4)


#
# Turn compiler warnings up to 11, at least with gcc.  I dont know how to
# do this with other compilers we might support and I'm leaving it up to
# the relevant platform maintainers...
# #include'ing the boost graph-library creates deprecated warnings
# thus we use now use -Wno-deprecated here.
#
if (CMAKE_BUILD_TOOL MATCHES "make")
    MESSAGE(STATUS "setting gcc options: -W -Wall -Wundef  
    -Wshadow -Winline -Wno-long-long -Wno-deprecated -Werror -pedantic-errors")
    
    # -pedantic -W -Wall -Wundef -Wfloat-equal -Wshadow -Winline -Wno-long-long
    add_definitions( -Werror -Wall -Wundef -Wno-error=uninitialized -Wshadow  -Wno-long-long -Wno-deprecated -pedantic -pedantic-errors)
    
    # some reasons why we don't enable certain errors:
    # this gives warning when comparing float/double with != or ==:  -Wfloat-equal
    # this gives warnings for BGL: -Winline
    
    
    add_definitions(-fPIC)
endif (CMAKE_BUILD_TOOL MATCHES "make")

# uncomment one of these to change build-type
 set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Profile)

MESSAGE(STATUS " CMAKE_BUILD_TYPE  = " ${CMAKE_BUILD_TYPE})

if (CMAKE_BUILD_TYPE MATCHES "Profile")
    set(CMAKE_CXX_FLAGS_PROFILE -pg)
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_PROFILE = " ${CMAKE_CXX_FLAGS_PROFILE})
endif (CMAKE_BUILD_TYPE MATCHES "Profile")


if (CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
endif (CMAKE_BUILD_TYPE MATCHES "Release")

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
endif (CMAKE_BUILD_TYPE MATCHES "Debug")

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
    MESSAGE(STATUS "Boost_LIB_VERSION: " ${Boost_LIB_VERSION})
    MESSAGE(STATUS "Boost_INCLUDE_DIRS are: " ${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "Boost_PYTHON_LIBRARY is: " ${Boost_PYTHON_LIBRARY})
    MESSAGE(STATUS "boost_LIBRARY_DIRS is: " ${Boost_LIBRARY_DIRS})
    MESSAGE(STATUS "Boost_LIBRARIES is: " ${Boost_LIBRARIES})    
endif()

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ) # Findqd.cmake is in this dir
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})

message( STATUS "build type = " ${CMAKE_BUILD_TYPE})

################ create version_string.hpp, http://stackoverflow.com/questions/3780667
# include the output directory, where the version_string.hpp file is generated
include_directories(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE(version_string.cmake)



# this defines the source-files
set(OVD_SRC
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram.cpp
    ${OpenVoronoi_SOURCE_DIR}/vertex.cpp
    ${OpenVoronoi_SOURCE_DIR}/edge.cpp
    ${OpenVoronoi_SOURCE_DIR}/checker.cpp
    ${OpenVoronoi_SOURCE_DIR}/facegrid.cpp
    ${OpenVoronoi_SOURCE_DIR}/common/point.cpp
    ${OpenVoronoi_SOURCE_DIR}/vertex_positioner.cpp
    ${OpenVoronoi_SOURCE_DIR}/common/numeric.cpp
)

set( OVD_INCLUDE_FILES
    ${OpenVoronoi_SOURCE_DIR}/graph.hpp
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram.hpp
    ${OpenVoronoi_SOURCE_DIR}/vertex.hpp
    ${OpenVoronoi_SOURCE_DIR}/edge.hpp
    ${OpenVoronoi_SOURCE_DIR}/site.hpp
    ${OpenVoronoi_SOURCE_DIR}/checker.hpp
    ${OpenVoronoi_SOURCE_DIR}/vertex_positioner.hpp
    ${OpenVoronoi_SOURCE_DIR}/facegrid.hpp
    
    ${OpenVoronoi_SOURCE_DIR}/offset.hpp
    ${OpenVoronoi_SOURCE_DIR}/offset2.hpp
    ${OpenVoronoi_SOURCE_DIR}/medial_axis.hpp
    ${OpenVoronoi_SOURCE_DIR}/polygon_interior.hpp
    
    ${CMAKE_CURRENT_BINARY_DIR}/version_string.hpp
)

set( OVD_INCLUDE_COMMON_FILES
    ${OpenVoronoi_SOURCE_DIR}/common/numeric.hpp  
    ${OpenVoronoi_SOURCE_DIR}/common/point.hpp
    ${OpenVoronoi_SOURCE_DIR}/common/halfedgediagram.hpp
)

set( OVD_INCLUDE_SOLVERS_FILES
    ${OpenVoronoi_SOURCE_DIR}/solvers/solution.hpp  
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver.hpp
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_lll.hpp
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_ppp.hpp
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_qll.hpp
)

# include directories
include_directories( ${OpenVoronoi_SOURCE_DIR} )

# the pure c++ shared library
add_library(
    libopenvoronoi
    SHARED
    ${OVD_SRC}
)
set_target_properties(libopenvoronoi PROPERTIES PREFIX "")

# the pure c++ shared library
add_library(
    openvoronoi_static
    STATIC
    ${OVD_SRC}
)

# this makes the Python module
add_library(
    openvoronoi
    MODULE
    py/open_voronoi_py.cpp
)
target_link_libraries(openvoronoi openvoronoi_static ${Boost_LIBRARIES} ${QD_LIBRARY}) 
set_target_properties(openvoronoi PROPERTIES PREFIX "") 


#
# this figures out where to install the Python modules
#
execute_process(
    COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib(0,0,\"/usr/local\")"
    OUTPUT_VARIABLE Python_site_packages
    OUTPUT_STRIP_TRAILING_WHITESPACE
) # on Ubuntu 11.10 this outputs: /usr/local/lib/python2.7/dist-packages


# strip away /usr/local/  because that is what CMAKE_INSTALL_PREFIX is set to
# also, since there is no leading "/", it makes ${Python_site_packages} a relative path.
STRING(REGEX REPLACE "/usr/local/(.*)$" "\\1" Python_site_packages "${Python_site_packages}" )

MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is : " ${CMAKE_INSTALL_PREFIX})
MESSAGE(STATUS "Python libraries will be installed to: " ${Python_site_packages})



    # this installs the python library
    install(
        TARGETS openvoronoi
        LIBRARY DESTINATION ${Python_site_packages}
    )

# this installs the c++ library
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
    install(
        FILES ${OVD_INCLUDE_COMMON_FILES}
        DESTINATION include/openvoronoi/common
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )
    install(
        FILES ${OVD_INCLUDE_SOLVERS_FILES}
        DESTINATION include/openvoronoi/solvers
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )


# "make package" stuff
 
execute_process(
    COMMAND dpkg --print-architecture
    OUTPUT_VARIABLE DEB_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)


#message( STATUS " git commit id: " ${GIT_COMMIT_ID})
# parse the commit id to get major, minor, patch
STRING(REGEX REPLACE "([0-9][0-9]).*" "\\1" GIT_MAJOR_VERSION "${GIT_COMMIT_ID}" )
#MESSAGE(STATUS "result= ${GIT_MAJOR_VERSION}")
STRING(REGEX REPLACE "[0-9][0-9].([0-9][0-9])-.*" "\\1" GIT_MINOR_VERSION "${GIT_COMMIT_ID}" )
#MESSAGE(STATUS "result= ${GIT_MINOR_VERSION}")
STRING(REGEX REPLACE "[0-9][0-9].[0-9][0-9]-(.*)-.*" "\\1" GIT_PATCH_VERSION "${GIT_COMMIT_ID}" )
#MESSAGE(STATUS "result= ${GIT_PATCH_VERSION}")

set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_CONTACT "Anders Wallin <anders.e.e.wallin@gmail.com>" CACHE STRING "email")
set(CPACK_PACKAGE_NAME "openvoronoi" CACHE STRING "name")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "2D voronoi diagram algorithm" CACHE STRING "name")
set(CPACK_PACKAGE_VENDOR https://github.com/aewallin/opencamlib)
set(CPACK_DEBIAN_PACKAGE_SECTION "science" CACHE STRING "name")
set(CPACK_DEBIAN_BUILD_DEPENDS debhelper python git cmake libboost-dev libqd0 libqd-dev libboost-python-dev libgomp1 CACHE STRINGS "name")
set(CPACK_DEBIAN_PACKAGE_DEPENDS python git cmake libboost-python libqd0 libgomp1 CACHE STRING "name")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${DEB_ARCHITECTURE} CACHE STRING "name")
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional CACHE STRING "name")
SET(CPACK_PACKAGE_VERSION "${GIT_MAJOR_VERSION}.${GIT_MINOR_VERSION}.${GIT_PATCH_VERSION}" CACHE STRING "name")

set_target_properties(libopenvoronoi PROPERTIES VERSION ${CPACK_PACKAGE_VERSION}) 
set_target_properties(openvoronoi PROPERTIES VERSION ${CPACK_PACKAGE_VERSION}) 

set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}_${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}" CACHE STRING "name")
set(CPACK_DEBIAN_DISTRIBUTION_NAME ubuntu CACHE STRING "name")
set(CPACK_DEBIAN_DISTRIBUTION_RELEASES oneiric CACHE STRING "name") 
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/debian/debian_package_description.txt CACHE STRING "${CPACK_PACKAGE_DESCRIPTION_FILE}")
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/debian/debian_copyright.txt CACHE STRING "name")
set(CPACK_DEBIAN_CHANGELOG "  * Latest development version." CACHE STRING "name")
set(DEB_SRC_DIR ${CMAKE_SOURCE_DIR} CACHE STRING "name" )

# run the source package generation as a target, so we don't have to do it every time
add_custom_target(spackage 
              ${CMAKE_COMMAND} -C ${CMAKE_BINARY_DIR}/CMakeCache.txt
              -P ${CMAKE_CURRENT_SOURCE_DIR}/debian/DebSourcePPA.cmake 
              # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            ) # get the changes file as an output from here?

include(CPack) # this should be last
