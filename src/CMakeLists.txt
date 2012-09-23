project(OpenVoronoi)

cmake_minimum_required(VERSION 2.4)

# figure out the gcc version
INCLUDE(gcc_version.cmake)


#
# Turn compiler warnings up to 11, at least with gcc.  I dont know how to
# do this with other compilers we might support and I'm leaving it up to
# the relevant platform maintainers...
# #include'ing the boost graph-library creates deprecated warnings
# thus we use now use -Wno-deprecated here.
#
if (CMAKE_BUILD_TOOL MATCHES "make")
    if(GCC_4_6)
        MESSAGE(STATUS "setting strict gcc options: -Wall -Werror etc...")
        add_definitions(-Werror -Wall -Wundef -Wno-error=uninitialized -Wshadow  -Wno-long-long -Wno-deprecated -pedantic -pedantic-errors)
    else(GCC_4_6)
        MESSAGE(STATUS "setting less strict gcc options.. (you have an old gcc version!)")
        add_definitions(-Wall -Wno-deprecated -pedantic-errors)
    endif(GCC_4_6)
    
    # some reasons why we don't enable certain errors:
    # -Wfloat-equal        gives warning when comparing float/double with != or ==:  
    # -Winline             gives warnings for BGL
    
    add_definitions(-fPIC)
endif (CMAKE_BUILD_TOOL MATCHES "make")


#
# Turn compiler warnings up to 11, at least with gcc.  I dont know how to
# do this with other compilers we might support and I'm leaving it up to
# the relevant platform maintainers...
# #include'ing the boost graph-library creates deprecated warnings
# thus we use now use -Wno-deprecated here.
#
#if (CMAKE_BUILD_TOOL MATCHES "make")
#    MESSAGE(STATUS "setting gcc options: -W -Wall -Wundef  -Wshadow -Winline -Wno-long-long -Wno-deprecated -Werror -pedantic-errors")
    # -pedantic -W -Wall -Wundef -Wfloat-equal -Wshadow -Winline -Wno-long-long
#    add_definitions( -Werror -Wall -Wundef -Wno-error=uninitialized -Wshadow  -Wno-long-long -Wno-deprecated -pedantic -pedantic-errors)
#endif (CMAKE_BUILD_TOOL MATCHES "make")


# uncomment one of these to change build-type

 set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Profile)
# set(CMAKE_BUILD_TYPE Coverage)

option(BUILD_PYTHON_MODULE "Build python-module? " ON) # this should make building the python-module optional
# when we build only the pure c++ library we should not need boost-ptyhon
# for a build without boost-python:
# cmake -DBUILD_PYTHON_MODULE=OFF ../src


option(BUILD_DOC "Build doxygen documentation? " ON)

#set(BUILD_TESTS 1)
option(BUILD_TESTS "Build tests?" OFF) # don't build tests

MESSAGE(STATUS " CMAKE_BUILD_TYPE  = " ${CMAKE_BUILD_TYPE})

if (CMAKE_BUILD_TYPE MATCHES "Profile")
    set(CMAKE_CXX_FLAGS_PROFILE "-p -g -DNDEBUG")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_PROFILE = " ${CMAKE_CXX_FLAGS_PROFILE})
endif (CMAKE_BUILD_TYPE MATCHES "Profile")


if (CMAKE_BUILD_TYPE MATCHES "Release")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_RELEASE = " ${CMAKE_CXX_FLAGS_RELEASE})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif (CMAKE_BUILD_TYPE MATCHES "Release")

if (CMAKE_BUILD_TYPE MATCHES "Debug")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS_DEBUG = " ${CMAKE_CXX_FLAGS_DEBUG})
    #MESSAGE(STATUS " CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
endif (CMAKE_BUILD_TYPE MATCHES "Debug")
MESSAGE(STATUS " CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})

if (CMAKE_BUILD_TYPE MATCHES "Coverage")
    set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage -fno-inline ${CMAKE_CXX_FLAGS_DEBUG}")
    MESSAGE(STATUS " CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
endif (CMAKE_BUILD_TYPE MATCHES "Coverage")



#
# this figures out the Python include directories and adds them to the
# header file search path
#

MESSAGE(STATUS "Your system is: ${CMAKE_SYSTEM_NAME}")
IF( ${BUILD_PYTHON_MODULE} MATCHES ON)
	find_package(PythonLibs REQUIRED)
	include_directories(${PYTHON_INCLUDE_DIRS})
	MESSAGE(STATUS "Python includes: " ${PYTHON_INCLUDE_DIRS})
ENDIF()

# find boost-python
IF( ${BUILD_PYTHON_MODULE} MATCHES ON)
	find_package( Boost COMPONENTS python REQUIRED)
	MESSAGE(STATUS "Boost_PYTHON_LIBRARY is: " ${Boost_PYTHON_LIBRARY})
ENDIF()

find_package( Boost COMPONENTS graph REQUIRED)
# alternative: find_package(Boost REQUIRED)
# find boost
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "Boost_LIB_VERSION: " ${Boost_LIB_VERSION})
    MESSAGE(STATUS "Boost_INCLUDE_DIRS are: " ${Boost_INCLUDE_DIRS})
    
    MESSAGE(STATUS "boost_LIBRARY_DIRS is: " ${Boost_LIBRARY_DIRS})
    MESSAGE(STATUS "Boost_LIBRARIES is: " ${Boost_LIBRARIES})    
endif()

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ) # Findqd.cmake is in this dir
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})
message( STATUS "libqd include dir = "${QD_INCLUDE_DIR})
message( STATUS "libqd library = "${QD_LIBRARY})

message( STATUS "build type = " ${CMAKE_BUILD_TYPE})

################ create version_string.hpp, http://stackoverflow.com/questions/3780667
INCLUDE(version_string.cmake)
set(MY_VERSION "${GIT_MAJOR_VERSION}.${GIT_MINOR_VERSION}.${GIT_PATCH_VERSION}" CACHE STRING "version")
# include the output directory, where version_string.hpp is generated
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# this defines the source-files
set(OVD_SRC
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram.cpp
    ${OpenVoronoi_SOURCE_DIR}/vertex.cpp
    ${OpenVoronoi_SOURCE_DIR}/edge.cpp
    ${OpenVoronoi_SOURCE_DIR}/checker.cpp
    ${OpenVoronoi_SOURCE_DIR}/common/point.cpp
    ${OpenVoronoi_SOURCE_DIR}/vertex_positioner.cpp
    ${OpenVoronoi_SOURCE_DIR}/common/numeric.cpp
    ${OpenVoronoi_SOURCE_DIR}/medial_axis_pocket.cpp
    ${OpenVoronoi_SOURCE_DIR}/medial_axis_walk.cpp
    ${OpenVoronoi_SOURCE_DIR}/offset.cpp
)

set( OVD_INCLUDE_FILES
    ${OpenVoronoi_SOURCE_DIR}/graph.hpp
    ${OpenVoronoi_SOURCE_DIR}/voronoidiagram.hpp
    ${OpenVoronoi_SOURCE_DIR}/vertex.hpp
    ${OpenVoronoi_SOURCE_DIR}/edge.hpp
    ${OpenVoronoi_SOURCE_DIR}/site.hpp
    ${OpenVoronoi_SOURCE_DIR}/checker.hpp
    ${OpenVoronoi_SOURCE_DIR}/vertex_positioner.hpp
    ${OpenVoronoi_SOURCE_DIR}/kdtree.hpp

    ${OpenVoronoi_SOURCE_DIR}/offset.hpp
    ${OpenVoronoi_SOURCE_DIR}/offset_sorter.hpp

    ${OpenVoronoi_SOURCE_DIR}/filter.hpp
    ${OpenVoronoi_SOURCE_DIR}/medial_axis_filter.hpp
    ${OpenVoronoi_SOURCE_DIR}/medial_axis_walk.hpp
    ${OpenVoronoi_SOURCE_DIR}/medial_axis_pocket.hpp
    ${OpenVoronoi_SOURCE_DIR}/polygon_interior_filter.hpp
    ${OpenVoronoi_SOURCE_DIR}/island_filter.hpp
    
    ${CMAKE_CURRENT_BINARY_DIR}/version_string.hpp
    ${CMAKE_SOURCE_DIR}/version.hpp
)

set( OVD_INCLUDE_COMMON_FILES
    ${OpenVoronoi_SOURCE_DIR}/common/numeric.hpp  
    ${OpenVoronoi_SOURCE_DIR}/common/point.hpp
    ${OpenVoronoi_SOURCE_DIR}/common/halfedgediagram.hpp
    
)

set( OVD_INCLUDE_UTIL_FILES
    ${OpenVoronoi_SOURCE_DIR}/utility/vd2svg.hpp    
    ${OpenVoronoi_SOURCE_DIR}/utility/simple_svg_1.0.0.hpp
)


set( OVD_INCLUDE_SOLVERS_FILES
    ${OpenVoronoi_SOURCE_DIR}/solvers/solution.hpp  
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver.hpp

    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_lll.hpp
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_lll_para.hpp

    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_ppp.hpp
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_qll.hpp
    
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_sep.hpp
    ${OpenVoronoi_SOURCE_DIR}/solvers/solver_alt_sep.hpp
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
set_target_properties(libopenvoronoi PROPERTIES VERSION ${MY_VERSION}) 
# uncomment to link against libqd here. an alternative is to link libqd when an application using openvoronoi is built.
target_link_libraries(libopenvoronoi ${Boost_LIBRARIES} ${QD_LIBRARY}) 

# c++lib for coverage testing
add_library(
    libopenvoronoi-coverage
    SHARED
    EXCLUDE_FROM_ALL # not built by default
    ${OVD_SRC}
)
set_target_properties(libopenvoronoi-coverage PROPERTIES PREFIX "")
set_target_properties( libopenvoronoi-coverage
    PROPERTIES
    CMAKE_BUILD_TYPE Coverage
    CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage -fno-inline ${CMAKE_CXX_FLAGS_DEBUG}"
    )

IF( ${BUILD_PYTHON_MODULE} MATCHES ON)

	# static library, linked with the python module below
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
	set_target_properties(openvoronoi PROPERTIES VERSION ${MY_VERSION}) 

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
        LIBRARY 
        DESTINATION ${Python_site_packages}
    )



    # this installs the python VTK-visualization helper library
    install(
        FILES ${OpenVoronoi_SOURCE_DIR}/py/ovdvtk.py
              ${OpenVoronoi_SOURCE_DIR}/py/ngc_writer.py
              ${OpenVoronoi_SOURCE_DIR}/py/ovdgenerators.py
        DESTINATION ${Python_site_packages}
        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )
    
ENDIF()
 
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
	FILES ${OVD_INCLUDE_UTIL_FILES}
	DESTINATION include/openvoronoi/utility
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
)
install(
	FILES ${OVD_INCLUDE_SOLVERS_FILES}
	DESTINATION include/openvoronoi/solvers
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
)


# tests
if( ${BUILD_TESTS} MATCHES ON)
    include(${CMAKE_SOURCE_DIR}/test/ovd_tests.cmake)
endif()

# doxygen documentation
include(doxygen.cmake)

# coverage-report target
message(STATUS "lcov cmake src dir ${CMAKE_SOURCE_DIR}")
add_custom_target(coverage-report
    COMMAND lcov --directory ./ --zerocounters
    # COMMAND make ExperimentalCoverage
    COMMAND ctest -R cpptest
    COMMAND lcov --directory ./ --capture  --output-file testcoverage.info # --base-directory ${CMAKE_SOURCE_DIR}
    COMMAND lcov --directory ./ --extract testcoverage.info \"*/openvoronoi*\" --output-file testcoverage_ext.info
    COMMAND genhtml --show-details --prefix \"${CMAKE_SOURCE_DIR}\" --output-directory coverage-report --title OpenVoronoi testcoverage_ext.info
)



# "make package" stuff
execute_process(
    COMMAND dpkg --print-architecture
    OUTPUT_VARIABLE DEB_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

include(package_details.cmake)

# run the source package generation as a target, so we don't have to do it every time
add_custom_target(spackage 
              ${CMAKE_COMMAND} -D SRC_DIR:STRING=${CMAKE_SOURCE_DIR} -D MY_VERSION:STRING=${MY_VERSION} -C ${CMAKE_SOURCE_DIR}/package_details.cmake
              -P ${CMAKE_CURRENT_SOURCE_DIR}/deb/DebSourcePPA.cmake 
              #ERROR_QUIET # WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            ) # get the changes file as an output from here?
# -C ${CMAKE_BINARY_DIR}/CMakeCache.txt
include(CPack) # this should be last
