project(OpenVoronoi)

cmake_minimum_required(VERSION 3.0)

# figure out the gcc version
INCLUDE(gcc_version.cmake)

# uninstall target
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# experimental use of clang, since it detects more warnings/errors than gcc
option(USE_CLANG "build application with clang" OFF) # OFF is the default
# Note: useful packages: clang, llvm, gcc-multilib
# see: http://stackoverflow.com/questions/7031126/switching-between-gcc-and-clang-llvm-using-cmake
if(USE_CLANG)
  SET (CMAKE_C_COMPILER             "/usr/bin/clang")
  SET (CMAKE_C_FLAGS                "-Wall -std=c99")
  SET (CMAKE_C_FLAGS_DEBUG          "-g")
  SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
  SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

  SET (CMAKE_CXX_COMPILER             "/usr/bin/clang++")
  SET (CMAKE_CXX_FLAGS                "-Wall")
  SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
  SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET (CMAKE_CXX_FLAGS_RELEASE        "-DNDEBUG") # -O3 gives "argument unused" ?
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

  add_definitions(-O3)
  # -Wno-error=uninitialized -Wno-deprecated -Wno-long-long
  add_definitions(-Werror -Wall -Wundef  -Wshadow -pedantic -pedantic-errors)
  add_definitions(-fPIC)
endif()

MESSAGE(STATUS " Your compiler is = " ${CMAKE_CXX_COMPILER})
MESSAGE(STATUS " Your linker is = " ${CMAKE_LINKER})


#
# Turn compiler warnings up to 11, at least with gcc.  I don't know how to
# do this with other compilers we might support and I'm leaving it up to
# the relevant platform maintainers...
# #include'ing the boost graph-library creates deprecated warnings
# thus we use now use -Wno-deprecated here.
#
if (UNIX AND NOT APPLE)
  MESSAGE(STATUS "setting strict gcc options: -Wall -Werror etc...")
  add_definitions(-Werror -Wall -Wundef -Wno-error=uninitialized -Wshadow  -Wno-long-long -Wno-deprecated -pedantic -pedantic-errors)

  # some reasons why we don't enable certain errors:
  # -Wfloat-equal        gives warning when comparing float/double with != or ==:  
  # -Winline             gives warnings for BGL

  add_definitions(-fPIC)
endif ()



# uncomment one of these to change build-type
# Release runs much faster than the others
# Debug/Profile/Coverage run slower but are useful for testing

# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_BUILD_TYPE Profile)
# set(CMAKE_BUILD_TYPE Coverage)
MESSAGE(STATUS " CMAKE_BUILD_TYPE  = " ${CMAKE_BUILD_TYPE})

option(BUILD_PYTHON_MODULE "Build python-module? " ON) # building the python-module optional, ON by default
# when we build only the pure c++ library we do not need boost-ptyhon
# command-line for a build without boost-python:
# cmake -DBUILD_PYTHON_MODULE=OFF ../src

option(BUILD_DOC "Build doxygen documentation? " ON)

option(BUILD_CPP_TESTS "Build c++ tests?" ON) 
option(BUILD_PY_TESTS "Build/configure Python tests?" ON) 


if (CMAKE_BUILD_TYPE MATCHES "Profile")
  set(CMAKE_CXX_FLAGS_PROFILE "-p -g -DNDEBUG")
  MESSAGE(STATUS " CMAKE_CXX_FLAGS_PROFILE = " ${CMAKE_CXX_FLAGS_PROFILE})
  set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_PROFILE})
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

# build the library for code-coverage testing
if (CMAKE_BUILD_TYPE MATCHES "Coverage")
  set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage -fno-inline ${CMAKE_CXX_FLAGS_DEBUG}")
  MESSAGE(STATUS " CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS})
endif (CMAKE_BUILD_TYPE MATCHES "Coverage")

# this figures out the Python include directories and adds them to the
# header file search path
MESSAGE(STATUS "Your system is: ${CMAKE_SYSTEM_NAME}")
IF( ${BUILD_PYTHON_MODULE} MATCHES ON)
  find_package (PythonLibs) # 2.7 EXACT REQUIRED 
  include_directories(${PYTHON_INCLUDE_DIRS})
  MESSAGE(STATUS "Python includes: " ${PYTHON_INCLUDE_DIRS})
ENDIF()



find_package( Boost REQUIRED )
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "Boost_LIB_VERSION: " ${Boost_LIB_VERSION})
  MESSAGE(STATUS "Boost_INCLUDE_DIRS are: " ${Boost_INCLUDE_DIRS})
  MESSAGE(STATUS "boost_LIBRARY_DIRS is: " ${Boost_LIBRARY_DIRS})
  MESSAGE(STATUS "Boost_LIBRARIES is: " ${Boost_LIBRARIES})
  
endif()

# find boost-python
IF( ${BUILD_PYTHON_MODULE} MATCHES ON)
  find_package( Boost COMPONENTS python REQUIRED)
  MESSAGE(STATUS "Boost_PYTHON_LIBRARY is: " ${Boost_PYTHON_LIBRARY})
ENDIF()

# find the qd library
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ) # Findqd.cmake is in this dir
find_package( qd REQUIRED ) 
include_directories(${QD_INCLUDE_DIR})

message( STATUS "libqd include dir = ${QD_INCLUDE_DIR}")
message( STATUS "libqd library = ${QD_LIBRARY}")
message( STATUS "build type = ${CMAKE_BUILD_TYPE}")

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
if (NOT APPLE)
  set_target_properties(libopenvoronoi PROPERTIES VERSION ${MY_VERSION}) 
endif (NOT APPLE)
# link against libqd here. 
# an alternative is to link libqd when an application using openvoronoi is built.
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


# the Python module
IF( ${BUILD_PYTHON_MODULE} MATCHES ON)

  message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
  message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
  message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
  # static c++ library, linked with the python module below
  add_library(
    openvoronoi_static
    STATIC
    ${OVD_SRC}
    )

  # the Python module itself
  add_library(
    openvoronoi
    MODULE
    py/open_voronoi_py.cpp
    )
  target_link_libraries(openvoronoi openvoronoi_static ${Boost_LIBRARIES} ${QD_LIBRARY} ${PYTHON_LIBRARIES}) 
  set_target_properties(openvoronoi PROPERTIES PREFIX "") 
  if (NOT APPLE)
    set_target_properties(openvoronoi PROPERTIES VERSION ${MY_VERSION}) 
  endif (NOT APPLE)
  # this figures out where to install the Python modules
  execute_process(
    COMMAND python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    OUTPUT_VARIABLE Python_site_packages
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ) # on Ubuntu 11.10 this outputs: /usr/local/lib/python2.7/dist-packages

  # strip away /usr/local/  because that is what CMAKE_INSTALL_PREFIX is set to
  # also, since there is no leading "/", it makes ${Python_site_packages} a relative path.
  # STRING(REGEX REPLACE "/usr/local/(.*)$" "\\1" Python_site_packages "${Python_site_packages}" )
  MESSAGE(STATUS "CMAKE_INSTALL_PREFIX is : " ${CMAKE_INSTALL_PREFIX})
  MESSAGE(STATUS "Python libraries will be installed to: " ${Python_site_packages})

  # this installs the python library
  install(
    TARGETS openvoronoi
    LIBRARY 
    DESTINATION ${Python_site_packages}
    )

  # this installs the python VTK-visualization helper py-files
  install(
    FILES ${OpenVoronoi_SOURCE_DIR}/py/ovdvtk.py
    ${OpenVoronoi_SOURCE_DIR}/py/ngc_writer.py
    ${OpenVoronoi_SOURCE_DIR}/py/ovdgenerators.py
    DESTINATION ${Python_site_packages}
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
    )
ENDIF() # end Python module

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
if( ${BUILD_CPP_TESTS} MATCHES ON)
  include(${CMAKE_SOURCE_DIR}/test/ovd_cpp_tests.cmake) # cmake file defines c++ tests
endif()
if( ${BUILD_PY_TESTS} MATCHES ON)
  include(${CMAKE_SOURCE_DIR}/test/ovd_py_tests.cmake) # cmake file defines Python tests
endif()

# doxygen documentation
include(doxygen.cmake)

# coverage-report target
# running "make coverage-report" produces a code-coverate report for the tests
# uses ctest, lcov, genhtml
message(STATUS "lcov cmake src dir ${CMAKE_SOURCE_DIR}")
add_custom_target(coverage-report
  COMMAND lcov --directory ./ --zerocounters
  # COMMAND make ExperimentalCoverage
  # COMMAND ctest -D Experimental Coverage -R cpptest
  COMMAND ctest -R cpptest -E "cpptest_arc\\|cpptest_ttt_glyph_big_7\\|cpptest_ttt_glyph_small_35" # exclude failing tests!
  COMMAND lcov --directory ./ --capture  --output-file testcoverage.info # --base-directory ${CMAKE_SOURCE_DIR}
  COMMAND lcov --directory ./ --extract testcoverage.info \"*/openvoronoi*\" --output-file testcoverage_ext.info
  COMMAND genhtml --show-details --prefix \"${CMAKE_SOURCE_DIR}\" --output-directory coverage-report --title OpenVoronoi testcoverage_ext.info
  )

# Debian package
# "make package" will build a binary package
# "make spackage" will build a source-package (for uploading to a PPA)
include(deb/package_details.cmake) 
add_custom_target(spackage 
  ${CMAKE_COMMAND} -D SRC_DIR:STRING=${CMAKE_SOURCE_DIR} 
  -D MY_VERSION:STRING=${MY_VERSION} 
  -C ${CMAKE_SOURCE_DIR}/deb/package_details.cmake # a lot of variables etc. defined here
  -P ${CMAKE_CURRENT_SOURCE_DIR}/deb/DebSourcePPA.cmake # the script that builds a debian package
  ) 
include(CPack) # this should be last
