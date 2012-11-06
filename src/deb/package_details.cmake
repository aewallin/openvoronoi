
# message(STATUS " called with SRC_DIR is = " ${SRC_DIR})
# message(STATUS " called with version  = " ${MY_VERSION})
message(STATUS "package_details.cmake: Configuring debian package details.")

set(CPACK_GENERATOR "DEB" CACHE string "generator" )
set(CPACK_PACKAGE_CONTACT "Anders Wallin <anders.e.e.wallin@gmail.com>" CACHE STRING "email")
set(CPACK_PACKAGE_NAME "openvoronoi" CACHE STRING "name")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "2D voronoi diagram algorithm" CACHE STRING "name")
set(CPACK_PACKAGE_VENDOR https://github.com/aewallin/openvoronoi CACHE STRING "web")
set(CPACK_DEBIAN_PACKAGE_SECTION "science" CACHE STRING "name3")

# the debian source-package builder wants dependencies as an array
set(DEBSRC_BUILD_DEPENDS debhelper python git cmake libboost-dev "libqd0 | libqd2c2a" libqd-dev libboost-python-dev libgomp1 CACHE STRINGS "name")

# we need to explicitly list the libboost-python versions here. why??
# quantal has 1.49.0
# precise has 1.48.0
# oneiric has 1.46.1
# natty/maverick has 1.42.0
# lucid has 1.40.0
set(DEBSRC_PACKAGE_DEPENDS python git cmake 
                # "libboost-python1.49.0 | libboost-python1.48.0 | libboost-python1.46.1 | libboost-python1.42.0 | libboost-python1.40.0"
                "libboost-python"
                "libqd0 | libqd2c2a" CACHE STRING "name")

# however CPack wants dependencies as a single comma separated string!
set(CPACK_DEBIAN_PACKAGE_DEPENDS)
foreach(DEP ${DEBSRC_PACKAGE_DEPENDS})
    set(CPACK_DEBIAN_PACKAGE_DEPEND "${CPACK_DEBIAN_PACKAGE_DEPENDS}, ${DEP}")
endforeach(DEP ${DEBSRC_PACKAGE_DEPENDS})  

set(CPACK_DEBIAN_BUILD_DEPENDS)
foreach(DEP ${DEBSRC_BUILD_DEPENDS})
    set(CPACK_DEBIAN_BUILD_DEPENDS "${CPACK_DEBIAN_BUILD_DEPENDS}, ${DEP}")
endforeach(DEP ${DEBSRC_BUILD_DEPENDS})  

# find out and set DEB_ARCHITECTURE
execute_process(
    COMMAND dpkg --print-architecture
    OUTPUT_VARIABLE DEB_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE
    
)
message(STATUS "package_details.cmake: your architecture is ${DEB_ARCHITECTURE}")

set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE ${DEB_ARCHITECTURE} CACHE STRING "name6")
set(CPACK_DEBIAN_PACKAGE_PRIORITY optional CACHE STRING "name7")
SET(CPACK_PACKAGE_VERSION ${MY_VERSION} CACHE STRING "name8")
set(CPACK_DEBIAN_DISTRIBUTION_NAME ubuntu CACHE STRING "name9")
set(CPACK_DEBIAN_DISTRIBUTION_RELEASES lucid maverick natty oneiric precise quantal CACHE STRING "name10") 
# message(STATUS " CMAKE_SOURCE_DIR is = " ${CMAKE_SOURCE_DIR})
if(${SRC_DIR} MATCHES "")
    set(CPACK_PACKAGE_DESCRIPTION_FILE ${SRC_DIR}/deb/debian_package_description.txt CACHE STRING "package description file")
    set(CPACK_RESOURCE_FILE_LICENSE ${SRC_DIR}/deb/debian_copyright.txt CACHE STRING "name11")
    set(DEB_SRC_DIR ${SRC_DIR} CACHE STRING "name13" )
else(${SRC_DIR} MATCHES "")
    set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_SOURCE_DIR}/deb/debian_package_description.txt CACHE STRING "package description file")
    set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/deb/debian_copyright.txt CACHE STRING "name11")
    set(DEB_SRC_DIR ${CMAKE_SOURCE_DIR} CACHE STRING "name13" )
endif(${SRC_DIR} MATCHES "")

#message(STATUS " descr file = ${CPACK_PACKAGE_DESCRIPTION_FILE}")


set(DEB_SRC_DIR ${CMAKE_SOURCE_DIR} CACHE STRING "name" )
