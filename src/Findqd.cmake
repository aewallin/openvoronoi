#  Try to find QD
#  see: http://crd.lbl.gov/~dhbailey/mpdist/
#
#  QD_FOUND - system has QD
#  QD_INCLUDE_DIR - the QGLViewer include directory
#  QD_LIBRARY - Link these to use QD
#

find_path(QD_INCLUDE_DIR
    NAMES qd/qd_real.h qd_real.h
    PATHS /usr/include
    /usr/local/include
)

find_library(QD_LIBRARY
    NAMES libqd qd
    PATHS /usr/lib
        /usr/local/lib
   ENV QDROOT
   ENV LD_LIBRARY_PATH
   ENV LIBRARY_PATH
)

    
if(QD_LIBRARY)
    set(QD_LIBRARY_ optimized ${QD_LIBRARY} debug ${QD_LIBRARY})
    set(QD_LIBRARY ${QD_LIBRARY_} CACHE FILEPATH "The QD library")
endif(QD_LIBRARY)

IF(QD_INCLUDE_DIR AND QD_LIBRARY)
    SET(QD_FOUND TRUE)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DQD_FOUND")
ENDIF(QD_INCLUDE_DIR AND QD_LIBRARY)

IF(QD_FOUND)
    #MESSAGE(STATUS "QD_INCLUDE_DIR is: ${QD_INCLUDE_DIR}")
    #MESSAGE(STATUS "QD_LIBRARY is: ${QD_LIBRARY}")
ELSE(QD_FOUND)
    MESSAGE(FATAL_ERROR "Could not find QD")
ENDIF(QD_FOUND)
