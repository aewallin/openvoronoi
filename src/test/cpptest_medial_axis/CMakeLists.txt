
# The next line tells CMake and CTest about "test_minimal".
SET(test_name "cpptest_medial_axis" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})

set(SOURCE_FILES medial_axis.cpp)
add_executable(${test_name} ${SOURCE_FILES})

find_library(OVD2_LIBRARY 
            NAMES openvoronoi
            PATHS /usr/local/lib/openvoronoi
            DOC "openvoronoi"
)
include_directories( /usr/local/include/openvoronoi )
    
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})

target_link_libraries(${test_name} ${OVD2_LIBRARY} ${QD_LIBRARY} )

ADD_TEST(${test_name} ${test_name})
