


SET(test_name "test_random_points" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})
set(SOURCE_FILES random_points.cpp)
add_executable(${test_name} ${SOURCE_FILES})

find_library(OVD2_LIBRARY 
            NAMES openvoronoi
            PATHS /usr/local/lib/openvoronoi
            DOC "openvoronoi"
            REQUIRED
)
#MESSAGE(STATUS "OVD_LIBRARY for c++ tests is now: " ${OVD2_LIBRARY})
include_directories( /usr/local/include/openvoronoi )
    
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})

find_package( Boost COMPONENTS program_options REQUIRED)

target_link_libraries(${test_name} ${OVD2_LIBRARY} ${QD_LIBRARY} ${Boost_LIBRARIES})


ADD_TEST(${test_name} ${test_name})
ADD_TEST(${test_name}_1000 ${test_name} --n 1000)
ADD_TEST(${test_name}_10000 ${test_name} --n 10000)
