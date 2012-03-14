


SET(test_name "cpptest_random_linesegs" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})
set(SOURCE_FILES random_linesegs.cpp)
add_executable(${test_name} ${SOURCE_FILES})

find_library(OVD2_LIBRARY 
            NAMES openvoronoi
            PATHS /usr/local/lib/openvoronoi
            DOC "openvoronoi"
            REQUIRED
)
#MESSAGE(STATUS "OVD_LIBRARY for c++ tests is now: " ${OVD2_LIBRARY})
#include_directories( /usr/local/include/openvoronoi )
    
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})

unset(Boost_LIBRARIES) # required because this contains boost-python when we come here
#MESSAGE(STATUS "linking ${test_name} with:  ${Boost_LIBRARIES}")
find_package( Boost COMPONENTS program_options REQUIRED)
#MESSAGE(STATUS "linking ${test_name} with:  ${Boost_LIBRARIES}")
target_link_libraries(${test_name} ${OVD2_LIBRARY} ${QD_LIBRARY} ${Boost_LIBRARIES})


ADD_TEST(${test_name} ${test_name}) # must provide --n option!
set_property(
    TEST ${test_name}
    PROPERTY WILL_FAIL TRUE
)
ADD_TEST(${test_name}_help ${test_name} --help)
set_property(
    TEST ${test_name}_help
    PROPERTY WILL_FAIL TRUE
)

#ADD_TEST(${test_name}_b ${test_name} --b 2)
ADD_TEST(${test_name}_200 ${test_name} --n 200)

# for coverage-testing this takes too long..
#ADD_TEST(${test_name}_10000 ${test_name} --n 10000)
