
SET(test_name "cpptest_arc" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})

set(SOURCE_FILES arc.cpp)
add_executable( ${test_name} ${SOURCE_FILES})
add_dependencies(${test_name}  libopenvoronoi)

#include_directories( ${CMAKE_SOURCE_DIR} )
    
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})

unset(Boost_LIBRARIES) # this may contain boost-python and other unwanted things
find_package( Boost COMPONENTS program_options REQUIRED)

target_link_libraries(${test_name} libopenvoronoi ${QD_LIBRARY} ${Boost_LIBRARIES})

ADD_TEST(${test_name} ${test_name})

ADD_TEST(${test_name}_help ${test_name} --help)
set_property(
    TEST ${test_name}_help
    PROPERTY WILL_FAIL TRUE
)

ADD_TEST(${test_name}_d ${test_name} --d ) # run in debug mode
ADD_TEST(${test_name}_r ${test_name} --r ) # reverse arc direction
