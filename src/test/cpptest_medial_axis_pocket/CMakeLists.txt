
# The next line tells CMake and CTest about "test_minimal".
SET(test_name "cpptest_medial_axis_pocket" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})

set(SOURCE_FILES medial_axis_pocket.cpp)
add_executable(${test_name} ${SOURCE_FILES})
add_dependencies(${test_name}  libopenvoronoi)

unset(Boost_LIBRARIES) # required because this contains boost-python when we come here
find_package( Boost COMPONENTS program_options REQUIRED)

target_link_libraries(${test_name} libopenvoronoi  ${Boost_LIBRARIES})

ADD_TEST(${test_name} ${test_name})

ADD_TEST(${test_name}_help ${test_name} --help)
set_property(
    TEST ${test_name}_help
    PROPERTY WILL_FAIL TRUE
)

ADD_TEST(${test_name}_debug ${test_name} --d )
