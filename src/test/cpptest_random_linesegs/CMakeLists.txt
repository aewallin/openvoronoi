SET(test_name "cpptest_random_linesegs" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})
set(SOURCE_FILES random_linesegs.cpp)
add_executable(${test_name} ${SOURCE_FILES})
add_dependencies(${test_name}  libopenvoronoi)

set(Boost_LIBRARIES) # required because this contains boost-python when we come here
find_package( Boost COMPONENTS program_options REQUIRED)
target_link_libraries(${test_name} libopenvoronoi ${Boost_LIBRARIES})

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

ADD_TEST(${test_name}_d ${test_name} --n 10 --d )

#ADD_TEST(${test_name}_b ${test_name} --b 2)
ADD_TEST(${test_name}_42 ${test_name} --n 42)

# for coverage-testing this takes too long..
#ADD_TEST(${test_name}_10000 ${test_name} --n 10000)
