
# The next line tells CMake and CTest about "test_minimal".
SET(test_name "cpptest_para_lines" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})

set(SOURCE_FILES para_lines.cpp)
add_executable(${test_name} ${SOURCE_FILES})
add_dependencies(${test_name}  libopenvoronoi)

unset(Boost_LIBRARIES) 
find_package( Boost COMPONENTS program_options REQUIRED)

target_link_libraries(${test_name} libopenvoronoi ${Boost_LIBRARIES})

ADD_TEST(${test_name} ${test_name})
ADD_TEST(${test_name}_r ${test_name} --r ) # reverse version

ADD_TEST(${test_name}_help ${test_name} --help)
set_property(
    TEST ${test_name}_help
    PROPERTY WILL_FAIL TRUE
)

ADD_TEST(${test_name}_d ${test_name} --d ) # run in debug mode
