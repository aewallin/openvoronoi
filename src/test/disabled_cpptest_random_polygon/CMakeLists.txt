SET(test_name "cpptest_random_polygon" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})

# the rpg lib
find_library(RPG_LIBRARY 
            NAMES randompolygon
            PATHS /usr/local/lib/randompolygon
            DOC "rpg"
            REQUIRED
)
include_directories( /usr/local/include/randompolygon )
if(RPG_LIBRARY_FOUND)
	set(SOURCE_FILES random_polygon.cpp)
	add_executable(${test_name} ${SOURCE_FILES})
	add_dependencies(${test_name}  libopenvoronoi)


	unset(Boost_LIBRARIES) # required because this contains boost-python when we come here
	find_package( Boost COMPONENTS program_options REQUIRED)

	target_link_libraries(${test_name} libopenvoronoi ${RPG_LIBRARY} ${Boost_LIBRARIES})

	ADD_TEST(${test_name} ${test_name}) # must provide --n option!

	ADD_TEST(${test_name}_help ${test_name} --help)
	set_property(
		TEST ${test_name}_help
		PROPERTY WILL_FAIL TRUE
	)

	ADD_TEST(${test_name}_10 ${test_name} --n 10 )
	ADD_TEST(${test_name}_n10_s5 ${test_name} --n 10 --s 5)
else()
	MESSAGE(STATUS "skipping c++ test: " ${test_name} ", dependency not found")
endif()