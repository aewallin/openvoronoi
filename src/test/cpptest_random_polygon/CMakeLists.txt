SET(test_name "cpptest_random_polygon" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})
set(SOURCE_FILES random_polygon.cpp)
add_executable(${test_name} ${SOURCE_FILES})

find_library(OVD2_LIBRARY 
            NAMES openvoronoi
            PATHS /usr/local/lib/openvoronoi
            DOC "openvoronoi"
            REQUIRED
)
# find the qd library
find_package( qd REQUIRED )
include_directories(${QD_INCLUDE_DIR})
unset(Boost_LIBRARIES) # required because this contains boost-python when we come here
find_package( Boost COMPONENTS program_options REQUIRED)
# the rpg lib
find_library(RPG_LIBRARY 
            NAMES librpg rpg
            PATHS /usr/local/lib/rpg
            DOC "rpg"
            REQUIRED
)
include_directories( /usr/local/include/rpg )
target_link_libraries(${test_name} ${OVD2_LIBRARY} ${RPG_LIBRARY} ${QD_LIBRARY} ${Boost_LIBRARIES})

ADD_TEST(${test_name} ${test_name}) # must provide --n option!

ADD_TEST(${test_name}_help ${test_name} --help)
set_property(
    TEST ${test_name}_help
    PROPERTY WILL_FAIL TRUE
)

ADD_TEST(${test_name}_10 ${test_name} --n 10 )
ADD_TEST(${test_name}_n10_s5 ${test_name} --n 10 --s 5)

