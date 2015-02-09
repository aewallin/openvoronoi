SET(test_name "pytest_2opt_rpg" )

MESSAGE(STATUS "configuring py test: " ${test_name})

set( 2OPT_RPG_CASES 5 10 15 20 30) # the number of vertices in the polygon
set( 2OPT_RPG_MAX_SEED 500) # run test for seeds 0,1,2,3,...,max-1
foreach( CASE ${2OPT_RPG_CASES} )
    ADD_TEST(${test_name}_pts${CASE}_reps${2OPT_RPG_MAX_SEED} python ${CMAKE_SOURCE_DIR}/test/${test_name}/2opt_random_polygon.py ${CASE} ${2OPT_RPG_MAX_SEED})
endforeach()

set( 2OPT_RPG_CASES 40 50 60 100 200) # 400 800 the number of vertices in the polygon
set( 2OPT_RPG_MAX_SEED 10) # run test for seeds 0,1,2,3,...,max-1
foreach( CASE ${2OPT_RPG_CASES} )
    ADD_TEST(${test_name}_pts${CASE}_reps${2OPT_RPG_MAX_SEED} python ${CMAKE_SOURCE_DIR}/test/${test_name}/2opt_random_polygon.py ${CASE} ${2OPT_RPG_MAX_SEED})
endforeach()

set( 2OPT_RPG_CASES 400 800) # 400 800 the number of vertices in the polygon
set( 2OPT_RPG_MAX_SEED 3) # run test for seeds 0,1,2,3,...,max-1
foreach( CASE ${2OPT_RPG_CASES} )
    ADD_TEST(${test_name}_pts${CASE}_reps${2OPT_RPG_MAX_SEED} python ${CMAKE_SOURCE_DIR}/test/${test_name}/2opt_random_polygon.py ${CASE} ${2OPT_RPG_MAX_SEED})
endforeach()
