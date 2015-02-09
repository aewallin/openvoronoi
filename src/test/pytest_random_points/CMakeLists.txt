SET(test_name "pytest_random_points" )

MESSAGE(STATUS "configuring py test: " ${test_name})
set( RAND_PT_CASES 128 256 512 1024) # 2048 4096 10000 20000
set( SEEDS 0 1 2 3 4 5 6 7 8 9 10)
foreach( CASE ${RAND_PT_CASES} )
    foreach( SEED ${SEEDS} )
        ADD_TEST(${test_name}_pts${CASE}_seed${SEED} python ${CMAKE_SOURCE_DIR}/test/${test_name}/random_points.py ${SEED} ${CASE})
    endforeach()
endforeach()

set( RAND_PT_CASES 2048 4096 10000 20000) # 2048 4096 10000 20000
set( SEEDS 0 1 )
foreach( CASE ${RAND_PT_CASES} )
    foreach( SEED ${SEEDS} )
        ADD_TEST(${test_name}_pts${CASE}_seed${SEED} python ${CMAKE_SOURCE_DIR}/test/${test_name}/random_points.py ${SEED} ${CASE})
    endforeach()
endforeach()
