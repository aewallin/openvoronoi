SET(test_name "pytest_random_points" )

MESSAGE(STATUS "configuring py test: " ${test_name})
set( RAND_PT_CASES 128 256 512 1024 ) # 2048 4096 10000 20000
foreach( CASE ${RAND_PT_CASES} )
    ADD_TEST(${test_name}_${CASE} python ${CMAKE_SOURCE_DIR}/test/${test_name}/random_points.py 0 ${CASE})
endforeach()
