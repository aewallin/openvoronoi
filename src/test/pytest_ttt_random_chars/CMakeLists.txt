SET(test_name "pytest_ttt_random_chars" )

MESSAGE(STATUS "configuring py test: " ${test_name})

foreach( CASE RANGE 25) # seeds 1..25. These take about 4s each to run.
    ADD_TEST(
        NAME ${test_name}_${CASE} 
        #CONFIGURATIONS Nightly
        COMMAND python ${CMAKE_SOURCE_DIR}/test/${test_name}/ttt_random_chars.py ${CASE})
endforeach()
