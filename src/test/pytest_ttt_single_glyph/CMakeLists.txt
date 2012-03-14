SET(test_name "pytest_ttt_single_glyph" )

MESSAGE(STATUS "configuring py test: " ${test_name})


foreach( CASE RANGE 25) # characters A..Z
    ADD_TEST(${test_name}_${CASE} python ${CMAKE_SOURCE_DIR}/test/${test_name}/ttt_single_glyph.py ${CASE})
endforeach()
