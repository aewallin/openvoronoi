SET(test_name "pytest_ttt_single_glyph" )

MESSAGE(STATUS "configuring py test: " ${test_name})

set( ROTS 0 45 90 135 180 225 270 315 )

foreach( CASE RANGE 25) # characters A..Z
    foreach( ROT ${ROTS} ) # rotation angles
        ADD_TEST(${test_name}_G${CASE}_R${ROT} python ${CMAKE_SOURCE_DIR}/test/${test_name}/ttt_single_glyph.py ${CASE} ${ROT})
    endforeach()
endforeach()
