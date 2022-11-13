SET(test_name "pytest_ttt_alphabet" )

MESSAGE(STATUS "configuring py test: " ${test_name})



set( TTT_ALPHA_CASES 200 100 50 25 )  # 12 6 3
# comment out the small values so that test doesn't take forever!
# cpu-times (t400s laptop):
#   25== takes 23 seconds 
#   12== takes ca 71 seconds
#    6== takes ca 384 seconds
#    3== takes ca ??  seconds
foreach( CASE ${TTT_ALPHA_CASES}) # characters A..Z
    ADD_TEST(
        NAME ${test_name}_${CASE} 
        #CONFIGURATIONS Nightly
        COMMAND python ${CMAKE_SOURCE_DIR}/test/${test_name}/ttt_alphabet.py ${CASE})
endforeach()
