SET(test_name "cpptest_ttt_glyph" )

MESSAGE(STATUS "configuring c++ test: " ${test_name})
# the ttt lib
find_library(TTT_LIBRARY 
            NAMES truetypetracer
            PATHS /usr/local/lib/truetypetracer
            DOC "ttt"
            REQUIRED
)
include_directories( /usr/local/include/truetypetracer )
	if( TTT_LIBRARY_FOUND )
	set(SOURCE_FILES ttt_glyph.cpp)
	add_executable(${test_name} ${SOURCE_FILES})
	add_dependencies(${test_name}  libopenvoronoi)

	unset(Boost_LIBRARIES) # required because this contains boost-python when we come here
	find_package( Boost COMPONENTS program_options REQUIRED)

	# freetype
	find_package(Freetype REQUIRED)
	if( FREETYPE_FOUND )
	   MESSAGE(STATUS "FREETYPE_LIBRARY is: " ${FREETYPE_LIBRARY})
	   MESSAGE(STATUS "FREETYPE_INCLUDE_DIRS is: " ${FREETYPE_INCLUDE_DIRS})
	endif()
	include_directories(${FREETYPE_INCLUDE_DIRS})

	target_link_libraries(${test_name} libopenvoronoi ${TTT_LIBRARY} ${Boost_LIBRARIES})

	ADD_TEST(${test_name} ${test_name})

	ADD_TEST(${test_name}_help ${test_name} --help)
	set_property(
		TEST ${test_name}_help
		PROPERTY WILL_FAIL TRUE
	)

	ADD_TEST(${test_name}_fail ${test_name} --c 53) # illegal char index
	set_property(
		TEST ${test_name}_fail
		PROPERTY WILL_FAIL TRUE
	)

	# --c options sets the glyph
	# 0,1,2,3...51  corresponds to "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
	set( CHARS   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25) #  
	set( CHARS2 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51) # 35

# --s option sets the scale. 3e-4 is a big character, 1e-6 is tiny
#set( SCALES 0 1 2 3 4 )
set( SCALES 0 ) # faster

# --d option sets subdivision of glyph, 50 is normal, 10 gives more segments
#set( SUBS 60 50 40 30 20 10 )
set( SUBS 50 ) # faster

	foreach( CASE ${CHARS} )
		foreach( SCALE ${SCALES} )
			foreach( SUB ${SUBS} )
				ADD_TEST(${test_name}_G${CASE}_s${SCALE}_sub${SUB} ${test_name} --c ${CASE} --s ${SCALE} --d ${SUB})
			endforeach()
		endforeach()
	endforeach()
	foreach( CASE ${CHARS2} )
		foreach( SCALE ${SCALES} )
			foreach( SUB ${SUBS} )
				ADD_TEST(${test_name}_g${CASE}_s${SCALE}_sub${SUB} ${test_name} --c ${CASE} --s ${SCALE} --d ${SUB})
			endforeach()
		endforeach()
	endforeach()
else()
	MESSAGE(STATUS "skipping c++ test: " ${test_name} ", dependency not found")
endif()
