include(CTest)
ENABLE_TESTING()


#
# C++ tests. Each test subdirectory has its own CMakeLists.txt file with
# build and test commands. 
file(GLOB TEST_SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/cpptest_*)
foreach(subdir ${TEST_SUBDIRS})
  add_subdirectory(${subdir})
endforeach()
