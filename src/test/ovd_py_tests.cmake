include(CTest)
ENABLE_TESTING()

# Python tests
# Each test subdirectory has its own CMakeLists.txt file with
# build and test commands. 

file(GLOB PYTEST_SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/pytest_*)
foreach(subdir ${PYTEST_SUBDIRS})
  add_subdirectory(${subdir})
endforeach()

