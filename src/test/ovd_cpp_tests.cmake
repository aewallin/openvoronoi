include(CTest)
ENABLE_TESTING()


#
# C++ tests. Each test subdirectory has its own CMakeLists.txt file with
# build and test commands. 
file(GLOB TEST_SUBDIRS test/cpptest_*)
subdirs(${TEST_SUBDIRS}) # calls the cmakelists.txt in each test_ subdirectory (?)
