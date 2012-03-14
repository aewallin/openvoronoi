ENABLE_TESTING()

# Python tests
#
#
file(GLOB PYTEST_SUBDIRS test/pytest_*)
subdirs(${PYTEST_SUBDIRS}) # calls the cmakelists.txt in each test_ subdirectory (?)


#
# C++ tests. Each test subdirectory has its own CMakeLists.txt file with
# build and test commands. 
file(GLOB TEST_SUBDIRS test/test_*)
subdirs(${TEST_SUBDIRS}) # calls the cmakelists.txt in each test_ subdirectory (?)
