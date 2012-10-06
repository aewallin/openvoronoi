include(CTest)
ENABLE_TESTING()

# Python tests
# Each test subdirectory has its own CMakeLists.txt file with
# build and test commands. 

file(GLOB PYTEST_SUBDIRS test/pytest_*)
subdirs(${PYTEST_SUBDIRS}) # calls the cmakelists.txt in each test_ subdirectory (?)

