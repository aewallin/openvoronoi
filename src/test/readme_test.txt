Running tests

The Python and C++ tests can be run with CTest. In the top-level
build-directory either "make test" or "ctest" will run all tests. You can
run only tests that have e.g. "ttt" in the test-name with
  $ ctest -R ttt
or
  $ make test ARGS="-R ttt"
for example, to run all python tests
$ ctest -R pytest
or tu run all c++ tests (these are usually much faster than the python tests)
$ ctest -R cpp

Notes

- Currently the tests do not produce any output (png or svg output could be
  an option?)
- The Python tests cannot be run without first calling "sudo make install".
  And if the Python tests are to respond to recoding in OVD, OVD must first
  be reinstalled. This means that developers will usually have an unstable
  version of OVD installed in "/usr/" or "/usr/local/". This seems to be on
  purpose, in order to verify that the packaging system works. Can we
  somehow test packaging separately, and only in that case use the installed
  version? Or can we stage and test the install in isolation using fakeroot,
  or something along those lines?


Adding new Python tests

Python test scripts are pretty simple to write: they call "exit(0)" to
indicate success, and "exit(-1)" to indicate failure. Adding a new Python
test script is also simple: place the new script beneath "src/test/pytest_<testname>", and
tell CTest about new test scripts by  in a CMakeLists.txt 
"src/test/pytest_<testname>/CMakeLists.txt" by  adding something like the following
  ADD_TEST(<test name> python ../src/test/pytest_<testname>/<test file>.py)

Adding new C++ tests

C++ tests are also easy to write: an "int main()" function returns "0" for
success, and "-1" for failure. Each C++ test has its own subdirectory starting
with "cpptest_", which should contain a CMakeLists.txt that calls 
  ADD_TEST( <test name> )

The easiest way to add a C++ test is to make a copy of the template
"cpptest_minimal" subdirectory

