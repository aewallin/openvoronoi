#use git for a pretty commit id
#uses 'git describe --tags', so tags are required in the repo
#create a tag with 'git tag <name>' and 'git push --tags'

#version_string.hpp will define VERSION_STRING to something like "test-1-g5e1fb47"
# where test is the name of the last tagged git revision, 1 is the number of commits since that tag,
# 'g' is ???, and 5e1fb47 is the first 7 chars of the git sha1 commit id.


find_package(Git)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags 
        RESULT_VARIABLE res_var 
        OUTPUT_VARIABLE GIT_COM_ID 
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
    if( NOT ${res_var} EQUAL 0 )
        message( WARNING "Git failed (not a repo, or no tags)." )
        file(READ "git-tag.txt" GIT_COMMIT_ID)
        message( STATUS "version_string.cmake read from file GIT_COMMIT_ID: " ${GIT_COMMIT_ID})
    else()
        string( REPLACE "\n" "" GIT_COMMIT_ID ${GIT_COM_ID} )
        message( STATUS "version_string.cmake git set GIT_COMMIT_ID: " ${GIT_COMMIT_ID})
    endif()
    
else()
    # if we don't have git, try to read git-tag from file instead
    file(READ "git-tag.txt" GIT_COMMIT_ID)
    
    #set( GIT_COMMIT_ID "unknown (git not found!)")
    message( STATUS "version_string.cmake read from file GIT_COMMIT_ID: " ${GIT_COMMIT_ID})
    #message( WARNING "Git not found. Reading tag from git-tag.txt instead: " ${GIT_COMMIT_ID})
endif()

# get the compiler and its version
exec_program(
  ${CMAKE_CXX_COMPILER}
  ARGS                    --version
  OUTPUT_VARIABLE _compiler_output
)
# we want only the first line of output
string(REGEX MATCH "^([^\n]+)" 
  compiler_version 
  "${_compiler_output}"
)

message(STATUS "version_string.cmake: C++ compiler: ${CMAKE_CXX_COMPILER}" )
message(STATUS "version_string.cmake: C++ compiler version: ${compiler_version}")
  

set( vstring "//version_string.hpp - written by cmake. changes will be lost!\n"
             "#ifndef VERSION_STRING\n"
             "#define VERSION_STRING \"${GIT_COMMIT_ID}\"\n"
             "#endif\n"
             " \n"
             "#ifndef BUILDTYPE_STRING\n"
             "#define BUILDTYPE_STRING \"${CMAKE_BUILD_TYPE}\"\n"
             "#endif\n"
             " \n"
             "#ifndef COMPILER_STRING\n"
             "#define COMPILER_STRING \"${compiler_version}\"\n"
             "#endif\n"
             " \n"
             "#ifndef SYSTEM_STRING\n"
             "#define SYSTEM_STRING \"${CMAKE_SYSTEM}\"\n"
             "#endif\n"
             " \n"
             "#ifndef PROCESSOR_STRING\n"
             "#define PROCESSOR_STRING \"${CMAKE_SYSTEM_PROCESSOR}\"\n"
             "#endif\n"
)

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/version_string.hpp ${vstring} )
set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/version_string.hpp
    PROPERTIES GENERATED TRUE
    HEADER_FILE_ONLY TRUE
)

STRING(REGEX REPLACE "([0-9][0-9]).*" "\\1" GIT_MAJOR_VERSION "${GIT_COMMIT_ID}" )
STRING(REGEX REPLACE "[0-9][0-9].([0-9][0-9])-.*" "\\1" GIT_MINOR_VERSION "${GIT_COMMIT_ID}" )
STRING(REGEX REPLACE "[0-9][0-9].[0-9][0-9]-(.*)-.*" "\\1" GIT_PATCH_VERSION "${GIT_COMMIT_ID}" )
