#include "version_string.hpp"

namespace ovd {

/// OpenVoronoi version as string, e.g. "12.01-102-gf516b17"
/// the format is: year.month-commits-commitID
std::string version() {
    return VERSION_STRING; // the git-commit-id, e.g. "12.01-102-gf516b17"
}
/// build type description string, e.g. "Release"  "Debug" or "Profile".
std::string build_type() {
    return BUILDTYPE_STRING;
}

/// Compiler description string. Set at build-time by version_string.cmake
std::string compiler() {
    return COMPILER_STRING;
}

/// OS description string. Set at build-time by version_string.cmake
std::string system() {
    return SYSTEM_STRING;
}

/// cpu description string. Set at build-time by version_string.cmake
std::string processor() {
    return PROCESSOR_STRING;
}

}
