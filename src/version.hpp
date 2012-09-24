#include "version_string.hpp"

namespace ovd {

/// return OpenVoronoi version as string, e.g. "12.01-102-gf516b17"
std::string version() {
    return VERSION_STRING; // the git-commit-id, e.g. "12.01-102-gf516b17"
}
/// return build type as string, e.g. "Release"  "Debug" or "Profile" 
std::string build_type() {
    return BUILDTYPE_STRING;
}

std::string compiler() {
    return COMPILER_STRING;
}

std::string system() {
    return SYSTEM_STRING;
}

std::string processor() {
    return PROCESSOR_STRING;
}

}
