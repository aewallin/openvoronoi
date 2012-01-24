#include "version_string.hpp"

namespace ovd {

std::string version() {
    return VERSION_STRING; // the git-commit-id, e.g. "12.01-102-gf516b17"
}
std::string build_type() {
    return BUILDTYPE_STRING; // "Release"  "Debug" or "Profile"
}

}
