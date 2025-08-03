# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

OpenVoronoi uses CMake as its build system. The standard build process:

```bash
mkdir bld
cd bld
cmake ../src
make
sudo make install
```

Key CMake options:
- `BUILD_PYTHON_MODULE=ON` (default) - Build Python bindings
- `BUILD_DOC=ON` (default) - Build Doxygen documentation  
- `BUILD_CPP_TESTS=ON` (default) - Build C++ tests
- `BUILD_PY_TESTS=ON` (default) - Build Python tests
- `USE_CLANG=OFF` (default) - Use Clang compiler

Build types: Release, Debug, Profile, Coverage

## Testing Commands

The project uses CTest for test management:

```bash
make test           # Run all tests
ctest               # Alternative test runner
ctest -R pattern    # Run specific tests (e.g., ctest -R ttt)
```

Test structure:
- C++ tests: `src/test/cpptest_*`
- Python tests: `src/test/pytest_*`

For code coverage:
```bash
make coverage-report
```

## Architecture Overview

OpenVoronoi is a C++ computational geometry library with Python bindings for computing 2D Voronoi diagrams. It uses an incremental topology-oriented algorithm optimized for point, line-segment, and circular-arc sites.

### Core Components

**Main Algorithm Classes:**
- `VoronoiDiagram` - Primary algorithm implementation in `src/`
- `VoronoiDiagramChecker` - Validation and sanity checking
- Graph structure: `Vertex`, `Edge` classes for topology representation

**Geometric Solvers (`src/solvers/`):**
- Multiple specialized solvers for different site type combinations
- Uses quad-double precision arithmetic (libqd) for numerical robustness
- Critical for handling edge cases in geometric computation

**Common Utilities (`src/common/`):**
- `Point` class - 2D point/vector operations
- Numeric utilities and precision handling

**Python Interface (`src/py/`):**
- Boost.Python bindings (maintained for compatibility)
- Comprehensive wrapper exposing C++ functionality
- Requires Python 3.8+ and Boost.Python

### Key Directory Structure

- `src/` - Core C++ algorithm implementation
- `src/solvers/` - Geometric computation solvers
- `src/utility/` - I/O utilities (SVG output, visualization)
- `cpp_examples/` - C++ usage demonstrations
- `python_examples/` - Extensive Python examples (60+ files)

### Application Focus

Primary use case is **2D offset generation for CNC machining**:
- Voronoi-based offset computation
- Medial axis filtering for pocket machining
- Toolpath calculation and optimization
- Support for complex manufacturing scenarios

The library emphasizes numerical robustness over speed, using quad-double precision arithmetic to handle challenging geometric edge cases that occur in real-world manufacturing applications.

### Dependencies

**Required:**
- cmake, libqd-dev (quad-double arithmetic)
- Boost graph library, graphviz

**Python Support:**
- Boost.Python (for Python bindings - required if BUILD_PYTHON_MODULE=ON)
- Python 3.8+ support (updated from Python 2.7)

**Optional:**
- VTK (Python visualization)
- Doxygen (documentation generation)

## Development Notes

- Compiler settings use strict warnings (`-Wall -Werror`)
- LGPL 2.1 licensed (changed from GPLv3 in 2018)
- Git-based version string generation
- Extensive test coverage across both C++ and Python interfaces

## Development Best Practices

- Whenever you make substantial changes, check that the test programs are at least mostly passing.