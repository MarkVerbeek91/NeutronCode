## Why

NeutronCode currently has no single verified build path: the Premake configuration is Unix-oriented, while the Code::Blocks project contains stale source paths. A CMake build provides one documented, reproducible entry point for the existing C99 program on both Linux and Windows without mixing build modernization with scientific-model changes.

## What Changes

- Add a root CMake project that builds the existing NeutronCode command-line executable from the current source tree.
- Define C99 requirements, useful compiler warnings, configuration behavior, include paths, and platform-appropriate dependencies.
- Support Debug and Release builds on Linux and Windows.
- Ignore generated CMake/build outputs.
- Leave Premake and Code::Blocks files unchanged in this change; their removal or migration is out of scope.

## Capabilities

### New Capabilities
- `cross-platform-build`: Configure and build the existing C99 command-line application using CMake on Linux and Windows.

### Modified Capabilities
- None.

## Impact

- Adds root `CMakeLists.txt` and a capability spec under `specs/cross-platform-build/`.
- Updates `.gitignore` and may add brief CMake build instructions to `README.md` if needed to satisfy the build contract.
- Includes existing C files under `main.c` and `src/`, and headers under `include/`; no application runtime behavior or scientific formulas are intentionally changed.
- Requires CMake and a supported C compiler (for example GCC/Clang on Linux and MSVC or MinGW on Windows). Links the math library only on platforms that provide it separately.
