## ADDED Requirements

### Requirement: Configure the application with CMake
The project SHALL provide a root CMake project that configures a `NeutronCode` command-line executable from the current application entry point and implementation sources, with the required project header include paths and C99 language support.

#### Scenario: Configure a clean checkout
- **WHEN** a developer runs `cmake -S . -B build` from a clean checkout with a supported CMake and C compiler
- **THEN** CMake generates build files for the `NeutronCode` executable without requiring generated files from another build system

#### Scenario: Build the current source set
- **WHEN** the configured build is built
- **THEN** the executable is linked from `main.c` and the required current `.c` implementation files under `src/`

### Requirement: Support Linux and Windows GCC-compatible toolchains
The project SHALL configure and build on Linux with GCC or Clang and on Windows with MinGW-w64. The project SHALL apply compiler options only when supported by the selected compiler and SHALL link a separate math library only on platforms that require it. Native MSVC support is not required by this capability.

#### Scenario: Build on Linux
- **WHEN** a developer configures and builds with a supported Linux GCC or Clang toolchain
- **THEN** the executable links successfully with required system libraries, including the math library where separately provided

#### Scenario: Build on Windows with MinGW-w64
- **WHEN** a developer configures and builds with CMake and a MinGW-w64 C compiler on Windows
- **THEN** the executable builds without attempting to link a Unix-only `m` library

### Requirement: Support Debug and Release configurations
The CMake project SHALL support Debug and Release configurations using the configuration mechanism of the selected generator, and SHALL enable useful compiler warnings without imposing unsupported compiler-specific flags.

#### Scenario: Build Debug configuration
- **WHEN** a developer selects the Debug configuration using the chosen CMake generator
- **THEN** CMake builds the executable with Debug configuration settings

#### Scenario: Build Release configuration
- **WHEN** a developer selects the Release configuration using the chosen CMake generator
- **THEN** CMake builds the executable with Release configuration settings

### Requirement: Keep generated build products out of version control
The repository SHALL ignore the documented out-of-source build directory and common generated CMake build files, while leaving source files and user-maintained inputs trackable.

#### Scenario: Build in the documented directory
- **WHEN** a developer configures the project using `-B build` and builds it
- **THEN** generated build files and binaries under that build directory do not appear as untracked source changes

### Requirement: Document the supported CMake workflow
The README SHALL provide the minimum CMake/compiler prerequisites and exact configure/build commands for Linux and Windows using the supported toolchains, including how to select Debug and Release configurations.

#### Scenario: Follow Linux build instructions
- **WHEN** a developer follows the documented Linux instructions on a supported toolchain
- **THEN** they can configure and build the application from a clean checkout

#### Scenario: Follow Windows build instructions
- **WHEN** a developer follows the documented Windows instructions with MinGW-w64 available to CMake
- **THEN** they can configure and build the application from a clean checkout
