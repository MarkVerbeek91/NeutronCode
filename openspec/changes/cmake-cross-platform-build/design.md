## Context

The repository is a C99 command-line scientific application with `main.c`, implementation files under `src/`, and headers under `include/`. `premake5.lua` currently uses a recursive file glob, sets Unix-style flags, and links `m`; the Code::Blocks project includes stale paths. There are no existing OpenSpec capability specs or automated tests.

The source currently includes POSIX-style command-line headers (`getopt.h` and `unistd.h`). Consequently, “Windows support” for this build increment means a Windows toolchain that supplies these interfaces, such as MinGW-w64. Supporting native MSVC would require a separate portability change and is not silently included here.

## Goals / Non-Goals

**Goals:**

- Provide one root CMake entry point for configuring and building the current executable.
- Build the current source set as C99 on Linux and Windows with GCC-compatible toolchains.
- Provide Debug and Release configurations, useful warnings, correct include paths, and platform-appropriate math-library linkage.
- Keep generated build products out of version control and document exact configure/build commands.

**Non-Goals:**

- Change scientific formulas, CLI behavior, configuration parsing, or runtime output.
- Add tests, CI, installation/package targets, or a new dependency manager.
- Support native MSVC in this increment; source portability blockers are outside scope.
- Delete or update Premake and Code::Blocks files, or clean up unrelated historical artifacts.

## Decisions

### CMake target and source selection

Use a root `CMakeLists.txt` and one executable target (`NeutronCode`). List `main.c` and the `.c` files currently under `src/` explicitly rather than using an unrestricted recursive glob. This makes the build input auditable and prevents accidental inclusion of generated or unrelated C files. Add the repository root and `include/` to include paths; existing headers include nested modules by their header basename, so add the required nested include directories or adjust only the build include directories to match current includes.

**Alternative considered:** recursive globbing mirrors Premake's current `**.c` behavior but can silently change the target as files are added; explicit sources are preferred for deterministic builds.

### C language and compiler configuration

Require C99 without compiler extensions. Use CMake's configuration model for Debug and Release rather than embedding configuration-specific flags in source settings. Enable a common warning set for GCC/Clang-compatible compilers (`-Wall -Wextra -Wpedantic`) and avoid applying GNU-only options to unsupported compilers. The supported Windows path is MinGW-w64; Linux GCC and Clang are supported. Do not advertise MSVC support until POSIX header/function portability is separately addressed.

**Alternative considered:** retaining Premake as a second supported build path contradicts the project decision that CMake is the sole supported build system and leaves behavior split across configurations.

### Math library linkage

Link the C math library (`m`) on Unix-like toolchains where it is a separate library. Do not request a library named `m` on Windows, where the math functions are provided by the runtime/toolchain.

### Generated output and build instructions

Use out-of-source CMake builds (`build/`) for documented workflows and ignore that directory and common CMake-generated files in `.gitignore`. Add concise configure/build commands to the README as part of this increment, without expanding into the full documentation increment.

## Risks / Trade-offs

- [Current source may fail on some Windows environments due to POSIX headers and APIs] → Specify MinGW-w64 as the supported Windows compiler for this increment and verify an actual Windows build; defer MSVC compatibility.
- [Warnings can expose pre-existing issues and vary by compiler] → Keep warnings enabled but do not expand this build-only change into unrelated warning cleanup; distinguish warnings from build failures.
- [An incomplete explicit source list can omit a translation unit] → Compare the CMake source list against tracked `.c` files and verify the final link on Linux and Windows.
- [CMake version or generator availability differs by environment] → State a minimum CMake version and provide generator-neutral commands where possible (`cmake -S . -B build`, `cmake --build build --config ...`).

## Migration Plan

1. Add and validate the root CMake project on Linux and Windows/MinGW-w64.
2. Add a concise build section and commands to the README and ignore generated build files.
3. Keep legacy build files in place during this change; they are not part of the supported workflow after CMake is verified, but removal is a separate maintenance decision.
4. Rollback consists of reverting the added CMake configuration, README build section, and `.gitignore` additions; no source or scientific behavior migration is required.

## Open Questions

- None for this increment. Native MSVC support can be considered later as a separate source-portability change.
