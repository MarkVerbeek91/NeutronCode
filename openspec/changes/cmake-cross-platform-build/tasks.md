## 1. Define the CMake target

- [ ] 1.1 Inventory tracked C implementation files and header include paths, then add a root `CMakeLists.txt` with an explicit source list for the `NeutronCode` executable.
- [ ] 1.2 Configure C99 without extensions, set the minimum supported CMake version, and enable compiler-appropriate warning options.
- [ ] 1.3 Add platform-conditional math-library linkage so Unix-like builds link `m` where required and Windows MinGW-w64 builds do not.
- [ ] 1.4 Configure Debug and Release builds through CMake generator configuration settings.

## 2. Ignore generated files and document usage

- [ ] 2.1 Update `.gitignore` for the documented out-of-source `build/` directory and common generated CMake files without ignoring source inputs.
- [ ] 2.2 Add concise README prerequisites and Linux configure/build commands, including Debug and Release selection.
- [ ] 2.3 Add concise README Windows/MinGW-w64 configure/build commands, including Debug and Release selection.

## 3. Verify supported builds

- [ ] 3.1 Configure and build Debug and Release on Linux with GCC or Clang; resolve build configuration defects without changing scientific behavior.
- [ ] 3.2 Configure and build Debug and Release on Windows with MinGW-w64; resolve build configuration defects without expanding scope to native MSVC.
- [ ] 3.3 Confirm the documented clean-checkout commands work and generated build outputs are ignored.
- [ ] 3.4 Run `git diff --check` and verify no Premake, Code::Blocks, scientific formula, or unrelated cleanup changes were introduced.
