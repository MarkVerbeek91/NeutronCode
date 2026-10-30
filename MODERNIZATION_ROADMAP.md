# NeutronCode Modernization Roadmap

This roadmap breaks modernization into small, independently specifiable OpenSpec changes. Implement each increment as a separate change with its own proposal, requirements, tasks, and review.

## Project decisions

- **Supported platforms:** Linux and Windows.
- **Build system:** CMake only. Legacy project files may be reviewed for removal in a later, separate increment.
- **Scientific status:** The model is theoretical at this point. There are no trusted sample outputs or expected neutron-production values to use as regression fixtures.
- **Scientific integrity:** Keep project scaffolding and scientific formula changes separate. Do not assert unverified model outputs as expected results.

## Current project observations

- The project is a modular C99 command-line application, with implementation files in `src/` and public headers in `include/`.
- `README.md` currently contains only a brief description; there is no evident automated test or CI setup.
- `premake5.lua` and `enviroment/NeutronCalculation.cbp` describe older build workflows. The Code::Blocks project references stale `scr/` paths.
- `main.c` appears to pass `&opt` to `fopen` and `arg_parameter_check` instead of the argument value supplied by `getopt` (`optarg`).
- `src/Auxiliary.c` parses configuration using fixed-order token reads, and uses unbounded `%s` conversions.
- `src/MathFunctions.c` compares a floating-point value to `NAN` with `==`, which does not detect NaN in C.
- Existing helper and generated files should be reviewed before deletion; some scripts reference old paths or external plotting/profiling tools.

## Increment 1 — CMake build for Linux and Windows

**Scope**

- Introduce a root `CMakeLists.txt` as the sole supported build entry point for the existing C99 executable.
- Include the current implementation sources and project headers, apply C99 settings and useful compiler warnings, and link required platform libraries (including `m` where applicable).
- Support Debug and Release configurations using the platform-appropriate CMake workflow.
- Add build directories and generated CMake artifacts to `.gitignore`.
- Leave Premake and Code::Blocks artifacts untouched in this increment; review their disposition separately after the CMake build is verified.

**Acceptance criteria**

- A clean checkout configures and builds the executable on Linux and Windows with CMake.
- Build instructions and exact configure/build commands are recorded for both platforms.
- Generated build outputs remain outside tracked source files.

## Increment 2 — README and example configuration

**Scope**

- Expand `README.md` to explain the application, its scientific reference, prerequisites, CMake build, invocation, configuration format, output behavior, and platform-specific notes.
- Rename `input_tempate.ini` to a correctly named template and update references.
- Explain the template's keys, booleans, and units where established by the code or source documentation; explicitly mark unknown or uncertain details rather than guessing.
- State that the model is theoretical and no trusted expected-output fixtures are currently available.

**Acceptance criteria**

- A new user can follow the README to build the project and locate the example input.
- The documented example matches the checked-in template.
- Documentation does not claim unverified scientific accuracy or output values.

## Increment 3 — Correct command-line argument handling

**Scope**

- Use the argument value supplied by `getopt` for `-f` and `-p` instead of the address of the option character.
- Define behavior for default input, missing option values, unknown options, and input-open failures.
- Report command-line and file errors to stderr and return a nonzero exit status on failure.
- Keep this change focused on the CLI; do not redesign configuration parsing in this increment.

**Acceptance criteria**

- `-f` opens the selected path and `-p` passes the selected parameter string to the existing parameter handler.
- Missing values, unknown options, and unreadable files fail with actionable diagnostics and a nonzero exit status.
- Default-input behavior remains documented and intentional.

## Increment 4 — Safe, predictable configuration parsing

**Scope**

- Replace fixed-order token scanning with line-oriented parsing for the documented configuration format.
- Define whether supported keys can appear in any order and how blank lines, comments, duplicate keys, unknown keys, and legacy spelling are handled.
- Bound string handling; validate numeric and boolean syntax; check read/open/write operations; report line- and key-specific errors.
- Specify physical-range validation only when justified by the model equations or cited source; do not invent scientific constraints.

**Acceptance criteria**

- The checked-in template parses successfully.
- Supported key-order behavior is covered by tests.
- Malformed, truncated, and unsupported entries fail safely with actionable diagnostics.
- Any retained legacy syntax is explicitly documented and tested.

## Increment 5 — Numerical utility correctness

**Scope**

- Add focused analytic tests for numerical integration and interpolation, including boundaries and invalid inputs where behavior is defined.
- Detect non-finite values using C math predicates such as `isnan` and `isfinite`, not equality comparisons against `NAN`.
- Define and implement behavior for invalid integration bounds, non-finite function results, and recursion limits.
- Keep this change limited to generic numerical utilities; do not alter scientific formulas.

**Acceptance criteria**

- Analytic cases with known mathematical results pass within documented floating-point tolerances.
- Edge cases follow explicit behavior and cannot cause unbounded recursion.
- Tests are runnable through CTest on Linux and Windows.

## Increment 6 — Testable model interface and theoretical validation

**Scope**

- Identify a narrow way to exercise core calculations without starting the CLI, plotting, or unrelated file output.
- Add configuration-to-calculation tests and model properties/invariants only when they follow from documented equations or cited sources.
- Record assumptions, units, and equation/source references alongside each scientific assertion.
- Do not add expected NPR/current numbers until authoritative, independently verified values become available.

**Acceptance criteria**

- Core calculations can be exercised in automated tests without invoking plotting.
- Every asserted model property has a documented theoretical basis.
- Test names and documentation distinguish mathematical properties from unverified expected model outputs.

## Increment 7 — Cross-platform CI and maintenance follow-up

**Scope**

- Add CI jobs to configure, build, and run CTest with CMake on Linux and Windows.
- Document incremental formatting and static-analysis expectations; address findings in scoped changes.
- Separately review old Code::Blocks/Premake files, generated dependency/index/backup artifacts, and profiling/gnuplot scripts.
- Remove obsolete items only after confirming they are not needed; update or retire broken helper scripts and document retained workflows.

**Acceptance criteria**

- CI builds and runs the available tests on both supported platforms.
- Retained helper scripts have a documented purpose and work with the current project layout, or are clearly marked as optional external-tool workflows.
- Artifact cleanup does not remove source data or user-maintained files.

## Suggested dependency order

1. CMake build.
2. README and example configuration.
3. CLI argument handling.
4. Configuration parsing.
5. Numerical utility tests and fixes.
6. Testable model interface and theoretically justified validation.
7. Cross-platform CI and maintenance follow-up.

OpenSpec changes may be reordered if dependencies dictate, but keep each change focused and preserve the boundary between generic project modernization and scientific-model changes.
