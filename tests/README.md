# Unit tests for Feltor


We use [catch2](https://github.com/catchorg/Catch2) v3 as a test framework for our unit tests.
(On linux can be installed through package manager).
We need catch version > 3.0.

## General

What unit tests do is show compilation and correct behaviour of individual functions or class members.
Tests are supposed to be written such that
 -# a test can run automatically and **does not need any user input**
 -# there is **no need to interpret some output** to figure out whether or not the test was passed or not

We use "Black box testing" i.e. we test only through the use of the public interface.

Tests should be small and follow the "AAA" design: "Arange", "Act", "Assert".

The messages logged in INFO are printed only when using the `-s` command line option or when a test fails.

All functions and global variables in `*t.cpp` programs need to be inline.

There is a `run-all-tests.sh` script that compiles and runs all catch2 tests for all 6 device combinations (cpu, omp, gpu, mpi+cpu, mpi+omp, mpi+gpu).

## Running the tests

How to run the test suite is described in the section "Running the test suite" of the main `README.adoc`.

A test that writes files with fixed names must not run concurrently with other such tests under `ctest -j`.
Register it with the `RESOURCE_LOCK` argument of `add_dg_test` (documented in `cmake/FeltorTests.cmake`), as the tests in `inc/file` do.

### Warnings as errors

A full test should also check that the code compiles without warnings. Configure the build directory once with
```bash
cmake --preset gpu -DCMAKE_COMPILE_WARNING_AS_ERROR=ON
```
to treat compiler warnings as errors. The setting stays in the cache of `build/gpu`, so later `cmake --preset gpu` calls and the `test-gpu` workflow keep it.
The presets themselves do not set it, because warnings depend on the compiler and its version, and a build should not fail only because a different compiler warns.

### How the presets in `CMakePresets.json` are organised

CMake has four kinds of presets. Each kind stands for one step on the command line, and each kind has its own namespace, so presets of different kinds may share a name:

| Kind | Step | Decides | Presets |
|---|---|---|---|
| configure | `cmake --preset <name>` | build directory, device and compiler flags, MPI, and which tests are registered | `cpu`, `omp`, `gpu`, `mpi-cpu`, `mpi-omp`, `mpi-gpu` (and `msvc-*` on Windows) |
| build | `cmake --build --preset <name>` | which targets to compile | `test-<backend>`: `dg_tests`, `dg_benchmarks`, `feltor_projects` |
| test | `ctest --preset <name>` | how registered tests run: output of failing tests, timeout per test, error if no tests are found | `test-<backend>` |
| workflow | `cmake --workflow --preset <name>` | only the order of steps; stops at the first failure | `test-<backend>`: configure `<backend>`, build `test-<backend>`, test `test-<backend>` |

 - Only configure presets define a build. Build and test presets name a configure preset and act on its build directory `build/<backend>`.
 - The test settings `FELTOR_TEST_MPI_RANKS`, `FELTOR_TEST_OMP_THREADS` and `MPIEXEC_PREFLAGS` are cache variables of the configure presets, because they decide which tests are registered; the test presets only decide how the registered tests run. Select tests with `ctest` options such as `-L`/`-LE` (see `README.adoc`).
 - `inherits` works only between presets of the same kind: `mpi-omp` takes all cache variables of `omp` and adds MPI.
 - The number of parallel jobs is not set in the presets, because it depends on the machine; how to set it is in `README.adoc`.
 - Personal presets belong in an untracked `CMakeUserPresets.json`, which may inherit from the presets in `CMakePresets.json`.
 - The `msvc-*` presets have no test presets yet.

## Convert a program to catch

 -# Include `#include "catch2/catch_test_macros.hpp"`
 -# Replace main with `TEST_CASE( "test name")` and remove the return 0;
 -# In MPI programs remove `MPI_Init` and  `MPI_Finalize` (catch calls it for us)
 -# Replace any `assert` with `CHECK`
 -# Replace any `std::cout` with `INFO`
 -# Make global functions in cpp file static
 -# Link to `-lCatch2` in Makefile

## Some notable features
 - `GENERATE` A "for" loop for tests
 - `TEST_CASE_TEMPLATE` A static for loop for tests for various types
