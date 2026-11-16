#!/bin/bash
# Full test of the library: run the CMake workflow preset test-<preset> for
# all 6 presets, i.e. configure, compile all targets and run all tests in
# build/<preset> (see README.adoc, "Running the test suite").
# Run before pull-requests, because the CI on github cannot run the GPU tests.
#
# Usage: tests/run-all-tests.sh [jobs]
#   jobs  number of parallel compile jobs and parallel tests
#         (default: number of cores, nproc)

jobs=${1:-$(nproc)}
if ! [[ "$jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "Usage: $0 [jobs], jobs must be a positive integer, got \"$jobs\"" >&2
    exit 2
fi
# cmake --workflow has no -j option
export CMAKE_BUILD_PARALLEL_LEVEL=$jobs
export CTEST_PARALLEL_LEVEL=$jobs

cd "$(dirname "$0")/.." || exit 1 # The presets live in the feltor directory
failed=""
for preset in cpu omp gpu mpi-cpu mpi-omp mpi-gpu
do
    cmake --workflow --preset test-$preset || failed="$failed $preset"
done
if [ -n "$failed" ]; then
    echo "Failed presets:$failed"
    exit 1
fi
echo "All presets passed"
