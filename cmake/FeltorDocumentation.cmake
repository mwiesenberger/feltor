# Generate and check the Doxygen documentation of the FELTOR libraries.
#
# Usage:
#   cmake -DDOXYGEN=<doxygen> -DFELTOR_SOURCE_DIR=<feltor> \
#         -DOUTPUT_DIR=<dir> [options] -P FeltorDocumentation.cmake
#
# Options:
#   -DFAIL_ON_WARNINGS=ON|OFF  ON (default): every run uses WARN_AS_ERROR =
#                              FAIL_ON_WARNINGS and the script fails at the end
#                              if any run warned. OFF: warnings are only printed.
#   -DHAVE_DOT=<bool>          Doxygen's HAVE_DOT: draw class and include graphs
#                              with graphviz. Default NO.
#   -DDOT_PATH=<dir>           Doxygen's DOT_PATH: directory of dot.
#                              Default empty: doxygen searches the PATH.
#   -DAWESOME_DIR=<dir>        doxygen-awesome-css directory: use the theme
#                              of the published documentation (doc/header.html,
#                              doc/footer.html, doc/menubar.css).
#                              Empty (default): plain Doxygen HTML.
#
# Runs doxygen on the Doxyfiles of exblas, dg, geometries, matrix and file in
# the order of doc/Makefile (tag files first, so that cross references
# resolve). Each library's HTML goes to OUTPUT_DIR/<library>/html, and
# OUTPUT_DIR/index.html links to the dg topics page. Every run includes the
# library's own Doxyfile and overrides only the settings above, the output
# location and the tag files. All runs execute even if one warns.

foreach(var DOXYGEN FELTOR_SOURCE_DIR OUTPUT_DIR)
    if(NOT DEFINED ${var})
        message(FATAL_ERROR "FeltorDocumentation.cmake: ${var} is not set")
    endif()
endforeach()
if(NOT DEFINED FAIL_ON_WARNINGS)
    set(FAIL_ON_WARNINGS ON)
endif()
if(NOT DEFINED HAVE_DOT)
    set(HAVE_DOT NO)
endif()

# Settings common to all runs
set(common_ "QUIET = YES" "WARN_LOGFILE ="
    "HAVE_DOT = ${HAVE_DOT}" "DOT_PATH = ${DOT_PATH}")
if(FAIL_ON_WARNINGS)
    list(APPEND common_ "WARN_AS_ERROR = FAIL_ON_WARNINGS")
else()
    list(APPEND common_ "WARN_AS_ERROR = NO")
endif()
if(AWESOME_DIR)
    if(NOT EXISTS "${AWESOME_DIR}/doxygen-awesome.css")
        message(FATAL_ERROR "FeltorDocumentation.cmake: no doxygen-awesome.css in AWESOME_DIR=${AWESOME_DIR}")
    endif()
    # Sidebar-only layout and fragment copy button as in the doxygen-awesome
    # README; doc/header.html loads and initialises the copy button
    set(doc_ "${FELTOR_SOURCE_DIR}/doc")
    list(APPEND common_
        "HTML_HEADER = ${doc_}/header.html"
        "HTML_FOOTER = ${doc_}/footer.html"
        "HTML_EXTRA_STYLESHEET = ${AWESOME_DIR}/doxygen-awesome.css"
        "HTML_EXTRA_STYLESHEET += ${AWESOME_DIR}/doxygen-awesome-sidebar-only.css"
        "HTML_EXTRA_STYLESHEET += ${doc_}/menubar.css"
        "HTML_EXTRA_FILES = ${AWESOME_DIR}/doxygen-awesome-fragment-copy-button.js"
        "HTML_COPY_CLIPBOARD = NO"
        "HTML_COLORSTYLE = LIGHT")
endif()

set(failed_runs_ "")

# run_doxygen(<name> <library dir relative to inc/> [<Doxyfile line>...])
function(run_doxygen name dir)
    set(out "${OUTPUT_DIR}/${name}")
    file(REMOVE_RECURSE "${out}")
    file(MAKE_DIRECTORY "${out}")
    string(JOIN "\n" lines
        "@INCLUDE = ${FELTOR_SOURCE_DIR}/inc/${dir}/Doxyfile"
        "OUTPUT_DIRECTORY = ${out}"
        ${common_} ${ARGN})
    file(WRITE "${out}/Doxyfile" "${lines}\n")
    message(STATUS "Doxygen: ${name}")
    # Relative paths in the Doxyfiles (INPUT, EXAMPLE_PATH) refer to the library directory
    execute_process(COMMAND "${DOXYGEN}" "${out}/Doxyfile"
        WORKING_DIRECTORY "${FELTOR_SOURCE_DIR}/inc/${dir}"
        RESULT_VARIABLE result)
    if(NOT result EQUAL 0)
        set(failed_runs_ ${failed_runs_} ${name} PARENT_SCOPE)
    endif()
endfunction()

set(T "${OUTPUT_DIR}")
run_doxygen(exblas dg/backend/exblas
    "GENERATE_TAGFILE = ${T}/exblas.tag")
run_doxygen(dg_tag dg
    "GENERATE_HTML = NO"
    "GENERATE_TAGFILE = ${T}/dg.tag")
run_doxygen(geometries geometries
    "TAGFILES = ${T}/dg.tag=../../dg/html"
    "GENERATE_TAGFILE = ${T}/geometries.tag")
run_doxygen(matrix matrix
    "TAGFILES = ${T}/dg.tag=../../dg/html"
    "GENERATE_TAGFILE = ${T}/matrix.tag")
run_doxygen(file file)
run_doxygen(dg dg
    "TAGFILES = ${T}/exblas.tag=../../exblas/html ${T}/geometries.tag=../../geometries/html ${T}/matrix.tag=../../matrix/html")

file(CREATE_LINK dg/html/topics.html "${OUTPUT_DIR}/index.html" SYMBOLIC)

if(failed_runs_)
    message(FATAL_ERROR "Doxygen warnings or errors (see above) in: ${failed_runs_}")
endif()
message(STATUS "Doxygen: open ${OUTPUT_DIR}/index.html")
