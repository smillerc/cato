# Use CMake to find a specific Python module
# https://cmake.org/pipermail/cmake/2011-January/041666.html
#
# Example Usage: FindPythonModule(PyQt4 REQUIRED)
#
# function(FindPythonModule module)
# 	string(TOUPPER ${module} module_upper)
# 	if(NOT PY_${module_upper})
# 		if(ARGC GREATER 1 AND ARGV1 STREQUAL "REQUIRED")
# 			set(${module}_FIND_REQUIRED TRUE)
# 		endif()
# 		# A module's location is usually a directory, but for binary modules
# 		# it's a .so file.
# 		execute_process(COMMAND "${PYTHON_EXEC}" "-c" "import ${module}"
# 			RESULT_VARIABLE _${module}_status
# 			OUTPUT_VARIABLE _${module}_location
#             ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
#         message(RESULT_VARIABLE $_${module}_status)
# 		if(NOT _${module}_status)
# 			set(PY_${module_upper} ${_${module}_location} CACHE STRING
# 				"Location of Python module ${module}")
# 		endif(NOT _${module}_status)
# 	endif(NOT PY_${module_upper})
#     find_package_handle_standard_args(PY_${module} DEFAULT_MSG PY_${module_upper})
# endfunction(FindPythonModule)
