# The Python module h5py (HDF5 files from Python), e.g. for the output files of LOKI/LokiOutput. Only the presence of
# the module is detected: no compilation against HDF5 (that is the HDF5 dependency of simplebuild). Without h5py, the
# packages depending on it are not built. Found via extdep_root in simplebuild.cfg.

EXECUTE_PROCESS(COMMAND "${Python_EXECUTABLE}" "-c" "import h5py;print (h5py.__version__)"
                OUTPUT_VARIABLE tmp RESULT_VARIABLE tmp_ec ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

if ("x${tmp_ec}" STREQUAL "x0")
  set(HAS_h5py 1)
  string(STRIP "${tmp}" ExtDep_h5py_VERSION)
  set(ExtDep_h5py_COMPILE_FLAGS "")
  set(ExtDep_h5py_LINK_FLAGS "")
else()
  set(HAS_h5py 0)
endif()
