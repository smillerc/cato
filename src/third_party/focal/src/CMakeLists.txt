add_library(
  focal
  Focal.f90
  Focal_Error.f90
  Focal_Memory.f90
  Focal_Profile.f90
  Focal_Setup.f90
  Focal_Debug.f90
  Focal_HostMemory.f90
  Focal_NoDebug.f90
  Focal_Query.f90
  Focal_Utils.f90)

target_link_libraries(focal clfortran fortranutils m_strings ${OpenCL_LIBRARIES})
