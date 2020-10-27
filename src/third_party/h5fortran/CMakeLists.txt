add_library(
  h5fortran
  src/interface.f90
  src/attributes.f90
  src/read.f90
  src/write.f90
  src/pathlib.f90
  src/reader.f90
  src/reader_lt.f90
  src/writer.f90
  src/writer_lt.f90)

target_include_directories(h5fortran PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(h5fortran ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
