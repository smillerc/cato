if(DEVICE_DETECT_DEBUG)
  add_definitions(-DDEVICE_DETECT_DEBUG=1)
endif(DEVICE_DETECT_DEBUG)

add_subdirectory(external)
add_subdirectory(src)
