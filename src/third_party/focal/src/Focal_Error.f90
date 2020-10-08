! -----------------------------------------------------------------------------
!  FOCAL
!
!   A modern Fortran abstraction layer for OpenCL
!   https://lkedward.github.io/focal-docs
!
! -----------------------------------------------------------------------------
!
! Copyright (c) 2020 Laurence Kedward
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! -----------------------------------------------------------------------------

submodule(Focal) Focal_Error
  !!  Implementation module for error handling

  !! @note This is an implementation submodule: it contains the code implementing the subroutines defined in the
  !!  corresponding header module file. See header module file (Focal.f90) for interface definitions. @endnote

use clfortran
implicit none

interface
    !! Interface to c function abort().
    !!  Used to print backtrace on error.
  subroutine c_abort() bind(C, name="abort")
  end subroutine
end interface

contains

module procedure fclDefaultErrorHandler

if(errcode /= CL_SUCCESS) then

  write(*, *) '(!) Fatal openCl error ', errcode, ' : ', trim(fclGetErrorString(errcode))
  write(*, *) '      at ', focalCall, ':', oclCall

  call c_abort()
end if

end procedure fclDefaultErrorHandler
! ---------------------------------------------------------------------------

module procedure fclHandleBuildError !(builderrcode,prog,ctx)

integer :: i

! Handle compilation error
if(builderrcode /= CL_SUCCESS) then

  write(*, *) '(!) Fatal openCl error while building kernel: ', builderrcode, ' : ', trim(fclGetErrorString(builderrcode))

  ! Iterate over context devices
  do i = 1, ctx%platform%numDevice

    call fclDumpBuildLog(ctx, prog, ctx%platform%devices(i))

  end do

  stop 1
end if

end procedure fclHandleBuildError
! ---------------------------------------------------------------------------

module procedure fclGetErrorString !(errcode)

select case(errcode)
case(CL_DEVICE_NOT_FOUND)
  errstr = 'CL_DEVICE_NOT_FOUND'

case(CL_DEVICE_NOT_AVAILABLE)
  errstr = 'CL_DEVICE_NOT_AVAILABLE'

case(CL_COMPILER_NOT_AVAILABLE)
  errstr = 'CL_COMPILER_NOT_AVAILABLE'

case(CL_MEM_OBJECT_ALLOCATION_FAILURE)
  errstr = 'CL_MEM_OBJECT_ALLOCATION_FAILURE'

case(CL_OUT_OF_RESOURCES)
  errstr = 'CL_OUT_OF_RESOURCES'

case(CL_OUT_OF_HOST_MEMORY)
  errstr = 'CL_OUT_OF_HOST_MEMORY'

case(CL_PROFILING_INFO_NOT_AVAILABLE)
  errstr = 'CL_PROFILING_INFO_NOT_AVAILABLE'

case(CL_MEM_COPY_OVERLAP)
  errstr = 'CL_MEM_COPY_OVERLAP'

case(CL_BUILD_PROGRAM_FAILURE)
  errstr = 'CL_BUILD_PROGRAM_FAILURE'

case(CL_MAP_FAILURE)
  errstr = 'CL_MAP_FAILURE'

case(CL_MISALIGNED_SUB_BUFFER_OFFSET)
  errstr = 'CL_MISALIGNED_SUB_BUFFER_OFFSET'

case(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
  errstr = 'CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST'

case(CL_COMPILE_PROGRAM_FAILURE)
  errstr = 'CL_COMPILE_PROGRAM_FAILURE'

case(CL_LINKER_NOT_AVAILABLE)
  errstr = 'CL_LINKER_NOT_AVAILABLE'

case(CL_LINK_PROGRAM_FAILURE)
  errstr = 'CL_LINK_PROGRAM_FAILURE'

case(CL_DEVICE_PARTITION_FAILED)
  errstr = 'CL_DEVICE_PARTITION_FAILED'

case(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
  errstr = 'CL_KERNEL_ARG_INFO_NOT_AVAILABLE'

case(CL_INVALID_VALUE)
  errstr = 'CL_INVALID_VALUE'

case(CL_INVALID_DEVICE_TYPE)
  errstr = 'CL_INVALID_DEVICE_TYPE'

case(CL_INVALID_PLATFORM)
  errstr = 'CL_INVALID_PLATFORM'

case(CL_INVALID_DEVICE)
  errstr = 'CL_INVALID_DEVICE'

case(CL_INVALID_CONTEXT)
  errstr = 'CL_INVALID_CONTEXT'

case(CL_INVALID_QUEUE_PROPERTIES)
  errstr = 'CL_INVALID_QUEUE_PROPERTIES'

case(CL_INVALID_COMMAND_QUEUE)
  errstr = 'CL_INVALID_COMMAND_QUEUE'

case(CL_INVALID_HOST_PTR)
  errstr = 'CL_INVALID_HOST_PTR'

case(CL_INVALID_MEM_OBJECT)
  errstr = 'CL_INVALID_MEM_OBJECT'

case(CL_INVALID_BINARY)
  errstr = 'CL_INVALID_BINARY'

case(CL_INVALID_BUILD_OPTIONS)
  errstr = 'CL_INVALID_BUILD_OPTIONS'

case(CL_INVALID_PROGRAM)
  errstr = 'CL_INVALID_PROGRAM'

case(CL_INVALID_PROGRAM_EXECUTABLE)
  errstr = 'CL_INVALID_PROGRAM_EXECUTABLE'

case(CL_INVALID_KERNEL_NAME)
  errstr = 'CL_INVALID_KERNEL_NAME'

case(CL_INVALID_KERNEL_DEFINITION)
  errstr = 'CL_INVALID_KERNEL_DEFINITION'

case(CL_INVALID_KERNEL)
  errstr = 'CL_INVALID_KERNEL'

case(CL_INVALID_ARG_INDEX)
  errstr = 'CL_INVALID_ARG_INDEX'

case(CL_INVALID_ARG_VALUE)
  errstr = 'CL_INVALID_ARG_VALUE'

case(CL_INVALID_ARG_SIZE)
  errstr = 'CL_INVALID_ARG_SIZE'

case(CL_INVALID_KERNEL_ARGS)
  errstr = 'CL_INVALID_KERNEL_ARGS'

case(CL_INVALID_WORK_DIMENSION)
  errstr = 'CL_INVALID_WORK_DIMENSION'

case(CL_INVALID_WORK_GROUP_SIZE)
  errstr = 'CL_INVALID_WORK_GROUP_SIZE'

case(CL_INVALID_WORK_ITEM_SIZE)
  errstr = 'CL_INVALID_WORK_ITEM_SIZE'

case(CL_INVALID_GLOBAL_OFFSET)
  errstr = 'CL_INVALID_GLOBAL_OFFSET'

case(CL_INVALID_EVENT_WAIT_LIST)
  errstr = 'CL_INVALID_EVENT_WAIT_LIST'

case(CL_INVALID_EVENT)
  errstr = 'CL_INVALID_EVENT'

case(CL_INVALID_OPERATION)
  errstr = 'CL_INVALID_OPERATION'

case(CL_INVALID_BUFFER_SIZE)
  errstr = 'CL_INVALID_BUFFER_SIZE'

case(CL_INVALID_GLOBAL_WORK_SIZE)
  errstr = 'CL_INVALID_GLOBAL_WORK_SIZE'

case(CL_INVALID_PROPERTY)
  errstr = 'CL_INVALID_PROPERTY'

case(CL_INVALID_COMPILER_OPTIONS)
  errstr = 'CL_INVALID_COMPILER_OPTIONS'

case(CL_INVALID_LINKER_OPTIONS)
  errstr = 'CL_INVALID_LINKER_OPTIONS'

case(CL_INVALID_DEVICE_PARTITION_COUNT)
  errstr = 'CL_INVALID_DEVICE_PARTITION_COUNT'

case(CL_PLATFORM_NOT_FOUND_KHR)
  errstr = 'CL_PLATFORM_NOT_FOUND_KHR'

case(NV_ILLEGAL_BUFFER_READ_WRITE)
  errstr = 'NVidia: Illegal read or write to a buffer'

case default
  errstr = 'UNKNOWN'

end select

end procedure fclGetErrorString
! ---------------------------------------------------------------------------

module procedure fclRuntimeError !(descrip)

write(*, *) '(!) Fatal runtime error: an incorrect Focal program has been written.'
if(present(descrip)) then
  write(*, *) '      at ', descrip
end if

call c_abort()

end procedure fclRuntimeError
! ---------------------------------------------------------------------------

end submodule Focal_Error
