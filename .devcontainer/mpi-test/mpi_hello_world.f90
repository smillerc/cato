program hello_world

  use mpi_f08
  use iso_fortran_env, only: int32, int64
  implicit none

  integer(int32) :: ierr, my_rank, size

  call mpi_init(ierr)

  call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call mpi_comm_size(MPI_COMM_WORLD, size, ierr)

  print *
  write(*, '(2(a,i2))') 'Hello Fortran World! I am rank ', my_rank, ' of size ', size
  print *

  call mpi_finalize(ierr)

end
