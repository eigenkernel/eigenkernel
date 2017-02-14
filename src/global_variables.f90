module ek_global_variables_m
  implicit none

  public :: g_block_size, g_version, g_mpi_wtime_init
  integer :: g_block_size = 64
  character(*), parameter :: g_version = '20160808'
  real(8) :: g_mpi_wtime_init = 0d0
end module ek_global_variables_m
