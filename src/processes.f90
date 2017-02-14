module ek_processes_m
  use mpi
  use ek_event_logger_m
  implicit none

  type ek_process_t
    integer :: my_rank, n_procs, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
  end type ek_process_t

  private
  public :: ek_process_t, setup_distribution, get_num_procs, layout_procs, &
       print_map_of_grid_to_processes, check_master, terminate

contains

  subroutine setup_distribution(proc)
    type(ek_process_t), intent(out) :: proc

    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call layout_procs(proc%n_procs, proc%n_procs_row, proc%n_procs_col)
    call blacs_get(-1, 0, proc%context)
    call blacs_gridinit(proc%context, 'R', proc%n_procs_row, proc%n_procs_col)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    if (proc%my_rank == 0) then
      print '("BLACS process grid: ", I0, " x ", I0, " (", I0, ")")', &
           proc%n_procs_row, proc%n_procs_col, proc%n_procs
    end if

    if (proc%my_proc_row >= proc%n_procs_row .or. proc%my_proc_col >= proc%n_procs_col) then
       call blacs_exit(0)
       stop '[Warning] setup_distribution: Out of process grid, process exit'
    end if
  end subroutine setup_distribution


  subroutine get_num_procs(num_mpi_procs, num_omp_procs)
    !$ use omp_lib

    integer, intent(out) :: num_mpi_procs, num_omp_procs

    integer :: ierr

    call mpi_comm_size(mpi_comm_world, num_mpi_procs, ierr)
    if (ierr /= 0) then
      call terminate('get_num_procs: mpi_comm_size failed', ierr)
    end if

    num_omp_procs = 1
    !$ num_omp_procs = omp_get_max_threads()
  end subroutine get_num_procs


  subroutine layout_procs(n_procs, n_procs_row, n_procs_col)
    integer, intent(in) :: n_procs
    integer, intent(out) :: n_procs_row, n_procs_col

    n_procs_row = int(sqrt(dble(n_procs + 1)))
    do while (mod(n_procs, n_procs_row) /= 0)
      n_procs_row = n_procs_row - 1
    end do
    n_procs_col = n_procs / n_procs_row
  end subroutine layout_procs


  subroutine map_grid_to_processes(context, num_procs_row, num_procs_col, map)
    integer, intent(in) :: context, num_procs_row, num_procs_col
    integer, intent(out) :: map(num_procs_row, num_procs_col)

    integer :: i, j
    integer :: blacs_pnum ! function

    do j = 0, num_procs_col - 1
      do i = 0, num_procs_row - 1
        map(i + 1, j + 1) = blacs_pnum(context, i, j)
      end do
    end do
  end subroutine map_grid_to_processes


  subroutine print_map_of_grid_to_processes()
    integer :: context, num_procs_row, num_procs_col, proc_row, proc_col
    integer :: my_rank, ierr
    integer, allocatable :: map(:, :)

    call blacs_get(-1, 0, context) ! Get default system context

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (ierr /= 0) then
      call terminate('print_map_of_grid_to_processes: mpi_comm_rank failed', ierr)
    end if
    if (my_rank /= 0) return

    call blacs_gridinfo(context, num_procs_row, num_procs_col, proc_row, proc_col)
    allocate(map(num_procs_row, num_procs_col))
    call map_grid_to_processes(context, num_procs_row, num_procs_col, map)

    print '("process numbers in BLACS grid is")'
    do proc_row = 1, num_procs_row
      do proc_col = 1, num_procs_col - 1
        write (*, '(i6, " ")', advance = 'no') map(proc_row, proc_col)
      end do
      print '(i6)', map(proc_row, num_procs_col)
    end do
  end subroutine print_map_of_grid_to_processes


  logical function check_master()
    integer :: my_rank, ierr

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    if (ierr /= 0) then
      call terminate('check_master: mpi_comm_rank failed', ierr)
    end if

    check_master = (my_rank == 0)
  end function check_master


  subroutine terminate(error_message, error_code)
    character(*), intent(in) :: error_message
    integer, intent(in) :: error_code

    integer :: ierr

    ! Print events before exit.
    if (check_master()) then
      call print_events()
    end if

    if (error_code == 0) then
      write (0, '("[Info] ", a)') error_message
    else
      write (0, '("[Error] ", a)') error_message
    end if
    call mpi_abort(mpi_comm_world, error_code, ierr)
  end subroutine terminate
end module ek_processes_m
