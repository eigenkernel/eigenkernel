program eigbench
  use mpi
  use fson
  use fson_value_m
  use fson_string_m
  use ek_distribute_matrix_m
  use ek_event_logger_m
  use ek_global_variables_m
  use ek_solver_main_m, only : eigen_solver
  use ek_command_argument_m, only : ek_argument_t, required_memory, &
       read_command_argument, validate_argument, print_command_argument, fson_setting_add
  use ek_matrix_io_m, only : ek_sparse_mat_t, read_matrix_file, print_eigenvectors
  use ek_processes_m, only : get_num_procs, check_master
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t
  use ek_verifier_m, only : eval_residual_norm, eval_orthogonality
  implicit none

  type(ek_argument_t) :: arg
  type(ek_sparse_mat_t) :: matrix_A, matrix_B
  type(ek_eigenpairs_types_union_t) :: eigenpairs
  double precision :: A_norm, rn_ave, rn_max, orthogonality
  double precision, allocatable :: ipratios(:)
  integer :: num_mpi_procs, num_omp_procs, j, ierr, ierr_mpi
  integer, parameter :: iunit = 10
  double precision :: time_start, time_start_part, time_end
  type(fson_value), pointer :: output
  type(ek_process_t) :: proc

  call mpi_init(ierr)
  if (ierr /= 0) then
    write (0, *) '[Error] eigen_test: mpi_init failed, error code is ', ierr
    stop
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  time_start = mpi_wtime()
  g_mpi_wtime_init = time_start
  time_start_part = time_start

  call read_command_argument(arg)

  if (check_master()) then
    print '("---------- Eigen Test start ----------")'
    print '("----- Configurations -----")'
    call print_command_argument(arg)
    print '("approximate required memory per process (Mbytes): ", f10.1)', &
         required_memory(arg) / real(2 ** 20)
    call get_num_procs(num_mpi_procs, num_omp_procs)
    print '("MPI processes: ", i0)', num_mpi_procs
    print '("OpenMP threads per process (may be inaccurate): ", i0)', num_omp_procs
  end if

  time_end = mpi_wtime()
  call add_event('main:read_command_argument', time_end - time_start_part)
  time_start_part = time_end

  call validate_argument(arg)
  output => fson_value_create()
  output%value_type = TYPE_OBJECT
  call fson_setting_add(arg, output)

  if (check_master()) then
    call read_matrix_file(arg%matrix_A_filename, arg%matrix_A_info, matrix_A, ierr)
  end if
  call mpi_bcast(ierr, 1, mpi_integer, 0, mpi_comm_world, ierr_mpi)
  if (ierr /= 0) then
    call terminate('read_matrix_file ' // trim(arg%matrix_A_filename) // ' failed',  ierr)
  end if

  if (arg%is_generalized_problem) then
    if (check_master()) then
      call read_matrix_file(arg%matrix_B_filename, arg%matrix_B_info, matrix_B, ierr)
    end if
    call mpi_bcast(ierr, 1, mpi_integer, 0, mpi_comm_world, ierr_mpi)
    if (ierr /= 0) then
      call terminate('read_matrix_file ' // trim(arg%matrix_B_filename) // ' failed',  ierr)
    end if
  end if

  time_end = mpi_wtime()
  call add_event('main:read_matrix_files', time_end - time_start_part)
  time_start_part = time_end

  call bcast_sparse_matrix(0, arg%matrix_A_info, matrix_A)
  if (arg%is_generalized_problem) then
    call bcast_sparse_matrix(0, arg%matrix_B_info, matrix_B)
  end if

  if (arg%is_dry_run) then
    if (check_master()) print '(/, "dry run mode, exit")'
    call mpi_finalize(ierr)
    stop
  end if

  time_end = mpi_wtime()
  call add_event('main:bcast_sparse_matrices', time_end - time_start_part)
  time_start_part = time_end

  if (check_master()) print '(/, "----- Solver Call -----")'
  if (arg%is_generalized_problem) then
    call eigen_solver(arg, matrix_A, eigenpairs, proc, matrix_B)
  else
    call eigen_solver(arg, matrix_A, eigenpairs, proc)
  end if

  time_end = mpi_wtime()
  call add_event('main:eigen_solver', time_end - time_start_part)
  time_start_part = time_end

  ! Print eigenvalues and eigenvectors if required
  if (check_master()) then
    open(iunit, file=arg%output_filename, status='unknown')
    do j=1,arg%n_vec
      if (eigenpairs%type_number == 1) then
        write (iunit, '(I8, " ", E26.16e3)') j, eigenpairs%local%values(j)
      else if (eigenpairs%type_number == 2) then
        write (iunit, '(I8, " ", E26.16e3)') j, eigenpairs%blacs%values(j)
      end if
    enddo
    close(iunit)
  end if

  if (arg%num_printed_vecs_ranges /= 0) then
    call print_eigenvectors(arg, eigenpairs)
  end if

  time_end = mpi_wtime()
  call add_event('main:print_eigenpairs', time_end - time_start_part)
  time_start_part = time_end

  allocate(ipratios(eigenpairs%blacs%desc(cols_)))
  if (arg%is_generalized_problem) then
    call get_ipratios(proc, eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, ipratios, matrix_B)
  else
    call get_ipratios(proc, eigenpairs%blacs%Vectors, eigenpairs%blacs%desc, ipratios)
  end if
  if (check_master()) then
    open(iunit, file=arg%ipratios_filename, status='unknown')
    do j = 1, eigenpairs%blacs%desc(cols_)
      write (iunit, '(I8, " ", E26.16e3)') j, ipratios(j)
    enddo
    close(iunit)
  end if

  time_end = mpi_wtime()
  call add_event('main:compute_and_print_ipratios', time_end - time_start_part)
  time_start_part = time_end

  if (arg%n_check_vec /= 0) then
    if (check_master()) print '(/, "----- Checker Call -----")'
    if (arg%is_generalized_problem) then
      call eval_residual_norm(arg, matrix_A, eigenpairs, &
           A_norm, rn_ave, rn_max, matrix_B)
    else
      call eval_residual_norm(arg, matrix_A, eigenpairs, &
           A_norm, rn_ave, rn_max)
    end if

    if (check_master()) then
      print '("A norm: ", e15.8)', A_norm
      print '("residual norm (average): ", e15.8)', rn_ave
      print '("residual norm (max):     ", e15.8)', rn_max
    end if
  end if

  time_end = mpi_wtime()
  call add_event('main:eval_residual_norm', time_end - time_start_part)
  time_start_part = time_end

  if (arg%ortho_check_index_start /= 0) then
    if (arg%is_generalized_problem) then
      call eval_orthogonality(arg, eigenpairs, orthogonality, matrix_B)
    else
      call eval_orthogonality(arg, eigenpairs, orthogonality)
    end if
    if (check_master()) then
      print '("orthogonality criterion: ", e15.8)', orthogonality
    end if
  end if

  time_end = mpi_wtime()
  call add_event('main:eval_orthogonality', time_end - time_start_part)
  call add_event('main', time_end - time_start)

  if (check_master()) then
    call fson_events_add(output)
    open(iunit, file=trim(arg%log_filename), status='unknown')
    call fson_print(iunit, output)
    close(iunit)
  end if

  call mpi_finalize(ierr)
end program eigbench
