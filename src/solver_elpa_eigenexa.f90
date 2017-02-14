module ek_solver_elpa_eigenexa_m
  use eigen_libs
  use ELPA1
  use ELPA2
  use mpi
  use ek_global_variables_m, only : g_block_size
  use ek_distribute_matrix_m, only : process, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use ek_descriptor_parameters_m
  use ek_eigenpairs_types_m, only : eigenpairs_types_union
  use ek_event_logger_m, only : add_event
  use ek_matrix_io_m, only : sparse_mat
  use ek_processes_m, only : check_master, setup_distribution, terminate
  use ek_solver_eigenexa_m

  implicit none
  private
  public :: solve_with_general_elpa_eigenexa, solve_with_general_elpa_eigenk

contains

  subroutine solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs
    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         desc_A_re(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success
    type(eigenpairs_types_union) :: eigenpairs_tmp

    double precision :: time_start, time_start_part, time_end
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    integer :: numroc

    time_start = mpi_wtime()
    time_start_part = time_start

    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)
#if	ELPA_VERSION >= 201502002
    ierr = get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)
#else
    call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)
#endif

    max_block_size = min(n / np_rows, n / np_cols)
    block_size = min(max_block_size, g_block_size)
    na_rows = numroc(n, block_size, my_prow, 0, np_rows)
    na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
    call descinit(sc_desc, n, n, block_size, block_size, 0, 0, proc%context, na_rows, info)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:init', time_end - time_start_part)
    time_start_part = time_end

    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:setup_distributed_matrices', time_end - time_start_part)
    time_start_part = time_end

    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:distribute_global_sparse_matrices', time_end - time_start_part)
    time_start_part = time_end

    ! Return of cholesky_real is stored in the upper triangle.
#if	ELPA_VERSION >= 201502001
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, .true., success)
#else
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif
    if (.not. success) then
      call terminate('solver_main, general_elpa_eigenexa: cholesky_real failed', 1)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:cholesky_real', time_end - time_start_part)
    time_start_part = time_end

#if	ELPA_VERSION >= 201502001
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, .true., success)
#else
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif
    ! invert_trm_real always returns fail
    !if (.not. success) then
    !  if (myid == 0) then
    !    print *, 'invert_trm_real failed'
    !  end if
    !  stop
    !end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:invert_trm_real', time_end - time_start_part)
    time_start_part = time_end

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    ! This operation can be done as below:
    ! call pdtrmm('Left', 'Upper', 'Trans', 'No_unit', na, na, 1.0d0, &
    !      b, 1, 1, sc_desc, a, 1, 1, sc_desc)
    ! but it is slow. Instead use mult_at_b_real.
    call mult_at_b_real('Upper', 'Full', n, n, &
         matrix_B_dist, na_rows, matrix_A2_dist, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols, matrix_A_dist, na_rows)
    deallocate(matrix_A2_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:mult_at_b_real', time_end - time_start_part)
    time_start_part = time_end

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:pdtrmm_right', time_end - time_start_part)
    time_start_part = time_end

    call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
    call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
    deallocate(matrix_A_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:pdgemr2d_1', time_end - time_start_part)
    time_start_part = time_end

    call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:eigen_solver_eigenexa', time_end - time_start_part)
    time_start_part = time_end

    deallocate(matrix_A_redist)
    eigenpairs%type_number = 2
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa_eigenexa: allocation failed', ierr)
    end if
    eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
    eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         eigenpairs_tmp%blacs%desc(context_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:pdgemr2d_2', time_end - time_start_part)
    time_start_part = time_end

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenexa:pdtrmm_EV', time_end - time_start_part)
    call add_event('solve_with_general_elpa_eigenexa', time_end - time_start)
  end subroutine solve_with_general_elpa_eigenexa


  subroutine solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs
    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         desc_A_re(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    logical :: success
    type(eigenpairs_types_union) :: eigenpairs_tmp
    double precision :: time_start, time_start_part, time_end
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    integer :: numroc

    time_start = mpi_wtime()
    time_start_part = time_start

    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)
#if	ELPA_VERSION >= 201502002
    ierr = get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)
#else
    call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
         mpi_comm_rows, mpi_comm_cols)
#endif

    max_block_size = min(n / np_rows, n / np_cols)
    block_size = min(max_block_size, g_block_size)
    na_rows = numroc(n, block_size, my_prow, 0, np_rows)
    na_cols = numroc(n, block_size, my_pcol, 0, np_cols)
    call descinit(sc_desc, n, n, block_size, block_size, 0, 0, proc%context, na_rows, info)
    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('A2', proc, n, n, desc_A2, matrix_A2_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:setup_matrices', time_end - time_start_part)
    time_start_part = time_end

    ! Return of cholesky_real is stored in the upper triangle.
#if	ELPA_VERSION >= 201502001
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, .true., success)
#else
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif
    if (.not. success) then
      call terminate('solver_main, general_elpa_eigenk: cholesky_real failed', 1)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:cholesky_real', time_end - time_start_part)
    time_start_part = time_end

#if	ELPA_VERSION >= 201502001
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, .true., success)
#else
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif
    ! invert_trm_real always returns fail
    !if (.not. success) then
    !  if (myid == 0) then
    !    print *, 'invert_trm_real failed'
    !  end if
    !  stop
    !end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:invert_trm_real', time_end - time_start_part)
    time_start_part = time_end

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    ! This operation can be done as below:
    ! call pdtrmm('Left', 'Upper', 'Trans', 'No_unit', na, na, 1.0d0, &
    !      b, 1, 1, sc_desc, a, 1, 1, sc_desc)
    ! but it is slow. Instead use mult_at_b_real.
    call mult_at_b_real('Upper', 'Full', n, n, &
         matrix_B_dist, na_rows, matrix_A2_dist, na_rows, &
         block_size, mpi_comm_rows, mpi_comm_cols, matrix_A_dist, na_rows)
    deallocate(matrix_A2_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:mult_at_b_real', time_end - time_start_part)
    time_start_part = time_end

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:pdtrmm_right', time_end - time_start_part)
    time_start_part = time_end

    call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
    call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
    deallocate(matrix_A_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:pdgemr2d_1', time_end - time_start_part)
    time_start_part = time_end

    call eigen_solver_eigenk(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:eigen_solver_eigenk', time_end - time_start_part)
    time_start_part = time_end

    deallocate(matrix_A_redist)
    eigenpairs%type_number = 2
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa_eigenk: allocation failed', ierr)
    end if
    eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
    eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         eigenpairs_tmp%blacs%desc(context_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:pdgemr2d_2', time_end - time_start_part)
    time_start_part = time_end

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_eigenk:pdtrmm_EV', time_end - time_start_part)
    call add_event('solve_with_general_elpa_eigenk', time_end - time_start)
  end subroutine solve_with_general_elpa_eigenk
end module ek_solver_elpa_eigenexa_m
