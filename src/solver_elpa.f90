module ek_solver_elpa_m
  use elpa
  use ELPA1
  use ELPA2
  use mpi
  use ek_global_variables_m, only : g_block_size
  use ek_distribute_matrix_m, only : ek_process_t, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use ek_descriptor_parameters_m
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t
  use ek_event_logger_m, only : add_event
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_processes_m, only : check_master, terminate

  implicit none
  private
  public :: solve_with_general_elpa_scalapack, solve_with_general_elpa1, solve_with_general_elpa2

contains

  subroutine solve_with_general_elpa_scalapack(n, proc, matrix_A, eigenpairs, matrix_B)
    class(elpa_t), pointer :: e
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         pdsyevd_lwork, pdsyevd_liwork, pdsyevd_trilwmin, &
         sc_desc(desc_size), ierr, info, mpierr

    integer, allocatable :: pdsyevd_iwork(:)
    double precision :: time_start, time_end
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :), pdsyevd_work(:)
    integer :: numroc, success

    time_start = mpi_wtime()

    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)

#if	ELPA_VERSION >= 201705003
    ierr = elpa_get_communicators(mpi_comm_world, my_prow, my_pcol, &
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
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa_scalapack: allocation failed', ierr)
    end if
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:setup_matrices', time_end - time_start)
    time_start = time_end

    if (elpa_init(20170403) /= elpa_ok) then
      print *, "ELPA API version not supported"
      stop
    endif
    e => elpa_allocate()

    call e%set("na", n, success) !Size of Matrix
    call e%set("nev", n, success) !Number of eigenvectors to be calculated
    call e%set("local_nrows", na_rows, success)
    call e%set("local_ncols", na_cols, success)
    call e%set("nblk", block_size, success)
    call e%set("mpi_comm_parent", mpi_comm_world, success)
    call e%set("process_row", my_prow, success)
    call e%set("process_col", my_pcol, success)
    success = e%setup()

    call e%set("solver", elpa_solver_1stage, success)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:setup_parameter', time_end - time_start)
    time_start = time_end

    ! Return of cholesky_real is stored in the upper triangle.
#if	ELPA_VERSION >= 201705003
    call e%cholesky(matrix_B_dist, success)
#else
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa_scalapack: cholesky_real failed', ierr)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:cholesky_real', time_end - time_start)
    time_start = time_end

#if	ELPA_VERSION >= 201705003
    call e%invert_triangular(matrix_B_dist, success)
#else
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa_scalapack: invert_trm_real failed', 1)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:invert_trm_real', time_end - time_start)
    time_start = time_end

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    call e%hermitian_multiply('U', 'Full', n, &
          matrix_B_dist, matrix_A2_dist, na_rows, na_cols, &
          matrix_A_dist, na_rows, na_cols, success)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:mult_at_b_real', time_end - time_start)
    time_start = time_end

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:pdtrmm_right', time_end - time_start)
    time_start = time_end

    pdsyevd_trilwmin = 3 * n + max(block_size * (na_rows + 1), 3 * block_size)
    pdsyevd_lwork = max(1 + 6 * n + 2 * na_rows * na_cols, pdsyevd_trilwmin) + 2 * n
    pdsyevd_liwork = 7 * n + 8 * np_cols + 2
    allocate(pdsyevd_work(pdsyevd_lwork), pdsyevd_iwork(pdsyevd_liwork), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa_scalapack: allocation failed', ierr)
    end if
    call pdsyevd('V', 'Upper', n, matrix_A_dist, 1, 1, desc_A, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         pdsyevd_work, pdsyevd_lwork, pdsyevd_iwork, pdsyevd_liwork, info)
    if (info /= 0) then
      call terminate('solver_main, general_elpa_scalapack: pdsyevd failed', 1)
    endif

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:pdsyevd', time_end - time_start)
    time_start = time_end

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    eigenpairs%type_number = 2

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa_scalapack:pdtrmm_EV', time_end - time_start)

    call elpa_deallocate(e)
    call elpa_uninit()
  end subroutine solve_with_general_elpa_scalapack


  subroutine solve_with_general_elpa1(n, proc, matrix_A, eigenpairs, matrix_B)
    class(elpa_t), pointer :: e
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr

    double precision :: time_start, time_end
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :)
    integer :: numroc, success

    time_start = mpi_wtime()

    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)

#if	ELPA_VERSION >= 201705003
    ierr = elpa_get_communicators(mpi_comm_world, my_prow, my_pcol, &
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
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa1: allocation failed', ierr)
    end if
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:setup_matrices', time_end - time_start)
    time_start = time_end

    if (elpa_init(20170403) /= elpa_ok) then
      print *, "ELPA API version not supported"
      stop
    endif
    e => elpa_allocate()

    call e%set("na", n, success) !Size of Matrix
    call e%set("nev", n, success) !Number of eigenvectors to be calculated
    call e%set("local_nrows", na_rows, success)
    call e%set("local_ncols", na_cols, success)
    call e%set("nblk", block_size, success)
    call e%set("mpi_comm_parent", mpi_comm_world, success)
    call e%set("process_row", my_prow, success)
    call e%set("process_col", my_pcol, success)
    success = e%setup()

    call e%set("solver", elpa_solver_1stage, success)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:setup_parameter', time_end - time_start)
    time_start = time_end

    ! Return of cholesky is stored in the upper triangle.
#if	ELPA_VERSION >= 201705003
    call e%cholesky(matrix_B_dist, success)
#else
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa1: cholesky_real failed', ierr)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:cholesky_real', time_end - time_start)
    time_start = time_end

#if	ELPA_VERSION >= 201705003
    call e%invert_triangular(matrix_B_dist, success)
#else
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa1: invert_trm_real failed', 1)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:invert_trm_real', time_end - time_start)
    time_start = time_end

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    call e%hermitian_multiply('U', 'Full', n, &
          matrix_B_dist, matrix_A2_dist, na_rows, na_cols, &
          matrix_A_dist, na_rows, na_cols, success)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:mult_at_b_real', time_end - time_start)
    time_start = time_end

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:pdtrmm_right', time_end - time_start)
    time_start = time_end

    call e%eigenvectors(matrix_A_dist, eigenpairs%blacs%values, eigenpairs%blacs%Vectors, success)

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa1: solve_evp_real failed', 1)
    endif

    call add_event('solve_evp_real:fwd', time_evp_fwd)
    call add_event('solve_evp_real:solve', time_evp_solve)
    call add_event('solve_evp_real:back', time_evp_back)
    call add_event('solve_evp_real', time_evp_fwd + time_evp_solve + time_evp_back)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:solve_evp_real', time_end - time_start)
    time_start = time_end

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    eigenpairs%type_number = 2

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa1:pdtrmm_EV', time_end - time_start)

    call elpa_deallocate(e)
    call elpa_uninit()
  end subroutine solve_with_general_elpa1


  subroutine solve_with_general_elpa2(n, proc, matrix_A, eigenpairs, matrix_B)
    class(elpa_t), pointer :: e
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs
    integer :: desc_A(desc_size), desc_A2(desc_size), desc_B(desc_size), &
         block_size, max_block_size, &
         myid, np_rows, np_cols, my_prow, my_pcol, &
         na_rows, na_cols, mpi_comm_rows, mpi_comm_cols, &
         sc_desc(desc_size), ierr, info, mpierr
    double precision :: time_start, time_end
    double precision, allocatable :: matrix_A_dist(:, :), matrix_A2_dist(:, :), matrix_B_dist(:, :)
    integer :: numroc, success

    time_start = mpi_wtime()

    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call blacs_gridinfo(proc%context, np_rows, np_cols, my_prow, my_pcol)

#if	ELPA_VERSION >= 201705003
ierr = elpa_get_communicators(mpi_comm_world, my_prow, my_pcol, &
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
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_elpa2: allocation failed', ierr)
    end if
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A2, matrix_A2_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:setup_matrices', time_end - time_start)
    time_start = time_end

    if (elpa_init(20170403) /= elpa_ok) then
      print *, "ELPA API version not supported"
      stop
    endif
    e => elpa_allocate()

    call e%set("na", n, success) !Size of Matrix
    call e%set("nev", n, success) !Number of eigenvectors to be calculated
    call e%set("local_nrows", na_rows, success)
    call e%set("local_ncols", na_cols, success)
    call e%set("nblk", block_size, success)
    call e%set("mpi_comm_parent", mpi_comm_world, success)
    call e%set("process_row", my_prow, success)
    call e%set("process_col", my_pcol, success)
    success = e%setup()

    call e%set("solver", elpa_solver_2stage, success)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:setup_parameter', time_end - time_start)
    time_start = time_end

    ! Return of cholesky is stored in the upper triangle.
#if	ELPA_VERSION >= 201705003
    call e%cholesky(matrix_B_dist, success)
#else
    call cholesky_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa2: cholesky_real failed', ierr)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:cholesky_real', time_end - time_start)
    time_start = time_end

#if	ELPA_VERSION >= 201705003
    call e%invert_triangular(matrix_B_dist, success)
#else
    call invert_trm_real(n, matrix_B_dist, na_rows, block_size, mpi_comm_rows, mpi_comm_cols, success)
#endif

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa2: invert_trm_real failed', 1)
    end if

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:invert_trm_real', time_end - time_start)
    time_start = time_end

    ! Reduce A as U^-T A U^-1
    ! A <- U^-T A
    call e%hermitian_multiply('U', 'Full', n, &
          matrix_B_dist, matrix_A2_dist, na_rows, na_cols, &
          matrix_A_dist, na_rows, na_cols, success)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:mult_at_b_real', time_end - time_start)
    time_start = time_end

    ! A <- A U^-1
    call pdtrmm('Right', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, matrix_A_dist, 1, 1, sc_desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:pdtrmm_right', time_end - time_start)
    time_start = time_end

    call e%eigenvectors(matrix_A_dist, eigenpairs%blacs%values, eigenpairs%blacs%Vectors, success)

    if (success == ELPA_ERROR) then
      call terminate('solver_main, general_elpa2: solve_evp_real failed', 1)
    endif

    call add_event('solve_evp_real_2stage:fwd', time_evp_fwd)
    call add_event('solve_evp_real_2stage:solve', time_evp_solve)
    call add_event('solve_evp_real_2stage:back', time_evp_back)
    call add_event('solve_evp_real_2stage', time_evp_fwd + time_evp_solve + time_evp_back)

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:solve_evp_real_2stage', time_end - time_start)
    time_start = time_end

    ! Z <- U^-1 Z
    call pdtrmm('Left', 'Upper', 'No_trans', 'No_unit', n, n, 1.0d0, &
         matrix_B_dist, 1, 1, sc_desc, eigenpairs%blacs%Vectors, 1, 1, sc_desc)

    eigenpairs%type_number = 2

    time_end = mpi_wtime()
    call add_event('solve_with_general_elpa2:pdtrmm_EV', time_end - time_start)

    call elpa_deallocate(e)
    call elpa_uninit()
  end subroutine solve_with_general_elpa2
end module ek_solver_elpa_m
