module ek_solver_scalapack_all_m
  use mpi
  use ek_descriptor_parameters_m
  use ek_distribute_matrix_m, only : &
       get_local_cols, gather_matrix, allgather_row_wise, setup_distributed_matrix, &
       distribute_global_sparse_matrix
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t
  use ek_event_logger_m, only : add_event
  use ek_generalized_to_standard_m, only : reduce_generalized, recovery_generalized
  use ek_processes_m, only : check_master, ek_process_t, terminate
  use ek_matrix_io_m, only : ek_sparse_mat_t
  implicit none

  private
  public :: eigen_solver_scalapack_all, solve_with_general_scalapack

contains

  subroutine eigen_solver_scalapack_all(proc, desc_A, A, eigenpairs)
    type(ek_process_t) :: proc
    integer, intent(in) :: desc_A(9)
    double precision, intent(in) :: A(:, :)
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    integer :: ierr, info
    integer :: dim, work_size, iwork_size, diag_size, subdiag_size
    integer :: eigenvectors_local_cols

    character(len = 1) :: uplo, side

    double precision, allocatable :: diag_local(:), subdiag_local(:)
    double precision, allocatable :: subdiag_global(:)
    double precision, allocatable :: tau(:), work(:), work_print(:)
    integer, allocatable :: iwork(:)
    double precision :: time_start, time_start_part, time_end
    integer :: numroc

    time_start = mpi_wtime()
    time_start_part = time_start

    eigenpairs%type_number = 2

    dim = desc_A(rows_)
    diag_size = numroc(dim, desc_A(block_col_), proc%my_proc_col, 0, proc%n_procs_col)
    subdiag_size = numroc(dim - 1, desc_A(block_col_), proc%my_proc_col, 0, proc%n_procs_col)
    work_size = max(desc_A(block_row_) * (desc_A(local_rows_) + 1), 3 * desc_A(block_row_))

    allocate(diag_local(diag_size), subdiag_local(subdiag_size), tau(diag_size), &
         work(work_size), work_print(desc_A(block_row_)), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver_scalapack_all: allocation failed', ierr)
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_scalapack_all:allocate', time_end - time_start_part)
    time_start_part = time_end

    uplo = 'l'
    call pdsytrd(uplo, dim, A, 1, 1, desc_A, diag_local, subdiag_local, tau, work, work_size, info)
    deallocate(work)
    if (proc%my_rank == 0) then
       print '("info(pdsytrd): ", i0)', info
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_scalapack_all:pdsytrd', time_end - time_start_part)
    time_start_part = time_end

    ! Diagonal elements of the tridiagonal matrix Initially
    allocate(eigenpairs%blacs%values(dim), subdiag_global(dim - 1), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver_scalapack_all: allocation failed', ierr)
    end if

    call allgather_row_wise(diag_local, proc%context, desc_A(block_col_), &
         eigenpairs%blacs%values)
    call allgather_row_wise(subdiag_local, proc%context, desc_A(block_col_), &
         subdiag_global)

    call setup_distributed_matrix('Eigenvectors', proc, dim, dim, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))
    eigenvectors_local_cols = get_local_cols(proc, eigenpairs%blacs%desc)

    work_size = 6 * dim + 2 * eigenpairs%blacs%desc(local_rows_) * &
         eigenvectors_local_cols
    iwork_size = 2 + 7 * dim + 8 * proc%n_procs_col
    allocate(work(work_size), iwork(iwork_size), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver_scalapack_all: allocation failed', ierr)
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_scalapack_all:gather1', time_end - time_start_part)
    time_start_part = time_end

    call pdstedc('i', dim, eigenpairs%blacs%values, subdiag_global, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         work, work_size, iwork, iwork_size, info)
    if (proc%my_rank == 0) then
       print '("info(pdstedc): ", i0)', info
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_scalapack_all:pdstedc', time_end - time_start_part)
    time_start_part = time_end

    side = 'l'
    work_size = work_size_for_pdormtr(side, uplo, dim, dim, 1, 1, 1, 1, desc_A, eigenpairs%blacs%desc)
    deallocate(work)
    allocate(work(work_size), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver_scalapack_all: allocation failed', ierr)
    end if

    call pdormtr(side, uplo, 'n', dim, dim, A, 1, 1, desc_A, tau, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, work, work_size, info)
    if (proc%my_rank == 0) then
       print '("info(pdormtr): ", i0)', info
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_scalapack_all:pdormtr', time_end - time_start_part)
    call add_event('eigen_solver_scalapack_all', time_end - time_start)
  end subroutine eigen_solver_scalapack_all


  subroutine solve_with_general_scalapack(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in) :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_B(desc_size)
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :)
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack:setup_matrices', time_end - time_start_part)
    time_start_part = time_end

    call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack:reduce_generalized', time_end - time_start_part)
    time_start_part = time_end

    call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack:eigen_solver_scalapack_all', time_end - time_start_part)
    time_start_part = time_end

    call recovery_generalized(n, n, matrix_B_dist, desc_B, &
         eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack:recovery_generalized', time_end - time_start_part)
    call add_event('solve_with_general_scalapack', time_end - time_start)
  end subroutine solve_with_general_scalapack


  integer function work_size_for_pdormtr(side, uplo, m, n, ia, ja, ic, jc, desc_A, desc_C) result(size)
    character(len = 1), intent(in) :: side, uplo
    integer, intent(in) :: m, n, ia, ja, ic, jc
    integer, intent(in) :: desc_A(9), desc_C(9)
    integer :: context, n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: mb_a, nb_a, mb_c, nb_c
    integer :: iaa, jaa, icc, jcc
    integer :: lcmq, npa0 , mpc0, nqc0
    integer :: iroffa, icoffa, iarow, iroffc, icoffc, icrow, iccol
    integer :: mi, ni
    logical :: is_upper, is_left

    integer :: numroc, indxg2p, ilcm
    logical :: lsame

    is_upper = lsame(uplo, 'U')
    is_left = lsame(side, 'L')

    mb_a = desc_A(5)
    nb_a = desc_A(6)
    mb_c = desc_C(5)
    nb_c = desc_C(6)

    context = desc_A(2)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (is_upper) then
      iaa = ia
      jaa = ja + 1
      icc = ic
      jcc = jc
    else
      iaa = ia + 1
      jaa = ja
      if (is_left) then
        icc = ic + 1
        jcc = jc
      else
        icc = ic
        jcc = jc + 1
      end if
    end if

    if (is_left) then
      mi = m - 1
      ni = n
    else
      mi = m
      ni = n - 1
    end if

    iroffc = mod(icc - 1, mb_c)
    icoffc = mod(jcc - 1, nb_c)
    icrow = indxg2p(icc, mb_c, my_proc_row, desc_C(7), n_procs_row)
    iccol = indxg2p(jcc, nb_c, my_proc_col, desc_C(8), n_procs_col)

    mpc0 = numroc(mi + iroffc, mb_c, my_proc_row, icrow, n_procs_row)
    nqc0 = numroc(ni + icoffc, nb_c, my_proc_col, iccol, n_procs_col)

    if (is_left) then
      size = max((nb_a * (nb_a - 1)) / 2, (nqc0 + mpc0) * nb_a) + nb_a * nb_a
    else
      iroffa = mod(iaa - 1, mb_a)
      icoffa = mod(jaa - 1, nb_a)
      iarow = indxg2p(iaa, mb_a, my_proc_row, desc_A(7), n_procs_row)
      npa0 = numroc(ni + iroffa, mb_a, my_proc_row, iarow, n_procs_row)
      lcmq = ilcm(n_procs_row, n_procs_col) / n_procs_col
      size = max((nb_a * (nb_a - 1)) / 2, &
           (nqc0 + max(npa0 + numroc(numroc(ni + icoffc, nb_a, 0, 0, n_procs_col), &
           nb_a, 0, 0, lcmq), mpc0)) * nb_a) + &
           nb_a * nb_a
    end if
  end function work_size_for_pdormtr
end module ek_solver_scalapack_all_m
