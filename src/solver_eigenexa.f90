module ek_solver_eigenexa_m
  use eigen_libs
  use mpi
  use ek_distribute_matrix_m, only : setup_distributed_matrix, distribute_global_sparse_matrix
  use ek_descriptor_parameters_m
  use ek_eigenpairs_types_m, only : eigenpairs_types_union
  use ek_event_logger_m, only : add_event
  use ek_generalized_to_standard_m, only : reduce_generalized, &
       reduce_generalized_new, recovery_generalized
  use ek_matrix_io_m, only : sparse_mat
  use ek_processes_m, only : check_master, terminate, process

  implicit none

  private
  public :: setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa, eigen_solver_eigenk, &
       solve_with_general_scalapack_eigenexa, solve_with_general_scalapack_eigenk, &
       solve_with_general_scalapacknew_eigenk

contains

  subroutine setup_distributed_matrix_for_eigenexa(dim, desc_A, matrix_A, eigenpairs)
    integer, intent(in) :: dim
    integer, intent(out) ::desc_A(desc_size)
    double precision, allocatable, intent(out) :: matrix_A(:, :)
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: nx, ny, context, info
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    if (check_master()) then
      print '( "Creating 2 distributed matrices for EigenExa with &
           &M, N, MB, NB: ", I0, ", ", I0, ", ", I0, ", ", I0 )', &
           dim, dim, 1, 1
    end if

    call eigen_init()
    call eigen_get_matdims(dim, nx, ny)
    context = eigen_get_blacs_context()

    call descinit(desc_A, dim, dim, 1, 1, 0, 0, context, nx, info)

    eigenpairs%type_number = 2
    call descinit(eigenpairs%blacs%desc, dim, dim, 1, 1, 0, 0, context, nx, info)
    if (info /= 0) then
      print '(a, i0)', 'info(descinit): ', info
      call terminate('setup_distributed_matrix_for_eigenexa: descinit failed', info)
    end if

    allocate(matrix_A(nx, ny), eigenpairs%blacs%Vectors(nx, ny), &
         eigenpairs%blacs%values(dim), stat = info)
    if (info /= 0) then
      call terminate('setup_distributed_matrix_for_eigenexa: allocation failed', info)
    end if

    matrix_A(:, :) = 0.0d0
    eigenpairs%blacs%Vectors(:, :) = 0.0d0

    time_end = mpi_wtime()
    call add_event('setup_distributed_matrix_for_eigenexa', time_end - time_start)
  end subroutine setup_distributed_matrix_for_eigenexa


  ! uplo: takes the value of 'L' or 'U' or (not present).
  !       Specifies how the entries of the input matrix is stored.
  !       If not present, it is assumed that both of upper and lower
  !       triangular parts are stored.
  subroutine eigen_solver_eigenexa(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    integer :: i, dim, nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128
    double precision :: time_start, time_start_part, time_end
#if	USE_EIGENEXA_WITH_TIMER
    double precision :: eigen_times(8)
#endif

    time_start = mpi_wtime()
    time_start_part = time_start

    dim = desc_mat(rows_)
    if (dim /= n_vec) then
      call terminate('eigen_solver_eigenexa: current version of EigenExa does not support partial eigenvector computation', 1)
    end if

    ! Unlike usual ScaLAPACK routines, EigenExa requires both of upper and lower
    ! triangular parts of the input matrix. Therefore copying from the lower
    ! triangular parts to upper one (or inverse) is needed.
    if (present(uplo)) then
      if (uplo == 'U') then ! Note: Not tested
        do i = 1, dim - 1
          call pdcopy(dim - i, mat, i, i + 1, desc_mat, dim, &
               mat, i + 1, i, desc_mat, 1)
        end do
      else if (uplo == 'L') then
        do i = 1, dim - 1
          call pdcopy(dim - i, mat, i + 1, i, desc_mat, 1, &
               mat, i, i + 1, desc_mat, dim)
        end do
      else
        call terminate("eigen_solver_eigenexa: uplo must be 'U' or 'L'", 1)
      end if
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenexa:transpose', time_end - time_start_part)
    time_start_part = time_end

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

#if	USE_EIGENEXA_WITH_TIMER
    call eigen_sx(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, eigen_times, &
         m_forward = m_forward, m_backward = m_backward)
    call add_event('eigen_sx:fwd', eigen_times(1))
    call add_event('!eigen_sx:fwd_Gflops', eigen_times(2))
    call add_event('eigen_sx:dc', eigen_times(3))
    call add_event('!eigen_sx:dc_Gflops', eigen_times(4))
    call add_event('eigen_sx:bak', eigen_times(5))
    call add_event('!eigen_sx:bak_Gflops', eigen_times(6))
    call add_event('eigen_sx', eigen_times(7))
    call add_event('!eigen_sx:total_Gflops', eigen_times(8))
#else
    call eigen_sx(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)
#endif

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenexa:eigen_sx', time_end - time_start_part)
    call add_event('eigen_solver_eigenexa', time_end - time_start)
  end subroutine eigen_solver_eigenexa


  subroutine eigen_solver_eigenk(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    integer :: i, dim, nx, ny, my_rank, ierr
    integer, parameter :: m_forward = 48, m_backward = 128
    double precision :: time_start, time_start_part, time_end
#if	USE_EIGENEXA_WITH_TIMER
    double precision :: eigen_times(8)
#endif

    time_start = mpi_wtime()
    time_start_part = time_start

    dim = desc_mat(rows_)
    if (dim /= n_vec) then
      call terminate('eigen_solver_eigenexa: current version of EigenExa does not support partial eigenvector computation', 1)
    end if

    ! Unlike usual ScaLAPACK routines, EigenExa requires both of upper and lower
    ! triangular parts of the input matrix. Therefore copying from the lower
    ! triangular parts to upper one (or inverse) is needed.
    if (present(uplo)) then
      if (uplo == 'U') then ! Note: Not tested
        do i = 1, dim - 1
          call pdcopy(dim - i, mat, i, i + 1, desc_mat, dim, &
               mat, i + 1, i, desc_mat, 1)
        end do
      else if (uplo == 'L') then
        do i = 1, dim - 1
          call pdcopy(dim - i, mat, i + 1, i, desc_mat, 1, &
               mat, i, i + 1, desc_mat, dim)
        end do
      else
        call terminate("eigen_solver_eigenk: uplo must be 'U' or 'L'", 1)
      end if
    end if

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenk:transpose', time_end - time_start_part)
    time_start_part = time_end

    call eigen_get_matdims(dim, nx, ny)

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)

#if	USE_EIGENEXA_WITH_TIMER
    call eigen_s(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, eigen_times, &
         m_forward = m_forward, m_backward = m_backward)
    call add_event('eigen_s:fwd', eigen_times(1))
    call add_event('!eigen_s:fwd_Gflops', eigen_times(2))
    call add_event('eigen_s:dc', eigen_times(3))
    call add_event('!eigen_s:dc_Gflops', eigen_times(4))
    call add_event('eigen_s:bak', eigen_times(5))
    call add_event('!eigen_s:bak_Gflops', eigen_times(6))
    call add_event('eigen_s', eigen_times(7))
    call add_event('!eigen_s:total_Gflops', eigen_times(8))
#else
    call eigen_s(dim, dim, mat, nx, &
         eigenpairs%blacs%values, eigenpairs%blacs%Vectors, nx, &
         m_forward = m_forward, m_backward = m_backward)
#endif

    time_end = mpi_wtime()
    call add_event('eigen_solver_eigenk:eigen_s', time_end - time_start_part)
    call add_event('eigen_solver_eigenk', time_end - time_start)
  end subroutine eigen_solver_eigenk


  subroutine solve_with_general_scalapack_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in) :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_B(desc_size), desc_A_re(desc_size)
    integer :: ierr
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs_tmp
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:setup_matrices', time_end - time_start_part)
    time_start_part = time_end

    call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:reduce_generalized', time_end - time_start_part)
    time_start_part = time_end

    call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
    deallocate(matrix_A_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:pdgemr2d_1', time_end - time_start_part)
    time_start_part = time_end

    call eigen_solver_eigenexa(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:eigen_solver_eigenexa', time_end - time_start_part)
    time_start_part = time_end

    deallocate(matrix_A_redist)
    eigenpairs%type_number = 2
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_scalapack_eigenexa: allocation failed', ierr)
    end if
    eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
    eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:setup_EV', time_end - time_start_part)
    time_start_part = time_end

    call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         eigenpairs_tmp%blacs%desc(context_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:pdgemr2d_2', time_end - time_start_part)
    time_start_part = time_end

    call recovery_generalized(n, n, matrix_B_dist, desc_B, &
         eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenexa:recovery_generalized', time_end - time_start_part)
    call add_event('solve_with_general_scalapack_eigenexa', time_end - time_start)
  end subroutine solve_with_general_scalapack_eigenexa


  subroutine solve_with_general_scalapack_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in) :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_B(desc_size), desc_A_re(desc_size), ierr
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs_tmp
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:setup_matrices', time_end - time_start_part)
    time_start_part = time_end

    call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:reduce_generalized', time_end - time_start_part)
    time_start_part = time_end

    call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
    deallocate(matrix_A_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:pdgemr2d_1', time_end - time_start_part)
    time_start_part = time_end

    call eigen_solver_eigenk(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:eigen_solver_eigenk', time_end - time_start_part)
    time_start_part = time_end

    deallocate(matrix_A_redist)
    eigenpairs%type_number = 2
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('eigen_solver, general_scalapack_eigenk: allocation failed', ierr)
    end if
    eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
    eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:setup_EV', time_end - time_start_part)
    time_start_part = time_end

    call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         eigenpairs_tmp%blacs%desc(context_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:pdgemr2d_2', time_end - time_start_part)
    time_start_part = time_end

    call recovery_generalized(n, n, matrix_B_dist, desc_B, &
         eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapack_eigenk:recovery_generalized', time_end - time_start_part)
    call add_event('solve_with_general_scalapack_eigenk', time_end - time_start)
  end subroutine solve_with_general_scalapack_eigenk


  subroutine solve_with_general_scalapacknew_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in) :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: desc_A(desc_size), desc_B(desc_size), desc_A_re(desc_size), ierr
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :), matrix_A_redist(:, :)
    type(eigenpairs_types_union) :: eigenpairs_tmp
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

    call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
    call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
    call setup_distributed_matrix_for_eigenexa(n, desc_A_re, matrix_A_redist, eigenpairs_tmp)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:setup_matrices', time_end - time_start_part)
    time_start_part = time_end

    call reduce_generalized_new(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:reduce_generalized_new', time_end - time_start_part)
    time_start_part = time_end

    call pdgemr2d(n, n, matrix_A_dist, 1, 1, desc_A, matrix_A_redist, 1, 1, desc_A_re, desc_A(context_))
    deallocate(matrix_A_dist)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:pdgemr2d_1', time_end - time_start_part)
    time_start_part = time_end

    call eigen_solver_eigenk(matrix_A_redist, desc_A_re, n, eigenpairs_tmp, 'L')

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:eigen_solver_eigenk', time_end - time_start_part)
    time_start_part = time_end

    deallocate(matrix_A_redist)
    eigenpairs%type_number = 2
    allocate(eigenpairs%blacs%values(n), stat = ierr)
    if (ierr /= 0) then
      call terminate('solve_with_general_scalapacknew_eigenk: allocation failed', ierr)
    end if
    eigenpairs%blacs%values(:) = eigenpairs_tmp%blacs%values(:)
    eigenpairs%blacs%desc(:) = eigenpairs_tmp%blacs%desc(:)
    call setup_distributed_matrix('Eigenvectors', proc, n, n, &
         eigenpairs%blacs%desc, eigenpairs%blacs%Vectors, desc_A(block_row_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:setup_EV', time_end - time_start_part)
    time_start_part = time_end

    call pdgemr2d(n, n, eigenpairs_tmp%blacs%Vectors, 1, 1, eigenpairs_tmp%blacs%desc, &
         eigenpairs%blacs%Vectors, 1, 1, eigenpairs%blacs%desc, &
         eigenpairs_tmp%blacs%desc(context_))

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:pdgemr2d_2', time_end - time_start_part)
    time_start_part = time_end

    call recovery_generalized(n, n, matrix_B_dist, desc_B, &
         eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)

    time_end = mpi_wtime()
    call add_event('solve_with_general_scalapacknew_eigenk:recovery_generalized', time_end - time_start_part)
    call add_event('solve_with_general_scalapacknew_eigenk', time_end - time_start)
  end subroutine solve_with_general_scalapacknew_eigenk
end module ek_solver_eigenexa_m
