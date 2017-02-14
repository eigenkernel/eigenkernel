module ek_verifier_m
  use mpi
  use ek_command_argument_m, only : ek_argument_t
  use ek_descriptor_parameters_m
  use ek_distribute_matrix_m, only : convert_sparse_matrix_to_dense, &
       setup_distributed_matrix, distribute_global_sparse_matrix
  use ek_eigenpairs_types_m, only: ek_eigenpairs_types_union_t, ek_eigenpairs_local_t, &
       ek_eigenpairs_blacs_t
  use ek_event_logger_m, only : add_event
  use ek_global_variables_m, only : g_block_size
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_processes_m, only : ek_process_t, terminate
  implicit none

  private
  public :: eval_residual_norm, eval_orthogonality

contains

  subroutine eval_residual_norm_local(arg, matrix_A, eigenpairs, &
       A_norm, res_norm_ave, res_norm_max, matrix_B)
    type(ek_argument_t), intent(in) :: arg
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(in) :: eigenpairs

    double precision, intent(out) :: A_norm, res_norm_ave, res_norm_max

    double precision, allocatable :: a(:,:)
    double precision, allocatable :: b(:,:)
    double precision, allocatable :: left(:), right(:) ! residual = left - right
    double precision, allocatable :: res_norm(:) ! residual norm

    integer :: j, dim, ierr
    double precision :: time_start, time_end
    double precision :: dlange  ! Function.

    time_start = mpi_wtime()

    dim = matrix_A%size

    allocate(left(dim), right(dim), res_norm(arg%n_check_vec), stat = ierr)
    if (ierr /= 0) then
      call terminate('eval_residual_norm_local: allocation failed', ierr)
    end if
    res_norm(:) = 0.0d0

    call convert_sparse_matrix_to_dense(matrix_A, a)
    if (arg%is_generalized_problem) then
      call convert_sparse_matrix_to_dense(matrix_B, b)
    end if

    do j = 1, arg%n_check_vec
      left(:) = matmul(a, eigenpairs%local%vectors(:, j))
      if (arg%is_generalized_problem) then
        right(:) = eigenpairs%local%values(j) * &
             matmul(b, eigenpairs%local%vectors(:, j))
      else
        right(:) = eigenpairs%local%values(j) * eigenpairs%local%vectors(:, j)
      end if
      res_norm(j) = sqrt(abs(dot_product(left - right, left - right) / &
           dot_product(eigenpairs%local%vectors(:, j), &
           eigenpairs%local%vectors(:, j))))
    enddo

    A_norm = dlange('F', dim, dim, a, dim, 0)
    res_norm_max = maxval(res_norm) / A_norm
    res_norm_ave = sum(res_norm) / A_norm / dble(arg%n_check_vec)

    time_end = mpi_wtime()
    call add_event('eval_residual_norm_local', time_end - time_start)
  end subroutine eval_residual_norm_local


  subroutine eval_residual_norm_blacs(arg, matrix_A, eigenpairs, &
       A_norm, res_norm_ave, res_norm_max, matrix_B)
    type(ek_argument_t), intent(in) :: arg
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_eigenpairs_blacs_t), intent(inout) :: eigenpairs
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    ! residual norm average, max
    double precision, intent(out) :: A_norm, res_norm_ave, res_norm_max

    type(ek_process_t) :: proc
    integer :: dim, desc_Residual(desc_size), desc_A(desc_size), desc_B(desc_size)
    double precision, allocatable :: Residual(:, :), matrix_A_dist(:, :), matrix_B_dist(:, :)
    ! ave_and_max is declared as array due to usage of bcast
    ! 3rd element is for the index of the max value (discarded currently)
    double precision :: res_norm, ave_and_max(3)
    integer :: j, block_size, owner_proc_col, ierr
    double precision :: time_start, time_start_part, time_end
    ! ScaLAPACK functions
    integer :: indxg2p
    double precision :: pdlange

    time_start = mpi_wtime()
    time_start_part = time_start

    if (arg%is_generalized_problem .and. .not. present(matrix_B)) then
      call terminate('eval_residual_norm_blacs: matrix_B is not provided', 1)
    end if

    if (trim(arg%solver_type) == 'eigenexa' .or. &
         trim(arg%solver_type) == 'general_eigenexa') then
      block_size = 1
    else if (arg%block_size > 0) then  ! Do not use the default block size
      block_size = arg%block_size
    else
      block_size = g_block_size
    end if

    ! call blacs_get(-1, 0, proc%context)
    ! Because context acquiring by blacs_get can fail in some environments,
    ! use context in the descriptor of eigenpairs instead.
    proc%context = eigenpairs%desc(context_)
    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    dim = arg%matrix_A_info%rows

    call setup_distributed_matrix('A', proc, dim, dim, desc_A, matrix_A_dist, &
         block_size = block_size)
    call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
    A_norm = pdlange('F', dim, dim, matrix_A_dist, 1, 1, desc_A, 0)

    if (arg%is_generalized_problem) then
      call setup_distributed_matrix('B', proc, dim, dim, &
           desc_B, matrix_B_dist, block_size = block_size)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
    end if

    call setup_distributed_matrix('Residual', proc, dim, arg%n_check_vec, &
         desc_Residual, Residual, block_size = block_size)

    time_end = mpi_wtime()
    call add_event('eval_residual_norm_blacs:setup_matrices', time_end - time_start_part)
    time_start_part = time_end

    if (arg%is_generalized_problem) then
      ! Residual <- B * Eigenvectors
      call pdsymm('L', 'L', dim, arg%n_check_vec, &
           1.0d0, matrix_B_dist, 1, 1, desc_B, &
           eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
           0.0d0, Residual, 1, 1, desc_Residual)
      time_end = mpi_wtime()
      call add_event('eval_residual_norm_blacs:pdsymm_B', time_end - time_start_part)
      time_start_part = time_end
    else
      ! Residual <- Eigenvectors
      call pdlacpy('A', dim, arg%n_check_vec, &
           eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
           Residual, 1, 1, desc_Residual)
      time_end = mpi_wtime()
      call add_event('eval_residual_norm_blacs:pdlacpy', time_end - time_start_part)
      time_start_part = time_end
    end if

    ! For each column j, Residual(*, j) <- Residual(*, j) * -eigenvalues(j)
    do j = 1, arg%n_check_vec
      call pdscal(dim, -1.0d0 * eigenpairs%values(j), &
           Residual, 1, j, desc_Residual, 1)
    end do
    time_end = mpi_wtime()
    call add_event('eval_residual_norm_blacs:pdscal', time_end - time_start_part)
    time_start_part = time_end

    ! Residual <- Residual + A * Eigenvectors
    call pdsymm('L', 'L', dim, arg%n_check_vec, &
         1.0d0, matrix_A_dist, 1, 1, desc_A, &
         eigenpairs%Vectors, 1, 1, eigenpairs%desc, &
         1.0d0, Residual, 1, 1, desc_Residual)
    time_end = mpi_wtime()
    call add_event('eval_residual_norm_blacs:pdsymm_R', time_end - time_start_part)
    time_start_part = time_end

    ! Store 2-norm of each residual column in the first row of the column
    do j = 1, arg%n_check_vec
      call pdnrm2(dim, res_norm, Residual, 1, j, desc_Residual, 1)
      owner_proc_col = indxg2p(j, desc_Residual(block_col_), 0, &
           desc_Residual(csrc_), proc%n_procs_col)
      if (proc%my_proc_row == 0 .and. proc%my_proc_col == owner_proc_col) then
        call pdelset(Residual, 1, j, desc_Residual, res_norm)
      end if
    end do
    time_end = mpi_wtime()
    call add_event('eval_residual_norm_blacs:pdnrm2', time_end - time_start_part)
    time_start_part = time_end

    ! Although these pd* routines return result values only in the scope of focused subvector,
    ! processor rank 0 must be in the scope because the subvector is the first row of Residual here
    call pdasum(arg%n_check_vec, ave_and_max(1), &
         Residual, 1, 1, desc_Residual, desc_Residual(rows_))
    call pdamax(arg%n_check_vec, ave_and_max(2), ave_and_max(3), &
         Residual, 1, 1, desc_Residual, desc_Residual(rows_))
    call mpi_bcast(ave_and_max(1), 3, mpi_double_precision, 0, mpi_comm_world, ierr)

    res_norm_ave = ave_and_max(1) / A_norm / dble(arg%n_check_vec)
    res_norm_max = ave_and_max(2) / A_norm

    time_end = mpi_wtime()
    call add_event('eval_residual_norm_blacs:finish', time_end - time_start_part)
    call add_event('eval_residual_norm_blacs', time_end - time_start)
  end subroutine eval_residual_norm_blacs


  subroutine eval_residual_norm(arg, matrix_A, eigenpairs, &
       A_norm, res_norm_ave, res_norm_max, matrix_B)
    type(ek_argument_t), intent(in) :: arg
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(inout) :: eigenpairs
    ! residual norm average, max
    double precision, intent(out) :: A_norm, res_norm_ave, res_norm_max

    if (eigenpairs%type_number == 1) then
      call eval_residual_norm_local(arg, matrix_A, eigenpairs, &
           A_norm, res_norm_ave, res_norm_max, matrix_B)
    else if (eigenpairs%type_number == 2) then
      if (arg%is_generalized_problem) then
        call eval_residual_norm_blacs(arg, matrix_A, eigenpairs%blacs, &
             A_norm, res_norm_ave, res_norm_max, matrix_B)
      else
        call eval_residual_norm_blacs(arg, matrix_A, eigenpairs%blacs, &
             A_norm, res_norm_ave, res_norm_max)
      end if
    else
      print '("[Warning] eval_residual_norm: output of this type is not supported yet")'
    end if
  end subroutine eval_residual_norm


  subroutine eval_orthogonality_blacs(index1, index2, eigenpairs, orthogonality, matrix_B)
    integer, intent(in) :: index1, index2
    type(ek_eigenpairs_blacs_t), intent(in) :: eigenpairs
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    double precision, intent(out) :: orthogonality

    type(ek_process_t) :: proc
    integer :: dim, n, block_size
    integer :: diag, i, j, owner_proc_row, owner_proc_col, info
    double precision :: scale
    double precision, allocatable :: InnerProducts(:, :), matrix_B_dist(:, :), BV(:, :)
    integer :: desc_InnerProducts(desc_size), desc_B(desc_size), desc_BV(desc_size)
    double precision :: time_start, time_start_part, time_end

    ! Functions.
    integer :: indxg2p, indxg2l, blacs_pnum
    double precision :: pdlange

    time_start = mpi_wtime()
    time_start_part = time_start

    dim = eigenpairs%desc(rows_)  ! Dimension of the original problem.
    n = index2 - index1 + 1  ! The number of vectors to be checked.
    block_size = eigenpairs%desc(block_row_)
    if (block_size /= eigenpairs%desc(block_col_)) then
      call terminate('eval_orthogonality_blacs: anisotropic block size not supported',1)
    end if
    proc%context = eigenpairs%desc(context_)
    call blacs_pinfo(proc%my_rank, proc%n_procs)
    call blacs_gridinfo(proc%context, proc%n_procs_row, proc%n_procs_col, &
         proc%my_proc_row, proc%my_proc_col)

    call setup_distributed_matrix('InnerProducts', proc, n, n, &
         desc_InnerProducts, InnerProducts, block_size = block_size)
    if (present(matrix_B)) then
      call setup_distributed_matrix('B', proc, dim, dim, &
           desc_B, matrix_B_dist, block_size = block_size)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call setup_distributed_matrix('BV', proc, dim, n, &
           desc_BV, BV, block_size = block_size)

      time_end = mpi_wtime()
      call add_event('eval_orthogonality_blacs:setup_matrices', time_end - time_start_part)
      time_start_part = time_end

      ! BV <- B * V
      call pdgemm('N', 'N', dim, n, dim, 1.0d0, &
           matrix_B_dist, 1, 1, desc_B, &
           eigenpairs%Vectors, 1, index1, eigenpairs%desc, &
           0.0d0, BV, 1, 1, desc_BV)

      time_end = mpi_wtime()
      call add_event('eval_orthogonality_blacs:pdgemm_R', time_end - time_start_part)
      time_start_part = time_end

      ! InnerProducts <- V' * BV
      call pdgemm('T', 'N', n, n, dim, 1.0d0, &
           eigenpairs%Vectors, 1, index1, eigenpairs%desc, &
           BV, 1, 1, desc_BV, &
           0.0d0, InnerProducts, 1, 1, desc_InnerProducts)

      time_end = mpi_wtime()
      call add_event('eval_orthogonality_blacs:pdgemm_L', time_end - time_start_part)
      time_start_part = time_end
    else
      ! InnerProducts <- V' * V
      call pdgemm('T', 'N', n, n, dim, 1.0d0, &
           eigenpairs%Vectors, 1, index1, eigenpairs%desc, &
           eigenpairs%Vectors, 1, index1, eigenpairs%desc, &
           0.0d0, InnerProducts, 1, 1, desc_InnerProducts)

      time_end = mpi_wtime()
      call add_event('eval_orthogonality_blacs:pdgemm', time_end - time_start_part)
      time_start_part = time_end
    end if

    ! Scale elements and subtract identity matrix from it.
    do diag = 1, n
      owner_proc_row = indxg2p(diag, desc_InnerProducts(block_row_), 0, 0, proc%n_procs_row)
      owner_proc_col = indxg2p(diag, desc_InnerProducts(block_col_), 0, 0, proc%n_procs_col)
      if (proc%my_proc_col == owner_proc_col .and. proc%my_proc_row == owner_proc_row) then
        i = indxg2l(diag, desc_InnerProducts(block_row_), 0, 0, proc%n_procs_row)
        j = indxg2l(diag, desc_InnerProducts(block_col_), 0, 0, proc%n_procs_col)
        scale = 1.0d0 / sqrt(InnerProducts(i, j))  ! Inverse of norm.
        InnerProducts(i, j) = 0.0d0  ! InnerProducts <- InnerProducts - I
      end if
      call mpi_bcast(scale, 1, mpi_double_precision, &
           blacs_pnum(proc%context, owner_proc_row, owner_proc_col), mpi_comm_world, info)
      call pdscal(n, scale, InnerProducts, 1, diag, desc_InnerProducts, 1)
      call pdscal(n, scale, InnerProducts, diag, 1, desc_InnerProducts, n)
    end do

    orthogonality = pdlange('F', n, n, InnerProducts, 1, 1, desc_InnerProducts, 0)

    time_end = mpi_wtime()
    call add_event('eval_orthogonality_blacs:finish', time_end - time_start_part)
    call add_event('eval_orthogonality_blacs', time_end - time_start)
  end subroutine eval_orthogonality_blacs


  subroutine eval_orthogonality(arg, eigenpairs, orthogonality, matrix_B)
    type(ek_argument_t), intent(in) :: arg
    type(ek_eigenpairs_types_union_t), intent(in) :: eigenpairs
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    double precision, intent(out) :: orthogonality

    if (eigenpairs%type_number == 2) then
      if (arg%is_generalized_problem .and. .not. present(matrix_B)) then
        call terminate('eval_orthogonality: matrix_B is not provided', 1)
      end if
      if (present(matrix_B)) then
        call eval_orthogonality_blacs(arg%ortho_check_index_start, &
             arg%ortho_check_index_end, eigenpairs%blacs, orthogonality, matrix_B)
      else
        call eval_orthogonality_blacs(arg%ortho_check_index_start, &
             arg%ortho_check_index_end, eigenpairs%blacs, orthogonality)
      end if
    else
      print '("[Warning] eval_orthogonality: orthogonality evaluator for output of this type is not implemeted yet")'
    end if
  end subroutine eval_orthogonality
end module ek_verifier_m
