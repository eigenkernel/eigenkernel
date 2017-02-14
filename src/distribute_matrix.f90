module ek_distribute_matrix_m
  use mpi
  use ek_command_argument_m, only : ek_matrix_info_t
  use ek_descriptor_parameters_m
  use ek_event_logger_m, only : add_event
  use ek_global_variables_m, only : g_block_size
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_processes_m, only : check_master, ek_process_t, layout_procs, terminate
  implicit none

  public :: get_local_cols, setup_distributed_matrix, &
       distribute_global_dense_matrix, distribute_global_sparse_matrix, &
       convert_sparse_matrix_to_dense, gather_matrix, gather_matrix_part, allgather_row_wise, &
       bcast_sparse_matrix

contains

  subroutine get_ipratios(proc, V, V_desc, ipratios, S_sparse)
    type(ek_process_t), intent(in) :: proc
    real(8), intent(in) :: V(:, :)
    integer, intent(in) :: V_desc(desc_size)
    real(8), intent(out) :: ipratios(V_desc(cols_))
    type(ek_sparse_mat_t), intent(in), optional :: S_sparse

    integer :: i, j, n_procs_row, n_procs_col, my_proc_row, my_proc_col
    real(8) :: elem, sum_power4(V_desc(cols_)), sum_power2(V_desc(cols_))
    integer :: indxg2p
    ! For overlap_mode.
    integer :: m, n
    real(8), allocatable :: S(:, :)
    integer :: S_desc(desc_size)
    real(8), allocatable :: SV(:, :)
    real(8) :: elem2
    logical :: is_overlap_mode
    character :: scope = 'A', topology = ' '

    is_overlap_mode = present(S_sparse)
    if (is_overlap_mode) then
      call setup_distributed_matrix('S', proc, S_sparse%size, S_sparse%size, S_desc, S)
      call distribute_global_sparse_matrix(S_sparse, S_desc, S)
      if (S_desc(rows_) /= S_desc(cols_) .or. S_desc(cols_) /= V_desc(rows_)) then
        call terminate('inconsistent matrix dimension', 1)
      end if
      m = size(V, 1)
      n = size(V, 2)
      allocate(SV(m, n))
      call pdgemm('No', 'No', S_desc(rows_), V_desc(cols_), S_desc(cols_), 1d0, S, 1, 1, S_desc, &
           V, 1, 1, V_desc, 0d0, SV, 1, 1, V_desc)
    end if

    call blacs_gridinfo(V_desc(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    sum_power4(:) = 0d0
    sum_power2(:) = 0d0
    do j = 1, V_desc(cols_)
      if (indxg2p(j, V_desc(block_col_), 0, 0, n_procs_col) == my_proc_col) then
        do i = 1, V_desc(rows_)
          if (indxg2p(i, V_desc(block_row_), 0, 0, n_procs_row) == my_proc_row) then
            call pdelget('', '', elem, V, i, j, V_desc)
            if (is_overlap_mode) then
              call pdelget('', '', elem2, SV, i, j, V_desc)
            end if
            sum_power4(j) = sum_power4(j) + elem ** 4d0
            if (is_overlap_mode) then
              sum_power2(j) = sum_power2(j) + elem * elem2
            else
              sum_power2(j) = sum_power2(j) + elem ** 2d0
            end if
          end if
        end do
      end if
    end do
    call dgsum2d(V_desc(context_), scope, topology, 1, V_desc(cols_), sum_power4, 1, -1, -1)
    call dgsum2d(V_desc(context_), scope, topology, 1, V_desc(cols_), sum_power2, 1, -1, -1)
    ipratios(:) = 0d0
    do j = 1, V_desc(cols_)  ! Redundant computation.
      ipratios(j) = sum_power4(j) / (sum_power2(j) ** 2d0)
    end do
  end subroutine get_ipratios


  integer function get_local_cols(proc, desc)
    type(ek_process_t), intent(in) :: proc
    integer, intent(in) :: desc(desc_size)

    integer :: numroc

    get_local_cols = max(1, numroc(desc(cols_), desc(block_col_), &
         proc%my_proc_col, 0, proc%n_procs_col))
  end function get_local_cols


  subroutine setup_distributed_matrix(name, proc, rows, cols, desc, mat, block_size)
    character(*), intent(in) :: name
    type(ek_process_t), intent(in) :: proc
    integer, intent(in) :: rows, cols
    integer, intent(out) :: desc(desc_size)
    double precision, intent(out), allocatable :: mat(:, :)
    integer, intent(in), optional :: block_size

    integer :: numroc
    integer :: max_block_size, actual_block_size, local_rows, info
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    if (present(block_size)) then
      actual_block_size = block_size
    else
      actual_block_size = g_block_size
    end if

    ! If there is a process which owns no entries in given block size
    ! configuration, diminish the block size and warn about in.
    max_block_size = max(min(rows / proc%n_procs_row, cols / proc%n_procs_col), 1)
    if (actual_block_size > max_block_size) then
      if (check_master()) then
        print '("[Warning] setup_distributed_matrix: size of matrix is very small relative to the number of processes")'
      end if
      actual_block_size = max_block_size
    end if

    if (check_master()) then
      print '( "Creating distributed matrix ", A, " with M, N, MB, NB: ", &
           &I0, ", ", I0, ", ", I0, ", ", I0 )', name, rows, cols, &
           actual_block_size, actual_block_size
    end if

    local_rows = max(1, numroc(rows, actual_block_size, &
         proc%my_proc_row, 0, proc%n_procs_row))

    call descinit(desc, rows, cols, actual_block_size, actual_block_size, &
         0, 0, proc%context, local_rows, info)
    if (info /= 0) then
      print *, 'info(descinit): ', info
      call terminate('setup_distributed_matrix: descinit failed', info)
    end if

    allocate(mat(1 : local_rows, 1 : get_local_cols(proc, desc)), stat = info)
    if (info /= 0) then
      print *, 'stat(allocate): ', info
      call terminate('setup_distributed_matrix: allocation failed', info)
    end if

    mat(:, :) = 0.0d0

    time_end = mpi_wtime()
    call add_event('setup_distributed_matrix', time_end - time_start)
  end subroutine setup_distributed_matrix


  subroutine convert_sparse_matrix_to_dense(mat_in, mat)
    use ek_matrix_io_m, only : ek_sparse_mat_t

    type(ek_sparse_mat_t), intent(in) :: mat_in
    double precision, intent(out), allocatable :: mat(:,:)

    integer :: k, i, j, n, ierr
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    n = mat_in%size

    allocate (mat(n, n), stat = ierr)
    if (ierr /= 0) then
      call terminate('convert_sparse_matrix_to_dense: allocation failed', ierr)
    end if

    mat(:, :) = 0.0d0

    do k = 1, mat_in%num_non_zeros
      i = mat_in%suffix(1, k)
      j = mat_in%suffix(2, k)
      mat(i, j) = mat_in%value(k)
      if (i /= j) then
        mat(j, i) = mat_in%value(k)
      endif
    enddo

    time_end = mpi_wtime()
    call add_event('convert_sparse_matrix_to_dense', time_end - time_start)
  end subroutine convert_sparse_matrix_to_dense


  subroutine gather_matrix(mat, desc, dest_row, dest_col, global_mat)
    double precision, intent(in) :: mat(:, :)
    integer, intent(in) :: desc(desc_size), dest_row, dest_col
    double precision, intent(out) :: global_mat(:, :)

    double precision, allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: m, n, mb, nb, m_local, n_local, context, m_recv, n_recv
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global
    integer :: ierr
    double precision :: time_start, time_end

    integer :: numroc, iceil

    time_start = mpi_wtime()

    m = desc(rows_)
    n = desc(cols_)
    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      if (m /= size(global_mat, 1) .or. n /= size(global_mat, 2)) then
        call terminate('gather_matrix: illegal matrix size', 1)
      end if
    end if
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(mat, 1)
    n_local = size(mat, 2)
    context = desc(context_)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      buf_size_row = numroc(m, mb, 0, 0, n_procs_row)
      buf_size_col = numroc(n, nb, 0, 0, n_procs_col)
      allocate(recv_buf(buf_size_row, buf_size_col), stat = ierr)
    else
      allocate(recv_buf(0, 0), stat = ierr) ! Only destination process uses receive buffer
    end if
    if (ierr /= 0) then
      call terminate('gather_matrix: allocation failed', ierr)
    end if

    do pr = 0, n_procs_row - 1
       do pc = 0, n_procs_col - 1
          if (my_proc_row == pr .and. my_proc_col == pc) then
             call dgesd2d(context, m_local, n_local, mat, m_local, dest_row, dest_col)
          end if
          if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
             m_recv = numroc(m, mb, pr, 0, n_procs_row)
             n_recv = numroc(n, nb, pc, 0, n_procs_col)
             call dgerv2d(context, m_recv, n_recv, recv_buf, buf_size_row, pr, pc)
             do br = 1, iceil(iceil(m, mb) - my_proc_row, n_procs_row)
                m1_global = 1 + mb * (n_procs_row * (br - 1) + pr)
                m2_global = min(m1_global + mb - 1, m)
                m1_local = 1 + mb * (br - 1)
                m2_local = m1_local + m2_global - m1_global
                do bc = 1, iceil(iceil(n, nb) - my_proc_col, n_procs_col)
                   n1_global = 1 + nb * (n_procs_col * (bc - 1) + pc)
                   n2_global = min(n1_global + nb - 1, n)
                   n1_local = 1 + nb * (bc - 1)
                   n2_local = n1_local + n2_global - n1_global
                   global_mat(m1_global : m2_global, n1_global : n2_global) = &
                        recv_buf(m1_local : m2_local, n1_local : n2_local)
                end do
             end do
          end if
       end do
    end do

    time_end = mpi_wtime()
    call add_event('gather_matrix', time_end - time_start)
  end subroutine gather_matrix


  subroutine gather_matrix_part(mat, desc, i, j, m, n, dest_row, dest_col, dest_mat)
    real(8), intent(in) :: mat(:, :)
    integer, intent(in) :: i, j, m, n, desc(desc_size), dest_row, dest_col
    real(8), intent(out) :: dest_mat(:, :)

    real(8), allocatable :: recv_buf(:, :)
    integer :: buf_size_row, buf_size_col
    integer :: rows, cols, mb, nb, m_local, n_local, context, m_send, n_send
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: pr, pc, br, bc
    integer :: m1_local, m2_local, m1_global, m2_global
    integer :: n1_local, n2_local, n1_global, n2_global
    integer :: ierr
    real(8) :: time_start, time_end

    integer :: numroc, iceil, indxg2l

    time_start = mpi_wtime()

    rows = desc(rows_)
    cols = desc(cols_)
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(mat, 1)
    n_local = size(mat, 2)
    context = desc(context_)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
      buf_size_row = numroc(rows, mb, 0, 0, n_procs_row)
      buf_size_col = numroc(cols, nb, 0, 0, n_procs_col)
      allocate(recv_buf(buf_size_row, buf_size_col), stat = ierr)
    else
      allocate(recv_buf(0, 0), stat = ierr) ! Only destination process uses receive buffer
    end if
    if (ierr /= 0) then
      call terminate('gather_matrix_part: allocation failed', ierr)
    end if

    do pr = 0, n_procs_row - 1
      do pc = 0, n_procs_col - 1
        m_send = numroc(rows, mb, pr, 0, n_procs_row)
        n_send = numroc(cols, nb, pc, 0, n_procs_col)
        if (my_proc_row == pr .and. my_proc_col == pc) then
          call dgesd2d(context, m_send, n_send, mat, desc(local_rows_), dest_row, dest_col)
        end if
        if (my_proc_row == dest_row .and. my_proc_col == dest_col) then
          call dgerv2d(context, m_send, n_send, recv_buf, buf_size_row, pr, pc)
          do br = 1, iceil(iceil(rows, mb) - pr, n_procs_row)
            m1_global = 1 + mb * (n_procs_row * (br - 1) + pr)
            m2_global = m1_global + mb - 1
            m1_global = max(m1_global, i)
            m2_global = min(m2_global, i + m - 1)
            if (m1_global > m2_global) then
              cycle
            end if
            m1_local = indxg2l(m1_global, mb, 0, 0, n_procs_row)
            m2_local = indxg2l(m2_global, mb, 0, 0, n_procs_row)
            do bc = 1, iceil(iceil(cols, nb) - pc, n_procs_col)
              n1_global = 1 + nb * (n_procs_col * (bc - 1) + pc)
              n2_global = n1_global + nb - 1
              n1_global = max(n1_global, j)
              n2_global = min(n2_global, j + n - 1)
              if (n1_global > n2_global) then
                cycle
              end if
              n1_local = indxg2l(n1_global, nb, 0, 0, n_procs_col)
              n2_local = indxg2l(n2_global, nb, 0, 0, n_procs_col)
              dest_mat(m1_global - i + 1 : m2_global - i + 1, n1_global - j + 1 : n2_global - j + 1) = &
                   recv_buf(m1_local : m2_local, n1_local : n2_local)
            end do
          end do
        end if
      end do
    end do

    time_end = mpi_wtime()
    call add_event('gather_matrix_part', time_end - time_start, .false.)
  end subroutine gather_matrix_part


  subroutine distribute_global_dense_matrix(global_mat, desc, local_mat)
    double precision, intent(in) :: global_mat(:, :) ! assumed to be same in all the procs
    integer, intent(in) :: desc(desc_size)
    double precision, intent(out) :: local_mat(:, :)

    integer :: m, n, mb, nb, m_local, n_local, context
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: m_start, n_start, m_end, n_end, m_start_global, n_start_global
    integer :: i, j
    double precision :: time_start, time_end

    integer :: iceil

    time_start = mpi_wtime()

    m = desc(rows_)
    n = desc(cols_)
    mb = desc(block_row_)
    nb = desc(block_col_)
    m_local = size(local_mat, 1)
    n_local = size(local_mat, 2)
    context = desc(context_)
    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)

    do j = 1, iceil(iceil(n, nb), n_procs_col)
       n_start = 1 + (j - 1) * nb
       if (n_start > n_local) then
          cycle
       end if
       n_end = min(j * nb, n_local)
       n_start_global = 1 + nb * (n_procs_col * (j - 1) + my_proc_col)
       do i = 1, iceil(iceil(m, mb), n_procs_row)
          m_start = 1 + mb * (i - 1)
          if (m_start > m_local) then
             cycle
          end if
          m_end = min(i * mb, m_local)
          m_start_global = 1 + mb * (n_procs_row * (i - 1) + my_proc_row)
          local_mat(m_start : m_end, n_start : n_end) = &
               global_mat(m_start_global : m_start_global + m_end - m_start, &
               n_start_global : n_start_global + n_end - n_start)
       end do
    end do

    time_end = mpi_wtime()
    call add_event('distribute_global_dense_matrix', time_end - time_start)
  end subroutine distribute_global_dense_matrix


  subroutine distribute_global_sparse_matrix(mat_in, desc, mat)
    type(ek_sparse_mat_t), intent(in) :: mat_in
    integer, intent(in) :: desc(desc_size)
    double precision, intent(out) :: mat(:,:)

    integer :: k, i, j
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    do k = 1, mat_in%num_non_zeros
      i = mat_in%suffix(1, k)
      j = mat_in%suffix(2, k)
      call pdelset(mat, i, j, desc, mat_in%value(k))
      if (i /= j) then
      call pdelset(mat, j, i, desc, mat_in%value(k))
      endif
    enddo

    time_end = mpi_wtime()
    call add_event('distribute_global_sparse_matrix', time_end - time_start)
  end subroutine distribute_global_sparse_matrix


  ! Supposing arrays of the same size are distributed by the block cyclic manner
  ! in each processor row, broadcast whole pieces of the array then share the
  ! array globally in the processor row.
  !
  ! array_local: local input, dimension LOCc(N) where N is the dimension of the
  ! global array.
  subroutine allgather_row_wise(array_local, context, block_size, array_global)
    double precision, intent(in) :: array_local(:)
    integer, intent(in) :: context, block_size
    double precision, intent(out) :: array_global(:)

    double precision, allocatable :: send_buf(:)
    integer :: n_global, n_local, max_buf_size, ierr
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col, sender_proc_col
    integer :: block_, head_global, tail_global, head_local
    double precision :: time_start, time_end

    integer :: numroc, iceil

    time_start = mpi_wtime()

    call blacs_gridinfo(context, n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    n_global = size(array_global)
    max_buf_size = numroc(n_global, block_size, 0, 0, n_procs_col)
    allocate(send_buf(max_buf_size), stat = ierr)
    if (ierr /= 0) then
      call terminate('allgather_row_wise: allocation failed', ierr)
    end if

    do sender_proc_col = 0, n_procs_col - 1
      n_local = numroc(n_global, block_size, sender_proc_col, 0, n_procs_col)
      if (my_proc_col == sender_proc_col) then
        if (n_local /= size(array_local)) then
          call terminate('distribute_matrix: wrong local array size', 1)
        end if
        send_buf(1 : n_local) = array_local(1 : n_local)
        call dgebs2d(context, 'Row', ' ', max_buf_size, 1, send_buf, max_buf_size)
      else
        call dgebr2d(context, 'Row', ' ', max_buf_size, 1, send_buf, max_buf_size, &
             my_proc_row, sender_proc_col)
      end if
      ! Slice entries in a process into blocks
      do block_ = 1, iceil(n_local, block_size)
        head_global = 1 + block_size * (n_procs_col * (block_ - 1) + sender_proc_col)
        tail_global = min(head_global + block_size - 1, n_global)
        head_local = 1 + block_size * (block_ - 1)
        array_global(head_global : tail_global) = &
             send_buf(head_local : head_local + tail_global - head_global)
      end do
    end do

  time_end = mpi_wtime()
  call add_event('allgather_row_wise', time_end - time_start)
  end subroutine allgather_row_wise


  subroutine bcast_sparse_matrix(root, mat_info, mat)
    integer, intent(in) :: root
    type(ek_matrix_info_t), intent(inout) :: mat_info
    type(ek_sparse_mat_t), intent(inout) :: mat

    integer :: my_rank, ierr
    double precision :: time_start, time_start_part, time_end

    time_start = mpi_wtime()
    time_start_part = time_start

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    call mpi_bcast(mat_info%rep, 10, mpi_character, root, mpi_comm_world, ierr)
    call mpi_bcast(mat_info%field, 7, mpi_character, root, mpi_comm_world, ierr)
    call mpi_bcast(mat_info%symm, 19, mpi_character, root, mpi_comm_world, ierr)
    call mpi_bcast(mat_info%rows, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(mat_info%cols, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(mat_info%entries, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(mat%size, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(mat%num_non_zeros, 1, mpi_integer, root, mpi_comm_world, ierr)
    if (my_rank /= root) then
      allocate(mat%suffix(2, mat_info%entries), mat%value(mat_info%entries), stat=ierr)
      if (ierr /= 0) then
        call terminate('bcast_sparse_matrix: allocation failed', ierr)
      end if
    end if

    time_end = mpi_wtime()
    call add_event('bcast_sparse_matrix:bcast_aux', time_end - time_start_part)
    time_start_part = time_end

    call mpi_bcast(mat%suffix, 2 * mat_info%entries, mpi_integer, root, mpi_comm_world, ierr)

    time_end = mpi_wtime()
    call add_event('bcast_sparse_matrix:bcast_suffix', time_end - time_start_part)
    time_start_part = time_end

    call mpi_bcast(mat%value, mat_info%entries, mpi_double_precision, root, mpi_comm_world, ierr)

    time_end = mpi_wtime()
    call add_event('bcast_sparse_matrix:bcast_value', time_end - time_start_part)
    call add_event('bcast_sparse_matrix', time_end - time_start)
  end subroutine bcast_sparse_matrix
end module ek_distribute_matrix_m
