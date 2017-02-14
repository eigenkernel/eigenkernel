module ek_matrix_io_m
  use mpi
  use ek_command_argument_m, only : ek_argument_t, ek_matrix_info_t
  use ek_descriptor_parameters_m
  use ek_event_logger_m, only : add_event
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t
  use ek_global_variables_m
  use ek_processes_m, only : check_master, terminate
  implicit none

  type ek_sparse_mat_t
    integer :: size, num_non_zeros
    double precision, allocatable :: value(:)
    integer, allocatable :: suffix(:, :)
  end type ek_sparse_mat_t

  private
  public :: read_matrix_file, print_matrix, ek_sparse_mat_t, print_eigenvectors

contains

  subroutine read_matrix_file(filename, info, matrix, ierr)
    character(*), intent(in) :: filename
    type(ek_matrix_info_t), intent(in) :: info
    type(ek_sparse_mat_t), intent(out) :: matrix
    integer, intent(out) :: ierr

    double precision :: time_start, time_start_part, time_end
    integer, parameter :: iunit = 8
    integer :: rows, cols, num_non_zeros

    time_start = mpi_wtime()
    time_start_part = time_start

    if (check_master()) then
      print '("start reading matrix file ", a)', trim(filename)
    end if

    matrix%size = info%rows
    matrix%num_non_zeros = info%entries

    allocate(matrix%suffix(2, info%entries), matrix%value(info%entries), stat = ierr)
    if (ierr /= 0) then
      return
    end if

    time_end = mpi_wtime()
    call add_event('read_matrix_file:allocate', time_end - time_start_part)
    time_start_part = time_end

    open(iunit, file=filename, status='old', iostat=ierr)
    if (ierr /= 0) then
      return
    end if
    ! read_matrix_file_header is added to skip comment lines
    call read_matrix_file_header(iunit, rows, cols, num_non_zeros)

    time_end = mpi_wtime()
    call add_event('read_matrix_file:header', time_end - time_start_part)
    time_start_part = time_end

    call read_matrix_file_value(0, iunit, info%rows, &
         info%entries, matrix%value, matrix%suffix)
    close(iunit)

    time_end = mpi_wtime()
    call add_event('read_matrix_file:value', time_end - time_start_part)
    call add_event('read_matrix_file', time_end - time_start)
  end subroutine read_matrix_file


  subroutine read_matrix_file_header(unit_num, rows, cols, num_non_zeros)
    integer, intent(in) :: unit_num

    integer :: ierr, rows, cols, num_non_zeros
    character(len=1024) :: chara_wrk

    ! Read the first line
    read(unit_num,'(a)') chara_wrk

    ! Plot the comment lines
    do while (.true.)
      read(unit_num,'(a)') chara_wrk
      if (index(chara_wrk, '%') /= 1) exit
    enddo

    read (chara_wrk, *, iostat = ierr) rows, cols, num_non_zeros
  end subroutine read_matrix_file_header


  subroutine read_matrix_file_value(verbose_level, unit_num, mat_size, num_non_zeros, mat_value, mat_suffix)
    integer, intent(in) :: verbose_level
    integer, intent(in) :: unit_num
    integer, intent(in) :: mat_size
    integer, intent(in) :: num_non_zeros
    double precision, allocatable, intent(inout) :: mat_value(:)
    integer, allocatable, intent(inout) :: mat_suffix(:,:)

    integer :: line_count, print_count
    logical :: debug_mode
    integer :: i, j, ierr
    double precision :: value_wrk

    if (verbose_level >= 100) then
      debug_mode = .true.
    else
      debug_mode = .false.
    endif

    if (debug_mode) write(*,'(a)')'@@ read_matrix_file_value'
    if (debug_mode) write(*,'(a)')' ONLY REAL-SYMMETRIC MATRIX is supported'

    if (debug_mode) write(*,*)'size(mat_value)=',size(mat_value)

    if (debug_mode) write(*,*)'size(mat_suffix,1)=',size(mat_suffix,1)
    if (debug_mode) write(*,*)'size(mat_suffix,2)=',size(mat_suffix,2)

    print_count = 1
    do line_count = 1, num_non_zeros
      if (line_count > num_non_zeros / 10 * print_count .and. check_master()) then
         write (0, '(A, F16.6, A, I0)') &
              '[Event', mpi_wtime() - g_mpi_wtime_init, '] read matrix element number ', line_count
         print_count = print_count + 1
      end if

      read(unit_num, *, iostat=ierr) i, j, value_wrk
      if (ierr /= 0) then
        call terminate('read_matrix_file_value: invalid format of matrix value', ierr)
      endif

      ! if (debug_mode) write(*,*)'i,j, matrix_value=',i, j, value_wrk

      if (i < 1 .or. i > mat_size .or. j < 1 .or. j > mat_size) then
        call terminate('read_matrix_file_value: index of matrix out of range', ierr)
      endif

      mat_value(line_count) = value_wrk
      mat_suffix(1,line_count) = i
      mat_suffix(2,line_count) = j
    enddo

    ! if (num_non_zeros /= size(mat_value,1)) then
    ! endif
  end subroutine read_matrix_file_value


  subroutine print_matrix(name, mat, m, n)
    character(*), intent(in) :: name
    double precision, intent(in) :: mat(:, :)

    integer :: i, j, m, n
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    if (m < 0) then
       m = size(mat, 1)
    end if
    if (n < 0) then
       n = size(mat, 2)
    end if
    do j = 1, n
       do i = 1, m
          print ' (A, "(", I6, ",", I6, ")=", D30.18) ', name, i, j, mat(i, j)
       end do
    end do

    time_end = mpi_wtime()
    call add_event('print_matrix', time_end - time_start)
  end subroutine print_matrix


  subroutine print_eigenvectors(arg, eigenpairs)
    type(ek_argument_t) :: arg
    type(ek_eigenpairs_types_union_t) :: eigenpairs

    double precision :: time_start, time_end
    integer :: len, i, j, digit, stat, desc(desc_size)
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col
    integer :: print_pcol, print_prow
    character(512) :: num_str, filename
    integer, parameter :: iunit = 10, max_num_digits = 8
    integer :: indxg2p  ! Functions.

    time_start = mpi_wtime()

    if (eigenpairs%type_number == 1) then
      stop '[Error] print_eigenvectors: printer for a local matrix not implemented yet'
    else if (eigenpairs%type_number == 2) then
      desc(:) = eigenpairs%blacs%desc(:)
      call blacs_gridinfo(desc(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
      do i = 1, arg%num_printed_vecs_ranges
        do j = arg%printed_vecs_ranges(1, i), arg%printed_vecs_ranges(2, i)
          print_pcol = indxg2p(j, desc(block_col_), 0, desc(csrc_), n_procs_col)
          print_prow = mod(mod(j - 1, desc(block_col_) * n_procs_col), n_procs_row)
          if (my_proc_col == print_pcol .and. my_proc_row == print_prow) then
            write (0, '(A, F16.6, A, I0, A, I0, A, I0, A)') &
                 '[Event', mpi_wtime() - g_mpi_wtime_init, '] print eigenvector ', j, &
                 ' on process (', my_proc_row, ', ', my_proc_col, ')'
            write (num_str, '(i0)') j
            len = len_trim(num_str)
            num_str(max_num_digits - len + 1 : max_num_digits) = num_str(1 : len)
            do digit = 1, max_num_digits - len
              num_str(digit:digit) = '0'
            end do
            filename = trim(arg%eigenvector_dir) // '/' // trim(num_str) // '.dat'
            if (arg%is_binary_output) then
              open(iunit, file=trim(filename), form="unformatted", access="sequential", &
                   status='replace', iostat=stat)
            else
              open(iunit, file=trim(filename), status='replace', iostat=stat)
            end if
            if (stat /= 0) then
              print *, 'iostat: ', stat
              call terminate('print_eigenvectors: cannot open ' // trim(filename), stat)
            end if
          end if

          call print_vector(eigenpairs%blacs%Vectors, j, desc, print_prow, print_pcol, iunit, arg%is_binary_output)

          if (my_proc_col == print_pcol .and. my_proc_row == print_prow) then
            close (iunit)
          end if
        end do
      end do
    end if

    time_end = mpi_wtime()
    call add_event('print_eigenvectors', time_end - time_start)
  end subroutine print_eigenvectors


  subroutine print_vector(A, j, desc_A, irprnt, icprnt, nout, is_binary)
    integer, intent(in) :: j, desc_A(desc_size), irprnt, icprnt, nout
    real(8), intent(in) :: A(:, :)
    logical, intent(in) :: is_binary
    integer :: n_procs_row, n_procs_col, my_proc_row, my_proc_col, own_proc_col
    integer :: m, mb, nb, block_num, local_j, i
    integer :: global_top, global_bot, local_top, local_bot
    real(8) :: work(desc_A(rows_))
    integer :: indxg2p, indxg2l, iceil

    m = desc_A(rows_)
    mb = desc_A(block_row_)
    nb = desc_A(block_col_)
    work(:) = 0d0
    call blacs_gridinfo(desc_A(context_), n_procs_row, n_procs_col, my_proc_row, my_proc_col)
    own_proc_col = indxg2p(j, nb, 0, desc_A(csrc_), n_procs_col)
    local_j = indxg2l(j, nb, 0, 0, n_procs_col)
    ! Exit if it is an unrelated node.
    if (my_proc_col /= own_proc_col .and. .not. &
         (my_proc_row == irprnt .and. my_proc_col == icprnt)) then
      return
    end if
    ! Copy local sub-elements of the vector to a global array to be summed.
    if (my_proc_col == own_proc_col) then
      do block_num = 0, iceil(m, mb) - 1
        global_top = mb * block_num + 1
        global_bot = min(global_top + mb - 1, m)
        if (my_proc_row == mod(block_num, n_procs_row)) then
          local_top = indxg2l(global_top, mb, 0, 0, n_procs_row)
          local_bot = indxg2l(global_bot, mb, 0, 0, n_procs_row)
          work(global_top : global_bot) = A(local_top : local_bot, local_j)
        end if
      end do
      call dgsum2d(desc_A(context_), 'Column', ' ', m, 1, work, m, irprnt, own_proc_col)
    end if
    ! Copy the complete vector to the printing node if needed.
    if (own_proc_col /= icprnt) then
      if (my_proc_row == irprnt .and. my_proc_col == own_proc_col) then
        call dgesd2d(desc_A(context_), m, 1, work, m, irprnt, icprnt)
      else if (my_proc_row == irprnt .and. my_proc_col == icprnt) then
        call dgerv2d(desc_A(context_), m, 1, work, m, irprnt, own_proc_col)
      end if
    end if
    if (my_proc_row == irprnt .and. my_proc_col == icprnt) then
       if (is_binary) then
          write(nout) work(1 : m)
       else
          do i = 1, m
             write(nout, "(I8, ' ', I8, ' ', E26.16e3)") i, j, work(i)
          end do
       end if
    end if
  end subroutine print_vector
end module ek_matrix_io_m
