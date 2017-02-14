module ek_command_argument_m
  use mpi
  use fson
  use fson_value_m
  use fson_string_m
  use ek_global_variables_m
  use ek_processes_m, only : check_master, terminate
  implicit none

  integer, parameter :: kMaxNumPrintedVecsRanges = 100

  ! Matrix Market format
  type ek_matrix_info_t
    character(len=10) :: rep
    character(len=7) :: field
    character(len=19) :: symm
    integer :: rows, cols, entries
  end type ek_matrix_info_t

  type ek_argument_t
    character(len=256) :: matrix_A_filename = ''
    character(len=256) :: matrix_B_filename = '' ! Empty means standard eigenvalue problem
    character(len=256) :: log_filename = 'log.json'
    type(ek_matrix_info_t) :: matrix_A_info, matrix_B_info
    character(len=256) :: solver_type
    character(len=256) :: output_filename = 'eigenvalues.dat'
    character(len=256) :: ipratios_filename = 'ipratios.dat'
    logical :: is_generalized_problem
    logical :: is_printing_grid_mapping = .false.
    logical :: is_dry_run = .false.
    logical :: is_binary_output = .false.
    integer :: block_size = 0  ! Zero means default block size
    integer :: n_vec = -1  ! The default -1 means 'all the vectors'.
    ! The default 0 means 'do not check any vector'.
    ! If -1 is specified, n_check_vec is set identical with n_vec.
    integer :: n_check_vec = 0
    ! When zero, orthogonality is not evaluated.
    integer :: ortho_check_index_start = 0
    integer :: ortho_check_index_end = 0
    character(len=256) :: eigenvector_dir = '.'
    integer :: num_printed_vecs_ranges = 0 ! Zero means do not print eigenvectors
    integer :: printed_vecs_ranges(2, kMaxNumPrintedVecsRanges)
    integer :: verbose_level = 0
  end type ek_argument_t

  private
  public :: ek_matrix_info_t, ek_argument_t, required_memory, &
       read_command_argument, validate_argument, print_command_argument, fson_setting_add

contains

  subroutine print_help()
    if (check_master()) then
      print *, 'Usage: eigen_test -s <solver_type> <options> <matrix_A> [<matrix_B>]'
      print *, 'Solver types are:'
      print *, '  scalapack (standard)'
      print *, '  scalapack_select (standard, selecting)'
      print *, '  general_scalapack (generalized)'
      print *, '  general_scalapack_select (generalized, selecting)'
      print *, '  eigensx (standard)'
      print *, '  general_scalapack_eigensx (generalized)'
      print *, '  general_scalapack_eigens (generalized)'
      print *, '  general_elpa_scalapack (generalized)'
      print *, '  general_elpa1 (generalized)'
      print *, '  general_elpa2 (generalized)'
      print *, '  general_elpa_eigensx (generalized)'
      print *, '  general_elpa_eigens (generalized)'
      print *, '  general_scalapacknew_eigens (generalized)'
      print *, 'Options are:'
      print *, '  -n <num>  (available with selecting solvers) Compute only &
           &<num> eigenpairs in ascending order of their eigenvalues'
      print *, '  -c <num>  Consider only <num> eigenvectors in residual norm &
           &checking. Default is 0. Set -1 to consider all the vectors'
      print *, '  -o <file>  Set output file name for eigenvalues to <file>'
      print *, '  -i <file>  Set output file name for ipratios to <file>'
      print *, '  -d <dir>  Set output files directory for eigenvectors to <dir>'
      print *, '  -p <num1>,<num2>  Specify range of the number of eigenvectors to be output'
      print *, '  -t <num1>,<num2>  Consider eigenvectors indexed <num1> to <num2>(included) in orthogonality checking'
      print *, '  -l <file>  Set output file name for elapse time log to <file>'
      print *, '  -h  Print this help and exit'
      print *, '  --block-size <n>  Change block size in block cyclic distribution'
      print *, '  --dry-run  Read command arguments and matrix files and instantly exit'
      print *, '  --print-grid-mapping  Print which process is assigned to each coordinate in BLACS grid'
      call flush(6)
    end if
  end subroutine print_help


  subroutine wrap_mminfo(filename, minfo, ierr)
    character(*), intent(in) :: filename
    type(ek_matrix_info_t), intent(out) :: minfo
    integer, intent(out) :: ierr

    integer, parameter :: iunit = 8

    open(unit=iunit, file=filename, status='old', iostat=ierr)
    if (ierr /= 0) then
      return
    end if
    call mminfo(iunit, minfo%rep, minfo%field, minfo%symm, &
         minfo%rows, minfo%cols, minfo%entries, ierr)
    close(iunit)
  end subroutine wrap_mminfo


  subroutine bcast_matrix_info(root, minfo)
    integer, intent(in) :: root
    type(ek_matrix_info_t), intent(inout) :: minfo

    integer :: ierr

    call mpi_bcast(minfo%rep, 10, mpi_character, root, mpi_comm_world, ierr)
    call mpi_bcast(minfo%field, 7, mpi_character, root, mpi_comm_world, ierr)
    call mpi_bcast(minfo%symm, 19, mpi_character, root, mpi_comm_world, ierr)
    call mpi_bcast(minfo%rows, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(minfo%cols, 1, mpi_integer, root, mpi_comm_world, ierr)
    call mpi_bcast(minfo%entries, 1, mpi_integer, root, mpi_comm_world, ierr)
  end subroutine bcast_matrix_info


  subroutine validate_argument(arg)
    type(ek_argument_t), intent(in) :: arg

    integer :: dim, i
    logical :: is_size_valid, is_solver_valid, is_n_vec_valid

    ! Is matrix size appropriate?
    dim = arg%matrix_A_info%rows
    is_size_valid = dim == arg%matrix_A_info%cols
    if (arg%is_generalized_problem) then
      is_size_valid = is_size_valid &
           .and. (dim == arg%matrix_B_info%rows) &
           .and. (dim == arg%matrix_B_info%cols)
    end if
    if (.not. is_size_valid) then
      call terminate('validate_argument: Matrix dimension mismatch', 1)
    end if

    ! Solver type and problem type matched?
    select case (trim(arg%solver_type))
    case ('lapack')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('scalapack')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('scalapack_select')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('general_scalapack')
      is_solver_valid = arg%is_generalized_problem
    case ('general_scalapack_select')
      is_solver_valid = arg%is_generalized_problem
    case ('eigensx')
      is_solver_valid = .not. arg%is_generalized_problem
    case ('general_scalapack_eigensx')
      is_solver_valid = arg%is_generalized_problem
    case ('general_scalapack_eigens')
      is_solver_valid = arg%is_generalized_problem
    case ('general_elpa_scalapack')
      is_solver_valid = arg%is_generalized_problem
    case ('general_elpa1')
      is_solver_valid = arg%is_generalized_problem
    case ('general_elpa2')
      is_solver_valid = arg%is_generalized_problem
    case ('general_elpa_eigensx')
      is_solver_valid = arg%is_generalized_problem
    case ('general_elpa_eigens')
      is_solver_valid = arg%is_generalized_problem
    case ('general_scalapacknew_eigens')
      is_solver_valid = arg%is_generalized_problem
    case default
      is_solver_valid = .false.
      call terminate("validate_argument: Unknown solver '" // &
           trim(arg%solver_type) // "'", 1)
    end select
    if (.not. is_solver_valid) then
      if (arg%is_generalized_problem) then
        call terminate("validate_argument: solver '" // &
             trim(arg%solver_type) // &
             "' is not for generalized eigenvalue problem", 1)
      else
        call terminate("validate_argument: solver '" // &
             trim(arg%solver_type) // &
             "' is not for standard eigenvalue problem", 1)
      end if
    end if

    ! For all eigenpair computation, n_vec = n?
    is_n_vec_valid = .false.
    select case (trim(arg%solver_type))
    case ('scalapack_select')
      is_n_vec_valid = .true.
    case ('general_scalapack_select')
      is_n_vec_valid = .true.
    case default
      is_n_vec_valid = arg%n_vec == dim
    end select
    if (.not. is_n_vec_valid) then
      call terminate("validate_argument: Solver '" // &
           trim(arg%solver_type) // &
           "' does not support partial eigenvalue computation", 1)
    end if

    do i = 1, arg%num_printed_vecs_ranges
      if (arg%printed_vecs_ranges(1, i) < 0 .or. arg%printed_vecs_ranges(2, i) < 0 .or. &
           arg%printed_vecs_ranges(2, i) > arg%n_vec .or. &
           arg%printed_vecs_ranges(1, i) > arg%printed_vecs_ranges(2, i)) then
        call terminate('validate_argument: Specified numbers with -p option are not valid', 1)
      end if
    end do

    if (arg%n_check_vec < 0 .or. arg%n_check_vec > arg%n_vec) then
      call terminate('validate_argument: Specified numbers with -c option are not valid', 1)
    end if

    if (arg%ortho_check_index_start < 0 .or. arg%ortho_check_index_end < 0 .or. &
         arg%ortho_check_index_end > arg%n_vec .or. &
         arg%ortho_check_index_start > arg%ortho_check_index_end) then
      call terminate('validate_argument: Specified numbers with -t option are not valid', 1)
    end if
  end subroutine validate_argument


  double precision function required_memory_lapack(arg)
    type(ek_argument_t), intent(in) :: arg

    double precision :: num_double, dim

    num_double = real(arg%matrix_A_info%entries)
    dim = real(arg%matrix_A_info%rows)
    num_double = num_double + dim * dim
    required_memory_lapack = 8.0d0 * num_double
  end function required_memory_lapack


  double precision function required_memory_parallel_standard(arg)
    ! This is just an approximation and partial eigenvector computation
    ! (means reduced columns of eigenvector storage) is not supported yet
    ! Generalized version below has the same problem
    type(ek_argument_t), intent(in) :: arg

    integer :: my_rank, n_procs
    double precision :: num_double, dim

    num_double = real(arg%matrix_A_info%entries)

    call blacs_pinfo(my_rank, n_procs)
    dim = real(arg%matrix_A_info%rows)
    ! 2 is for the input matrix and eigenvectors
    num_double = num_double + dim * dim * 2.0d0 / real(n_procs)

    required_memory_parallel_standard = 8.0d0 * num_double
  end function required_memory_parallel_standard


  double precision function required_memory_parallel_generalized(arg)
    type(ek_argument_t), intent(in) :: arg

    integer :: my_rank, n_procs
    double precision :: num_double, dim

    num_double = real(arg%matrix_A_info%entries + arg%matrix_B_info%entries)

    call blacs_pinfo(my_rank, n_procs)
    dim = real(arg%matrix_A_info%rows)
    ! 3 is for the input matrices (A and B) and eigenvectors
    num_double = num_double + dim * dim * 3.0d0 / real(n_procs)

    required_memory_parallel_generalized = 8.0d0 * num_double
  end function required_memory_parallel_generalized


  subroutine arg_str_to_printed_vecs_ranges(arg_str, num_printed_vecs_ranges, printed_vecs_ranges)
    character(len=*), intent(in) :: arg_str
    integer, intent(out) :: num_printed_vecs_ranges, printed_vecs_ranges(2, kMaxNumPrintedVecsRanges)

    integer :: i, j, k1, k2
    num_printed_vecs_ranges = 1
    k1 = 1
    do
      i = index(arg_str(k1 :), ',')
      if (i == 0) then
        k2 = len_trim(arg_str)
      elseif (i == 1) then
        call terminate('arg_str_to_printed_vecs_ranges: invalid comma placement', 1)
      else
        k2 = k1 + i - 2
      end if

      j = index(arg_str(k1 : k2), '-')
      if (j == 0) then
        read (arg_str(k1 : k2), *) printed_vecs_ranges(1, num_printed_vecs_ranges)
        printed_vecs_ranges(2, num_printed_vecs_ranges) = printed_vecs_ranges(1, num_printed_vecs_ranges)
      elseif (j == 1) then
        call terminate('arg_str_to_printed_vecs_ranges: invalid hyphen placement', 1)
      else
        read (arg_str(k1 : k1 + j - 2), *) printed_vecs_ranges(1, num_printed_vecs_ranges)
        read (arg_str(k1 + j : k2), *) printed_vecs_ranges(2, num_printed_vecs_ranges)
      end if

      k1 = k2 + 2
      if (k1 >= len_trim(arg_str)) then
        exit
      end if
      num_printed_vecs_ranges = num_printed_vecs_ranges + 1
      if (num_printed_vecs_ranges > kMaxNumPrintedVecsRanges) then
        print *, 'arg_str_to_printed_vecs_ranges: too many ranges ', num_printed_vecs_ranges, &
             ' (> ', kMaxNumPrintedVecsRanges, ')'
        call terminate('arg_str_to_printed_vecs_ranges: too many ranges', 1)
      end if
    end do
    print *, 'arg_str_to_printed_vecs_ranges: num_printed_vecs_ranges ', num_printed_vecs_ranges
    do i = 1, num_printed_vecs_ranges
      print *, 'arg_str_to_printed_vecs_ranges: range ', i, ' : ', &
           printed_vecs_ranges(1, i), ' - ', printed_vecs_ranges(2, i)
    end do
  end subroutine arg_str_to_printed_vecs_ranges


  double precision function required_memory(arg)
    type(ek_argument_t), intent(in) :: arg

    select case (trim(arg%solver_type))
    case ('lapack')
      required_memory = required_memory_lapack(arg)
    case ('scalapack')
      required_memory = required_memory_parallel_standard(arg)
    case ('scalapack_select')
      required_memory = required_memory_parallel_standard(arg)
    case ('general_scalapack')
      required_memory = required_memory_parallel_generalized(arg)
    case ('general_scalapack_select')
      required_memory = required_memory_parallel_generalized(arg)
    case default
      required_memory = -1.0d0 ! Required memory unknown for this solver
    end select
  end function required_memory


  subroutine read_command_argument(arg)
    type(ek_argument_t), intent(out) :: arg

    integer :: argi = 1
    character(len=256) :: arg_str
    integer :: i, ierr, ierr_mpi

    do while (argi <= command_argument_count())
      call get_command_argument(argi, arg_str)

      if (arg_str(1:1) == '-') then
        select case (trim(arg_str(2:)))
        case ('s')
          call get_command_argument(argi + 1, arg_str)
          arg%solver_type = trim(arg_str)
          argi = argi + 1
        case ('n')
          call get_command_argument(argi + 1, arg_str)
          read (arg_str, *) arg%n_vec
          argi = argi + 1
        case ('c')
          call get_command_argument(argi + 1, arg_str)
          read (arg_str, *) arg%n_check_vec
          argi = argi + 1
        case ('o')
          call get_command_argument(argi + 1, arg_str)
          arg%output_filename = trim(arg_str)
          argi = argi + 1
        case ('i')
          call get_command_argument(argi + 1, arg_str)
          arg%ipratios_filename = trim(arg_str)
          argi = argi + 1
        case ('d')
          call get_command_argument(argi + 1, arg_str)
          arg%eigenvector_dir = trim(arg_str)
          argi = argi + 1
        case ('p')
          call get_command_argument(argi + 1, arg_str)
          call arg_str_to_printed_vecs_ranges(arg_str, arg%num_printed_vecs_ranges, arg%printed_vecs_ranges)
          argi = argi + 1
        case ('t')
          call get_command_argument(argi + 1, arg_str)
          i = index(arg_str, ',')
          if (i == 0) then
            call terminate('read_command_argument: wrong format for -t option', 1)
          else
            read (arg_str(1 : i - 1), *) arg%ortho_check_index_start
            read (arg_str(i + 1 :), *) arg%ortho_check_index_end
          end if
          argi = argi + 1
        case ('v')
          arg%verbose_level = 1
        case ('l')
          call get_command_argument(argi + 1, arg_str)
          arg%log_filename = trim(arg_str)
          argi = argi + 1
        case ('h')
          call print_help()
          call terminate('read_command_argument: help printed', 0)
        case ('-block-size')
          call get_command_argument(argi + 1, arg_str)
          read (arg_str, *) arg%block_size
          argi = argi + 1
        case ('-dry-run')
          arg%is_dry_run = .true.
        case ('-print-grid-mapping')
          arg%is_printing_grid_mapping = .true.
        case ('-binary')
          arg%is_binary_output = .true.
        case default
          call print_help()
          call terminate('read_command_argument: unknown option' // &
               trim(arg_str), 1)
        end select
      else if (len_trim(arg%matrix_A_filename) == 0) then
        ! The first non-option argument specifies the (left) input matrix
        arg%matrix_A_filename = trim(arg_str)
      else
        arg%matrix_B_filename = trim(arg_str)
      end if
      argi = argi + 1
    enddo

    if (len_trim(arg%matrix_A_filename) == 0) then
      call terminate('read_command_argument: Matrix A file not specified', 1)
    end if
    arg%is_generalized_problem = (len_trim(arg%matrix_B_filename) /= 0)

    if (check_master()) then
      call wrap_mminfo(arg%matrix_A_filename, arg%matrix_A_info, ierr)
    end if
    call mpi_bcast(ierr, 1, mpi_integer, 0, mpi_comm_world, ierr_mpi)
    if (ierr /= 0) then
      call terminate('mminfo ' // trim(arg%matrix_A_filename) // ' failed',  ierr)
    end if
    call bcast_matrix_info(0, arg%matrix_A_info)

    if (arg%is_generalized_problem) then
      if (check_master()) then
        call wrap_mminfo(arg%matrix_B_filename, arg%matrix_B_info, ierr)
      end if
      call mpi_bcast(ierr, 1, mpi_integer, 0, mpi_comm_world, ierr_mpi)
      if (ierr /= 0) then
        call terminate('mminfo ' // trim(arg%matrix_B_filename) // ' failed',  ierr)
      end if
      call bcast_matrix_info(0, arg%matrix_B_info)
    end if

    if (arg%n_vec == -1) then ! unspecified in command line arguments
      arg%n_vec = arg%matrix_A_info%rows
    end if

    if (arg%n_check_vec == -1) then
      arg%n_check_vec = arg%n_vec
    end if
  end subroutine read_command_argument


  subroutine print_matrix_info(name, info)
    character(*), intent(in) :: name
    type(ek_matrix_info_t), intent(in) :: info

    print "('matrix ', A, ' field: ', A)", name, trim(info%field)
    print "('matrix ', A, ' symm: ', A)", name, trim(info%symm)
    print "('matrix ', A, ' rows: ', I0)", name, info%rows
    print "('matrix ', A, ' cols: ', I0)", name, info%cols
    print "('matrix ', A, ' entries: ', I0)", name, info%entries
  end subroutine print_matrix_info


  subroutine print_command_argument(arg)
    type(ek_argument_t), intent(in) :: arg

    if (arg%is_generalized_problem) then
      print '("problem type: generalized")'
    else
      print '("problem type: standard")'
    end if

    print '("matrix A file: ", a)', trim(arg%matrix_A_filename)
    call print_matrix_info('A', arg%matrix_A_info)

    if (arg%is_generalized_problem) then
      print '("matrix B file: ", a)', trim(arg%matrix_B_filename)
      call print_matrix_info('B', arg%matrix_B_info)
    end if

    print '("solver: ", a)', trim(arg%solver_type)
    print '("eigenvalues output file: ", a)', trim(arg%output_filename)
    print '("ipratios output file: ", a)', trim(arg%ipratios_filename)
    print '("required eigenpairs: ", i0)', arg%n_vec
    print '("verified eigenpairs: ", i0)', arg%n_check_vec
    print '("log output file: ", a)', trim(arg%log_filename)
  end subroutine print_command_argument


  subroutine fson_setting_add(arg, output)
    type(ek_argument_t), intent(in) :: arg
    type(fson_value), pointer, intent(inout) :: output

    type(fson_value), pointer :: setting_in_fson, setting_elem
    type(fson_string), pointer :: str
    character(len=1024) :: argv
    integer :: index_arg

    setting_in_fson => fson_value_create()
    setting_in_fson%value_type = TYPE_OBJECT
    call fson_set_name('setting', setting_in_fson)

    ! Set version.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_OBJECT
    call fson_set_name('version', setting_elem)
    call fson_set_as_string(g_version, setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set command.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_STRING
    call fson_set_name('command', setting_elem)
    str => fson_string_create()
    do index_arg = 0, command_argument_count()
      call getarg(index_arg, argv)
      call fson_string_append(str, trim(argv))
      if (index_arg < command_argument_count()) then
        call fson_string_append(str, ' ')
      end if
    end do
    setting_elem%value_string => str
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set file names.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_OBJECT
    call fson_set_name('matrix_A_filename', setting_elem)
    call fson_set_as_string(trim(arg%matrix_A_filename), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_OBJECT
    call fson_set_name('matrix_B_filename', setting_elem)
    call fson_set_as_string(trim(arg%matrix_B_filename), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_OBJECT
    call fson_set_name('log_filename', setting_elem)
    call fson_set_as_string(trim(arg%log_filename), setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set dimension.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('dimension', setting_elem)
    setting_elem%value_integer = arg%matrix_A_info%rows
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set solver.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_OBJECT
    call fson_set_name('solver', setting_elem)
    call fson_set_as_string(arg%solver_type, setting_elem)
    call fson_value_add(setting_in_fson, setting_elem)

    ! Set blocksize.
    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('g_block_size', setting_elem)
    setting_elem%value_integer = g_block_size
    call fson_value_add(setting_in_fson, setting_elem)

    setting_elem => fson_value_create()
    setting_elem%value_type = TYPE_INTEGER
    call fson_set_name('block_size', setting_elem)
    setting_elem%value_integer = arg%block_size
    call fson_value_add(setting_in_fson, setting_elem)

    call fson_value_add(output, setting_in_fson)
  end subroutine fson_setting_add
end module ek_command_argument_m
