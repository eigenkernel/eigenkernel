module solver_main
  use mpi
  use descriptor_parameters
  use command_argument, only : argument
  use matrix_io, only : sparse_mat
  use distribute_matrix, only : setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use eigenpairs_types, only: eigenpairs_types_union, eigenpairs_blacs
  use event_logger_m, only : add_event
  use generalized_to_standard, only : reduce_generalized, recovery_generalized
  use global_variables, only : g_block_size
  use processes, only : process, setup_distribution, print_map_of_grid_to_processes, &
       check_master, terminate

  implicit none

  private
  public :: eigen_solver

contains

  subroutine eigen_solver(arg, matrix_A, eigenpairs, matrix_B)
    use solver_lapack, only : eigen_solver_lapack
    use solver_scalapack_all, only : eigen_solver_scalapack_all, solve_with_general_scalapack
    use solver_scalapack_select, only : eigen_solver_scalapack_select
    use solver_eigenexa, only : setup_distributed_matrix_for_eigenexa, eigen_solver_eigenexa, eigen_solver_eigenk, &
         solve_with_general_scalapack_eigenexa, solve_with_general_scalapack_eigenk, &
         solve_with_general_scalapacknew_eigenk
    use solver_elpa
    use solver_elpa_eigenexa

    type(argument) :: arg
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    integer :: n, desc_A(desc_size), desc_B(desc_size)
    type(process) :: proc
    double precision, allocatable :: matrix_A_dist(:, :), matrix_B_dist(:, :)

    n = arg%matrix_A_info%rows
    if (arg%is_printing_grid_mapping) call print_map_of_grid_to_processes()
    if (arg%block_size > 0) then  ! Do not use the default block size.
      g_block_size = arg%block_size
    end if

    if (trim(arg%solver_type) /= 'lapack') then
      call setup_distribution(proc)
    end if

    select case (trim(arg%solver_type))
    case ('lapack')
      call eigen_solver_lapack(matrix_A, eigenpairs)
    case ('scalapack')
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_all(proc, desc_A, matrix_A_dist, eigenpairs)
    case ('scalapack_select')
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
    case ('general_scalapack')
      call solve_with_general_scalapack(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_scalapack_select')
      call setup_distributed_matrix('A', proc, n, n, desc_A, matrix_A_dist)
      call setup_distributed_matrix('B', proc, n, n, desc_B, matrix_B_dist)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call distribute_global_sparse_matrix(matrix_B, desc_B, matrix_B_dist)
      call reduce_generalized(n, matrix_A_dist, desc_A, matrix_B_dist, desc_B)
      call eigen_solver_scalapack_select(proc, desc_A, matrix_A_dist, &
           arg%n_vec, eigenpairs)
      call recovery_generalized(n, arg%n_vec, matrix_B_dist, desc_B, &
           eigenpairs%blacs%Vectors, eigenpairs%blacs%desc)
    case ('eigensx')
      call setup_distribution(proc)
      call setup_distributed_matrix_for_eigenexa(n, desc_A, matrix_A_dist, eigenpairs)
      call distribute_global_sparse_matrix(matrix_A, desc_A, matrix_A_dist)
      call eigen_solver_eigenexa(matrix_A_dist, desc_A, arg%n_vec, eigenpairs)
    case ('general_scalapack_eigensx')
      call solve_with_general_scalapack_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_scalapack_eigens')
      call solve_with_general_scalapack_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa_scalapack')
      call solve_with_general_elpa_scalapack(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa1')
      call solve_with_general_elpa1(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa2')
      call solve_with_general_elpa2(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa_eigensx')
      call solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_elpa_eigens')
      call solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    case ('general_scalapacknew_eigens')
      call solve_with_general_scalapacknew_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    case default
      call terminate('eigen_solver: Unknown solver', 1)
    end select
  end subroutine eigen_solver
end module solver_main
