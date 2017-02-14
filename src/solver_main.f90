module ek_solver_main_m
  use mpi
  use ek_descriptor_parameters_m
  use ek_command_argument_m, only : ek_argument_t
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_distribute_matrix_m, only : setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use ek_eigenpairs_types_m, only: ek_eigenpairs_types_union_t, ek_eigenpairs_blacs_t
  use ek_event_logger_m, only : add_event
  use ek_generalized_to_standard_m, only : reduce_generalized, recovery_generalized
  use ek_global_variables_m, only : g_block_size
  use ek_processes_m, only : ek_process_t, setup_distribution, print_map_of_grid_to_processes, &
       check_master, terminate

  implicit none

  private
  public :: eigen_solver

contains

  subroutine eigen_solver(arg, matrix_A, eigenpairs, proc, matrix_B)
    use ek_solver_lapack_m, only : eigen_solver_lapack
    use ek_solver_scalapack_all_m, only : eigen_solver_scalapack_all, solve_with_general_scalapack
    use ek_solver_scalapack_select_m, only : eigen_solver_scalapack_select
    use ek_solver_eigenexa_m, only : setup_distributed_matrix_for_eigenexa, &
         eigen_solver_eigenexa, eigen_solver_eigenk, &
         solve_with_general_scalapack_eigenexa, solve_with_general_scalapack_eigenk, &
         solve_with_general_scalapacknew_eigenk
    use ek_solver_elpa_m
    use ek_solver_elpa_eigenexa_m

    type(ek_argument_t) :: arg
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs
    type(ek_process_t), intent(out) :: proc

    integer :: n, desc_A(desc_size), desc_B(desc_size)
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
end module ek_solver_main_m
