module ek_solver_elpa_eigenexa_m
  use ek_global_variables_m, only : g_block_size
  use ek_distribute_matrix_m, only : ek_process_t, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use ek_descriptor_parameters_m
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_processes_m, only : check_master, terminate

  implicit none
  private
  public :: solve_with_general_elpa_eigenexa, solve_with_general_elpa_eigenk

contains

  subroutine solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solve_elpa_eigenexa: ELPA is not supported in this build', 1)
    end subroutine solve_with_general_elpa_eigenexa


    subroutine solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solve_elpa_eigenexa: ELPA is not supported in this build', 1)
  end subroutine solve_with_general_elpa_eigenk
end module ek_solver_elpa_eigenexa_m
