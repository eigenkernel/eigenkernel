module ek_solver_elpa_m
  use ek_distribute_matrix_m, only : ek_process_t
  use ek_descriptor_parameters_m
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_processes_m, only : check_master, terminate

  implicit none
  private
  public :: solve_with_general_elpa_scalapack, solve_with_general_elpa1, solve_with_general_elpa2

contains

  subroutine solve_with_general_elpa_scalapack(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solver_elpa: ELPA is not supported in this build', 1)
  end subroutine solve_with_general_elpa_scalapack


  subroutine solve_with_general_elpa1(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solver_elpa: ELPA is not supported in this build', 1)
  end subroutine solve_with_general_elpa1


  subroutine solve_with_general_elpa2(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in), optional :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solver_elpa: ELPA is not supported in this build', 1)
  end subroutine solve_with_general_elpa2
end module ek_solver_elpa_m
