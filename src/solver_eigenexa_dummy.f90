module ek_solver_eigenexa_m
  use ek_distribute_matrix_m, only : desc_size
  use ek_processes_m, only : ek_process_t, terminate
  use ek_matrix_io_m, only : ek_sparse_mat_t
  use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t

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
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('setup_distributed_matrix_for_eigenexa: EigenExa is not supported in this build', 1)
  end subroutine setup_distributed_matrix_for_eigenexa


  subroutine eigen_solver_eigenexa(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(ek_eigenpairs_types_union_t), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    call terminate('eigen_solver_eigenexa: EigenExa is not supported in this build', 1)
  end subroutine eigen_solver_eigenexa


  subroutine eigen_solver_eigenk(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(ek_eigenpairs_types_union_t), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    call terminate('eigen_solver_eigenk: EigenExa is not supported in this build', 1)
  end subroutine eigen_solver_eigenk


  subroutine solve_with_general_scalapack_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in) :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solve_with_general_scalapack_eigenexa: EigenExa is not supported in this build', 1)
  end subroutine solve_with_general_scalapack_eigenexa


  subroutine solve_with_general_scalapack_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in) :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solve_with_general_scalapack_eigenk: EigenExa is not supported in this build', 1)
  end subroutine solve_with_general_scalapack_eigenk


  subroutine solve_with_general_scalapacknew_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(ek_process_t), intent(in) :: proc
    type(ek_sparse_mat_t), intent(in) :: matrix_A
    type(ek_sparse_mat_t), intent(in) :: matrix_B
    type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

    call terminate('solve_with_general_scalapacknew_eigenk: EigenExa is not supported in this build', 1)
  end subroutine solve_with_general_scalapacknew_eigenk
end module ek_solver_eigenexa_m
