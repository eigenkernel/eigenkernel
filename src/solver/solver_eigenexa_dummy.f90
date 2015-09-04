module solver_eigenexa
  use distribute_matrix, only : desc_size
  use processes, only : process, terminate
  use matrix_io, only : sparse_mat
  use eigenpairs_types, only : eigenpairs_types_union

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

    call terminate('setup_distributed_matrix_for_eigenexa: EigenExa is not supported in this build', 1)
  end subroutine setup_distributed_matrix_for_eigenexa


  subroutine eigen_solver_eigenexa(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    call terminate('eigen_solver_eigenexa: EigenExa is not supported in this build', 1)
  end subroutine eigen_solver_eigenexa


  subroutine eigen_solver_eigenk(mat, desc_mat, n_vec, eigenpairs, uplo)
    double precision, intent(inout) :: mat(:, :)
    integer, intent(in) :: desc_mat(desc_size), n_vec
    type(eigenpairs_types_union), intent(inout) :: eigenpairs
    character, intent(in), optional :: uplo

    call terminate('eigen_solver_eigenk: EigenExa is not supported in this build', 1)
  end subroutine eigen_solver_eigenk


  subroutine solve_with_general_scalapack_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in) :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('solve_with_general_scalapack_eigenexa: EigenExa is not supported in this build', 1)
  end subroutine solve_with_general_scalapack_eigenexa


  subroutine solve_with_general_scalapack_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in) :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('solve_with_general_scalapack_eigenk: EigenExa is not supported in this build', 1)
  end subroutine solve_with_general_scalapack_eigenk


  subroutine solve_with_general_scalapacknew_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in) :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('solve_with_general_scalapacknew_eigenk: EigenExa is not supported in this build', 1)
  end subroutine solve_with_general_scalapacknew_eigenk
end module solver_eigenexa
