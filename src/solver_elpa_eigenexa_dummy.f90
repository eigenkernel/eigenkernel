module solver_elpa_eigenexa
  use global_variables, only : g_block_size
  use distribute_matrix, only : process, setup_distributed_matrix, &
       gather_matrix, distribute_global_sparse_matrix
  use descriptor_parameters
  use eigenpairs_types, only : eigenpairs_types_union
  use matrix_io, only : sparse_mat
  use processes, only : check_master, terminate

  implicit none
  private
  public :: solve_with_general_elpa_eigenexa, solve_with_general_elpa_eigenk

contains

  subroutine solve_with_general_elpa_eigenexa(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('solve_elpa_eigenexa: ELPA is not supported in this build', 1)
    end subroutine solve_with_general_elpa_eigenexa


    subroutine solve_with_general_elpa_eigenk(n, proc, matrix_A, eigenpairs, matrix_B)
    integer, intent(in) :: n
    type(process), intent(in) :: proc
    type(sparse_mat), intent(in) :: matrix_A
    type(sparse_mat), intent(in), optional :: matrix_B
    type(eigenpairs_types_union), intent(out) :: eigenpairs

    call terminate('solve_elpa_eigenexa: ELPA is not supported in this build', 1)
  end subroutine solve_with_general_elpa_eigenk
end module solver_elpa_eigenexa
