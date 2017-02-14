module ek_solver_lapack_m
  implicit none

  private
  public :: eigen_solver_lapack

contains

  subroutine eigen_solver_lapack(mat, eigenpairs)
   use ek_matrix_io_m, only : ek_sparse_mat_t
   use ek_distribute_matrix_m, only : convert_sparse_matrix_to_dense
   use ek_eigenpairs_types_m, only : ek_eigenpairs_types_union_t

   type(ek_sparse_mat_t), target, intent(in) :: mat
   type(ek_eigenpairs_types_union_t), intent(out) :: eigenpairs

   integer :: n, lda, lwork, info

   double precision, allocatable :: work(:)

   eigenpairs%type_number = 1

   call convert_sparse_matrix_to_dense(mat, eigenpairs%local%vectors)

   n = mat%size
   lda = n
   lwork = n * n  ! Note: (lwork > 3*n-1 ) should be satisfied.

   allocate(eigenpairs%local%values(n), work(lwork))

   call dsyev("V", "U", n, eigenpairs%local%vectors, lda, &
        eigenpairs%local%values, work, lwork, info)
  end subroutine eigen_solver_lapack
end module ek_solver_lapack_m
