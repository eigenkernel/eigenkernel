module ek_eigenpairs_types_m
  type ek_eigenpairs_local_t ! type number = 1
    double precision, allocatable :: values(:)
    double precision, allocatable :: vectors(:, :)
  end type ek_eigenpairs_local_t

  type ek_eigenpairs_blacs_t ! type number = 2
    double precision, allocatable :: values(:)
    integer :: desc(9)
    double precision, allocatable :: Vectors(:, :)
  end type ek_eigenpairs_blacs_t

  type ek_eigenpairs_types_union_t
    integer :: type_number
    type(ek_eigenpairs_local_t) :: local
    type(ek_eigenpairs_blacs_t) :: blacs
  end type ek_eigenpairs_types_union_t
end module ek_eigenpairs_types_m
