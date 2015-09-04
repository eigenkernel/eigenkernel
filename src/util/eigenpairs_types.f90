module eigenpairs_types
  type eigenpairs_local ! type number = 1
    double precision, allocatable :: values(:)
    double precision, allocatable :: vectors(:, :)
  end type eigenpairs_local

  type eigenpairs_blacs ! type number = 2
    double precision, allocatable :: values(:)
    integer :: desc(9)
    double precision, allocatable :: Vectors(:, :)
  end type eigenpairs_blacs

  type eigenpairs_types_union
    integer :: type_number
    type(eigenpairs_local) :: local
    type(eigenpairs_blacs) :: blacs
  end type eigenpairs_types_union
end module eigenpairs_types
