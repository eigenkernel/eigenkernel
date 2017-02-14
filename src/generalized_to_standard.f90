module ek_generalized_to_standard_m
  use mpi
  use ek_descriptor_parameters_m
  use ek_event_logger_m, only : add_event
  use ek_processes_m, only : check_master, terminate
  implicit none

  private
  public :: reduce_generalized, reduce_generalized_new, recovery_generalized

contains

  subroutine reduce_generalized(dim, A, desc_A, B, desc_B)
    integer, intent(in) :: dim, desc_A(9), desc_B(9)
    double precision, intent(inout) :: A(:, :), B(:, :)

    integer :: info
    double precision :: scale
    double precision :: time_start, time_end

    time_start = mpi_wtime()

    ! B = LL', overwritten to B
    call pdpotrf('L', dim, B, 1, 1, desc_B, info)
    if (info /= 0) then
      if (check_master()) then
        print '("info(pdpotrf): ", i0)', info
      end if
      call terminate('reduce_generalized: pdpotrf failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('reduce_generalized:pdpotrf', time_end - time_start)
    time_start = time_end

    ! Reduction to standard problem by A <- L^(-1) * A * L'^(-1)
    call pdsygst(1, 'L', dim, A, 1, 1, desc_A, B, 1, 1, desc_B, scale, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdsygst): ", i0)', info
      call terminate('reduce_generalized: pdsygst failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('reduce_generalized:pdsygst', time_end - time_start)
  end subroutine reduce_generalized


  subroutine reduce_generalized_new(dim, A, desc_A, B, desc_B)
    integer, intent(in) :: dim, desc_A(9), desc_B(9)
    double precision, intent(inout) :: A(:, :), B(:, :)

    integer :: nprow, npcol, myrow, mycol, info
    integer :: nb, np0, nq0, lwork
    double precision :: scale
    double precision, allocatable :: work(:)
    double precision :: time_start, time_end
    integer :: numroc

    time_start = mpi_wtime()

    ! B = LL', overwritten to B
    call pdpotrf('L', dim, B, 1, 1, desc_B, info)
    if (info /= 0) then
      if (check_master()) then
        print '("info(pdpotrf): ", i0)', info
      end if
      call terminate('reduce_generalized_new: pdpotrf failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('reduce_generalized_new:pdpotrf', time_end - time_start)
    time_start = time_end

    ! Reduction to standard problem by A <- L^(-1) * A * L'^(-1)
    call blacs_gridinfo(desc_A(context_), nprow, npcol, myrow, mycol)
    nb = desc_A(block_row_)
    np0 = numroc(dim, nb, 0, 0, nprow)
    nq0 = numroc(dim, nb, 0, 0, npcol)
    lwork = 2 * np0 * nb + nq0 * nb + nb * nb
    allocate(work(lwork))
    call pdsyngst(1, 'L', dim, A, 1, 1, desc_A, B, 1, 1, desc_B, scale, work, lwork, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdsyngst): ", i0)', info
      call terminate('reduce_generalized_new: pdsyngst failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('reduce_generalized_new:pdsyngst', time_end - time_start)
  end subroutine reduce_generalized_new


  subroutine recovery_generalized(dim, n_vec, B, desc_B, Vectors, desc_Vectors)
    integer, intent(in) :: dim, n_vec, desc_B(9), desc_Vectors(9)
    double precision, intent(in) :: B(:, :)
    double precision, intent(inout) :: Vectors(:, :)

    double precision :: time_start, time_end
    integer :: info

    time_start = mpi_wtime()

    ! Recovery eigenvectors by V <- L'^(-1) * V
    call pdtrtrs('L', 'T', 'N', dim, n_vec, B, 1, 1, desc_B, &
         Vectors, 1, 1, desc_Vectors, info)
    if (info /= 0) then
      if (check_master()) print '("info(pdtrtrs): ", i0)', info
      call terminate('reduce_generalized: pdtrtrs failed', info)
    end if

    time_end = mpi_wtime()
    call add_event('recovery_generalized', time_end - time_start)
  end subroutine recovery_generalized
end module ek_generalized_to_standard_m
