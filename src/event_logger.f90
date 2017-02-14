module ek_event_logger_m
  use mpi
  use fson
  use fson_value_m
  use fson_string_m
  use ek_global_variables_m
  implicit none

  type event_t
    type(event_t), pointer :: next_node => null()
    character(len=128) :: name
    integer :: num_repeated
    real(8) :: val
  end type event_t

  private
  type(event_t), pointer :: events => null()

  public :: add_event, num_events, print_events, fson_events_add

contains

  subroutine add_event(name, val, to_print)
    character(*), intent(in) :: name
    real(8), intent(in) :: val
    logical, optional, intent(in) :: to_print  ! Print log to stderr in default.
    type(event_t), pointer :: new_event, p

    integer :: my_rank, ierr
    double precision :: t
    logical :: is_found, to_print_actual

    if (present(to_print)) then
      to_print_actual = to_print
    else
      to_print_actual = .true.
    end if

    call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
    t = mpi_wtime() - g_mpi_wtime_init
    if (to_print_actual .and. my_rank == 0) then
      write (0, '(A, F16.6, 3A, E24.16e3)') '[Event', t, '] ', name, ',', val
    end if

    is_found = .false.
    p => events
    do while (associated(p))
      if (trim(name) == trim(p%name)) then
        p%num_repeated = p%num_repeated + 1
        p%val = p%val + val
        is_found = .true.
        exit
      end if
      p => p%next_node
    end do

    if (.not. is_found) then
      allocate(new_event)
      new_event%next_node => events
      new_event%name = name
      new_event%num_repeated = 1
      new_event%val = val
      events => new_event
    end if
  end subroutine add_event


  integer function num_events()
    type(event_t), pointer :: p

    num_events = 0
    p => events
    do while (associated(p))
      num_events = num_events + 1
      p => p%next_node
    end do
  end function num_events


  subroutine print_events()
    integer :: i, n
    character(len=128), allocatable :: names(:)
    integer, allocatable :: nums_repeated(:)
    double precision, allocatable :: vals(:)
    type(event_t), pointer :: p

    n = num_events()
    allocate(names(n), nums_repeated(n), vals(n))

    p => events
    do i = 1, n
      names(n - i + 1) = p%name
      nums_repeated(n - i + 1) = p%num_repeated
      vals(n - i + 1) = p%val
      p => p%next_node
    end do

    do i = 1, n
      print *, trim(names(i)), nums_repeated(i), vals(i)
    end do
  end subroutine print_events


  subroutine fson_events_add(output)
    type(fson_value), pointer, intent(in) :: output

    type(fson_value), pointer :: events_in_fson, event, event_elem
    type(event_t), pointer :: p

    events_in_fson => fson_value_create()
    events_in_fson%value_type = TYPE_ARRAY
    call fson_set_name('events', events_in_fson)

    p => events
    do while (associated(p))
      event => fson_value_create()
      event%value_type = TYPE_OBJECT

      event_elem => fson_value_create()
      call fson_set_name('name', event_elem)
      call fson_set_as_string(trim(p%name), event_elem)
      call fson_value_add(event, event_elem)

      event_elem => fson_value_create()
      event_elem%value_type = TYPE_INTEGER
      call fson_set_name('num_repeated', event_elem)
      event_elem%value_integer = p%num_repeated
      call fson_value_add(event, event_elem)

      event_elem => fson_value_create()
      event_elem%value_type = TYPE_REAL
      call fson_set_name('val', event_elem)
      event_elem%value_real = p%val
      call fson_value_add(event, event_elem)

      call fson_value_add(events_in_fson, event)
      p => p%next_node
    end do

    call fson_value_add(output, events_in_fson)
  end subroutine fson_events_add
end module ek_event_logger_m
