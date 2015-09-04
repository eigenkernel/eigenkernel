! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

!
! File:   string.f95
! Author: josephalevin
!
! Created on March 7, 2012, 7:40 PM
!
! Modified on April 23, 2015, EigenKernel developer team

module fson_string_m
    implicit none
    private

    public :: fson_string, fson_string_create, fson_string_destroy
    public :: fson_string_length, fson_string_append, fson_string_clear
    public :: fson_string_equals, fson_string_copy

    integer, parameter :: BLOCK_SIZE = 32

    type fson_string
        character (len = BLOCK_SIZE) :: chars
        integer :: index = 0
        type(fson_string), pointer :: next => null()
    end type fson_string

    interface fson_string_append
        module procedure append_chars, append_string
    end interface fson_string_append

    interface fson_string_copy
        module procedure copy_chars
    end interface fson_string_copy

    interface fson_string_equals
        module procedure equals_string
    end interface fson_string_equals

    interface fson_string_length
        module procedure string_length
    end interface fson_string_length

contains

    !
    ! FSON STRING CREATE
    !
    function fson_string_create(chars) result(new)
        character(len=*), optional :: chars
        type(fson_string), pointer :: new

        allocate(new)

        ! append chars if available
        if(present(chars)) then
            call append_chars(new, chars)
        end if
    end function fson_string_create

    !
    ! FSON STRING CREATE
    !
    recursive subroutine fson_string_destroy(this)
        type(fson_string), pointer :: this

        if(associated(this % next)) then
            call fson_string_destroy(this % next)
        end if

        nullify (this % next)
        nullify (this)
    end subroutine fson_string_destroy

    !
    ! ALLOCATE BLOCK
    !
    subroutine allocate_block(this)
        type(fson_string), pointer :: this
        type(fson_string), pointer :: new

        if (.not.associated(this % next)) then
            allocate(new)
            this % next => new
        end if
    end subroutine allocate_block

    !
    ! APPEND_STRING
    !
    subroutine append_string(str1, str2)
        type(fson_string), pointer :: str1, str2
        integer length, i

        length = string_length(str2)

        do i = 1, length
            call append_char(str1, get_char_at(str2, i))
        end do
    end subroutine append_string

    !
    ! APPEND_CHARS
    !
    subroutine append_chars(str, c)
        type(fson_string), pointer :: str
        character (len = *), intent(in) :: c
        integer length, i

        length = len(c)

        do i = 1, length
            call append_char(str, c(i:i))
        end do
    end subroutine append_chars

    !
    ! APPEND_CHAR
    !
    recursive subroutine append_char(str, c)
        type(fson_string), pointer :: str
        character, intent(in) :: c

        if (str % index .GE. BLOCK_SIZE) then
            !set down the chain
            call allocate_block(str)
            call append_char(str % next, c)
        else
            ! set local
            str % index = str % index + 1
            str % chars(str % index:str % index) = c
        end if
    end subroutine append_char

    !
    ! COPY CHARS
    !
    subroutine copy_chars(this, to)
        type(fson_string), pointer :: this
        character(len = *), intent(inout) :: to
        integer :: length, i

        length = min(string_length(this), len(to))

        do i = 1, length
            to(i:i) = get_char_at(this, i)
        end do

        ! pad with nothing
        do i = length + 1, len(to)
            to(i:i) = ""
        end do
    end subroutine copy_chars

    !
    ! CLEAR
    !
    recursive subroutine fson_string_clear(this)
        type(fson_string), pointer :: this

        if (associated(this % next)) then
            call fson_string_clear(this % next)
            deallocate(this % next)
            nullify (this % next)
        end if

        this % index = 0
    end subroutine fson_string_clear

    !
    ! SIZE
    !
    recursive integer function string_length(str) result(count)
        type(fson_string), pointer :: str

        count = str % index

        if (str % index == BLOCK_SIZE .AND. associated(str % next)) then
            count = count + string_length(str % next)
        end if
    end function string_length

    !
    ! GET CHAR AT
    !
    recursive character function get_char_at(this, i) result(c)
        type(fson_string), pointer :: this
        integer, intent(in) :: i

        if (i .LE. this % index) then
            c = this % chars(i:i)
        else
            c = get_char_at(this % next, i - this % index)
        end if
    end function get_char_at

    !
    ! EQUALS STRING
    !
    logical function equals_string(this, other) result(equals)
        type(fson_string), pointer :: this, other
        integer :: i
        equals = .false.

        if(fson_string_length(this) .ne. fson_string_length(other)) then
            equals = .false.
            return
        else if(fson_string_length(this) == 0) then
            equals = .true.
            return
        end if

        do i=1, fson_string_length(this)
            if(get_char_at(this, i) .ne. get_char_at(other, i)) then
                equals = .false.
                return
            end if
        end do

        equals = .true.
    end function equals_string
end module fson_string_m


! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

!
! File:   value_m.f95
! Author: josephalevin
!
! Created on March 7, 2012, 10:14 PM
!
! Modified on April 23, 2015, EigenKernel developer team

module fson_value_m
    use fson_string_m

    implicit none

    private

    public :: fson_value, fson_value_create, fson_value_destroy
    public :: fson_value_add, fson_value_get, fson_value_count, fson_value_print

    !constants for the value types
    integer, public, parameter :: TYPE_UNKNOWN = -1
    integer, public, parameter :: TYPE_NULL = 0
    integer, public, parameter :: TYPE_OBJECT = 1
    integer, public, parameter :: TYPE_ARRAY = 2
    integer, public, parameter :: TYPE_STRING = 3
    integer, public, parameter :: TYPE_INTEGER = 4
    integer, public, parameter :: TYPE_REAL = 5
    integer, public, parameter :: TYPE_LOGICAL = 6
    integer, public, parameter :: TYPE_REAL_ARRAY = 7

    !
    ! FSON VALUE
    !
    type fson_value
        type(fson_string), pointer :: name => null()
        integer :: value_type = TYPE_UNKNOWN
        logical :: value_logical
        integer :: value_integer
        double precision :: value_real
        double precision, allocatable :: value_real_array(:)
        type(fson_string), pointer :: value_string => null()
        type(fson_value), pointer :: next => null()
        type(fson_value), pointer :: parent => null()
        type(fson_value), pointer :: children => null()
    end type fson_value

    !
    ! FSON VALUE GET
    !
    ! Use either a 1 based index or member name to get the value.
    interface fson_value_get
        module procedure get_by_index
        module procedure get_by_name_chars
        module procedure get_by_name_string
    end interface fson_value_get

contains

    !
    ! FSON VALUE CREATE
    !
    function fson_value_create() result(new)
        type(fson_value), pointer :: new

        allocate(new)
    end function fson_value_create

    !
    ! FSON VALUE DESTROY
    !
    recursive subroutine fson_value_destroy(this)
        type(fson_value), pointer :: this

        if(allocated(this % value_real_array)) then
            deallocate(this % value_real_array)
        end if

        if(associated(this % children)) then
            call fson_value_destroy(this % children)
            nullify(this % children)
        end if

        if(associated(this % next)) then
            call fson_value_destroy(this % next)
            nullify (this % next)
        end if

        if(associated(this % name)) then
            call fson_string_destroy(this % name)
            nullify (this % name)
        end if

        if(associated(this % value_string)) then
            call fson_string_destroy(this % value_string)
            nullify (this % value_string)
        end if

        nullify(this)
    end subroutine fson_value_destroy

    !
    ! FSON VALUE ADD
    !
    ! Adds the memeber to the linked list
    subroutine fson_value_add(this, member)
        type(fson_value), pointer :: this, member, p

        ! associate the parent
        member % parent => this

        ! add to linked list
        if (associated(this % children)) then
            ! get to the tail of the linked list
            p => this % children
            do while (associated(p % next))
                p => p % next
            end do

            p % next => member
        else
            this % children => member
        end if
    end subroutine

    !
    ! FSON_VALUE_COUNT
    !
    integer function fson_value_count(this) result(count)
        type(fson_value), pointer :: this, p

        count = 0

        p => this % children

        do while (associated(p))
            count = count + 1
            p => p % next
        end do
    end function

    !
    ! GET BY INDEX
    !
    function get_by_index(this, index) result(p)
        type(fson_value), pointer :: this, p
        integer, intent(in) :: index
        integer :: i

        p => this % children

        do i = 1, index - 1
            p => p % next
        end do
    end function get_by_index

    !
    ! GET BY NAME CHARS
    !
    function get_by_name_chars(this, name) result(p)
        type(fson_value), pointer :: this, p
        character(len=*), intent(in) :: name

        type(fson_string), pointer :: string

        ! convert the char array into a string
        string => fson_string_create(name)

        p => get_by_name_string(this, string)
    end function get_by_name_chars

    !
    ! GET BY NAME STRING
    !
    function get_by_name_string(this, name) result(p)
        type(fson_value), pointer :: this, p
        type(fson_string), pointer :: name
        integer :: i

        if(this % value_type .ne. TYPE_OBJECT) then
            nullify(p)
            return
        end if

        do i=1, fson_value_count(this)
            p => fson_value_get(this, i)
            if (fson_string_equals(p%name, name)) then
                return
            end if
        end do

        ! didn't find anything
        nullify(p)
    end function get_by_name_string

    !
    ! FSON VALUE PRINT
    !
    recursive subroutine fson_value_print(iunit, this, is_compact_, is_bol_, indent_)
        integer, intent(in) :: iunit
        type(fson_value), pointer :: this, element
        logical, optional, intent(in) :: is_compact_, is_bol_
        integer, optional, intent(in) :: indent_
        character (len = 1024) :: tmp_chars
        integer :: indent, i, count, spaces, spaces1
        logical :: is_compact, is_bol

        if (present(is_compact_)) then
          is_compact = is_compact_
        else
          is_compact = .false.
        end if

        if (present(indent_) .and. .not. is_compact) then
            indent = indent_
        else
            indent = 0
        end if
        spaces = indent * 2
        if (is_compact) then
          spaces1 = 0
        else
          spaces1 = spaces + 2
        end if

        if (present(is_bol_)) then
            is_bol = is_bol_
        else
            is_bol = .true.
        end if

        if (is_bol) then
            write(iunit, '(A)', advance='no') repeat(" ", spaces)
        end if

        select case (this % value_type)
        case(TYPE_OBJECT)
            write(iunit, '(A)') "{"
            count = fson_value_count(this)
            do i = 1, count
                element => fson_value_get(this, i)
                call fson_string_copy(element % name, tmp_chars)
                write(iunit, '(4A)', advance='no') &
                     repeat(" ", spaces1), '"', trim(tmp_chars), '": '
                ! recursive print of the element
                call fson_value_print(iunit, element, is_compact, .false., indent + 1)
                if (i < count) then
                    write(iunit, '(A)') ","
                else
                    write(iunit, '()')
                end if
            end do
            write(iunit, '(2A)', advance='no') repeat(" ", spaces), "}"
        case (TYPE_ARRAY)
            write(iunit, '(A)') "["
            count = fson_value_count(this)
            do i = 1, count
                element => fson_value_get(this, i)
                ! recursive print of the element
                call fson_value_print(iunit, element, is_compact, .true., indent + 1)
                if (i < count) then
                    write(iunit, '(A)') ","
                else
                    write(iunit, '()')
                end if
            end do
            write(iunit, '(2A)', advance='no') repeat(" ", spaces), "]"
        case (TYPE_NULL)
            write(iunit, '(A)', advance='no') "null"
        case (TYPE_STRING)
            call fson_string_copy(this % value_string, tmp_chars)
            write(iunit, '(3A)', advance='no') '"', trim(tmp_chars), '"'
        case (TYPE_LOGICAL)
            if (this % value_logical) then
                write(iunit, '(A)', advance='no') "true"
            else
                write(iunit, '(A)', advance='no') "false"
            end if
        case (TYPE_INTEGER)
            write(iunit, '(I0)', advance='no') this % value_integer
        case (TYPE_REAL)
            write(iunit, '(E24.16e3)', advance='no') this % value_real
        case (TYPE_REAL_ARRAY)
            count = size(this % value_real_array)
            write(iunit, '(A)') "["
            do i = 1, count
                write(iunit, '(A, E24.16e3)', advance='no') &
                     repeat(" ", spaces1), &
                     this % value_real_array(i)
                if (i < count) then
                    write(iunit, '(A)') ","
                else
                    write(iunit, '()')
                end if
            end do
            write(iunit, '(2A)', advance='no') repeat(" ", spaces), "]"
        end select
    end subroutine fson_value_print
end module fson_value_m


!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.

!
! File:   fson_path_m.f95
! Author: Joseph A. Levin
!
! Created on March 10, 2012, 11:01 PM
!
! Modified on April 23, 2015, EigenKernel developer team

module fson_path_m
    use fson_value_m
    use fson_string_m

    implicit none

    private

    public :: fson_path_get

    interface fson_path_get
        module procedure get_by_path
        module procedure get_integer
        module procedure get_real
        module procedure get_double
        module procedure get_logical
        module procedure get_chars
        module procedure get_array
    end interface fson_path_get

contains

    !
    ! GET BY PATH
    !
    ! $     = root
    ! @     = this
    ! .     = child object member
    ! []    = child array element
    !
    recursive subroutine get_by_path(this, path, p)
        type(fson_value), pointer :: this, p
        character(len=*), intent(inout) :: path
        integer :: i, length, child_i
        character :: c
        logical :: array

        ! default to assuming relative to this
        p => this

        child_i = 1

        array = .false.

        length = len_trim(path)

        do i=1, length
            c = path(i:i)
            select case (c)
                case ("$")
                    ! root
                    do while (associated (p % parent))
                        p => p % parent
                    end do
                    child_i = i + 1
                case ("@")
                    ! this
                    p => this
                    child_i = i + 1
                case (".")
                    ! get child member from p
                    if (child_i < i) then
                        p => fson_value_get(p, path(child_i:i-1))
                    else
                        child_i = i + 1
                        cycle
                    end if

                    if(.not.associated(p)) then
                        return
                    end if

                    child_i = i+1
                case ("[")
                    ! start looking for the array element index
                    array = .true.
                    child_i = i + 1
                case ("]")
                    if (.not.array) then
                        print *, "ERROR: Unexpected ], not missing preceding ["
                        call exit(1)
                    end if
                    array = .false.
                    child_i = parse_integer(path(child_i:i-1))
                    p => fson_value_get(p, child_i)

                    child_i= i + 1
            end select
        end do

        ! grab the last child if present in the path
        if (child_i <= length) then
            p => fson_value_get(p, path(child_i:i-1))
            if(.not.associated(p)) then
                return
            else
            end if
        end if
    end subroutine get_by_path

    !
    ! PARSE INTEGER
    !
    integer function parse_integer(chars) result(integral)
        character(len=*) :: chars
        character :: c
        integer :: tmp, i

        integral = 0
        do i=1, len_trim(chars)
            c = chars(i:i)
            select case(c)
                case ("0":"9")
                    ! digit
                    read (c, '(i1)') tmp

                    ! shift
                    if(i > 1) then
                        integral = integral * 10
                    end if
                    ! add
                    integral = integral + tmp

                case default
                    return
            end select
        end do
    end function parse_integer

    !
    ! GET INTEGER
    !
    subroutine get_integer(this, path, value)
        type(fson_value), pointer :: this, p
        character(len=*), optional :: path
        integer :: value

        nullify(p)
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if

        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            call exit(1)
        end if

        if(p % value_type == TYPE_INTEGER) then
            value = p % value_integer
        else if (p % value_type == TYPE_REAL) then
            value = int(p % value_real)
        else if (p % value_type == TYPE_LOGICAL) then
            if (p % value_logical) then
                value = 1
            else
                value = 0
            end if
        else
            print *, "Unable to resolve value to integer: ", path
            call exit(1)
        end if
    end subroutine get_integer

    !
    ! GET REAL
    !
    subroutine get_real(this, path, value)
        type(fson_value), pointer :: this, p
        character(len=*), optional :: path
        real :: value

        nullify(p)

        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if

        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            call exit(1)
        end if

        if(p % value_type == TYPE_INTEGER) then
            value = p % value_integer
        else if (p % value_type == TYPE_REAL) then
            value = real(p % value_real)
        else if (p % value_type == TYPE_LOGICAL) then
            if (p % value_logical) then
                value = 1
            else
                value = 0
            end if
        else
            print *, "Unable to resolve value to real: ", path
            call exit(1)
        end if
    end subroutine get_real

    !
    ! GET DOUBLE
    !
    subroutine get_double(this, path, value)
        type(fson_value), pointer :: this, p
        character(len=*), optional :: path
        double precision :: value

        nullify(p)

        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if

        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            call exit(1)
        end if

        if(p % value_type == TYPE_INTEGER) then
            value = p % value_integer
        else if (p % value_type == TYPE_REAL) then
            value = p % value_real
        else if (p % value_type == TYPE_LOGICAL) then
            if (p % value_logical) then
                value = 1
            else
                value = 0
            end if
        else
            print *, "Unable to resolve value to double: ", path
            call exit(1)
        end if
    end subroutine get_double

    !
    ! GET LOGICAL
    !
    subroutine get_logical(this, path, value)
        type(fson_value), pointer :: this, p
        character(len=*), optional :: path
        logical :: value

        nullify(p)

        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if

        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            call exit(1)
        end if

        if(p % value_type == TYPE_INTEGER) then
            value = (p % value_integer > 0)
        else if (p % value_type == TYPE_LOGICAL) then
            value = p % value_logical
        else
            print *, "Unable to resolve value to real: ", path
            call exit(1)
        end if
    end subroutine get_logical

    !
    ! GET CHARS
    !
    subroutine get_chars(this, path, value)
        type(fson_value), pointer :: this, p
        character(len=*), optional :: path
        character(len=*) :: value

        nullify(p)

        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if

        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            call exit(1)
        end if

        if(p % value_type == TYPE_STRING) then
            call fson_string_copy(p % value_string, value)
        else
            print *, "Unable to resolve value to characters: ", path
            call exit(1)
        end if
    end subroutine get_chars

    !
    ! GET ARRAY
    !
    subroutine get_array(this, path, array_callback)
        type(fson_value), pointer :: this, p, element
        character(len=*), optional :: path
        integer :: index, count

        ! ELEMENT CALLBACK
        interface
            subroutine array_callback(element, index, count)
                use fson_value_m
                type(fson_value), pointer :: element
                integer :: index, count
            end subroutine array_callback
        end interface

        nullify(p)

        ! resolve the path to the value
        if(present(path)) then
            call get_by_path(this=this, path=path, p=p)
        else
            p => this
        end if

        if(.not.associated(p)) then
            print *, "Unable to resolve path: ", path
            call exit(1)
        end if

        if(p % value_type == TYPE_ARRAY) then
            count = fson_value_count(p)
            do index=1, count
                element => fson_value_get(p, index)
                call array_callback(element, index, count)
            end do
        else
            print *, "Resolved value is not an array. ", path
            call exit(1)
        end if
    end subroutine get_array
end module fson_path_m


! Copyright (c) 2012 Joseph A. Levin
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or
! substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
! LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
! DEALINGS IN THE SOFTWARE.


! FSON MODULE
!
! File:   fson.f95
! Author: Joseph A. Levin
!
! Created on March 6, 2012, 7:48 PM
!
! Modified on April 23, 2015, EigenKernel developer team

module fson
    use fson_value_m, fson_print => fson_value_print, fson_destroy => fson_value_destroy
    use fson_string_m
    use fson_path_m, fson_get => fson_path_get

    implicit none

    private

    public :: fson_parse, fson_value, fson_get, fson_print, fson_destroy
    public :: fson_set_as_string, fson_set_name, fson_set_as_real_array, fson_set_as_cmplx_array

    ! FILE IOSTAT CODES
    integer, parameter :: end_of_file = -1
    integer, parameter :: end_of_record = -2

    ! PARSING STATES
    integer, parameter :: STATE_LOOKING_FOR_VALUE = 1
    integer, parameter :: STATE_IN_OBJECT = 2
    integer, parameter :: STATE_IN_PAIR_NAME = 3
    integer, parameter :: STATE_IN_PAIR_VALUE = 4

    ! POP/PUSH CHARACTER
    integer :: pushed_index = 0
    character (len = 10) :: pushed_char

contains

    !
    ! FSON PARSE
    !
    function fson_parse(file, unit) result(p)
        type(fson_value), pointer :: p
        integer, optional, intent(inout) :: unit
        character(len = *), intent(in) :: file
        logical :: unit_available
        integer :: u
        ! init the pointer to null
        nullify(p)

        ! select the file unit to use
        if (present(unit)) then
            u = unit
        else
            ! find the first available unit
            unit_available = .false.
            u = 20

            do while (.not.unit_available)
                inquire(unit = u, exist = unit_available)
                u = u + 1
            end do
        end if

        ! open the file
        open (unit = u, file = file, status = "old", action = "read", form = "formatted", position = "rewind")

        ! create the value and associate the pointer
        p => fson_value_create()

        ! parse as a value
        call parse_value(unit = u, value = p)

        ! close the file
        if( .not. present(unit)) then
            close (u)
        end if
    end function fson_parse

    !
    ! PARSE_VALUE
    !
    recursive subroutine parse_value(unit, value)
        integer, intent(inout) :: unit
        type(fson_value), pointer :: value
        logical :: eof
        character :: c

        ! for some unknown reason the next pointer is getting messed with the pop
        type(fson_value), pointer :: hack

        ! start the hack
        hack => value % next

        ! pop the next non whitespace character off the file
        c = pop_char(unit, eof = eof, skip_ws = .true.)

        ! finish the hack; set the next pointer to whatever it was before the pop
        value % next => hack

        if (eof) then
            return
        else
            select case (c)
            case ("{")
                ! start object
                value % value_type = TYPE_OBJECT
                call parse_object(unit, value)
            case ("[")
                ! start array
                value % value_type = TYPE_ARRAY
                call parse_array(unit, value)
            case ("]")
                ! end an empty array
                nullify(value)
            case ('"')
                ! string
                value % value_type = TYPE_STRING
                value % value_string => parse_string(unit)
            case ("t")
                !true
                value % value_type = TYPE_LOGICAL
                call parse_for_chars(unit, "rue")
                value % value_logical = .true.
            case ("f")
                !false
                value % value_type = TYPE_LOGICAL
                value % value_logical = .false.
                call parse_for_chars(unit, "alse")
            case ("n")
                value % value_type = TYPE_NULL
                call parse_for_chars(unit, "ull")
            case("-", "0": "9")
                call push_char(c)
                call parse_number(unit, value)
            case default
                print *, "ERROR: Unexpected character while parsing value. '", c, "' ASCII=", iachar(c)
                call exit (1)
            end select
        end if
    end subroutine parse_value

    !
    ! PARSE OBJECT
    !
    recursive subroutine parse_object(unit, parent)
        integer, intent(inout) :: unit
        type(fson_value), pointer :: parent, pair

        logical :: eof
        character :: c

        ! pair name
        c = pop_char(unit, eof = eof, skip_ws = .true.)
        if (eof) then
            print *, "ERROR: Unexpected end of file while parsing start of object."
            call exit (1)
        else if ("}" == c) then
            ! end of an empty object
            return
        else if ('"' == c) then
            pair => fson_value_create()
            pair % name => parse_string(unit)
        else
            print *, "ERROR: Expecting string: '", c, "'"
            call exit (1)
        end if

        ! pair value
        c = pop_char(unit, eof = eof, skip_ws = .true.)
        if (eof) then
            print *, "ERROR: Unexpected end of file while parsing object member. 1"
            call exit (1)
        else if (":" == c) then
            ! parse the value
            call parse_value(unit, pair)
            call fson_value_add(parent, pair)
        else
            print *, "ERROR: Expecting : and then a value. ", c
            call exit (1)
        end if

        ! another possible pair
        c = pop_char(unit, eof = eof, skip_ws = .true.)
        if (eof) then
            return
        else if ("," == c) then
            ! read the next member
            call parse_object(unit = unit, parent = parent)
        else if ("}" == c) then
            return
        else
            print *, "ERROR: Expecting end of object.", c
            call exit (1)
        end if
    end subroutine parse_object

    !
    ! PARSE ARRAY
    !
    recursive subroutine parse_array(unit, array)
        integer, intent(inout) :: unit
        type(fson_value), pointer :: array, element

        logical :: eof
        character :: c

        ! try to parse an element value
        element => fson_value_create()
        call parse_value(unit, element)

        ! parse value will disassociate an empty array value
        if (associated(element)) then
            call fson_value_add(array, element)
        end if

        ! popped the next character
        c = pop_char(unit, eof = eof, skip_ws = .true.)

        if (eof) then
            return
        else if ("," == c) then
            ! parse the next element
            call parse_array(unit, array)
        else if ("]" == c) then
            ! end of array
            return
        end if
    end subroutine parse_array

    !
    ! PARSE STRING
    !
    function parse_string(unit) result(string)
        integer, intent(inout) :: unit
        type(fson_string), pointer :: string

        logical :: eof
        character :: c, last

        string => fson_string_create()

        do
            c = pop_char(unit, eof = eof, skip_ws = .false.)
            if (eof) then
                print *, "Expecting end of string"
                call exit(1)!
            else if ('"' == c .and. last .ne. '\\') then
                exit
            else
                last = c
                call fson_string_append(string, c)
            end if
        end do
    end function parse_string

    !
    ! PARSE FOR CHARACTERS
    !
    subroutine parse_for_chars(unit, chars)
        integer, intent(in) :: unit
        character(len = *), intent(in) :: chars
        integer :: i, length
        logical :: eof
        character :: c

        length = len_trim(chars)

        do i = 1, length
            c = pop_char(unit, eof = eof, skip_ws = .true.)
            if (eof) then
                print *, "ERROR: Unexpected end of file while parsing array."
                call exit (1)
            else if (c .ne. chars(i:i)) then
                print *, "ERROR: Unexpected character.'", c,"'", chars(i:i)
                call exit (1)
            end if
        end do
    end subroutine parse_for_chars

    !
    ! PARSE NUMBER
    !
    subroutine parse_number(unit, value)
        integer, intent(inout) :: unit
        type(fson_value), pointer :: value
        logical :: eof  !, negative, scientific, decimal
        character :: c
        !integer :: integral, exp, digit_count
        !double precision :: frac

        integer, parameter :: max_num_length = 256
        character(len=max_num_length) :: buf
        !logical :: is_real
        integer :: i  !, n
        do i = 1, max_num_length
          c = pop_char(unit, eof = eof, skip_ws = .true.)
          if (eof) then
            print *, "ERROR: Unexpected end of file while parsing number."
            call exit(1)
          end if
          if (c == ',' .or. c == ']' .or. c == '}') then
            call push_char(c)
            exit
          end if
          buf(i:i) = c
        end do

        if (i == max_num_length) then
          print *, "ERROR: too long number found, ", buf
          call exit(1)
        else if (index(buf(1 : i - 1), '.') > 0) then
          value%value_type = TYPE_REAL
          read(buf(1 : i - 1), *) value%value_real
        else
          value%value_type = TYPE_INTEGER
          read(buf(1 : i - 1), *) value%value_integer
        end if
    end subroutine

    !
    ! POP CHAR
    !
    recursive character function pop_char(unit, eof, skip_ws) result(popped)
        integer, intent(in) :: unit
        logical, intent(out) :: eof
        logical, intent(in), optional :: skip_ws

        integer :: ios
        character :: c
        logical :: ignore

        eof = .false.
        if (.not.present(skip_ws)) then
            ignore = .false.
        else
            ignore = skip_ws
        end if

        do
            if (pushed_index > 0) then
                ! there is a character pushed back on, most likely from the number parsing
                c = pushed_char(pushed_index:pushed_index)
                pushed_index = pushed_index - 1
            else
                read (unit = unit, fmt = "(a)", advance = "no", iostat = ios) c
            end if
            if (ios == end_of_record) then
                cycle
            else if (ios == end_of_file) then
                eof = .true.
                exit
            else if (iachar(c) <= 32) then
                ! non printing ascii characters
                cycle
            else if (ignore .and. c == " ") then
                cycle
            else
                popped = c
                exit
            end if
        end do
    end function pop_char

    !
    ! PUSH CHAR
    !
    subroutine push_char(c)
        character, intent(inout) :: c
        pushed_index = pushed_index + 1
        pushed_char(pushed_index:pushed_index) = c

    end subroutine push_char

    subroutine fson_set_as_string(str, val)
        character(len=*), intent(in) :: str
        type(fson_value), pointer, intent(inout) :: val

        type(fson_string), pointer :: str_fson

        val%value_type = TYPE_STRING
        str_fson => fson_string_create()
        call fson_string_append(str_fson, str)
        val%value_string => str_fson
    end subroutine fson_set_as_string

    subroutine fson_set_name(name, val)
        character(len=*), intent(in) :: name
        type(fson_value), pointer, intent(inout) :: val

        type(fson_string), pointer :: name_fson

        name_fson => fson_string_create()
        call fson_string_append(name_fson, name)
        val%name => name_fson
    end subroutine fson_set_name

    subroutine fson_set_as_real_array(length, xs, val)
        integer, intent(in) :: length
        double precision, intent(in) :: xs(:)
        type(fson_value), pointer, intent(inout) :: val

        val%value_type = TYPE_REAL_ARRAY

        if (allocated(val%value_real_array)) then
          deallocate(val%value_real_array)
        end if

        allocate(val%value_real_array(length))
        val%value_real_array(1:length) = xs(1:length)
    end subroutine fson_set_as_real_array


    subroutine fson_set_as_cmplx_array(length, xs, val)
        integer, intent(in) :: length
        complex(kind(0d0)), intent(in) :: xs(:)
        type(fson_value), pointer, intent(inout) :: val

        type(fson_value), pointer :: real_part, imag_part

        val%value_type = TYPE_OBJECT
        allocate(real_part)
        allocate(imag_part)
        call fson_set_name('real', real_part)
        call fson_set_name('imag', imag_part)
        call fson_set_as_real_array(length, dble(xs(1:length)), real_part)
        call fson_set_as_real_array(length, dimag(xs(1:length)), imag_part)
        call fson_value_add(val, real_part)
        call fson_value_add(val, imag_part)
    end subroutine fson_set_as_cmplx_array
end module fson
