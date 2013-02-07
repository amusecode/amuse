module FortranSocketsInterface
    integer HEADER_FLAGS, HEADER_CALL_ID, HEADER_FUNCTION_ID, HEADER_CALL_COUNT, &
        HEADER_INTEGER_COUNT, HEADER_LONG_COUNT, HEADER_FLOAT_COUNT, &
        HEADER_DOUBLE_COUNT, HEADER_BOOLEAN_COUNT, HEADER_STRING_COUNT, &
        HEADER_SIZE

    parameter (HEADER_FLAGS=1, HEADER_CALL_ID=2, HEADER_FUNCTION_ID=3, &
        HEADER_CALL_COUNT=4, HEADER_INTEGER_COUNT=5, HEADER_LONG_COUNT=6, &
        HEADER_FLOAT_COUNT=7, HEADER_DOUBLE_COUNT=8, &
        HEADER_BOOLEAN_COUNT=9, HEADER_STRING_COUNT=10, &
        HEADER_SIZE=10)

    interface
        subroutine receive_integers &
            (ints, length) &
            bind(c, name='forsockets_receive_integers')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: ints
            integer (c_int32_t), value :: length
        end subroutine receive_integers

        subroutine receive_longs &
            (longs, length) &
            bind(c, name='forsockets_receive_longs')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: longs
            integer (c_int32_t), value :: length
        end subroutine receive_longs

        subroutine receive_floats &
            (floats, length) &
            bind(c, name='forsockets_receive_floats')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: floats
            integer (c_int32_t), value :: length
        end subroutine receive_floats

        subroutine receive_doubles &
            (doubles, length) &
            bind(c, name='forsockets_receive_doubles')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: doubles
            integer (c_int32_t), value :: length
        end subroutine receive_doubles

        subroutine receive_booleans &
            (booleans, length) &
            bind(c, name='forsockets_receive_booleans')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: booleans
            integer (c_int32_t), value :: length
        end subroutine receive_booleans

        subroutine receive_string &
            (string, length) &
            bind(c, name='forsockets_receive_string')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: string
            integer (c_int32_t), value :: length
        end subroutine receive_string

        subroutine send_integers &
            (ints, length) &
            bind(c, name='forsockets_send_integers')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: ints
            integer (c_int32_t), value :: length
        end subroutine send_integers

        subroutine send_longs &
            (longs, length) &
            bind(c, name='forsockets_send_longs')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: longs
            integer (c_int32_t), value :: length
        end subroutine send_longs

        subroutine send_floats &
            (floats, length) &
            bind(c, name='forsockets_send_floats')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: floats
            integer (c_int32_t), value :: length
        end subroutine send_floats

        subroutine send_doubles &
            (doubles, length) &
            bind(c, name='forsockets_send_doubles')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: doubles
            integer (c_int32_t), value :: length
        end subroutine send_doubles

        subroutine send_booleans &
            (booleans, length) &
            bind(c, name='forsockets_send_booleans')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: booleans
            integer (c_int32_t), value :: length
        end subroutine send_booleans

        subroutine send_string &
            (string, length) &
            bind(c, name='forsockets_send_string')
            use iso_c_binding
            implicit none
            type (c_ptr), value :: string
            integer (c_int32_t), value :: length
        end subroutine send_string

        subroutine forsockets_init &
            (port) &
            bind(c, name='forsockets_init')
            use iso_c_binding
            implicit none
            integer (c_int32_t), value :: port
        end subroutine forsockets_init

        subroutine forsockets_close &
            () &
            bind(c, name='forsockets_close')
            use iso_c_binding
            implicit none
        end subroutine forsockets_close

    end interface
end module FortranSocketsInterface
