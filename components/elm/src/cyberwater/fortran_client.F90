module c_interface_combined
    use, intrinsic :: iso_c_binding

    implicit none

    interface
        function send_data_to_server(arr, n) bind(C, name="send_data_to_server")
            import :: c_int, c_double
            real(c_double), intent(in) :: arr(*)
            integer(c_int), value, intent(in) :: n
            integer(c_int) :: send_data_to_server
        end function send_data_to_server

        function fetch_data_from_server(arr, n) bind(C, name="fetch_data_from_server")
            import :: c_int, c_double
            real(c_double), intent(out) :: arr(*)
            integer(c_int), value, intent(in) :: n
            integer(c_int) :: fetch_data_from_server
        end function fetch_data_from_server
    end interface

end module c_interface_combined
