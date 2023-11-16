module c_interface_combined
    use, intrinsic :: iso_c_binding

    implicit none

    interface
        function send_data_to_server(arr, n) bind(C, name="send_data_to_server")
            import :: c_int, c_float
            real(c_float), intent(in) :: arr(*)
            integer(c_int), value, intent(in) :: n
            integer(c_int) :: send_data_to_server
        end function send_data_to_server

        function fetch_data_from_server(arr, n) bind(C, name="fetch_data_from_server")
            import :: c_int, c_float
            real(c_float), intent(out) :: arr(*)
            integer(c_int), value, intent(in) :: n
            integer(c_int) :: fetch_data_from_server
        end function fetch_data_from_server
    end interface

end module c_interface_combined

! The main program starts here
program combined_client
    use c_interface_combined
    implicit none

    integer :: status_send, status_fetch
    real, dimension(5) :: arr_send = [1.2, 3.4, 5.6, 7.8, 9.10]
    real, dimension(5) :: arr_fetch

    ! Send float array to server
    status_send = send_data_to_server(arr_send, size(arr_send))
    if (status_send /= 0) then
        print *, "Failed to send data to server"
    else
        print *, "Data sent successfully!"
    end if

    ! Fetch float array from server
    status_fetch = fetch_data_from_server(arr_fetch, size(arr_fetch))
    if (status_fetch /= 0) then
        print *, "Failed to fetch data from server"
    else
        print *, "Data received:", arr_fetch
    end if

end program combined_client

