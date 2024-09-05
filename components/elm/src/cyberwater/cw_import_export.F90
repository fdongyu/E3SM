module cw_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use abortutils   , only: endrun
  use spmdMod      , only : mpicom, masterproc
  use decompmod    , only : bounds_type, ldecomp
  use elm_cpl_indices
  use mpi
  use cw_cpl_indices
  use high_level_api

  implicit none

  !--------------------------------------------------------------------
  !
  ! Private module variables
  !
  !--------------------------------------------------------------------

  ! MPI variables
  integer :: nCellsGlobal
  integer, dimension(:), allocatable :: nCellsDisplacement, indexToCellIDGathered
  integer, dimension(:), allocatable :: nCellsPerProc

contains


  !-----------------------------------------------------------------------
  subroutine cw_import_mct( bounds, x2l)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model
    use elm_varctl     , only: iulog
    use domainMod      , only : ldomain
    use shr_kind_mod   , only: r8 => shr_kind_r8

    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    !
    ! !LOCAL VARIABLES:
    integer  :: g, i, n
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code

    ! MPI variables
    integer :: iProc, nProcs
    real(r8),allocatable :: glat(:), glon(:)  ! global
    real(r8),allocatable :: Sa_z(:)  ! bottom atm level height    m
    real(r8),allocatable :: Sa_u(:)  ! bottom atm level zon wind  m/s
    real(r8),allocatable :: Sa_v(:)  ! bottom atm level mer wind  m/s

    real(r8),allocatable :: lat_recv(:), lon_recv(:)  ! local array for testing scatter
    real(r8),allocatable :: arr_receive_lat(:) !, arr_receive_lon, arr_receive_Sa_z
    !real(c_double), dimension(:), allocatable :: arr_receive_lat, arr_receive_lon, arr_receive_Sa_z
    integer :: send_status, var_receive_size, receive_status

    character(len=*), parameter :: sub = 'cw_import_mct'

    if (masterproc) then
      write(6,*) 'Importing variables from drv'
    end if


    ! local grid order 
    allocate(gindex(bounds%begg:bounds%endg),stat=ier) 
    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = bounds%endg - bounds%begg + 1        ! local grid size
    gsize = ldomain%ni * ldomain%nj              ! global grid size

    call MPI_COMM_SIZE( mpicom, nProcs, ier)     ! need to update mpicom to mpicom_lnd


    if (masterproc) then
       allocate(nCellsPerProc(nProcs))
       allocate(nCellsDisplacement(nProcs))
    end if

    ! Gather local grid size
    call MPI_GATHER( lsize, 1, MPI_INTEGER, nCellsPerProc, 1, MPI_INTEGER, &
                     0, mpicom, ier)

    ! Set Displacement variable for GATHERV command
    if (masterproc) then
       nCellsGlobal = sum(nCellsPerProc)
       if (nCellsGlobal .NE. gsize) then
          call endrun('global size not consistent!')
       end if
       allocate(indexToCellIDGathered(nCellsGlobal))
       nCellsDisplacement(1) = 0
       if (nProcs > 1) then
          do iProc=2,nProcs
              nCellsDisplacement(iProc) = nCellsDisplacement(iProc-1) + nCellsPerProc(iProc-1)
          enddo
       endif
    endif
         
    ! Gather indexToCellID
    call MPI_GATHERV( gindex, lsize, MPI_INTEGER, indexToCellIDGathered, &
                    nCellsPerProc, nCellsDisplacement, MPI_INTEGER, 0, mpicom, ier) 
 
    call MPI_Barrier(mpicom, ier)

    if (masterproc) then
       write (6,*) nCellsPerProc
    end if


    ! Allocate global array only on master processor
    if (masterproc) then
       allocate(glat(gsize))
       allocate(glon(gsize))
       allocate(Sa_z(gsize))
       allocate(Sa_u(gsize))
       allocate(Sa_v(gsize))
    end if 


    ! Gather lat and lon first
    call MPI_GATHERV(ldomain%latc, lsize, MPI_DOUBLE, glat, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(ldomain%lonc, lsize, MPI_DOUBLE, glon, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) 

    ! Gather variables before sending to CW
    call MPI_GATHERV(x2l(index_x2l_Sa_z,:), lsize, MPI_DOUBLE, Sa_z, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Sa_u,:), lsize, MPI_DOUBLE, Sa_u, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Sa_v,:), lsize, MPI_DOUBLE, Sa_v, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)


    call MPI_Barrier(mpicom, ier)


    ! send data to server 
    if (masterproc) then

      write(iulog,*) "Send data to server !"

      send_status = send_data_with_retries(index_a2l_latitude, glat, 25, 10)
      print *, "------ Sleeping for 10 seconds ------"
      call sleep(5)

      send_status = send_data_with_retries(index_a2l_longitude, glon, 25, 10)
      print *, "------ Sleeping for 10 seconds ------"
      call sleep(5)

      send_status = send_data_with_retries(index_Sa_z, Sa_z, 25, 10)
      print *, "------ Sleeping for 10 seconds ------"
      call sleep(5)
    

      write(iulog,*) "Receive data from server !"


      ! receive varid=4 from CyberWater
      if (check_data_availability_with_retries(4, 25, 10) == 1) then
        allocate(arr_receive_lat(gsize))
        receive_status = receive_data_with_retries(4, arr_receive_lat, 25 ,10)

        print *, "Received data:"
        do i = 1, gsize
            print *, "arr_receive(", i, ") = ", arr_receive_lat(i)
        end do

        print *, "------ Sleeping for 10 seconds ------"
        call sleep(10)
      else
        ! Handle the case where data is not available
        print *, "Data is not available after retries."
        call endrun('Data is not available after retries!')
      endif



      ! Clean up allocated resources
      deallocate(arr_receive_lat)
!      deallocate(arr_receive_lon)
!      deallocate(arr_receive_Sa_z)

    end if
    call MPI_Barrier(mpicom, ier)


    !! test code for Scatter back after receiving from CyberWater
!    allocate(lat_recv(lsize))
!    allocate(lon_recv(lsize))

!    call MPI_SCATTERV(glon, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    lon_recv, lsize, MPI_DOUBLE, 0, mpicom, ier)

!    call MPI_Barrier(mpicom, ier)

!    write(6,*) "Checking MPI scatter!"
!    write(6,*) "ldomain%lonc", ldomain%lonc
!    write(6,*) "lon_recv", lon_recv

!    deallocate(lat_recv)
!    deallocate(lon_recv)



    ! free memory
    deallocate(gindex)


    if (masterproc) then
       deallocate(glat)
       deallocate(glon)
       deallocate(Sa_z)
       deallocate(Sa_u)
       deallocate(Sa_v)
       deallocate(nCellsPerProc)
       deallocate(nCellsDisplacement)
       deallocate(indexToCellIDGathered)
    end if



  end subroutine cw_import_mct

  !===============================================================================

!  subroutine cw_export_mct( bounds, l2x)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the Cyberwater model to the coupler 
    ! 
    ! !USES:
!    use elm_varctl     , only : iulog
!    use domainMod      , only : ldomain
    
    ! !ARGUMENTS:
!    implicit none
!    type(bounds_type) , intent(in)    :: bounds  ! bounds
!    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    !
    ! !LOCAL VARIABLES:
!    integer  :: g,i   ! indices
!    integer  :: dtime ! time step
!    real(r8),allocatable :: lat_recv(:), lon_recv(:)  ! local


    ! receive data from CyberWater
!    allocate(lon_recv(lsize))

    ! glon is the 1D longitude received from CyberWater 
!    call MPI_SCATTERV(glon, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    lon_recv, lsize, MPI_DOUBLE, 0, mpicom, ier)

!    call MPI_Barrier(mpicom, ier)
       
!    l2x(:,:) = 0.0_r8

!    do g = bounds%begg,bounds%endg
!       i = 1 + (g-bounds%begg)
!       l2x(index_l2x_Sl_t,i)        =  lnd2atm_vars%t_rad_grc(g)
!       l2x(index_l2x_Sl_snowh,i)    =  lnd2atm_vars%h2osno_grc(g)
!       l2x(index_l2x_Sl_avsdr,i)    =  lnd2atm_vars%albd_grc(g,1)
!       l2x(index_l2x_Sl_anidr,i)    =  lnd2atm_vars%albd_grc(g,2)
!       l2x(index_l2x_Sl_avsdf,i)    =  lnd2atm_vars%albi_grc(g,1)
!       l2x(index_l2x_Sl_anidf,i)    =  lnd2atm_vars%albi_grc(g,2)
!       l2x(index_l2x_Sl_tref,i)     =  lnd2atm_vars%t_ref2m_grc(g)
!       l2x(index_l2x_Sl_qref,i)     =  lnd2atm_vars%q_ref2m_grc(g)
!       l2x(index_l2x_Sl_u10,i)      =  lnd2atm_vars%u_ref10m_grc(g)
!       l2x(index_l2x_Fall_taux,i)   = -lnd2atm_vars%taux_grc(g)
!       l2x(index_l2x_Fall_tauy,i)   = -lnd2atm_vars%tauy_grc(g)
!       l2x(index_l2x_Fall_lat,i)    = -lnd2atm_vars%eflx_lh_tot_grc(g)
!       l2x(index_l2x_Fall_sen,i)    = -lnd2atm_vars%eflx_sh_tot_grc(g)
!       l2x(index_l2x_Fall_lwup,i)   = -lnd2atm_vars%eflx_lwrad_out_grc(g)
!       l2x(index_l2x_Fall_evap,i)   = -lnd2atm_vars%qflx_evap_tot_grc(g)
!       l2x(index_l2x_Fall_swnet,i)  =  lnd2atm_vars%fsa_grc(g)
!    end do

!    deallocate(lon_recv)
    
!  end subroutine cw_export_mct


  function get_value() result(value)
     ! Implement the logic to get the value
     ! Return the value
     integer :: value
     ! Example: Return a random value between 0 and 20
     value = randint(0, 20) ! randint is a placeholder, replace it with your actual function
  end function get_value

  ! Placeholder for the randint function
  function randint(lower, upper) result(value)
    integer, intent(in) :: lower, upper
    integer :: value
    real :: harvest  ! Declare harvest as real
    call random_seed()
    call random_number(harvest)
    value = lower + int(real(upper - lower + 1) * harvest)
  end function randint


end module cw_import_export
