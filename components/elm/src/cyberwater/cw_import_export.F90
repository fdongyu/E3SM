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

  ! Sending Params
  integer :: send_retry_time  = 10
  integer  :: send_retry_count = 25
  real(r8) :: recv_delay       = 5.0_r8

  ! Receiving Params
  integer :: recv_retry_time  = 10
  integer  :: recv_retry_count = 25
  logical  :: use_delay = .true.


contains


  !-----------------------------------------------------------------------
  subroutine cw_import_mct(bounds, x2l, id)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model
    use elm_varctl     , only: iulog
    use domainMod      , only : ldomain
    use shr_kind_mod   , only: r8 => shr_kind_r8

    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    type(SessionID)    , intent(in)    :: id
    !
    ! !LOCAL VARIABLES:
    integer  :: g, i, j, n
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code

    ! MPI variables
    integer :: iProc, nProcs
    ! Global 
    real(r8),allocatable :: glat(:), glon(:)  ! latitude, longitude
    real(r8),allocatable :: Sa_z(:)  ! bottom atm level height    m
    real(r8),allocatable :: Sa_vel(:)  ! bottom atm level zon wind velocity sqrt(Sa_u*Sa_u+Sa_v*Sa_v)  m/s
    real(r8),allocatable :: Sa_u(:)
    real(r8),allocatable :: Sa_v(:)
    real(r8),allocatable :: Sa_shum(:) ! bottom atm level spec hum Pa
    real(r8),allocatable :: Sa_pbot(:) ! bottom atm level pressue  Pa
    real(r8),allocatable :: Sa_tbot(:) ! bottom atm level temp     degree
    real(r8),allocatable :: Faxa_lwdn(:) ! downward longwave heat flux W m-2
    real(r8),allocatable :: Faxa_precip(:) ! total precipitation rainc+rainl+snowc+snowl mm/s
    real(r8),allocatable :: Faxa_rainc(:) ! convective precipitation rate mm/s
    real(r8),allocatable :: Faxa_rainl(:) ! large-scale precipitation rate mm/s
    real(r8),allocatable :: Faxa_snowc(:) ! convective snow rate (water equivalent) mm/s
    real(r8),allocatable :: Faxa_snowl(:) ! large-scale snow rate (water equivalent) mm/s
    real(r8),allocatable :: Faxa_sw(:)    ! total solar short wave radiation W m-2
    real(r8),allocatable :: Faxa_swndr(:) ! direct near-infrared incident solar radiation W m-2 
    real(r8),allocatable :: Faxa_swvdr(:) ! direct visible indicent solar radiation W m-2
    real(r8),allocatable :: Faxa_swndf(:) ! diffuse near-infrared incident solar radiation W m-2 
    real(r8),allocatable :: Faxa_swvdf(:) ! diffuse visible incident solar radiation W m-2
    real(r8),allocatable :: Sa_co2prog(:) ! prognostic CO2 at the lowest model level 1e-6 mol/mol
    !real(r8),allocatable :: Sa_co2diag(:) ! diagnostic CO2 at the lowest model level 1e-6 mol/mol

    real(r8) :: forc_rainc           ! rainxy Atm flux mm/s
    real(r8) :: forc_rainl           ! rainxy Atm flux mm/s
    real(r8) :: forc_snowc           ! snowfxy Atm flux  mm/s
    real(r8) :: forc_snowl           ! snowfxl Atm flux  mm/s
    real(r8) :: swndf, swndr, swvdf, swvdr
    real(r8) :: ubot, vbot

     
    integer :: send_status

    real(r8), save :: last_exit_time = 0.0  ! Time when subroutine last exited
    real(r8) :: entry_time, outside_time
    logical, save :: first_call = .true.  ! To handle the first call scenario


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
       allocate(Sa_vel(gsize))
       allocate(Sa_u(gsize))
       allocate(Sa_v(gsize))
       allocate(Sa_shum(gsize))
       allocate(Sa_pbot(gsize))
       allocate(Sa_tbot(gsize))
       allocate(Faxa_lwdn(gsize))
       allocate(Faxa_precip(gsize))
       allocate(Faxa_rainc(gsize))
       allocate(Faxa_rainl(gsize))
       allocate(Faxa_snowc(gsize))
       allocate(Faxa_snowl(gsize))
       allocate(Faxa_sw(gsize))
       allocate(Faxa_swndr(gsize))
       allocate(Faxa_swvdr(gsize))
       allocate(Faxa_swndf(gsize))
       allocate(Faxa_swvdf(gsize))
       allocate(Sa_co2prog(gsize))
       !allocate(Sa_co2diag(gsize))
    end if 


    ! Gather lat and lon first
    call MPI_GATHERV(ldomain%latc, lsize, MPI_DOUBLE, glat, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(ldomain%lonc, lsize, MPI_DOUBLE, glon, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) 

    ! Gather variables before sending to CW
    call MPI_GATHERV(x2l(index_x2l_Sa_z,:), lsize, MPI_DOUBLE, Sa_z, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) ! Atm state m
    call MPI_GATHERV(x2l(index_x2l_Sa_u,:), lsize, MPI_DOUBLE, Sa_u, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) ! Atm state m/s
    call MPI_GATHERV(x2l(index_x2l_Sa_v,:), lsize, MPI_DOUBLE, Sa_v, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) ! Atm state m/s
    call MPI_GATHERV(x2l(index_x2l_Sa_shum,:), lsize, MPI_DOUBLE, Sa_shum, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Sa_pbot,:), lsize, MPI_DOUBLE, Sa_pbot, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Sa_tbot,:), lsize, MPI_DOUBLE, Sa_tbot, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_lwdn,:), lsize, MPI_DOUBLE, Faxa_lwdn, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_rainc,:), lsize, MPI_DOUBLE, Faxa_rainc, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_rainl,:), lsize, MPI_DOUBLE, Faxa_rainl, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_snowc,:), lsize, MPI_DOUBLE, Faxa_snowc, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_snowl,:), lsize, MPI_DOUBLE, Faxa_snowl, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swndr,:), lsize, MPI_DOUBLE, Faxa_swndr, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swvdr,:), lsize, MPI_DOUBLE, Faxa_swvdr, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swndf,:), lsize, MPI_DOUBLE, Faxa_swndf, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swvdf,:), lsize, MPI_DOUBLE, Faxa_swvdf, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Sa_co2prog,:), lsize, MPI_DOUBLE, Sa_co2prog, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    !call MPI_GATHERV(x2l(index_x2l_Sa_co2diag,:), lsize, MPI_DOUBLE, Sa_co2diag, nCellsPerProc, &
    !                    nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)

    call MPI_Barrier(mpicom, ier)



    ! send data to server 
    if (masterproc) then

      ! Prepare data 
      ! This calculation depends on the hydrological model in CyberWater
      do i = 0, gsize
        ubot       = Sa_u(i)         ! m/s
        vbot       = Sa_v(i)         ! m/s
        forc_rainc = Faxa_rainc(i)   ! mm/s
        forc_rainl = Faxa_rainl(i)   ! mm/s
        forc_snowc = Faxa_snowc(i)   ! mm/s
        forc_snowl = Faxa_snowl(i)   ! mm/s
        swndr      = Faxa_swndr(i)   ! forc_solldxy Atm flux  W/m^2
        swvdr      = Faxa_swvdr(i)   ! forc_solsxy  Atm flux  W/m^2
        swndf      = Faxa_swndf(i)   ! forc_solldxy Atm flux  W/m^2
        swvdf      = Faxa_swvdf(i)   ! forc_solsdxy Atm flux  W/m^2

        Sa_vel(i)       = sqrt(ubot*ubot + vbot*vbot)
        Faxa_precip(i)  = forc_rainc + forc_rainl + forc_snowc + forc_snowl
        Faxa_sw(i)      = swndr + swvdr + swndf + swvdf ! solar radiation
      end do


      call cpu_time(entry_time)

      ! If not the first call, calculate the time spent outside this subroutine
      if (.not. first_call) then
        outside_time = entry_time - last_exit_time
        print *, 'Time spent since last execution: ', outside_time, ' seconds'
      else
        first_call = .false.
      endif

      write(iulog,*) "Send data to server !"

      send_status = send_data_with_retries(index_a2l_latitude, glat, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_a2l_longitude, glon, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_z, Sa_z, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_vel, Sa_vel, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_shum, Sa_shum, send_retry_count, send_retry_time, recv_delay)
      
      send_status = send_data_with_retries(index_Sa_pbot, Sa_pbot, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_tbot, Sa_tbot, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Faxa_lwdn, Faxa_lwdn, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Faxa_precip, Faxa_precip, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Faxa_sw, Faxa_sw, send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_co2prog, Sa_co2prog, send_retry_count, send_retry_time, recv_delay)
    
      ! Record the exit time
      call cpu_time(last_exit_time)

    end if
    call MPI_Barrier(mpicom, ier)


    ! free memory
    deallocate(gindex)


    if (masterproc) then
      deallocate(glat)
      deallocate(glon)
      deallocate(Sa_z)
      deallocate(Sa_vel)
      deallocate(Sa_u)
      deallocate(Sa_v)
      deallocate(Sa_shum)
      deallocate(Sa_pbot)
      deallocate(Sa_tbot)
      deallocate(Faxa_lwdn)
      deallocate(Faxa_precip)
      deallocate(Faxa_rainc)
      deallocate(Faxa_rainl)
      deallocate(Faxa_snowc)
      deallocate(Faxa_snowl)
      deallocate(Faxa_sw)
      deallocate(Faxa_swndr)
      deallocate(Faxa_swvdr)
      deallocate(Faxa_swndf)
      deallocate(Faxa_swvdf)
      deallocate(Sa_co2prog)
      deallocate(nCellsPerProc)
      deallocate(nCellsDisplacement)
      deallocate(indexToCellIDGathered)
    end if



  end subroutine cw_import_mct

  !===============================================================================

  subroutine cw_export_mct(bounds, l2x, id)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the Cyberwater model to the coupler 
    ! 
    ! !USES:
    use elm_varctl     , only : iulog
    use domainMod      , only : ldomain
    use shr_kind_mod   , only : r8 => shr_kind_r8
    
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds  ! bounds
    real(r8)          , intent(out)   :: l2x(:,:)! land to coupler export state on land grid
    type(SessionID)   , intent(in)    :: id
    
    !
    ! !LOCAL VARIABLES:
    integer  :: g, i, j, n
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code

    ! MPI variables
    integer :: iProc, nProcs

    ! DATA VARIABLES
    integer :: receive_var_id
    real(r8),allocatable :: arr_receive_lat(:)
    real(r8),allocatable :: lat_recv(:)

    character(len=*), parameter :: sub = 'cw_export_mct'

    ! local grid order
    allocate(gindex(bounds%begg:bounds%endg),stat=ier)
    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = bounds%endg - bounds%begg + 1        ! local grid size
    gsize = ldomain%ni * ldomain%nj              ! global grid size

    call MPI_COMM_SIZE( mpicom, nProcs, ier)

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


    allocate(lat_recv(lsize)) ! scattered local array after receiving from CyberWater

    if (masterproc) then
      write(iulog,*) "Receive data from server !"
      
      receive_var_id = 4
      allocate(arr_receive_lat(gsize))
      call retrieve_variable_data(id, receive_var_id, arr_receive_lat, gsize)

      ! Scatter
      call MPI_SCATTERV(arr_receive_lat, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
                    lat_recv, lsize, MPI_DOUBLE, 0, mpicom, ier)

      ! Clean up after scattering
      deallocate(arr_receive_lat)

    end if
    call MPI_Barrier(mpicom, ier)

    write(6,*) "Checking MPI scatter!"
    write(6,*) "ldomain%lonc", ldomain%latc
    write(6,*) "lat_recv", lat_recv

    deallocate(lat_recv)
       
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

    ! Clean up allocated resources
    if (masterproc) then
      deallocate(nCellsPerProc)
      deallocate(nCellsDisplacement)
      deallocate(indexToCellIDGathered)
    end if

  end subroutine cw_export_mct


  subroutine retrieve_variable_data(session_id, var_id, arr_receive, gsize)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Receive data for individual variable
    !
    ! !USES:
    use elm_varctl     , only : iulog
    use shr_kind_mod   , only : r8 => shr_kind_r8

    implicit none

    ! !ARGUMENTS:
    type(SessionID)   , intent(in)    :: session_id
    integer, intent(in)    :: var_id, gsize
    real(r8),intent(inout) :: arr_receive(:)

    !
    ! !LOCAL VARIABLES:
    integer :: receive_var_size, receive_status

    ! Check data availability with retries
    if (check_data_availability_with_retries(var_id, recv_retry_count, recv_retry_time) == 1) then
      ! Retrieve the variable size
      receive_var_size = retrieve_variable_size(session_id, var_id)
      print *, "Variable size for ID", var_id, "is", receive_var_size

      ! Check if the received size matches the expected global size
      if (receive_var_size .NE. gsize) then
        call endrun('Variable receive size not equal to global mesh size!')
      endif

      ! Receive data with retries
      receive_status = receive_data_with_retries(var_id, arr_receive, recv_retry_count, recv_retry_time, use_delay=use_delay)
      if (receive_status == 1) then
        print *, "Received data ID=", var_id
      else
        print *, "Failed to receive data for variable ID:", var_id
      endif

      ! Sleep for a defined period
      print *, "------ Sleeping for 1 seconds ------"
      call sleep(1)

    else
      ! Handle the case where data is not available after retries
      print *, "Data is not available after retries."
      call endrun('Data is not available after retries!')
    endif

  end subroutine retrieve_variable_data

end module cw_import_export
