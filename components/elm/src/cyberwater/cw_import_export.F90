module cw_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use abortutils   , only: endrun
  use spmdMod      , only : mpicom, masterproc
  use decompmod    , only : bounds_type, ldecomp
  use elm_instMod  , only : lnd2atm_vars, atm2lnd_vars
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
  !real(r8) :: recv_delay       = 5.0_r8
  real(r8) :: recv_delay       = 0.0_r8

  ! Receiving Params
  integer :: recv_retry_time  = 10
  integer  :: recv_retry_count = 25
  !logical  :: use_delay = .true.
  logical  :: use_delay = .false.


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
    integer :: g, i, j, n
    integer :: i1, i2, subdomain_size
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
    real(r8) :: t_start, t_end, t_elapsed


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
!    gsize = ldomain%ni * ldomain%nj              ! global grid size
    call MPI_Allreduce(lsize, gsize, 1, MPI_INTEGER, MPI_SUM, mpicom, ier) ! active global grid size

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

    call MPI_Barrier(mpicom, ier)



    ! send data to server 
    if (masterproc) then

      ! indexing array
      subdomain_size = size(subdomain_ind)
      i1 = 1
      i2 = min(subdomain_size, 10)   ! first 10
      write(iulog,*) 'subdomain_ind length = ', subdomain_size
      write(iulog,*) 'First ', i2, ' entries:'
      write(iulog,'(A)') '  idx      lat            lon'
      write(iulog,'(I6,2X,ES16.8,2X,ES16.8)') ( subdomain_ind(i), glat(subdomain_ind(i)), &
                                            glon(subdomain_ind(i)), i=i1,i2 )



      ! Prepare data 
      ! This calculation depends on the hydrological model in CyberWater
      do i = 1, gsize
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

      t_start = MPI_WTIME()

      send_status = send_data_with_retries(index_a2l_latitude, glat(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_a2l_longitude, glon(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_z, Sa_z(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_vel, Sa_vel(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_shum, Sa_shum(subdomain_ind), send_retry_count, send_retry_time, recv_delay)
      
      send_status = send_data_with_retries(index_Sa_pbot, Sa_pbot(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_tbot, Sa_tbot(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Faxa_lwdn, Faxa_lwdn(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Faxa_precip, Faxa_precip(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Faxa_sw, Faxa_sw(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      send_status = send_data_with_retries(index_Sa_co2prog, Sa_co2prog(subdomain_ind), send_retry_count, send_retry_time, recv_delay)

      t_end = MPI_WTIME()
      t_elapsed = t_end - t_start
      write(iulog,*) "Time used for sending data from server: ", t_elapsed, " seconds"
    
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

  subroutine cw_export_mct(bounds, x2l, l2x, id)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the data to be sent from the Cyberwater model to the coupler 
    ! 
    ! !USES:
    use elm_varctl     , only : iulog
    use domainMod      , only : ldomain
    use shr_kind_mod   , only : r8 => shr_kind_r8
    
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds   ! bounds
    real(r8)          , intent(in)    :: x2l(:,:) ! driver import state to land model
    real(r8)          , intent(out)   :: l2x(:,:) ! land to coupler export state on land grid
    type(SessionID)   , intent(in)    :: id
    
    !
    ! !LOCAL VARIABLES:
    integer  :: g, i, j, n
    integer  :: subdomain_size
    integer,allocatable :: gindex(:)  ! Number the local grid points
    integer :: lsize,gsize            ! GS Map size
    integer :: ier                    ! Error code

    ! MPI variables
    integer :: iProc, nProcs

    ! DATA VARIABLES
    integer :: receive_var_id
    real(r8),allocatable :: arr_receive_lat(:)
    real(r8),allocatable :: lat_recv(:)
    ! Global, variables read from x2l
    real(r8),allocatable :: Sa_u(:)    ! zonal wind velocity
    real(r8),allocatable :: Sa_v(:)    ! meridional wind velocity
    real(r8),allocatable :: Sa_shum(:) ! bottom atm level spec hum Pa
    real(r8),allocatable :: Sa_tbot(:) ! bottom atm level temp     degree
    real(r8),allocatable :: Faxa_swndr(:) ! direct near-infrared incident solar radiation W m-2
    real(r8),allocatable :: Faxa_swvdr(:) ! direct visible indicent solar radiation W m-2
    real(r8),allocatable :: Faxa_swndf(:) ! diffuse near-infrared incident solar radiation W m-2
    real(r8),allocatable :: Faxa_swvdf(:) ! diffuse visible incident solar radiation W m-2

    ! Global, variable received from cyberwater
    real(r8),allocatable :: Sl_t_recv(:)    ! Surface temperature K
    real(r8),allocatable :: Sl_snowh_recv(:)! Surface snow water equivalent m
    real(r8),allocatable :: Sl_albd_recv(:) ! Surface albedo average 
    real(r8),allocatable :: Fall_lat_recv(:)! Latent heat flux W m-2
    real(r8),allocatable :: Fall_sen_recv(:)! Sensible heat flux W m-2
    real(r8),allocatable :: Fall_lwup_recv(:) ! Upward longwave heat flux W m-2
    real(r8),allocatable :: Fall_evap_recv(:) ! Evaporation water flux  kg m-2 s-1
    real(r8),allocatable :: Fall_swnet_recv(:)! Heat flux shortwave net W m-2
    real(r8),allocatable :: Sl_ram1_recv(:)   ! Aerodynamic resistance  s/m

    ! Derived, variables
    real(r8),allocatable :: Sa_vel(:)  ! bottom atm level zon wind velocity sqrt(Sa_u*Sa_u+Sa_v*Sa_v)  m/s
    real(r8),allocatable :: Faxa_swndr_cw(:) ! direct near-infrared incident solar radiation W m-2
    real(r8),allocatable :: Faxa_swvdr_cw(:) ! direct visible indicent solar radiation W m-2
    real(r8),allocatable :: Faxa_swndf_cw(:) ! diffuse near-infrared incident solar radiation W m-2 
    real(r8),allocatable :: Faxa_swvdf_cw(:) ! diffuse visible incident solar radiation W m-2

    ! Local for scatter
    real(r8),allocatable :: Sl_t_local(:)
    real(r8),allocatable :: Sl_snowh_local(:)
    real(r8),allocatable :: Faxa_swndr_cw_local(:)
    real(r8),allocatable :: Faxa_swvdr_cw_local(:)
    real(r8),allocatable :: Faxa_swndf_cw_local(:)
    real(r8),allocatable :: Faxa_swvdf_cw_local(:)
    real(r8),allocatable :: Sa_shum_local(:)
    real(r8),allocatable :: Sa_tbot_local(:)
    real(r8),allocatable :: Sa_vel_local(:)
    real(r8),allocatable :: Sa_u_local(:)
    real(r8),allocatable :: Sa_v_local(:)
    real(r8),allocatable :: Fall_lat_local(:)
    real(r8),allocatable :: Fall_sen_local(:)
    real(r8),allocatable :: Fall_lwup_local(:)
    real(r8),allocatable :: Fall_evap_local(:)
    real(r8),allocatable :: Fall_swnet_local(:)
    real(r8),allocatable :: Sl_ram1_local(:)

    ! Local
    real(r8) :: ubot, vbot
    real(r8) :: swndf, swndr, swvdf, swvdr, sw
    real(r8) :: fswndf, fswndr, fswvdf, fswvdr
    real(r8) :: t_start, t_end, t_elapsed
    

    character(len=*), parameter :: sub = 'cw_export_mct'

    ! local grid order
    allocate(gindex(bounds%begg:bounds%endg),stat=ier)
    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = bounds%endg - bounds%begg + 1        ! local grid size
    !gsize = ldomain%ni * ldomain%nj              ! global grid size
    call MPI_Allreduce(lsize, gsize, 1, MPI_INTEGER, MPI_SUM, mpicom, ier) ! active global grid size

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

    ! Allocate global array for x2l variables
    if (masterproc) then
      allocate(Sa_vel(gsize))
      allocate(Sa_u(gsize))
      allocate(Sa_v(gsize))
      allocate(Sa_shum(gsize))
      allocate(Sa_tbot(gsize))
      allocate(Faxa_swndr(gsize))
      allocate(Faxa_swvdr(gsize))
      allocate(Faxa_swndf(gsize))
      allocate(Faxa_swvdf(gsize))
      allocate(Faxa_swndr_cw(gsize))
      allocate(Faxa_swvdr_cw(gsize))
      allocate(Faxa_swndf_cw(gsize))
      allocate(Faxa_swvdf_cw(gsize))

      ! Global variables, received from CyberWater
      subdomain_size = size(subdomain_ind)
      allocate(Sl_t_recv(subdomain_size)) 
      allocate(Sl_snowh_recv(subdomain_size))
      allocate(Sl_albd_recv(subdomain_size))
      allocate(Fall_lat_recv(subdomain_size))
      allocate(Fall_sen_recv(subdomain_size))
      allocate(Fall_lwup_recv(subdomain_size))
      allocate(Fall_evap_recv(subdomain_size))
      allocate(Fall_swnet_recv(subdomain_size))
      allocate(Sl_ram1_recv(subdomain_size))
    end if


    call MPI_GATHERV(x2l(index_x2l_Sa_u,:), lsize, MPI_DOUBLE, Sa_u, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) ! Atm state m/s

    call MPI_GATHERV(x2l(index_x2l_Sa_u,:), lsize, MPI_DOUBLE, Sa_u, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) ! Atm state m/s
    call MPI_GATHERV(x2l(index_x2l_Sa_v,:), lsize, MPI_DOUBLE, Sa_v, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) ! Atm state m/s
    call MPI_GATHERV(x2l(index_x2l_Sa_shum,:), lsize, MPI_DOUBLE, Sa_shum, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Sa_tbot,:), lsize, MPI_DOUBLE, Sa_tbot, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    ! Surface albedo
    call MPI_GATHERV(x2l(index_x2l_Faxa_swndr,:), lsize, MPI_DOUBLE, Faxa_swndr, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swvdr,:), lsize, MPI_DOUBLE, Faxa_swvdr, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swndf,:), lsize, MPI_DOUBLE, Faxa_swndf, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(x2l(index_x2l_Faxa_swvdf,:), lsize, MPI_DOUBLE, Faxa_swvdf, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)

    call MPI_Barrier(mpicom, ier)
    

    ! scattered local array after receiving from CyberWater
    allocate(Sl_t_local(lsize))
    allocate(Sl_snowh_local(lsize))
    allocate(Faxa_swndr_cw_local(lsize))
    allocate(Faxa_swvdr_cw_local(lsize))
    allocate(Faxa_swndf_cw_local(lsize))
    allocate(Faxa_swvdf_cw_local(lsize))
    allocate(Sa_shum_local(lsize))
    allocate(Sa_tbot_local(lsize))
    allocate(Sa_vel_local(lsize))
    allocate(Sa_u_local(lsize))
    allocate(Sa_v_local(lsize))
    allocate(Fall_lat_local(lsize))
    allocate(Fall_sen_local(lsize))
    allocate(Fall_lwup_local(lsize))
    allocate(Fall_evap_local(lsize))
    allocate(Fall_swnet_local(lsize))
    allocate(Sl_ram1_local(lsize))


    if (masterproc) then
      write(iulog,*) "Receive data from server !"
      
      ! Receive data from server
      t_start = MPI_WTIME()
      call retrieve_variable_data(id, index_Sl_t, Sl_t_recv, subdomain_size)
      call retrieve_variable_data(id, index_Sl_snowh, Sl_snowh_recv, subdomain_size)
      call retrieve_variable_data(id, index_Sl_albd, Sl_albd_recv, subdomain_size)
      call retrieve_variable_data(id, index_Fall_lat, Fall_lat_recv, subdomain_size)
      call retrieve_variable_data(id, index_Fall_sen, Fall_sen_recv, subdomain_size)
      call retrieve_variable_data(id, index_Fall_lwup, Fall_lwup_recv, subdomain_size)
      call retrieve_variable_data(id, index_Fall_evap, Fall_evap_recv, subdomain_size)
      call retrieve_variable_data(id, index_Fall_swnet, Fall_swnet_recv, subdomain_size)
      call retrieve_variable_data(id, index_Sl_ram1, Sl_ram1_recv, subdomain_size)
      t_end = MPI_WTIME()
      t_elapsed = t_end - t_start
      write(iulog,*) "Time used for receiving data from server: ", t_elapsed, " seconds"

      ! Process the receive data on master proc
      ! four albedo, zonal and meridional wind stress 
      do i = 1, subdomain_size
        !swndr      = Faxa_swndr(i)   ! forc_solldxy Atm flux  W/m^2
        !swvdr      = Faxa_swvdr(i)   ! forc_solsxy  Atm flux  W/m^2
        !swndf      = Faxa_swndf(i)   ! forc_solldxy Atm flux  W/m^2
        !swvdf      = Faxa_swvdf(i)   ! forc_solsdxy Atm flux  W/m^2
        !sw         = swndr + swvdr + swndf + swvdf
        !fswndr     = swndr / sw
        !fswvdr     = swvdr / sw
        !fswndf     = swndf / sw
        !fswvdf     = swvdf / sw
        Faxa_swndr_cw(i) = 0.25 * Sl_albd_recv(i)  
        Faxa_swvdr_cw(i) = 0.25 * Sl_albd_recv(i)
        Faxa_swndf_cw(i) = 0.25 * Sl_albd_recv(i)
        Faxa_swvdf_cw(i) = 0.25 * Sl_albd_recv(i)

!        ubot       = Sa_u(i)         ! m/s
!        vbot       = Sa_v(i)         ! m/s
!        Sa_vel(i)  = sqrt(ubot*ubot + vbot*vbot)
      end do

!      ! Scatter received data
!      call MPI_SCATTERV(Sl_t, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sl_t_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sl_snowh, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sl_snowh_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Faxa_swndr_cw, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Faxa_swndr_cw_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Faxa_swvdr_cw, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Faxa_swvdr_cw_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Faxa_swndf_cw, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Faxa_swndf_cw_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Faxa_swvdf_cw, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Faxa_swvdf_cw_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sa_shum, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sa_shum_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sa_tbot, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sa_tbot_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sa_vel, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sa_vel_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sa_u, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sa_u_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sa_v, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sa_v_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Fall_lat, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Fall_lat_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Fall_sen, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Fall_sen_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Fall_lwup, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Fall_lwup_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Fall_evap, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Fall_evap_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Fall_swnet, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Fall_swnet_local, lsize, MPI_DOUBLE, 0, mpicom, ier)
!      call MPI_SCATTERV(Sl_ram1, nCellsPerProc, nCellsDisplacement, MPI_DOUBLE, &
!                    Sl_ram1_local, lsize, MPI_DOUBLE, 0, mpicom, ier)

    end if
    call MPI_Barrier(mpicom, ier)

!    write(6,*) "Checking MPI scatter!"
!    write(6,*) "ldomain%lonc", ldomain%latc


!    ! Update lnd2atm_vars using CyberWater data
!    do g = bounds%begg,bounds%endg
!       i = 1 + (g-bounds%begg)
!       lnd2atm_vars%t_rad_grc(g)    = Sl_t_local(g)
!       lnd2atm_vars%h2osno_grc(g)   = Sl_snowh_local(g)
!       lnd2atm_vars%albd_grc(g,1)   = Faxa_swvdr_cw_local(g) ! Direct albedo (visible radiation)
!       lnd2atm_vars%albd_grc(g,2)   = Faxa_swndr_cw_local(g) ! Direct albedo (near-infrared radiation)
!       lnd2atm_vars%albi_grc(g,1)   = Faxa_swvdf_cw_local(g) ! Diffuse albedo (visible radiation)
!       lnd2atm_vars%albi_grc(g,2)   = Faxa_swndf_cw_local(g) ! Diffuse albedo (near-infrared radiation)
!       lnd2atm_vars%t_ref2m_grc(g)  = Sa_tbot_local(g)
!       lnd2atm_vars%q_ref2m_grc(g)  = Sa_shum_local(g)
!       lnd2atm_vars%u_ref10m_grc(g) = Sa_vel_local(g)
!       lnd2atm_vars%u_ref10m_with_gusts_grc(g) = Sa_vel_local(g) ! need to revisit wind gust
!       lnd2atm_vars%taux_grc(g)     = (-1) * atm2lnd_vars%forc_rho_not_downscaled_grc(g) * Sa_u_local(g) / Sl_ram1_local(g)
!       lnd2atm_vars%tauy_grc(g)     = (-1) * atm2lnd_vars%forc_rho_not_downscaled_grc(g) * Sa_v_local(g) / Sl_ram1_local(g)
!       lnd2atm_vars%eflx_lh_tot_grc(g)    = Fall_lat_local(g)
!       lnd2atm_vars%eflx_sh_tot_grc(g)    = Fall_sen_local(g)
!       lnd2atm_vars%eflx_lwrad_out_grc(g) = Fall_lwup_local(g)
!       lnd2atm_vars%qflx_evap_tot_grc(g)  = Fall_evap_local(g)
!       lnd2atm_vars%fsa_grc(g)            = Fall_swnet_local(g)
!    end do
       

!    This section is commented permantently
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

    deallocate(Sl_t_local)
    deallocate(Sl_snowh_local)
    deallocate(Faxa_swndr_cw_local)
    deallocate(Faxa_swvdr_cw_local)
    deallocate(Faxa_swndf_cw_local)
    deallocate(Faxa_swvdf_cw_local)
    deallocate(Sa_shum_local)
    deallocate(Sa_tbot_local)
    deallocate(Sa_vel_local)
    deallocate(Sa_u_local)
    deallocate(Sa_v_local)
    deallocate(Fall_lat_local)
    deallocate(Fall_sen_local)
    deallocate(Fall_lwup_local)
    deallocate(Fall_evap_local)
    deallocate(Fall_swnet_local)
    deallocate(Sl_ram1_local)

    ! Clean up allocated resources
    if (masterproc) then
      deallocate(nCellsPerProc)
      deallocate(nCellsDisplacement)
      deallocate(indexToCellIDGathered)
      deallocate(Sa_vel)
      deallocate(Sa_u)
      deallocate(Sa_v)
      deallocate(Sa_shum)
      deallocate(Sa_tbot)
      deallocate(Faxa_swndr)
      deallocate(Faxa_swvdr)
      deallocate(Faxa_swndf)
      deallocate(Faxa_swvdf)
      deallocate(Faxa_swndr_cw)
      deallocate(Faxa_swvdr_cw)
      deallocate(Faxa_swndf_cw)
      deallocate(Faxa_swvdf_cw)
      ! Global variables, received from CyberWater
      deallocate(Sl_t_recv)
      deallocate(Sl_snowh_recv)
      deallocate(Sl_albd_recv)
      deallocate(Fall_lat_recv)
      deallocate(Fall_sen_recv)
      deallocate(Fall_lwup_recv)
      deallocate(Fall_evap_recv)
      deallocate(Fall_swnet_recv)
      deallocate(Sl_ram1_recv)
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
!      print *, "------ Sleeping for 1 seconds ------"
!      call sleep(1)

    else
      ! Handle the case where data is not available after retries
      print *, "Data is not available after retries."
      call endrun('Data is not available after retries!')
    endif

  end subroutine retrieve_variable_data

end module cw_import_export
