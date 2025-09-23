module cyberwaterMod

  use spmdMod      , only : mpicom, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use abortutils   , only: endrun
  use shr_sys_mod  , only : shr_sys_flush, shr_sys_abort
  use fileutils    , only : getavu, relavu
  use elm_nlUtilsMod   , only : find_nlgroup_name
  use controlMod       , only: NLFilename
  use mct_mod      , only : mct_avect
  use decompmod    , only : bounds_type, ldecomp
  use seq_cdata_mod, only : seq_cdata
  use domainMod    , only : ldomain
  use mpi
  use cw_import_export, only : cw_import_mct, cw_export_mct 
  use cw_cpl_indices
  use high_level_api

  implicit none

  private


  ! MPI variables
  integer :: nCellsGlobal
  integer, dimension(:), allocatable :: nCellsDisplacement, indexToCellIDGathered
  integer, dimension(:), allocatable :: nCellsPerProc

  integer, public    :: iulog = 6
  type(session_data) :: sd
  type(SessionID) :: id


  public :: cyberwater_init
  public :: cyberwater_run
  public :: cyberwater_final

contains


  !-----------------------------------------------------------------------
  subroutine cyberwater_init(bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cyberwater
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds    ! bounds
    ! !LOCAL VARIABLES:
    integer :: n
    integer :: lsize, gsize           ! GS Map size
    integer :: iProc, nProcs
    integer :: ier, nml_error         ! Error code
    integer :: nu_nml                               ! unit for namelist file
    integer :: session_status
    integer :: initiator_id, invitee_id
    integer :: subdomain_size
    real(r8):: gmin
    integer,allocatable  :: gindex(:)  ! Number the local grid points
    real(r8),allocatable :: glat(:), glon(:)  ! latitude, longitude
    real(r8),allocatable :: lats_in(:), lons_in(:)
    character(len=256):: url               = ' '
    character(len=256):: stream_cyberwater = ' '     ! cyberwater parameter input file
    character(len=256):: subdomain_latlon_file = ' '  ! cyberwater index input file
    logical  :: lexist                               ! File exists
    logical  :: is_glon_0360
    character(len=*),parameter :: subname = '(cyberwater_init)'

    ! local grid order 
    allocate(gindex(bounds%begg:bounds%endg),stat=ier)
    do n = bounds%begg, bounds%endg
       gindex(n) = ldecomp%gdc2glo(n)
    end do
    lsize = bounds%endg - bounds%begg + 1                                  ! local grid size
    !write (6,*) "lsize is ", lsize

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
       write(6,*) 'nCellsGlobal is ', nCellsGlobal
       write(6,*) 'gsize is ', gsize
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

    !if (masterproc) then
    !   write (6,*) nCellsPerProc
    !end if


    ! Allocate global array only on master processor
    if (masterproc) then
       allocate(glat(gsize))
       allocate(glon(gsize))
    end if

    call MPI_GATHERV(ldomain%latc, lsize, MPI_DOUBLE, glat, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(ldomain%lonc, lsize, MPI_DOUBLE, glon, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)


    namelist /cyberwater_streams / stream_cyberwater 
    inquire (file = trim(NLFilename), exist = lexist)
    if ( .not. lexist ) then
      write(iulog,*) subname // ' ERROR: NLFilename_dnstrm does NOT exist:'&
           //trim(NLFilename)
      call shr_sys_abort(trim(subname)//' ERROR NLFilename_dnstrm does not exist')
    end if

    if (masterproc) then

      write(iulog,*) 'Initializing coupling with cyberwater'

      nu_nml = getavu()
      write(iulog,*) 'Read in cyberwater_streams namelist from: ', trim(NLFilename)
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error)
      call find_nlgroup_name(nu_nml, 'cyberwater_streams', status=nml_error)
      if (nml_error == 0) then
          read(nu_nml, nml=cyberwater_streams,iostat=nml_error)
          if (nml_error /= 0) then
              call endrun(msg='ERROR reading cyberwater namelist')
          end if 
      end if     
      close(nu_nml)
      call relavu( nu_nml )
      
      write(iulog,*) 'stream cyberwater is:', stream_cyberwater


      ! query cyberwater parameter namelist
      namelist /cyberwater_inparm / &
              url, initiator_id, invitee_id, subdomain_latlon_file

      write(iulog,*) 'Read in cyberwater parameter namelist from: ', trim(stream_cyberwater)
      open( nu_nml, file=trim(stream_cyberwater), status='old', iostat=nml_error)
      call find_nlgroup_name(nu_nml, 'cyberwater_inparm', status=nml_error)
      if (nml_error == 0) then
          read(nu_nml, nml=cyberwater_inparm,iostat=nml_error)
          if (nml_error /= 0) then
              call endrun(msg='ERROR reading cyberwater namelist')
          end if
      end if
      close(nu_nml)
      call relavu( nu_nml )
      
      write(iulog,*) 'url is:', trim(url)
      write(iulog,*) 'initiator_id is:', initiator_id
      write(iulog,*) 'invitee_id is:',   invitee_id 
      write(iulog,*) 'subdomain_latlon_file: ', trim(subdomain_latlon_file)

      ! 1) read pairs from ASCII file
      if (len_trim(subdomain_latlon_file) == 0) call endrun('subdomain_latlon_file not set')
      call read_latlon_file(subdomain_latlon_file, lats_in, lons_in)
      subdomain_size = size(lats_in)
      write(iulog,*) 'Read subdomain size: ', subdomain_size
      if (subdomain_size == 0) call endrun('subdomain_index not set')

      if (subdomain_size > 0) then
        write(iulog,*) 'First 10: ', lats_in(1:min(10,subdomain_size)), lons_in(1:min(10,subdomain_size))
        write(iulog,*) 'Last 10 : ', lats_in(max(1,subdomain_size-9):subdomain_size), lons_in(max(1,subdomain_size-9):subdomain_size)
      end if

      ! 2) (optional) normalize longitudes to match glon convention
      gmin = minval(glon)
      if (gmin >= 0.0) then
          is_glon_0360 = .true. 
      else
          is_glon_0360 = .false.
      end if
      call normalize_lon_array(lons_in, want_0360 = is_glon_0360)  ! set is_glon_0360 accordingly

      ! 3) build subdomain_ind by nearest-neighbor on the global 1D lat/lon
      write(iulog,*) '--- Debug lat/lon arrays ---'
      write(iulog,*) 'glat size =', size(glat), ' glon size =', size(glon)
      write(iulog,*) 'glat(1:5) = ', glat(1:min(5,size(glat)))
      write(iulog,*) 'glon(1:5) = ', glon(1:min(5,size(glon)))
      write(iulog,*) 'glat(end-4:end) = ', glat(max(1,size(glat)-4):size(glat))
      write(iulog,*) 'glon(end-4:end) = ', glon(max(1,size(glon)-4):size(glon))

      write(iulog,*) 'lats_in size =', size(lats_in), ' lons_in size =', size(lons_in)
      write(iulog,*) 'lats_in(1:min(5,size(lats_in))) = ', lats_in(1:min(5,size(lats_in)))
      write(iulog,*) 'lons_in(1:min(5,size(lons_in))) = ', lons_in(1:min(5,size(lons_in)))

      call latlon_to_indices(glat, glon, lats_in, lons_in, subdomain_ind)
      write(iulog,*) 'Index First 10 ', subdomain_ind(1:10)


      write(iulog,*) 'Start session'

      ! Set the server URL at runtime
      call set_server_url(trim(url))

      ! User sets values directly
      sd%source_model_ID = index_EAM ! 2001
      sd%destination_model_ID = index_VIC5 ! 2005
      sd%initiator_id = initiator_id
      sd%invitee_id = invitee_id
      sd%input_variables_ID = [index_a2l_latitude, index_a2l_longitude, index_Sa_z, &
                               index_Sa_vel, index_Sa_shum, index_Sa_pbot, &
                               index_Sa_tbot, index_Faxa_lwdn, index_Faxa_precip, &
                               index_Faxa_sw, index_Sa_co2prog]
      sd%input_variables_size = [subdomain_size, subdomain_size, subdomain_size, &
                                 subdomain_size, subdomain_size, subdomain_size, &
                                 subdomain_size, subdomain_size, subdomain_size, &
                                 subdomain_size, subdomain_size]
      sd%output_variables_ID = [index_Sl_t, index_Sl_snowh, index_Sl_albd, index_Fall_lat, index_Fall_sen, index_Fall_lwup, index_Fall_evap, index_Fall_swnet, index_Sl_ram1]
      sd%output_variables_size = [subdomain_size, subdomain_size, subdomain_size, &
                                  subdomain_size, subdomain_size, subdomain_size, &
                                  subdomain_size, subdomain_size, subdomain_size]

      ! Write the session_ID for the whole program
      id = start_session(sd)
      call set_session_id(id)

      print *, "------ Sleeping for 10 seconds ------"
      call sleep(10)

      ! Join the session
      ! Check the session status first
      session_status = retrieve_session_status(id)
      if (session_status == 1) then
        print *, "Session status is 'created'"
      else if (session_status == 2) then
        print *, "Session status is 'active'"
      else if (session_status == 3) then
        print *, "Session status is 'partial end'"
      else
        print *, "Current status of session: ", session_status
      endif 

      call stream_updates()

       deallocate(nCellsPerProc)
       deallocate(nCellsDisplacement)
       deallocate(indexToCellIDGathered)
       deallocate(gindex)
       deallocate(glat)
       deallocate(glon)

    else
      subdomain_size = 0
    end if ! end of if-masterproc if-block 

    !call MPI_Barrier(mpicom, ier)

    call MPI_Bcast(subdomain_size, 1, MPI_INTEGER, 0, mpicom, ier)

    if (.not. masterproc) then
      if (allocated(subdomain_ind)) deallocate(subdomain_ind)
      allocate(subdomain_ind(subdomain_size))
    end if

    call mpi_bcast(subdomain_ind, subdomain_size, MPI_INTEGER, 0, mpicom, ier)

  end subroutine cyberwater_init


  !-----------------------------------------------------------------------
  subroutine cyberwater_run(EClock, bounds, cdata_l, x2l_l, l2x_l)
    !
    ! !DESCRIPTION:
    ! Run CyberWater
    !
    ! !USES:
    use shr_kind_mod    ,  only : r8 => shr_kind_r8
    !use clm_time_manager,  only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use shr_file_mod    ,  only : shr_file_setLogUnit, shr_file_setLogLevel
    use seq_timemgr_mod ,  only : seq_timemgr_EClockGetData
    use perf_mod         , only : t_startf, t_stopf, t_barrierf
    use ESMF
    !
    ! !ARGUMENTS:
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
    type(bounds_type), intent(in)    :: bounds    ! bounds
    type(seq_cdata)  , intent(inout) :: cdata_l   ! Input driver data for land model
    type(mct_aVect)  , intent(inout) :: x2l_l     ! Import state to land model
    type(mct_aVect)  , intent(inout) :: l2x_l     ! Export state from land model
    !
    ! !LOCAL VARIABLES
    integer      :: ymd_sync             ! Sync date (YYYYMMDD)
    integer      :: yr_sync              ! Sync current year
    integer      :: mon_sync             ! Sync current month
    integer      :: day_sync             ! Sync current day
    integer      :: tod_sync             ! Sync current time of day (sec)
    integer      :: ymd                  ! ELM current date (YYYYMMDD)
    integer      :: yr                   ! ELM current year
    integer      :: mon                  ! ELM current month
    integer      :: day                  ! ELM current day
    integer      :: tod                  ! ELM current time of day (sec)
    integer      :: dtime                ! time step increment (sec)
    integer      :: nstep                ! time step index
    integer      :: shrlogunit,shrloglev ! old values for share log unit and log level
    integer      :: ier                    ! Error code
    character(len=32), parameter    :: sub = "cyberwater_run"


    call seq_timemgr_EClockGetData(EClock, &
         curr_ymd=ymd, curr_tod=tod_sync,  &
         curr_yr=yr_sync, curr_mon=mon_sync, curr_day=day_sync)

    if (masterproc) then
      write(iulog, '(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)' ) yr_sync,mon_sync,day_sync,tod_sync
    end if

    ! update the clock in cyberwater

    ! import data from coupler
    call t_startf ('cyberwater_import')
    call cw_import_mct(bounds, x2l_l%rattr, id)
    call t_stopf ('cyberwater_import')

    call MPI_Barrier(mpicom, ier)


    ! run cyberwater
    call t_barrierf('sync_cyberwater_run', mpicom)
    call t_startf ('cyberwater_run')
    if (masterproc) then
      write(iulog,*) 'Running CyberWater'
    end if
    call t_stopf ('cyberwater_run')

    call MPI_Barrier(mpicom, ier)


    ! export data to coupler
    call t_startf ('cyberwater_export')
    call cw_export_mct(bounds, x2l_l%rattr, l2x_l%rattr, id)
    call t_stopf ('cyberwater_export')

    call MPI_Barrier(mpicom, ier)

    ! Reset shr logging to my original values
    !call shr_file_setLogUnit (shrlogunit)
    !call shr_file_setLogLevel(shrloglev)

  end subroutine cyberwater_run


  !-----------------------------------------------------------------------
  subroutine cyberwater_final()
    !
    ! !DESCRIPTION:
    ! Destroy RDy object
    !
    ! !USES:
    !
    implicit none
    !
    ! !LOCAL VARIABLES
    character(len=32), parameter    :: sub = "cyberwater_final"
    integer :: result

    sd%initiator_id = 35

    if (masterproc) then
      write(iulog,*) 'Finishing CyberWater'

      result = join_sse_thread_f()
      write(iulog,*) 'End session'
      call end_session_now(sd%initiator_id)

      if (allocated(subdomain_ind)) deallocate(subdomain_ind)

    end if

  end subroutine cyberwater_final


  subroutine read_latlon_file(fname, lats, lons)
    implicit none
    character(len=*), intent(in) :: fname
    real(kind=8),     allocatable, intent(out) :: lats(:), lons(:)

    integer :: u, ios, sz, i, n, p, q
    character(:), allocatable :: buf, tok
    real(kind=8) :: val
    real(kind=8), allocatable :: tmp(:)
    integer :: count

    inquire(file=trim(fname), size=sz, iostat=ios)
    if (ios /= 0 .or. sz <= 0) call endrun('Cannot stat file: '//trim(fname))

    open(newunit=u, file=trim(fname), access='stream', form='unformatted', status='old', iostat=ios)
    if (ios /= 0) call endrun('Cannot open file: '//trim(fname))

    allocate(character(len=sz) :: buf)
    read(u, iostat=ios) buf
    close(u)
    if (ios /= 0) call endrun('Stream read failed: '//trim(fname))

    ! normalize separators (comma/tab/CR/LF -> space)
    do i=1, len(buf)
      if (buf(i:i)==',' .or. buf(i:i)==achar(9) .or. buf(i:i)==achar(10) .or. buf(i:i)==achar(13)) buf(i:i)=' '
    end do

    ! pass 1: count tokens
    n = 0; p = 1
    do
      call next_token(buf, p, q)
      if (q < p) exit
      n = n + 1
      p = q + 1
    end do
    if (n < 2 .or. mod(n,2) /= 0) call endrun('latlon file must contain pairs: '//trim(fname))

    allocate(lats(n/2), lons(n/2))
    ! pass 2: read pairs
    p = 1; count = 0
    do i = 1, n/2
      call next_token(buf, p, q); read(buf(p:q),*,iostat=ios) lats(i); if (ios/=0) call endrun('bad lat token')
      p = q + 1
      call next_token(buf, p, q); read(buf(p:q),*,iostat=ios) lons(i); if (ios/=0) call endrun('bad lon token')
      p = q + 1
    end do

  contains
    subroutine next_token(s, p, q)
      character(len=*), intent(in) :: s
      integer, intent(inout) :: p
      integer, intent(out)   :: q
      integer :: L
      L = len(s)
      do while (p <= L .and. s(p:p) == ' ')
        p = p + 1
      end do
      if (p > L) then
        q = p - 1
        return
      end if
      q = p
      do while (q < L .and. s(q+1:q+1) /= ' ')
        q = q + 1
      end do
    end subroutine
  end subroutine

  subroutine normalize_lon_array(lon, want_0360)
    implicit none
    real(kind=8), intent(inout) :: lon(:)
    logical, intent(in)         :: want_0360
    integer :: i
    if (want_0360) then
      do i=1,size(lon)
        lon(i) = modulo(lon(i), 360.0d0)
        if (lon(i) < 0.d0) lon(i) = lon(i) + 360.d0
      end do
    else
      do i=1,size(lon)
        lon(i) = modulo(lon(i)+180.0d0, 360.0d0)
        if (lon(i) < 0.d0) lon(i) = lon(i) + 360.0d0
        lon(i) = lon(i) - 180.0d0
      end do
    end if
  end subroutine


subroutine latlon_to_indices(glat1d, glon1d, lats_in, lons_in, ind)
  implicit none
  real(kind=8), intent(in)  :: glat1d(:), glon1d(:)   ! global grid, length = NG
  real(kind=8), intent(in)  :: lats_in(:), lons_in(:)  ! query points, length = M
  integer,     allocatable, intent(out) :: ind(:)      ! nearest 1-based indices

  integer :: NG, M, i
  real(kind=8) :: gmin, gmax, lonq, dlat_scale
  real(kind=8), allocatable :: d2(:), dlon(:)
  integer :: jmin(1)
  logical :: grid0360
  real(kind=8), parameter :: deg2rad = 3.141592653589793238462643d0/180.0d0

  NG = size(glat1d);  M = size(lats_in)
  if (size(glon1d) /= NG) call endrun('glat1d/glon1d size mismatch')
  if (size(lons_in) /= M) call endrun('lats_in/lons_in size mismatch')

  gmin = minval(glon1d);  gmax = maxval(glon1d)
  grid0360 = (gmin >= 0.0d0 .and. gmax <= 360.0d0)

  allocate(ind(M))
  allocate(d2(NG), dlon(NG))

  do i = 1, M
     ! Normalize the query longitude to the same domain as the grid
     lonq = lons_in(i)
     if (grid0360) then
        lonq = modulo(lonq, 360.0d0); if (lonq < 0.d0) lonq = lonq + 360.d0
     else
        lonq = modulo(lonq + 180.0d0, 360.0d0); if (lonq < 0.d0) lonq = lonq + 360.d0
        lonq = lonq - 180.0d0
     end if

     ! Shortest wrapped Δlon in degrees (handles the 0/360 seam cleanly)
     dlon = abs(glon1d - lonq)
     if (grid0360) dlon = min(dlon, 360.0d0 - dlon)
     ! For -180..180 grids the above still works (period is 360)

     ! Scale longitude differences by cos(latitude) so degrees ~ km
     dlat_scale = cos(lats_in(i)*deg2rad)
     ! Guard near poles
     if (dlat_scale < 1.0d-12) dlat_scale = 1.0d-12

     ! Squared “distance” in degree-space with wrap-aware lon
     d2 = (glat1d - lats_in(i))**2 + (dlon * dlat_scale)**2

     ! Nearest neighbor
     jmin = minloc(d2)
     ind(i) = jmin(1)
  end do

  deallocate(d2, dlon)
end subroutine latlon_to_indices
 

end module cyberwaterMod
