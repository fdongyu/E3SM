module cyberwaterMod

  use spmdMod      , only : mpicom, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use abortutils   , only: endrun
  use shr_sys_mod  , only : shr_sys_flush, shr_sys_abort
  use fileutils    , only : getavu, relavu
  use elm_nlUtilsMod   , only : find_nlgroup_name
  use controlMod       , only: NLFilename
  use mct_mod      , only : mct_avect
  use decompmod    , only : bounds_type
  use seq_cdata_mod, only : seq_cdata
  use domainMod    , only : ldomain
  use mpi
  use cw_import_export
  use cw_cpl_indices
  use high_level_api

  implicit none

  private

  integer, public    :: iulog = 6
  type(session_data) :: sd
  type(SessionID) :: id


  public :: cyberwater_init
  public :: cyberwater_run
  public :: cyberwater_final

contains


  !-----------------------------------------------------------------------
  subroutine cyberwater_init()
    !
    ! !DESCRIPTION:
    ! Initialize cyberwater
    !
    ! !USES:
    !
    implicit none
    !
    ! !LOCAL VARIABLES:
    integer :: gsize                  ! GS Map size
    integer :: ier, nml_error         ! Error code
    integer  :: nu_nml                               ! unit for namelist file
    integer :: session_status
    integer :: initiator_id, invitee_id
    character(len=256):: url               = ' '
    character(len=256):: stream_cyberwater = ' '     ! cyberwater parameter input file
    logical  :: lexist                               ! File exists
    character(len=*),parameter :: subname = '(cyberwater_init)'


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
              url, initiator_id, invitee_id

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


      write(iulog,*) 'Start session'

      gsize = ldomain%ni * ldomain%nj                          ! global grid size
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
      sd%input_variables_size = [gsize, gsize, gsize, &
                                 gsize, gsize, gsize, &
                                 gsize, gsize, gsize, &
                                 gsize, gsize]
      sd%output_variables_ID = [index_Sl_t, index_Sl_snowh, index_Sl_albd, index_Fall_lat, index_Fall_sen, index_Fall_lwup, index_Fall_evap, index_Fall_swnet, index_Sl_ram1]
      sd%output_variables_size = [gsize, gsize, gsize, &
                                  gsize, gsize, gsize, &
                                  gsize, gsize, gsize]

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

    end if

    call MPI_Barrier(mpicom, ier)

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
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)

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

    sd%initiator_id = 35

    if (masterproc) then
      write(iulog,*) 'Finishing CyberWater'
      write(iulog,*) 'End session'
      call end_session_now(sd%initiator_id)
    end if

  end subroutine cyberwater_final


end module cyberwaterMod
