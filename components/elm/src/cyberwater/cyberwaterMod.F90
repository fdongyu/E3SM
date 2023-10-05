module cyberwaterMod

  use spmdMod      , only : mpicom, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush

  implicit none

  private

  integer, public    :: iulog = 6


  public :: cyberwater_init
  public :: cyberwater_run
  public :: cyberwater_final

contains

  !-----------------------------------------------------------------------
  subroutine cyberwater_init()
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    !
    implicit none
    !
    ! !LOCAL VARIABLES:

    if (masterproc) then
      write(iulog,*) 'Hello CyberWater'
    end if


  end subroutine cyberwater_init


  !-----------------------------------------------------------------------
  subroutine cyberwater_run(EClock)
    !
    ! !DESCRIPTION:
    ! Initialize RDycore
    !
    ! !USES:
    !use clm_time_manager,  only : get_curr_date, get_nstep, get_curr_calday, get_step_size
    use shr_file_mod    ,  only : shr_file_setLogUnit, shr_file_setLogLevel
    use seq_timemgr_mod ,  only : seq_timemgr_EClockGetData
    use perf_mod         , only : t_startf, t_stopf, t_barrierf
    use ESMF
    !
    ! !ARGUMENTS:
    type(ESMF_Clock) , intent(inout) :: EClock    ! Input synchronization clock from driver
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
    !lnd_import_mct(x2l)
    call t_stopf ('cyberwater_import')
    
    ! run cyberwater
    call t_barrierf('sync_cyberwater_run', mpicom)
    call t_startf ('cyberwater_run')
    if (masterproc) then
      write(iulog,*) 'Running CyberWater'
    end if
    call t_stopf ('cyberwater_run')

    ! export data to coupler
    call t_startf ('cyberwater_export')
    !call lnd_export_mct(l2x)
    call t_stopf ('cyberwater_export')

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

    if (masterproc) then
      write(iulog,*) 'Finishing CyberWater'
    end if



  end subroutine cyberwater_final


end module cyberwaterMod
