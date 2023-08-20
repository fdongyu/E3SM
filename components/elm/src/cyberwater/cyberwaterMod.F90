module cyberwaterMod

  use spmdMod      , only : mpicom, masterproc
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush

  implicit none

  private

  integer, public    :: iulog = 6


  public :: cyberwater_init

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


end module cyberwaterMod
