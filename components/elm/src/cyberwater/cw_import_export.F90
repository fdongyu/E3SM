module cw_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8
  use spmdMod      , only : mpicom, masterproc
  use decompmod    , only : bounds_type
  use elm_cpl_indices

  implicit none

contains


  !-----------------------------------------------------------------------
  subroutine cw_import_mct( bounds, x2l)

    !---------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Convert the input data from the coupler to the land model
    use elm_varctl     , only: iulog
    use shr_kind_mod   , only: r8 => shr_kind_r8

    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds   ! bounds
    real(r8)           , intent(in)    :: x2l(:,:) ! driver import state to land model
    !
    ! !LOCAL VARIABLES:
    integer  :: g, i

    if (masterproc) then
      write(6,*) 'Importing variables from drv'
    end if

    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       !write(6,*) x2l(1,i)
       write (6,*) x2l(index_x2l_Sa_z,i)         ! bottom atm level height   m
       write (6,*) x2l(index_x2l_Sa_u,i)         ! bottom atm level zon wind m/s
       write (6,*) x2l(index_x2l_Sa_v,i)         ! bottom atm level mer wind m/s
       write (6,*) x2l(index_x2l_Sa_ptem,i)      ! bottom atm level pot temp Atm State K
       write (6,*) x2l(index_x2l_Sa_shum,i)      ! bottom atm level spec hum Atm state kg/kg
       write (6,*) x2l(index_x2l_Sa_pbot,i)      ! bottom atm level pressure Atm state Pa
       write (6,*) x2l(index_x2l_Sa_tbot,i)      ! bottom atm level temp     Atm state K
       write (6,*) x2l(index_x2l_Faxa_lwdn,i)    ! downward lw heat flux     Atm flux  W/m^2
       write (6,*) x2l(index_x2l_Faxa_rainc,i)   ! prec: liquid "convective" mm/s
       write (6,*) x2l(index_x2l_Faxa_rainl,i)   ! prec: liquid "large scale"  mm/s
       write (6,*) x2l(index_x2l_Faxa_snowc,i)   ! prec: frozen "convective" mm/s
       write (6,*) x2l(index_x2l_Faxa_snowl,i)   ! prec: frozen "large scale"  mm/s
       write (6,*) x2l(index_x2l_Faxa_swndr,i)   ! sw: nir direct  downward Atm flux  W/m^2
       write (6,*) x2l(index_x2l_Faxa_swvdr,i)   ! sw: vis direct  downward Atm flux  W/m^2 
       write (6,*) x2l(index_x2l_Faxa_swndf,i)   ! sw: nir diffuse downward Atm flux  W/m^2
       write (6,*) x2l(index_x2l_Faxa_swvdf,i)   ! sw: vis diffuse downward Atm flux  W/m^2

       write (6,*) x2l(index_x2l_Sa_co2prog,i)   ! bottom atm level prognostic co2 
       write (6,*) x2l(index_x2l_Sa_co2diag,i)   ! bottom atm level diagnostic co2

       ! atmosphere coupling, for prognostic/prescribed aerosols
       write (6,*) x2l(index_x2l_Faxa_bcphidry,i) ! flux: Black Carbon hydrophilic dry deposition
       write (6,*) x2l(index_x2l_Faxa_bcphodry,i) ! flux: Black Carbon hydrophobic dry deposition
       write (6,*) x2l(index_x2l_Faxa_bcphiwet,i) ! flux: Black Carbon hydrophilic wet deposition
       write (6,*) x2l(index_x2l_Faxa_ocphidry,i) ! flux: Organic Carbon hydrophilic dry deposition
       write (6,*) x2l(index_x2l_Faxa_ocphodry,i) ! flux: Organic Carbon hydrophobic dry deposition
       write (6,*) x2l(index_x2l_Faxa_ocphiwet,i) ! flux: Organic Carbon hydrophilic dry deposition
       write (6,*) x2l(index_x2l_Faxa_dstwet1,i)  ! flux: Size 1 dust -- wet deposition
       write (6,*) x2l(index_x2l_Faxa_dstdry1,i)  ! flux: Size 1 dust -- dry deposition
       write (6,*) x2l(index_x2l_Faxa_dstwet2,i)  ! flux: Size 2 dust -- wet deposition
       write (6,*) x2l(index_x2l_Faxa_dstdry2,i)  ! flux: Size 2 dust -- dry deposition
       write (6,*) x2l(index_x2l_Faxa_dstwet3,i)  ! flux: Size 3 dust -- wet deposition
       write (6,*) x2l(index_x2l_Faxa_dstdry3,i)  ! flux: Size 3 dust -- dry deposition
       write (6,*) x2l(index_x2l_Faxa_dstwet4,i)  ! flux: Size 4 dust -- wet deposition
       write (6,*) x2l(index_x2l_Faxa_dstdry4,i)  ! flux: Size 4 dust -- dry deposition

       
    end do

  end subroutine cw_import_mct

end module cw_import_export
