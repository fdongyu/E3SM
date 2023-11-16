module cw_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8
  use abortutils   , only: endrun
  use spmdMod      , only : mpicom, masterproc
  use decompmod    , only : bounds_type, ldecomp
  use elm_cpl_indices
  use mpi
  use c_interface_combined

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
    integer :: Nvar=3
    real(r8),allocatable :: globalArray(:,:)  ! global array to pack all data
    real(r8),allocatable :: globalArray1D(:)  ! global array to pack all data
    real(r8),allocatable :: llat(:), llon(:)  ! local
    real(r8),allocatable :: glat(:), glon(:)  ! global
    real(r8),allocatable :: Sa_z(:)  ! bottom atm level height (m),  Sa_z_l(:) 

    ! communicate variables
    integer :: status_send, status_fetch
    real, dimension(5) :: arr_send = [1.2, 3.4, 5.6, 7.8, 9.10]
    real, dimension(5) :: arr_fetch

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


    ! Allocate local array
    !allocate(llon(bounds%begg:bounds%endg),stat=ier)
    !allocate(Sa_z_l(lsize))

    ! Allocate global array only on master processor
    if (masterproc) then
       allocate(glat(gsize))
       allocate(glon(gsize))
       allocate(Sa_z(gsize))
       allocate(globalArray(Nvar,gsize))
       allocate(globalArray1D(Nvar*gsize))
    end if 

    ! Gather lat and lon first
    call MPI_GATHERV(ldomain%latc, lsize, MPI_DOUBLE, glat, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)
    call MPI_GATHERV(ldomain%lonc, lsize, MPI_DOUBLE, glon, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier) 

    ! Gather variables before sending to CW
    call MPI_GATHERV(x2l(index_x2l_Sa_z,:), lsize, MPI_DOUBLE, Sa_z, nCellsPerProc, &
                        nCellsDisplacement, MPI_DOUBLE, 0, mpicom, ier)

    call MPI_Barrier(mpicom, ier)

    if (masterproc) then
       globalArray(1,:) = glat(:)
       globalArray(2,:) = glon(:)
       globalArray(3,:) = Sa_z(:)
       globalArray1D(1:gsize) = glat(:)
       globalArray1D(gsize+1:gsize*2) = glon(:)
       globalArray1D(gsize*2+1:gsize*3) = Sa_z(:)
!       write (6,*) glat
!       write (6,*) glon
!       write (6,*) Sa_z
!       write (6,*) globalArray
        write (6,*) globalArray1D
    end if


    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
!       Sa_z_l(i) = x2l(index_x2l_Sa_z,i)
!       write (6,*) ldomain%lonc(g), ldomain%latc(g)
!       write (6,*) x2l(index_x2l_Sa_z,i)         ! bottom atm level height   m
!       write (6,*) x2l(index_x2l_Sa_u,i)         ! bottom atm level zon wind m/s
!       write (6,*) x2l(index_x2l_Sa_v,i)         ! bottom atm level mer wind m/s
!       write (6,*) x2l(index_x2l_Sa_ptem,i)      ! bottom atm level pot temp Atm State K
!       write (6,*) x2l(index_x2l_Sa_shum,i)      ! bottom atm level spec hum Atm state kg/kg
!       write (6,*) x2l(index_x2l_Sa_pbot,i)      ! bottom atm level pressure Atm state Pa
!       write (6,*) x2l(index_x2l_Sa_tbot,i)      ! bottom atm level temp     Atm state K
!       write (6,*) x2l(index_x2l_Faxa_lwdn,i)    ! downward lw heat flux     Atm flux  W/m^2
!       write (6,*) x2l(index_x2l_Faxa_rainc,i)   ! prec: liquid "convective" mm/s
!       write (6,*) x2l(index_x2l_Faxa_rainl,i)   ! prec: liquid "large scale"  mm/s
!       write (6,*) x2l(index_x2l_Faxa_snowc,i)   ! prec: frozen "convective" mm/s
!       write (6,*) x2l(index_x2l_Faxa_snowl,i)   ! prec: frozen "large scale"  mm/s
!       write (6,*) x2l(index_x2l_Faxa_swndr,i)   ! sw: nir direct  downward Atm flux  W/m^2
!       write (6,*) x2l(index_x2l_Faxa_swvdr,i)   ! sw: vis direct  downward Atm flux  W/m^2 
!       write (6,*) x2l(index_x2l_Faxa_swndf,i)   ! sw: nir diffuse downward Atm flux  W/m^2
!       write (6,*) x2l(index_x2l_Faxa_swvdf,i)   ! sw: vis diffuse downward Atm flux  W/m^2

!       write (6,*) x2l(index_x2l_Sa_co2prog,i)   ! bottom atm level prognostic co2 
!       write (6,*) x2l(index_x2l_Sa_co2diag,i)   ! bottom atm level diagnostic co2

       ! atmosphere coupling, for prognostic/prescribed aerosols
!       write (6,*) x2l(index_x2l_Faxa_bcphidry,i) ! flux: Black Carbon hydrophilic dry deposition
!       write (6,*) x2l(index_x2l_Faxa_bcphodry,i) ! flux: Black Carbon hydrophobic dry deposition
!       write (6,*) x2l(index_x2l_Faxa_bcphiwet,i) ! flux: Black Carbon hydrophilic wet deposition
!       write (6,*) x2l(index_x2l_Faxa_ocphidry,i) ! flux: Organic Carbon hydrophilic dry deposition
!       write (6,*) x2l(index_x2l_Faxa_ocphodry,i) ! flux: Organic Carbon hydrophobic dry deposition
!       write (6,*) x2l(index_x2l_Faxa_ocphiwet,i) ! flux: Organic Carbon hydrophilic dry deposition
!       write (6,*) x2l(index_x2l_Faxa_dstwet1,i)  ! flux: Size 1 dust -- wet deposition
!       write (6,*) x2l(index_x2l_Faxa_dstdry1,i)  ! flux: Size 1 dust -- dry deposition
!       write (6,*) x2l(index_x2l_Faxa_dstwet2,i)  ! flux: Size 2 dust -- wet deposition
!       write (6,*) x2l(index_x2l_Faxa_dstdry2,i)  ! flux: Size 2 dust -- dry deposition
!       write (6,*) x2l(index_x2l_Faxa_dstwet3,i)  ! flux: Size 3 dust -- wet deposition
!       write (6,*) x2l(index_x2l_Faxa_dstdry3,i)  ! flux: Size 3 dust -- dry deposition
!       write (6,*) x2l(index_x2l_Faxa_dstwet4,i)  ! flux: Size 4 dust -- wet deposition
!       write (6,*) x2l(index_x2l_Faxa_dstdry4,i)  ! flux: Size 4 dust -- dry deposition
    end do


    if (masterproc) then
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
    end if



    ! free memory
    deallocate(gindex)
    !deallocate(Sa_z_l)


    if (masterproc) then
       deallocate(glat)
       deallocate(glon)
       deallocate(Sa_z)
       deallocate(nCellsPerProc)
       deallocate(nCellsDisplacement)
       deallocate(indexToCellIDGathered)
    end if


  end subroutine cw_import_mct

end module cw_import_export
