module cw_import_export

  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_sys_mod  , only : shr_sys_flush
  use abortutils   , only: endrun
  use spmdMod      , only : mpicom, masterproc
  use decompmod    , only : bounds_type, ldecomp
  use elm_cpl_indices
  use mpi
  use cw_cpl_indices
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

  type, public :: e3sm2cw_type

     integer :: id_source_model  ! identifier for source model
     integer :: id_dest_model  ! identifier for destination model 
     integer :: dataSize     ! data size
     real(r8),allocatable :: dataArray1D(:)  ! 1D data array
  
  end type e3sm2cw_type


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
    !integer :: Nvar=3
    !real(r8),allocatable :: globalArray1D(:)  ! global array to pack all data
    real(r8),allocatable :: glat(:), glon(:)  ! global
    real(r8),allocatable :: Sa_z(:)  ! bottom atm level height    m
    real(r8),allocatable :: Sa_u(:)  ! bottom atm level zon wind  m/s
    real(r8),allocatable :: Sa_v(:)  ! bottom atm level mer wind  m/s

    type(e3sm2cw_type) :: lat_sent, lon_sent

    ! communicate variables
    integer :: status_send, status_fetch
    real(r8),allocatable :: arr_fetch(:)
    integer :: Narr_fetch=2000

    ! periodically check if data is sent/received
    integer :: counter, threshold, value
    logical :: data_received
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
       !allocate(globalArray1D(Nvar*gsize))
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

    ! send_data_to_DEserver(int source_model_ID, int destination_model_ID, int variable_ID, int variable_unit,  double data_array, int num_elements, int error_val)
    if (masterproc) then
       write (6,*) 'source_model_ID, destination_model_ID, variable_ID, variable_unit, num_elements, data_array', index_EAM, index_VIC5, index_a2l_latitude, index_degree, gsize, glat
       write (6,*) 'source_model_ID, destination_model_ID, variable_ID, variable_unit, num_elements, data_array', index_EAM, index_VIC5, index_a2l_longitude, index_degree, gsize, glon
       write (6,*) 'source_model_ID, destination_model_ID, variable_ID, variable_unit, num_elements, data_array', index_EAM, index_VIC5, index_Sa_z, index_meter, gsize, Sa_z
       write (6,*) 'source_model_ID, destination_model_ID, variable_ID, variable_unit, num_elements, data_array', index_EAM, index_VIC5, index_Sa_u, index_meter, gsize, Sa_u
       write (6,*) 'source_model_ID, destination_model_ID, variable_ID, variable_unit, num_elements, data_array', index_EAM, index_VIC5, index_Sa_v, index_meter, gsize, Sa_v
    end if

    ! check if data is sent
    if (masterproc) then

       ! Initialize counter and threshold
       counter = 0
       threshold = 100
       data_received = .false.

       ! Loop until data is received or counter exceeds the threshold
       do while (.not. data_received .and. counter < threshold)
          counter = counter + 1
      
          ! Check if data is received 
          value = get_value()   ! Replace get_value() with check flag to obtain the value
          write(6,*) "random value=", value
          if (value <= 10) then ! Data is received
             data_received = .true.
          end if

          if (.not. data_received) then
             write (6,*) "Data not received. Sleeping for 30 seconds..."
             call sleep(30) 
          end if
       end do

       ! Check if data is received
       if (.not. data_received) then
          call endrun( sub//' ERROR: Maximum threshold reached. Data sent not successful' )
       else
          write(6, *) "Data is received"
       end if

    end if

!    if (masterproc) then
!       globalArray1D(1:gsize) = glat(:)
!       globalArray1D(gsize+1:gsize*2) = glon(:)
!       globalArray1D(gsize*2+1:gsize*3) = Sa_z(:)
!       write (6,*) glat
!       write (6,*) glon
!       write (6,*) Sa_z
!        write (6,*) globalArray1D
!    end if


!    if (masterproc) then
!       lon_sent%id_var_e3sm = 1
!       lon_sent%id_model_cw = 201
!       lon_sent%dataSize = gsize

!       ! Allocate and assign values to array_attribute
!       allocate(lon_sent%dataArray1D(gsize))
!       lon_sent%dataArray1D = glon

!       ! print for testing
!       write (6,*) 'lon_sent%id_var_e3sm=', lon_sent%id_var_e3sm
!       write (6,*) 'lon_sent%id_model_cw=', lon_sent%id_model_cw
!       write (6,*) 'lon_sent%dataSize=', lon_sent%dataSize
!       write (6,*) 'lon_sent%dataArray1D=', lon_sent%dataArray1D

!    end if


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
       !deallocate(globalArray1D)
       !deallocate(lon_sent%dataArray1D)
    end if


  end subroutine cw_import_mct



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
