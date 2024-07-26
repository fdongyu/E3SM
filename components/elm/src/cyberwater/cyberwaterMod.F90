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
!  use cw_import_export
!  use cw_cpl_indices
!  use data_exchange

  implicit none

  private

  integer, public    :: iulog = 6
!  integer(c_int), dimension(5) :: id
!  type(session_data) :: sd


  public :: cyberwater_init
!  public :: cyberwater_run
!  public :: cyberwater_final

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
      

    end if

  end subroutine cyberwater_init

end module cyberwaterMod
