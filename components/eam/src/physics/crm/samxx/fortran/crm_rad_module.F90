module crm_rad_module
   ! Module to encapsulate data and methods specific to radiation. This exists
   ! as a separate module with separate derived types because radiation is
   ! handled in a special way when interacting between the GCM and CRM.
   ! The radiation calculations may be on a separate grid than the CRM, and the
   ! quantities may be aggregated and updated between calls in a different way than
   ! the other diagnostic quantities. This module should also contain methods to
   ! update the radiation in the future, should we choose to put the radiation
   ! calculations on the CRM.
   use params,       only: crm_rknd

   implicit none

   public crm_rad_type

   type crm_rad_type
      ! Radiative heating
      real(crm_rknd), allocatable :: qrad(:,:,:,:)

      ! Quantities used by the radiation code. Note that these are strange in that they are 
      ! time-averages, but spatially-resolved.
      real(crm_rknd), allocatable :: temperature(:,:,:,:) ! rad temperature
      real(crm_rknd), allocatable :: qv (:,:,:,:) ! rad vapor
      real(crm_rknd), allocatable :: qc (:,:,:,:) ! rad cloud water
      real(crm_rknd), allocatable :: qi (:,:,:,:) ! rad cloud ice
      real(crm_rknd), allocatable :: cld(:,:,:,:) ! rad cloud fraction

      ! Only relevant when using 2-moment microphysics
      real(crm_rknd), allocatable :: nc(:,:,:,:) ! rad cloud droplet number (#/kg)
      real(crm_rknd), allocatable :: ni(:,:,:,:) ! rad cloud ice crystal number (#/kg)
      real(crm_rknd), allocatable :: qs(:,:,:,:) ! rad cloud snow (kg/kg)
      real(crm_rknd), allocatable :: ns(:,:,:,:) ! rad cloud snow crystal number (#/kg)
   end type crm_rad_type

contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_rad_type
   subroutine crm_rad_initialize(rad)
      class(crm_rad_type), intent(inout) :: rad

      ! Nullify pointers
      ! rad%qrad => null()
      ! rad%temperature => null()
      ! rad%qv => null()
      ! rad%qi => null()
      ! rad%cld => null()

      ! rad%nc => null()
      ! rad%ni => null()
      ! rad%qs => null()
      ! rad%ns => null()

   end subroutine crm_rad_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_rad_finalize(rad)
      class(crm_rad_type), intent(inout) :: rad

      ! Nullify pointers
      ! rad%qrad => null()
      ! rad%temperature => null()
      ! rad%qv => null()
      ! rad%qi => null()
      ! rad%cld => null()

      ! rad%nc => null()
      ! rad%ni => null()
      ! rad%qs => null()
      ! rad%ns => null()

   end subroutine crm_rad_finalize
   !------------------------------------------------------------------------------------------------

end module crm_rad_module
