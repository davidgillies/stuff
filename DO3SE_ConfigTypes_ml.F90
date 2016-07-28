!-------------------------------------------------------------------------------
!> Configuration structures for the DO3SE model, including subroutines for
!! sanity-checking values and applying various pre-procesisng steps.
!-------------------------------------------------------------------------------
module DO3SE_ConfigTypes_ml

  use DO3SE_ModelConstants_ml, only: UNDEF, IUNDEF, MAX_LAYERS, MAX_LC
  use DO3SE_Util_ml
#include "interop_types.h"
#include "DO3SE_Util_ml.h"

  implicit none
  private

  integer :: i

  !> Location properties.
  TYPE, public :: Location_t
    REAL :: lat = UNDEF     !< Latitude (degrees North)
    REAL :: lon = UNDEF     !< Longitude (degrees East)
    REAL :: elev = UNDEF    !< Elevation (m above sea level)
    REAL :: albedo = UNDEF  !< Surface albedo (fraction)
    REAL :: Rsoil = 100     !< Soil resistance (s m-1)

    LOGICAL :: OTC = .false.  !< Is open top chamber experiment?
    REAL :: z_u = UNDEF       !< Measurement height for windspeed (m)
    REAL :: z_O3 = UNDEF      !< Measurement height for O3 (m)
    REAL :: h_u = UNDEF       !< Canopy height for windspeed measurement, default = target canopy (m)
    REAL :: h_O3 = UNDEF      !< Canopy height for O3 measurement, default = target canopy (m)
  end type Location_t
  public :: check_Location

  !> Configuration for meteorological data.
  TYPE, public :: MetConfig_t
    !> Method for supplying CO2 concentration:
    !!    - "constant": Use CO2_constant value
    !!    - "input":    CO2 concentration supplied
    CHARACTER(len=16) :: CO2_method = "constant"
    !> Constant CO2 concentration (ppm)
    REAL :: CO2_constant = 391.0
    !> Method for supplying O3 concentration:
    !!    - "constant": Use O3_constant value
    !!    - "offset":   Add O3_constant value to supplied O3 concentration
    !!    - "input":    O3 concentration supplied
    CHARACTER(len=16) :: O3_method = "input"
    !> Constant or offset amount for O3 concenttration (ppb)
    REAL :: O3_constant = UNDEF
  end type MetConfig_t

  !> Land cover season definition.
  TYPE, public :: Season_t
    !> Growing season method:
    !!    - "constant":         SGS and EGS supplied
    !!    - "all year":         SGS = 1 and EGS = 365
    !!    - "forest latitude":  Use EMEP forest latitude model
    !!    - "wheat latitude":   Use wheat latitude model
    CHARACTER(len=16) :: growing_season_method = "constant"
    ! Start and end of growing season
    INTEGER :: SGS = IUNDEF
    INTEGER :: EGS = IUNDEF

    !> Accumulation period method:
    !!    - "constant":         Astart and Aend supplied
    !!    - "growing season":   Astart = SGS and Aend = EGS
    !!    - "wheat latitude":   Use wheat latitude model
    CHARACTER(len=16) :: accumulation_period_method = "constant"
    ! Start and end of accumulation period
    INTEGER :: Astart = IUNDEF
    INTEGER :: Aend = IUNDEF

    ! LAI piecewise linear function parameters
    REAL :: LAI_a = UNDEF       !< LAI value at SGS (m2 m-2)
    REAL :: LAI_b = UNDEF       !< LAI value at SGS + LAI_1 (m2 m-2)
    REAL :: LAI_c = UNDEF       !< LAI value at EGS - LAI_2 (m2 m-2)
    REAL :: LAI_d = UNDEF       !< LAI value at EGS (m2 m-2)
    INTEGER :: LAI_1 = IUNDEF   !< Time from LAI_a to LAI_b (days)
    INTEGER :: LAI_2 = IUNDEF   !< Time from LAI_c to LAI_d (days)

    !> SAI method:
    !!    - "LAI":          SAI = LAI
    !!    - "forest":       Forest method: SAI = LAI + 1
    !!    - "wheat":        Based on wheat lifecycle, requires LAI PLF parameters
    CHARACTER(len=16) :: SAI_method = "LAI"
  end type Season_t
  public :: check_Season

  !> Configuration for multiplicative stomatal conductance functions.
  TYPE, public :: GstoConfig_t
    REAL :: fmin = UNDEF    !< Minimum stomatal conductance (fraction)
    REAL :: gmax = UNDEF    !< Maximum stomatal conductance (mmol O3 m-2 PLA s-1)
    REAL :: gmorph = 1.0    !< Sun/shade morphology factor (fraction)

    !> Stomatal conductance method:
    !!    - "multiplicative":   Use DO3SE multiplicative model
    !!    - "photosynthesis":   Use Farquar-based photosynthesis model (hybrid
    !!                          with some multiplicative components)
    CHARACTER(len=16) :: method = "multiplicative"

    !> f_phen method:
    !!    - "disabled":         f_phen supplied (or left at default value of 1.0)
    !!    - "simple day PLF":   "single hump" method using 3 values (_a, _c and
    !!                          _e) and 2 slope periods (_1 and _4)
    !!    - "complex day PLF":  "double hump" method using 5 values (_a to _e)
    !!                          and 4 slopes: _1 and _4 for start and end, _2
    !!                          and _3 for the middle section, starting at _limA
    !!                          and ending at _limB
    CHARACTER(len=16) :: f_phen_method = "simple day PLF"
    INTEGER :: f_phen_limA = IUNDEF   !< Start of soil water limitation
    INTEGER :: f_phen_limB = IUNDEF   !< End of soil water limitation
    REAL :: f_phen_a = UNDEF          !< f_phen at SGS
    REAL :: f_phen_b = UNDEF
    REAL :: f_phen_c = UNDEF
    REAL :: f_phen_d = UNDEF
    REAL :: f_phen_e = UNDEF          !< f_phen at EGS
    INTEGER :: f_phen_1 = IUNDEF      !< Time from f_phen_a to f_phen_b (days)
    INTEGER :: f_phen_2 = IUNDEF      !< Time from f_phen_b to f_phen_c (days)
    INTEGER :: f_phen_3 = IUNDEF      !< Time from f_phen_c to f_phen_d (days)
    INTEGER :: f_phen_4 = IUNDEF      !< Time from f_phen_d to f_phen_e (days)

    !> leaf_f_phen method:
    !!    - "disabled": leaf_f_phen supplied (or left at default value of 1.0)
    !!    - "f_phen":   leaf_f_phen = f_phen
    !!    - "day PLF":  "single hump" PLF between Astart and Aend
    CHARACTER(len=16) :: leaf_f_phen_method = "f_phen"
    REAL :: leaf_f_phen_a = UNDEF       !< f_phen at Astart
    REAL :: leaf_f_phen_b = UNDEF       !< f_phen at mid-season peak
    REAL :: leaf_f_phen_c = UNDEF       !< f_phen at Aend
    INTEGER :: leaf_f_phen_1 = IUNDEF   !< Time from _a to _b (days)
    INTEGER :: leaf_f_phen_2 = IUNDEF   !< Time from _b to _c (days)

    !> f_light method:
    !!    - "disabled":   f_light and leaf_f_light supplied (or left at default
    !!                    value of 1.0)
    !!    - "enabled":    f_light and leaf_f_light calculated
    CHARACTER(len=16) :: f_light_method = "enabled"
    REAL :: f_lightfac = 0.006  !< Single leaf f_light coefficient
    REAL :: cosA = 0.5          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)

    !> f_temp method:
    !!    - "disabled":       f_temp supplied (or left at default value of 1.0)
    !!    - "default":        Normal bell-shaped function over T_min -> T_opt -> T_max
    !!    - "square high":    Same as "default", but straight lines from
    !!                        (T_opt, 1.0) -> (T_max, 1.0) -> (T_max, 0.0)
    CHARACTER(len=16) :: f_temp_method = "default"
    REAL :: T_min = UNDEF     !< Minimum temperature (degrees C)
    REAL :: T_opt = UNDEF     !< Optimum temperature, for max. gsto (degrees C)
    REAL :: T_max = UNDEF     !< Maximum temperature (degrees C)

    !> f_VPD method:
    !!    - "disabled":   f_VPD supplied (or left at default value of 1.0)
    !!    - "linear":     Linear f_VPD relationship between VPD_max and VPD_min
    !!    - "log":        Simple, unparameterised, logarithmic relationship
    CHARACTER(len=16) :: f_VPD_method = "linear"
    REAL :: VPD_min = UNDEF     !< VPD for minimum gsto (kPa)
    REAL :: VPD_max = UNDEF     !< VPD for maximum gsto (kPa)

    !> f_SW method:
    !!    - "disabled":     f_SW supplied (or left at default value of 1.0)
    !!    - "fSWP exp":     Use fSWP exponential curve (see fSWP_exp_curve)
    !!    - "fSWP linear":  Use linear fSWP function (see SWP_min and SWP_max)
    !!    - "fLWP exp":     Use fSWP exponential curve, but with LWP instead of SWP
    !!    - "fPAW":         Use fPAW relationship
    CHARACTER(len=16) :: f_SW_method = "disabled"

    ! fSWP linear parameters:
    REAL :: SWP_min = UNDEF   !< SWP for minimum gsto (MPa)
    REAL :: SWP_max = UNDEF   !< SWP for maximum gsto (MPa)

    !> fSWP exponential curve:
    !!    - "custom":         fSWP_exp_a and fSWP_exp_b supplied
    !!    - "temperate":      a = 0.355, b = -0.706
    !!    - "mediterranean":  a = 0.619, b = -1.024
    CHARACTER(len=16) :: fSWP_exp_curve = "temperate"
    REAL :: fSWP_exp_a = UNDEF
    REAL :: fSWP_exp_b = UNDEF

    !> f_O3 method:
    !!    - "disabled":   f_O3 supplied (or left at default value of 1.0)
    !!    - "wheat":      Wheat f_O3 method
    !!    - "potato":     Potato f_O3 method
    CHARACTER(len=16) :: f_O3_method = "disabled"

    !> Critical daily VPD threshold above which stomatal conductance will stop
    !! increasing (kPa).
    REAL :: VPD_crit = 1000.0
  end type GstoConfig_t
  public :: check_GstoConfig

  !> Parameters for photosynthesis-based stomatal conductance.
  TYPE, public :: PnGstoConfig_t
    REAL :: g_sto_0 = UNDEF       !< Closed stomata conductance (umol m-2 s-1)
    REAL :: m = UNDEF             !< Species-specific sensitivity to An (dimensionless)

    !> V/J max method:
    !!    - "input":      V_cmax_25 and J_max_25 supplied as hourly inputs
    !!    - "constant":   Use constant V_cmax_25 and J_max_25 values (below)
    CHARACTER(len=16) :: V_J_method = "constant"

    REAL :: V_cmax_25 = UNDEF     !< Maximum catalytic rate at 25 degrees (umol m-2 s-1)
    REAL :: J_max_25 = UNDEF      !< Maximum rate of electron transport at 25 degrees (umol m-2 s-1)

    !> D_0 method:
    !!    - "constant":   Use constant D_0 value
    !!    - "f_VPD":      Determine D_0 from multiplicative f_VPD method (at f_VPD = 0.5)
    !!
    CHARACTER(len=16) :: D_0_method = "f_VPD"
    REAL :: D_0 = UNDEF   !< "The VPD at which g_sto is reduced by a factor of 2" (kPa) (Leuning et al. 1998)

    !> Tleaf method:
    !!    - "input":      Tleaf_C supplied
    !!    - "ambient":    Use ambient air temperature
    !!    - "Nikolov":    Use estimation method based on from Nikolov (1995)
    !!    - "EB":         Use estimation method based on "An Introduction to
    !!                    Environmental Biophysics" (Campbell & Norman, 1998)
    !!    - "de Boeck":   Use estimation method based on "Leaf temperatures in
    !!                    glasshouses and open-top chambers" (de Boeck, 2012)
    CHARACTER(len=16) :: Tleaf_method = "ambient"
    !> Threshold (from 0) to consider leaf energy balance equation as "balanced"
    REAL :: Tleaf_balance_threshold = 0.001
    !> Factor to apply to energy balance when adjusting leaf temperature
    REAL :: Tleaf_adjustment_factor = 0.02
    !> Maximum number of iterations to find Tleaf solution
    INTEGER :: Tleaf_max_iterations = 50

    !> O3 effect method:
    !!    - "disabled":     No effect
    !!    - "martin2000":   Use method from Martin (2000) paper, adapted to DO3SE flux
    CHARACTER(len=16) :: O3_method = "disabled"

    ! Parameters for martin2000 O3 effect method
    REAL :: K_z = 24.0      !< Coefficient for ozone damage (dimensionless)
    REAL :: F_0 = 1.0       !< Threshold flux for ozone damage (nmol O3 m-2 PLA s-1)

    !> Phenology method:
    !!    - "disabled":     No effect
    !!    - "leaf_f_phen":  Apply multiplicative leaf_f_phen to V_cmax_25 and J_max_25
    CHARACTER(len=16) :: phenology_method = "disabled"
  end type PnGstoConfig_t
  public :: check_PnGstoConfig

  !> Land cover properties.
  TYPE, public :: LandCover_t
    CHARACTER(len=128) :: name = ""
    type(Season_t) :: season = Season_t()
    type(GstoConfig_t) :: gsto = GstoConfig_t()
    type(PnGstoConfig_t) :: pn_gsto = PnGstoConfig_t()
    REAL :: height = UNDEF  !< Constant canopy height (m)
    REAL :: Lm = UNDEF      !< Leaf dimension (m)
    REAL :: Y = UNDEF       !< POD threshold (nmol O3 m-2 s-1)
  end type
  public :: check_LandCover
  public :: read_LandCover
  public :: read_all_LandCover

  !> Land cover configuration
  TYPE, public :: LandCoverConfig_t
    INTEGER :: nL = 1           !< Number of layers
    INTEGER :: nLC = 1          !< Number of land covers configured
    INTEGER :: primary_LC = 1   !< Primary land cover for height, LAI, etc.

    !> Land cover parameters
    type(LandCover_t), dimension(MAX_LC) :: LCs = LandCover_t()

    !> Height method:
    !!    - "input":      Input canopy height
    !!    - "constant":   Use constant height from primary land cover
    CHARACTER(len=16) :: height_method = "constant"

    !> Distribution of layers heights within canopy; height of top of layer as
    !! a fraction of the canopy height.
    REAL, dimension(MAX_LAYERS) :: layer_height = 1

    !> LAI method:
    !!    - "input":          Input all LAI values
    !!    - "input total":    Input total LAI, split according to fLAI
    !!    - "estimate total": Estimate total LAI from primary land cover,
    !!                        split according to fLAI
    CHARACTER(len=16) :: LAI_method = "estimate total"

    !> SAI method:
    !!    - "input":          Input all SAI values
    !!    - "input total":    Input total SAI, split according to fLAI
    !!    - "estimate total": Estimate total SAI from primary land cover and
    !!                        total LAI, split according to fLAI
    CHARACTER(len=16) :: SAI_method = "estimate total"

    !> Distribution of LAI/SAI in canopy (default: uniform distribution)
    REAL, dimension(MAX_LAYERS,MAX_LC) :: fLAI = 1
  end type
  public :: check_LandCoverConfig

contains

  subroutine check_Location(location)
    type(Location_t), intent(inout) :: location

    ASSERT_DEFINED(location%lat)
    ASSERT_DEFINED(location%lon)
    ASSERT_DEFINED(location%elev)
    ASSERT_DEFINED(location%albedo)
    ASSERT_DEFINED(location%Rsoil)

    if (.not. location%OTC) then
      ASSERT_DEFINED(location%z_u)
      ASSERT_DEFINED(location%z_O3)
    end if
  end subroutine check_Location

  subroutine check_Season(season, loc)
    type(Season_t), intent(inout) :: season
    type(Location_t), intent(in) :: loc

    select case (season%growing_season_method)
    case ("constant")
      ! Nothing to do, SGS and EGS supplied
    case ("all year")
      season%SGS = 1
      season%EGS = 366
    case ("forest latitude")
      season%SGS = nint(105 + ((loc%lat - 50) * 1.5) + ((loc%elev / 1000) * 10))
      season%EGS = nint(297 - ((loc%lat - 50) * 2.0) - ((loc%elev / 1000) * 10))
    case ("wheat latitude")
      season%SGS = nint((2.57 * loc%lat + 40) - 50)
      season%EGS = nint((2.57 * loc%lat + 40) + 42)
    case default
      UNKNOWN_STRING(season%growing_season_method)
    end select
    ASSERT_DEFINED(season%SGS)
    ASSERT_DEFINED(season%EGS)

    select case (season%accumulation_period_method)
    case ("constant")
      ! Nothing to do, Astart and Aend supplied
    case ("growing season")
      season%Astart = season%SGS
      season%Aend = season%EGS
    case ("wheat latitude")
      season%Astart = nint((2.57 * loc%lat + 40) - 15)
      season%Aend   = nint((2.57 * loc%lat + 40) + 40)
    case default
      UNKNOWN_STRING(season%accumulation_period_method)
    end select
    ASSERT_DEFINED(season%Astart)
    ASSERT_DEFINED(season%Aend)

    ASSERT_DEFINED(season%LAI_a)
    ASSERT_DEFINED(season%LAI_b)
    ASSERT_DEFINED(season%LAI_c)
    ASSERT_DEFINED(season%LAI_d)
    ASSERT_DEFINED(season%LAI_1)
    ASSERT_DEFINED(season%LAI_2)

    select case (season%SAI_method)
    case ("LAI", "forest", "wheat")
      ! Nothing to do
    case default
      UNKNOWN_STRING(season%SAI_method)
    end select
  end subroutine check_Season

  subroutine check_GstoConfig(gc)
    type(GstoConfig_t), intent(inout) :: gc

    ASSERT_DEFINED(gc%fmin)
    ASSERT_DEFINED(gc%gmax)

    ! Disable multiplicative gsto components that photosynthesis gsto replaces
    if (gc%method == "photosynthesis") then
      gc%f_light_method = "disabled"
      gc%f_temp_method = "disabled"
    end if

    ! Ensure necessary f_phen parameters are present
    if (gc%f_phen_method /= "disabled") then
      ASSERT_DEFINED(gc%f_phen_a)
      ASSERT_DEFINED(gc%f_phen_c)
      ASSERT_DEFINED(gc%f_phen_e)
      ASSERT_DEFINED(gc%f_phen_1)
      ASSERT_DEFINED(gc%f_phen_4)
      if (gc%f_phen_method == "complex day PLF") then
        ASSERT_DEFINED(gc%f_phen_b)
        ASSERT_DEFINED(gc%f_phen_d)
        ASSERT_DEFINED(gc%f_phen_limA)
        ASSERT_DEFINED(gc%f_phen_limB)
        ASSERT_DEFINED(gc%f_phen_2)
        ASSERT_DEFINED(gc%f_phen_3)
      end if
    end if

    select case (gc%leaf_f_phen_method)
    case ("disabled")
      ! Nothing to do
    case ("f_phen")
      ! *Don't* check the f_phen_method, because a disabled f_phen will instead
      ! nicely propagate to a disabled leaf_f_phen.
      !call assert(gc%f_phen_method /= "disabled", &
      !            "gsto%leaf_f_phen_method=f_phen but gsto%f_phen_method=disabled")
    case ("day PLF")
      ! Ensure necessary leaf_f_phen parameters are present
      ASSERT_DEFINED(gc%leaf_f_phen_a)
      ASSERT_DEFINED(gc%leaf_f_phen_b)
      ASSERT_DEFINED(gc%leaf_f_phen_c)
      ASSERT_DEFINED(gc%leaf_f_phen_1)
      ASSERT_DEFINED(gc%leaf_f_phen_2)
    case default
      UNKNOWN_STRING(gc%leaf_f_phen_method)
    end select

    ! Ensure necessary f_temp parameters are present
    if (gc%f_temp_method /= "disabled") then
      ASSERT_DEFINED(gc%T_min)
      ASSERT_DEFINED(gc%T_opt)
      ASSERT_DEFINED(gc%T_max)
    end if

    select case (gc%f_VPD_method)
    case ("disabled")
      ! Nothing to do
    case ("linear")
      ! Ensure necessary f_VPD parameters are present
      ASSERT_DEFINED(gc%VPD_max)
      ASSERT_DEFINED(gc%VPD_min)
    case ("log")
      ! Nothing to do
    case default
      UNKNOWN_STRING(gc%f_VPD_method)
    end select

    ! Parameterise fSWP exponential curve
    select case (gc%fSWP_exp_curve)
    case ("custom")
      ! Nothing to do
    case ("temperate")
      gc%fSWP_exp_a = 0.355
      gc%fSWP_exp_b = -0.706
    case ("mediterranean")
      gc%fSWP_exp_a = 0.619
      gc%fSWP_exp_b = -1.024
    case default
      UNKNOWN_STRING(gc%fSWP_exp_curve)
    end select

    ! Check that appropriate f_SW parameters are defined
    select case (gc%f_SW_method)
    case ("fSWP exp", "fLWP exp")
      ASSERT_DEFINED(gc%fSWP_exp_a)
      ASSERT_DEFINED(gc%fSWP_exp_b)
    case ("fSWP linear")
      ASSERT_DEFINED(gc%SWP_min)
      ASSERT_DEFINED(gc%SWP_max)
    end select
  end subroutine check_GstoConfig

  subroutine check_PnGstoConfig(pn_gsto, gsto)
    type(PnGstoConfig_t), intent(inout) :: pn_gsto
    type(GstoConfig_t), intent(in) :: gsto

    ASSERT_DEFINED(pn_gsto%g_sto_0)
    ASSERT_DEFINED(pn_gsto%m)

    select case (pn_gsto%V_J_method)
    case ("input")
      ! Nothing to do
    case ("constant")
      ASSERT_DEFINED(pn_gsto%V_cmax_25)
      ASSERT_DEFINED(pn_gsto%J_max_25)
    case default
      UNKNOWN_STRING(pn_gsto%V_J_method)
    end select

    !select case (pn_gsto%D_0_method)
    !case ("constant")
    !  ASSERT_DEFINED(pn_gsto%D_0)
    !case ("f_VPD")
      ! Check that f_VPD is enabled
    !  call assert(gsto%f_VPD_method /= "disabled", &
    !              "pn_gsto%D_0_method=f_VPD but gsto%f_VPD_method=disabled")
    !case default
    !  UNKNOWN_STRING(pn_gsto%D_0_method)
    !end select

    select case (pn_gsto%Tleaf_method)
    case ("input", "ambient", "Nikolov", "EB", "de Boeck")
      ! Nothing to do
    case default
      UNKNOWN_STRING(pn_gsto%Tleaf_method)
    end select

    select case (pn_gsto%O3_method)
    case ("disabled", "martin2000")
      ! Nothing to do
    case default
      UNKNOWN_STRING(pn_gsto%O3_method)
    end select

    select case (pn_gsto%phenology_method)
    case ("disabled")
      ! Nothing to do
    case ("leaf_f_phen")
      ! *Don't* check the leaf_f_phen_method, because a disabled leaf_f_phen
      ! will instead nicely propagate to disabling photosynthesis phenology effect.
      !call assert(gsto%leaf_f_phen_method /= "disabled", &
      !            "pn_gsto%phenology_method=leaf_f_phen but gsto%leaf_f_phen_method=disabled")
    case default
      UNKNOWN_STRING(pn_gsto%phenology_method)
    end select
  end subroutine check_PnGstoConfig

  subroutine check_LandCover(LC, loc)
    type(LandCover_t), intent(inout) :: LC  ! Land cover to check
    type(Location_t), intent(in) :: loc     ! Location for model run

    call check_Season(LC%season, loc)
    call check_GstoConfig(LC%gsto)
    if (LC%gsto%method == "photosynthesis") then
      call check_PnGstoConfig(LC%pn_gsto, LC%gsto)
    end if
    ASSERT_DEFINED(LC%height)
    ASSERT_DEFINED(LC%Lm)
    ASSERT_DEFINED(LC%Y)
  end subroutine check_LandCover

  subroutine check_LandCoverConfig(LCC, loc)
    type(LandCoverConfig_t), intent(inout) :: LCC
    type(Location_t), intent(in) :: loc

    integer :: i

    ASSERT(LCC%nLC >= 1)
    ASSERT(LCC%nLC <= MAX_LC)
    ASSERT(LCC%primary_LC >= 1)
    ASSERT(LCC%primary_LC <= LCC%nLC)
    ASSERT(LCC%nL >= 1)
    ASSERT(LCC%nL <= MAX_LAYERS)

    ! Check every land cover
    do i = 1, LCC%nLC
      call check_LandCover(LCC%LCs(i), loc)
    end do

    select case (LCC%height_method)
    case ("input", "constant")
    case default
      UNKNOWN_STRING(LCC%height_method)
    end select

    ! Remove irrelevant parts of layer_height
    LCC%layer_height(LCC%nL+1:) = 0
    ASSERT(all(LCC%layer_height(:LCC%nL) > 0.0))
    ASSERT(all(LCC%layer_height(:LCC%nL) <= 1.0))
    call assert(all(LCC%layer_height(:LCC%nL-1) > (LCC%layer_height(2:LCC%nL))), &
                "layer_height must decrease through the layers")

    select case (LCC%LAI_method)
    case ("input", "input total", "estimate total")
      ! Nothing to do
    case default
      UNKNOWN_STRING(LCC%LAI_method)
    end select

    select case (LCC%SAI_method)
    case ("input", "input total", "estimate total")
      ! Nothing to do
    case default
      UNKNOWN_STRING(LCC%SAI_method)
    end select

    ! Remove irrelevant parts of fLAI
    LCC%fLAI(:,LCC%nLC+1:) = 0
    LCC%fLAI(LCC%nL+1:,:) = 0
    ! Check there are no negatives, and that at least one value exists
    ASSERT(all(LCC%fLAI >= 0) .and. any(LCC%fLAI > 0))
    ! Normalise fLAI to sum to 1
    LCC%fLAI = LCC%fLAI / sum(LCC%fLAI)
  end subroutine check_LandCoverConfig

  !> Read a land cover definition from a namelist file.  Returns .true. or
  !! .false. depending on whether or not the read succeeded.
  logical function read_LandCover(unit, LC)
    integer, intent(in) :: unit
    type(LandCover_t), intent(out), target :: LC

    character(len=128), pointer :: name
    type(Season_t), pointer :: season
    type(GstoConfig_t), pointer :: gsto
    type(PnGstoConfig_t), pointer :: pn_gsto
    real, pointer :: Lm
    real, pointer :: Y
    namelist /DO3SE_LandCover/ LC, name, season, gsto, pn_gsto, Lm, Y

    LC = LandCover_t()
    name => LC%name
    season => LC%season
    gsto => LC%gsto
    pn_gsto => LC%pn_gsto
    Lm => LC%Lm
    Y => LC%Y

    read (unit=unit, nml=DO3SE_LandCover, end=100)
    read_LandCover = .true.
    return
100 read_LandCover = .false.
    return
  end function read_LandCover

  !> Read all land cover definitions from a namelist file.  Returns an
  !! allocatable array of the exact size to hold all land covers.
  function read_all_LandCover(filename) result(LCs)
    character(len=*), intent(in) :: filename
    type(LandCover_t), dimension(:), allocatable :: LCs

    integer :: unit
    type(LandCover_t), dimension(MAX_LC) :: buffer
    integer :: nLC
    logical :: success

    nLC = 0
    open (newunit=unit, file=filename, status="old", action="read", position="rewind")
    do while (nLC < MAX_LC)
      nLC = nLC + 1
      success = read_LandCover(unit, buffer(nLC))
      if (.not. success) then
        nLC = nLC - 1
        exit
      end if
    end do
    close (unit)

    LCs = buffer(1:nLC)
  end function read_all_LandCover

end module DO3SE_ConfigTypes_ml
