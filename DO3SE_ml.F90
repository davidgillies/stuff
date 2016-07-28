module DO3SE_ml

  use DO3SE_ModelConstants_ml
  use DO3SE_PhysicalConstants_ml
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"
  use DO3SE_ConfigTypes_ml
  use DO3SE_Met_ml
  use DO3SE_Phenology_ml
  use DO3SE_Gsto_ml
  use DO3SE_Photosynthesis_ml
  use DO3SE_Resistance_ml
  use DO3SE_SMD_ml
#include "interop_types.h"

  implicit none
  public


  !> Variables that do not vary throughout the canopy.
  TYPE, public :: V_t
    INTEGER :: dd = IUNDEF                  !< Day of year
    INTEGER :: hr = IUNDEF                  !< Hour of day (0--23)

    REAL :: canopy_height = UNDEF           !< Overall canopy height

    type(MetData_t) :: met = MetData_t()    !< Meteorological data
    REAL :: VPD_dd = UNDEF                  !< Daily VPD sum during daylight hours (kPa)

    !> Multi-layer canopy O3 resistance model
    type(ResistanceModel_t) :: rmodel_O3 = ResistanceModel_t()
    REAL :: Vd = UNDEF                      !< Velocity of O3 deposition to top of canopy (m s-1)

    !> Single-layer canopy water vapour resistance model
    type(ResistanceModel_t) :: rmodel_H2O = ResistanceModel_t()
    type(PM_State_t) :: PM = PM_State_t()   !< Penman-Monteith results and accumulators
    type(SMDData_t) :: SMD = SMDData_t()    !< Results of SMD calculations
    LOGICAL :: Es_blocked = .false.         !< Should soil evaporation be blocked?
  end type V_t

  !> Variables that vary by canopy layer.
  TYPE, public :: ML_t
    REAL :: layer_height = UNDEF            !< Height of top of layer (m)

    !> Meteorological data
    type(MicroMetData_t) :: met = MicroMetData_t()
  end type ML_t

  !> Variables that vary by canopy component (land cover).
  TYPE, public :: MC_t
    REAL :: LC_dist = UNDEF
  end type MC_t

  !> Varibales that vary by canopy layer and component.
  TYPE, public :: MLMC_t
    REAL :: LAI = UNDEF                     !< Leaf area index (m2 m-2)
    REAL :: SAI = UNDEF                     !< Stand area index (m2 m-2)
    REAL :: LAIsunfrac = UNDEF              !< Fraction of LAI that is sunlit

    !> Multiplicative stomatal conductance parameters
    type(GstoParams_t) :: gsto_params = GstoParams_t()

    ! Photosynthetic stomatal conductance parameters/results
    REAL :: V_cmax_25 = UNDEF               !< Maximum catalytic rate at 25 degrees (umol m-2 s-1)
    REAL :: J_max_25 = UNDEF                !< Maximum rate of electron transport at 25 degrees (umol m-2 s-1)
    REAL :: A_n = UNDEF                     !< Net CO2 assimilation rate (umol CO2 m-2 PLA s-1)
    REAL :: FO3_eff = 0.0                   !< (Accumulated) effective ozone dose (nmol O3 m-2 PLA)
    REAL :: g_sv = UNDEF
    REAL :: g_bv = UNDEF

    REAL :: leaf_gsto = UNDEF               !< Leaf stomatal conductance (mmol O3 m-2 PLA s-1)
    REAL :: mean_gsto = UNDEF               !< Canopy mean stomatal conductance (mmol O3 m-2 PLA s-1)
    REAL :: bulk_gsto = UNDEF               !< Canopy total stomatal conductance (mmol O3 m-2 PLA s-1)

    !> Leaf-level O3 resistance model
    type(LeafResistanceModel_t) :: leaf_rmodel_O3 = LeafResistanceModel_t()

    REAL :: Fst = UNDEF                     !< Stomatal ozone flux (nmol O3 m-2 PLA s-1)
    REAL :: POD_0 = 0.0                     !< Phytotoxic Ozone Dose, no threshold (mmol m-2 PLA)
    REAL :: POD_Y = 0.0                     !< Phytotoxic Ozone Dose above threshold Y (mmol m-2 PLA)
    REAL :: OT_0 = UNDEF                    !< (ppm)
    REAL :: OT_40 = UNDEF                   !< (ppm)
    REAL :: AOT_0 = 0.0                     !< (ppm h)
    REAL :: AOT_40 = 0.0                    !< (ppm h)
  end type MLMC_t


  type, public :: DO3SE_State_t

    !
    ! Configuration
    !

    !> Location properties
    type(Location_t), pointer :: location
    !> Land cover configuration
    type(LandCoverConfig_t), pointer :: LC_conf
    type(LandCover_t), dimension(:), pointer :: LCs
    integer, pointer :: nL
    integer, pointer :: nLC
    real, dimension(:,:), pointer :: fLAI
    !> Meteorology configuration
    type(MetConfig_t), pointer :: met_conf
    !> Configuration for SMD module
    type(SMDConfig_t), pointer :: SMD_conf

    !
    ! Hourly data
    !

    type(V_t) :: V
    type(ML_t), dimension(:), allocatable :: ML
    type(MC_t), dimension(:), allocatable :: MC
    type(MLMC_t), dimension(:,:), allocatable :: MLMC

  contains

    procedure :: start_of_day
    procedure :: end_of_day
    procedure :: is_daylight

    procedure :: calc_phenology
    procedure :: calc_met_data
    procedure :: calc_micromet_data
    procedure :: calc_gsto_parameters
    procedure :: calc_gsto
    procedure :: calc_ozone_deposition
    procedure :: calc_ozone_dose
    procedure :: calc_soil_moisture
    procedure :: calc_soil_moisture_changes

    procedure :: init
    procedure :: run_hour

  end type DO3SE_State_t

contains

  !> Check configuration, allocate model variables, and initialise accumulators.
  subroutine init(this)
    class(DO3SE_State_t), intent(inout) :: this

    call check_Location(this%location)
    call check_LandCoverConfig(this%LC_conf, this%location)
    call check_SMDConfig(this%SMD_conf)

    this%nL => this%LC_conf%nL
    this%nLC => this%LC_conf%nLC
    this%LCs => this%LC_conf%LCs(1:this%nLC)
    this%fLAI => this%LC_conf%fLAI(1:this%nL,1:this%nLC)

    SAFE_ALLOC_1D(this%ML, this%nL)
    SAFE_ALLOC_1D(this%MC, this%nLC)
    SAFE_ALLOC_2D(this%MLMC, this%nL, this%nLC)

    call init_ResistanceModel(this%V%rmodel_O3, this%nL)

    if (this%SMD_conf%source == "P-M") then
      this%V%PM = PM_state_t()
      this%V%SMD = soil_moisture_from_SWC(this%SMD_conf, this%SMD_conf%initial_SWC)
    end if
  end subroutine init

  !> Is it currently the first hour of the day?
  logical function start_of_day(this)
    class(DO3SE_State_t), intent(in) :: this

    start_of_day = this%V%hr == 0
  end function start_of_day

  !> Is it currently the last hour of the day?
  logical function end_of_day(this)
    class(DO3SE_State_t), intent(in) :: this

    end_of_day = this%V%hr == 23
  end function end_of_day

  !> Is it currently daylight?  Uses the accepted criteria for accumulating
  !! OT40, i.e. when global radiation > 50.0 W m-2.
  logical function is_daylight(this)
    class(DO3SE_State_t), intent(in) :: this

    is_daylight = this%V%met%R > 50.0
  end function is_daylight

  !> Calculate canopy height, LAI, SAI, etc.
  !!
  !! Supplying LAI and SAI for multi-layer multi-component models can be
  !! awkward.  To generate a nL*nLC array for LAI/SAI, three schemes are
  !! possible:
  !!
  !!   - All LAI/SAI values are known and supplied as input.
  !!   - Total LAI/SAI is known and supplied as input in the top-left cell
  !!     (e.g. LAI(1,1)).  Values are divided among layers and land covers
  !!     according to a (normalised) nL*nLC fLAI array.
  !!   - Total LAI/SAI is estimated from the properties of the "primary land
  !!     cover" and divided according to fLAI.
  !!
  !! Canopy height is either an input or taken from the "primary land cover".
  subroutine calc_phenology(this)
    class(DO3SE_State_t), intent(inout) :: this

    associate (LC => this%LCs(this%LC_conf%primary_LC))
      select case (this%LC_conf%height_method)
      case ("input")
        ! Nothing to do
      case ("constant")
        ! Use height of primary land cover
        this%V%canopy_height = LC%height
      case default
        UNKNOWN_STRING(this%LC_conf%height_method)
      end select

      ! Calculate height of top of each canopy
      this%ML(:)%layer_height = this%V%canopy_height * this%LC_conf%layer_height(:this%nL)

      select case (this%LC_conf%LAI_method)
      case ("input")
        ! Nothing to do
      case ("input total")
        ! Spread single LAI value to layers and LCs
        this%MLMC(:,:)%LAI = this%MLMC(1,1)%LAI * this%fLAI(:,:)
      case ("estimate total")
        ! Use primary land cover's estimate of total LAI
        this%MLMC(1,1)%LAI = LAI_day_PLF(LC%season, this%V%dd)
        ! Spread single LAI value to layers and LCs
        this%MLMC(:,:)%LAI = this%MLMC(1,1)%LAI * this%fLAI(:,:)
      case default
        UNKNOWN_STRING(this%LC_conf%LAI_method)
      end select

      select case (this%LC_conf%SAI_method)
      case ("input")
        ! Nothing to do
      case ("input total")
        ! Spread single SAI value to layers and LCs
        this%MLMC(:,:)%SAI = this%MLMC(1,1)%SAI * this%fLAI(:,:)
      case ("estimate total")
        ! Use primary land cover's estimate of total SAI
        select case (LC%season%SAI_method)
        case ("LAI")
          this%MLMC(1,1)%SAI = sum(this%MLMC(:,:)%LAI)
        case ("forest")
          this%MLMC(1,1)%SAI = sum(this%MLMC(:,:)%LAI) + 1.0
        case ("wheat")
          this%MLMC(1,1)%SAI = SAI_wheat(LC%season, this%V%dd, LAI_day_PLF(LC%season, this%V%dd))
        case default
          UNKNOWN_STRING(LC%season%SAI_method)
        end select
        ! Spread single SAI value to layers and LCs
        this%MLMC(:,:)%SAI = this%MLMC(1,1)%SAI * this%fLAI(:,:)
      case default
        UNKNOWN_STRING(this%LC_conf%SAI_method)
      end select

      ! Calculate the distribution of LAI between land covers
      if (sum(this%MLMC(:,:)%LAI) <= 0.0) then
        this%MC(:)%LC_dist = 1.0 / (this%nL * this%nLC)
      else
        this%MC(:)%LC_dist = sum(this%MLMC(:,:)%LAI, 1) / sum(this%MLMC(:,:)%LAI)
      end if
    end associate
  end subroutine calc_phenology

  !> Calculate parts of meteorological data that need to be calculated.
  !!
  !! This includes:
  !!   * applying conversions between various kinds of similar input (e.g.
  !!     R vs. PAR vs. PPFD)
  !!   * estimates for missing data (e.g. net radiation)
  !!   * accumulators (VPD sum, thermal time)
  subroutine calc_met_data(this)
    class(DO3SE_State_t), intent(inout) :: this

    ! Fixup met data
    call met_humidity(this%V%met)
    call met_radiation(this%V%met, this%location, this%V%dd, this%V%hr)
    call met_CO2(this%met_conf, this%V%met)

    ! Reset VPD sum at start of day
    if (this%start_of_day()) then
      this%V%VPD_dd = 0
    end if
    ! Only accumulate VPD sum during daylight hours
    if (this%is_daylight()) then
      this%V%VPD_dd = this%V%VPD_dd + this%V%met%VPD
    end if
  end subroutine calc_met_data

  !> Calculate meterological data for the top of the canopy and within
  !! the canopy.
  subroutine calc_micromet_data(this)
    class(DO3SE_State_t), intent(inout) :: this

    call PAR_sun_shade(this%V%met%Idrctt, this%V%met%Idfuse, this%V%met%sinB, &
                       sum(this%LCs(1:this%nLC)%gsto%cosA)/this%nLC, sum(this%MLMC(:,:)%LAI), &
                       this%ML(1)%met%PARsun, this%ML(1)%met%PARshade)
    this%ML(:)%met%PARsun = this%ML(1)%met%PARsun
    this%ML(:)%met%PARshade = this%ML(1)%met%PARshade
    ! TODO: multi-layer PAR

    call met_windspeed(this%V%met, this%ML(1)%met, this%location, this%V%canopy_height)
    if (this%location%OTC) then
      ! In OTC, assume uniform windspeed being driven by a fan.
      this%ML(:)%met%u = this%ML(1)%met%u
    else
      ! Estimate windspeed within the canopy.  Uses total SAI, and an LAI-weighted mean leaf width.
      call multi_layer_windspeed(this%V%canopy_height, sum(this%MC(:)%LC_dist * this%LCs(:)%Lm), &
                                 sum(this%MLMC(:,:)%SAI), this%ML(1)%met%u, this%ML(:)%layer_height, &
                                 this%ML(:)%met%u)
    end if

    ! Estimate sunlit LAI fractions (used later in f_light)
    ! TODO: remove the per-land-cover component of this?
    this%MLMC(:,:)%LAIsunfrac = MLMC_sunlit_LAI(this%MLMC(:,:)%LAI, this%V%met%sinB)
  end subroutine calc_micromet_data

  !> Calculate the various effects on stomatal conductance for a specified layer
  !! and land cover.
  subroutine calc_gsto_parameters(this, LC, V, ll, xx)
    class(DO3SE_State_t), intent(inout) :: this
    type(LandCover_t), intent(in) :: LC
    type(V_t), intent(in) :: V
    type(ML_t), intent(in) :: ll
    type(MLMC_t), intent(inout) :: xx

    real :: delta_V_cmax, mult_V_cmax

    associate (season => LC%season, &
               gc => LC%gsto, &
               pgc => LC%pn_gsto, &
               gp => xx%gsto_params)
      ! Initialise gsto parameters from configuration
      gp = GstoParams_t(fmin=gc%fmin, gmax=gc%gmax, gmorph=gc%gmorph)

      ! Calculate f_phen
      select case (gc%f_phen_method)
      case ("disabled")
        ! Nothing to do
      case ("simple day PLF")
        gp%f_phen = f_phen_simple_PLF(gc, season%SGS, season%EGS, V%dd)
      case ("complex day PLF")
        gp%f_phen = f_phen_complex_PLF(gc, season%SGS, season%EGS, V%dd)
      case default
        UNKNOWN_STRING(gc%f_phen_method)
      end select

      ! Calculate leaf_f_phen
      select case (gc%leaf_f_phen_method)
      case ("disabled")
        ! Nothing to do
      case ("f_phen")
        gp%leaf_f_phen = gp%f_phen
      case ("day PLF")
        gp%leaf_f_phen = leaf_f_phen_PLF(gc, season%Astart, season%Aend, V%dd)
      case default
        UNKNOWN_STRING(gc%leaf_f_phen_method)
      end select

      select case (gc%f_light_method)
      case ("disabled")
        ! Nothing to do
      case ("enabled")
        if (xx%LAI > 0 .and. V%met%sinB > 0) then
          ! Calculate f_light and leaf_f_light
          ! TODO: attenuate PAR properly through the canopy
          gp%f_light = f_light(gc%f_lightfac, ll%met%PARsun, ll%met%PARshade, xx%LAIsunfrac)
          ! TODO: "grassland multilayer" model used leaf_flight = Flightsun, i.e.
          !       leaf_f_light(gc%f_lightfac, ll%met%PARsun) - which version is right?
          gp%leaf_f_light = leaf_f_light(gc%f_lightfac, V%met%PAR)
        else
          gp%f_light = 0.0
          gp%leaf_f_light = 0.0
        end if
      case default
        UNKNOWN_STRING(gc%f_light_method)
      end select

      ! Calculate f_temp
      select case (gc%f_temp_method)
      case ("disabled")
        ! Nothing to do
      case ("default")
        gp%f_temp = f_temp(V%met%Ts_C, gc%T_min, gc%T_opt, gc%T_max, gc%fmin)
      case ("square high")
        gp%f_temp = f_temp_square_high(V%met%Ts_C, gc%T_min, gc%T_opt, gc%T_max, gc%fmin)
      case default
        UNKNOWN_STRING(gc%f_temp_method)
      end select

      ! Calculate f_VPD
      select case (gc%f_VPD_method)
      case ("disabled")
        ! Nothing to do
      case ("linear")
        gp%f_VPD = f_VPD_linear(V%met%VPD, gc%VPD_max, gc%VPD_min, gc%fmin)
      case ("log")
        gp%f_VPD = f_VPD_log(V%met%VPD, gc%fmin)
      case default
        UNKNOWN_STRING(gc%f_VPD_method)
      end select

      ! Calculate f_SW
      select case (gc%f_SW_method)
      case ("disabled")
        ! Nothing to do
      case ("fSWP exp")
        gp%f_SW = f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, gc%fmin, V%SMD%SWP)
      case ("fSWP linear")
        gp%f_SW = f_SWP_linear(gc%SWP_min, gc%SWP_max, gc%fmin, V%SMD%SWP)
      ! TODO: implement LWP
      !case ("fLWP exp")
      !  gp%f_SW = f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, gc%fmin, V%SMD%LWP)
      case ("fPAW")
        gp%f_SW = f_PAW(this%SMD_conf%ASW_FC, gc%fmin, V%SMD%ASW)
      case default
        UNKNOWN_STRING(gc%f_SW_method)
      end select

      ! Calculate f_O3
      select case (gc%f_O3_method)
      case ("disabled")
        ! Nothing to do
      case ("wheat")
        gp%f_O3 = ((1+(xx%POD_0/11.5)**10)**(-1))
      case ("potato")
        gp%f_O3 = ((1+(xx%AOT_0/40)**5)**(-1))
      case default
        UNKNOWN_STRING(gc%f_O3_method)
      end select

      ! Photosynthetic gsto: V_cmax_25 and J_max_25
      select case (pgc%V_J_method)
      case ("input")
        ! Nothing to do
      case ("constant")
        xx%V_cmax_25 = pgc%V_cmax_25
        xx%J_max_25 = pgc%J_max_25
      case default
        UNKNOWN_STRING(pgc%V_J_method)
      end select

      ! Photosynthetic gsto: O3 effect on V_cmax_25
      select case (pgc%O3_method)
      case ("disabled")
        ! Nothing to do
      case ("martin2000")
        ! Percentage reduction in V_cmax (converting FO3_eff from nmol to mmol)
        delta_V_cmax = pgc%K_z * (xx%FO3_eff / 1e6)
        ! Convert to multiplier
        mult_V_cmax = 1.0 - (delta_V_cmax / 100)
        ! Reduce V_cmax and J_max
        xx%V_cmax_25 = max(0.1 * xx%V_cmax_25, min(xx%V_cmax_25 * mult_V_cmax, xx%V_cmax_25))
        xx%J_max_25 = max(0.1 * xx%J_max_25, min(xx%J_max_25 * mult_V_cmax, xx%J_max_25))
      case default
        UNKNOWN_STRING(pgc%V_J_method)
      end select

      ! Photosynthetic gsto: phenology effect
      select case (pgc%phenology_method)
      case ("disabled")
        ! Nothing to do
      case ("leaf_f_phen")
        xx%V_cmax_25 = xx%V_cmax_25 * gp%leaf_f_phen
        xx%J_max_25 = xx%J_max_25 * gp%leaf_f_phen
      case default
        UNKNOWN_STRING(pgc%phenology_method)
      end select
    end associate
  end subroutine calc_gsto_parameters

  !> Calculate stomatal conductances (leaf, mean and bulk) for a specified
  !! layer and land cover.  Must be run after calc_gsto_parameters.
  subroutine calc_gsto(this, loc, LC, V, ll, xx)
    class(DO3SE_State_t), intent(inout) :: this
    type(Location_t), intent(in) :: loc
    type(LandCover_t), intent(in) :: LC
    type(V_t), intent(in) :: V
    type(ML_t), intent(inout) :: ll
    type(MLMC_t), intent(inout) :: xx

    real :: D_0

    associate (gc => LC%gsto, &
               pgc => LC%pn_gsto, &
               gp => xx%gsto_params)
      ! Calculate gsto
      select case (gc%method)
      case ("multiplicative")
        xx%leaf_gsto = apply_VPD_crit(gc%VPD_crit, V%VPD_dd, xx%leaf_gsto, gsto_leaf(gp))
        xx%mean_gsto = apply_VPD_crit(gc%VPD_crit, V%VPD_dd, xx%mean_gsto, gsto_mean(gp))
      case ("photosynthesis")
        select case (pgc%D_0_method)
        case ("constant")
          D_0 = pgc%D_0
        case ("f_VPD")
          select case (gc%f_VPD_method)
          case ("linear")
            D_0 = inverse_f_VPD_linear(0.5, gc%VPD_max, gc%VPD_min, gc%fmin)
          case ("log")
            D_0 = inverse_f_VPD_log(0.5, gc%fmin)
          end select
        end select

        ! TODO: make this multi-layer aware
        ! TODO: Tleaf_C isn't multi-component aware, and will get overwritten
        !       repeatedly!  Fixing this should be able to make `ll` intent(in)
        call gsto_pn(LC%pn_gsto, xx%V_cmax_25, xx%J_max_25, D_0, &
                     LC%Lm, V%met%Ts_C, ll%met%u, V%met%CO2, V%met%PPFD, &
                     V%met%Rn, V%met%R, loc%albedo, V%met%P, V%met%eact, &
                     ll%met%Tleaf_C, xx%A_n, xx%leaf_gsto, xx%g_sv, xx%g_bv)
        xx%mean_gsto = xx%leaf_gsto
      case default
        UNKNOWN_STRING(gc%method)
      end select
      ! Scale mean gsto up to bulk gsto
      xx%bulk_gsto = xx%mean_gsto * xx%LAI
    end associate
  end subroutine

  !> Calculate multi-layer ozone resistance model, deposition velocity, and
  !! multi-layer ozone concentration.
  subroutine calc_ozone_deposition(this)
    class(DO3SE_State_t), intent(inout) :: this

    integer :: iL

    ! Calculate resistance model for O3 over the target canopy
    call init_ResistanceModel(this%V%rmodel_O3, this%nL)
    this%V%rmodel_O3%Ra_c = Ra_simple(this%V%met%ustar, (this%V%canopy_height * (CANOPY_D + CANOPY_Z0)), &
                               50.0, this%V%canopy_height * CANOPY_D)
    this%V%rmodel_O3%Ra = Ra_simple(this%V%met%ustar, this%V%canopy_height, 50.0, this%V%canopy_height * CANOPY_D)
    this%V%rmodel_O3%Rb = Rb(this%V%met%ustar, DIFF_O3)
    do iL = 1, this%nL
      this%V%rmodel_O3%Rinc(iL) = Rinc(sum(this%MLMC(iL,:)%SAI), this%V%canopy_height, this%V%met%ustar)
      !this%V%rmodel_O3%Rinc(iL) = Rinc_prototype(sum(this%MLMC(iL,:)SAI), this%V%met%ustar)
      this%V%rmodel_O3%Rext(iL) = Rext(sum(this%MLMC(iL,:)%SAI))
      this%V%rmodel_O3%Rsto(iL) = Rsto(sum(this%MLMC(iL,:)%bulk_gsto))
    end do
    this%V%rmodel_O3%Rgs = this%location%Rsoil
    this%V%rmodel_O3%Rsur(1:this%nL) = Rsur(this%nL, &
                                            this%V%rmodel_O3%Rb, &
                                            this%V%rmodel_O3%Rsto(1:this%nL), &
                                            this%V%rmodel_O3%Rext(1:this%nL), &
                                            sum(this%MLMC(:,:)%LAI, dim=2), &
                                            sum(this%MLMC(:,:)%SAI, dim=2))
    this%V%rmodel_O3%Rtotal(1:this%nL) = Rtotal(this%nL, &
                                                this%V%rmodel_O3%Rsur(1:this%nL), &
                                                this%V%rmodel_O3%Rinc(1:this%nL), &
                                                this%V%rmodel_O3%Rgs)

    ! Vd calculation duplicated here just for comparison purposes
    ! TODO: better way to expose Vd
    this%V%Vd = deposition_velocity(this%V%rmodel_O3)

    ! Calculate O3 concentration at 50m and canopy height
    call met_O3(this%met_conf, this%V%met, this%ML(1)%met, this%location, this%V%rmodel_O3, this%V%canopy_height)
    ! Calculate per-layer O3 concentrations
    if (this%location%OTC) then
      ! In OTC, assume uniform O3 driven by external factors.
      this%ML(:)%met%O3 = this%ML(1)%met%O3
    else
      ! Use multi-layer resistance model for per-layer O3 concentration.
      call multi_layer_O3(this%V%rmodel_O3, this%V%met%O3_50, this%ML(:)%met%O3)
    end if
  end subroutine calc_ozone_deposition

  !> Calculate the leaf-level ozone resistance model per layer and land cover,
  !! and calculate the various measures of ozone dose.
  subroutine calc_ozone_dose(this)
    class(DO3SE_State_t), intent(inout) :: this

    integer :: iL, iLC
    real :: O3_ppb_to_nmol

    ! TODO: move this to its own subroutine?
    this%MLMC(:,:)%leaf_rmodel_O3 = LeafResistanceModel_t()
    do iL = 1, this%nL
      do iLC = 1, this%nLC
      associate (LC => this%LCs(iLC), ll => this%ML(iL), xx => this%MLMC(iL,iLC))
        xx%leaf_rmodel_O3%Rb = leaf_rb(leaf_gb(LEAF_G_O3, LC%Lm, ll%met%u))
      end associate
      end do
    end do
    this%MLMC(:,:)%leaf_rmodel_O3%Rext = Rext(1.0)
    this%MLMC(:,:)%leaf_rmodel_O3%Rsto = Rsto(this%MLMC(:,:)%leaf_gsto)

    O3_ppb_to_nmol = O3_ppb_to_nmol_factor(this%V%met%Ts_C, this%V%met%P)
    do iL = 1, this%nL
      do iLC = 1, this%nLC
      associate (LC => this%LCs(iLC), xx => this%MLMC(iL,iLC))
        ! TODO: fix reliance on MAX_RSTO...
        if (xx%leaf_gsto > 0) then
          xx%Fst = this%ML(iL)%met%O3 * O3_ppb_to_nmol * stomatal_flux_rate(xx%leaf_rmodel_O3)
        else
          xx%Fst = 0
        end if

        xx%POD_0 = xx%POD_0 + ((xx%Fst*DT)/1000000)
        xx%POD_Y = xx%POD_Y + ((max(0.0, xx%Fst - LC%Y)*DT)/1000000)

        ! Default OT0/OT40 to 0
        xx%OT_0 = 0
        xx%OT_40 = 0

        ! Only accumulate OT when global radiation > 50 W m-2
        if (this%is_daylight()) then
          ! Only accumulate OT0 when leaf_fphen > 0
          if (xx%gsto_params%leaf_f_phen > 0) then
            xx%OT_0 = this%ML(iL)%met%O3 / 1000
          end if

          ! Only accumulate OT40 when fphen > 0
          if (xx%gsto_params%f_phen > 0) then
            xx%OT_40 = max(0.0, this%ML(iL)%met%O3 - 40) / 1000
          end if
        end if

        ! Accumulate OT0/OT40
        xx%AOT_0 = xx%AOT_0 + xx%OT_0
        xx%AOT_40 = xx%AOT_40 + xx%OT_40

        if (LC%pn_gsto%O3_method == "martin2000") then
          ! Effective ozone dose for effect on photosynthesis
          xx%FO3_eff = xx%FO3_eff + (xx%Fst - LC%pn_gsto%F_0) * DT
          xx%FO3_eff = max(0.0, xx%FO3_eff)
        end if
      end associate
      end do
    end do
  end subroutine calc_ozone_dose

  !> Calculate SMD values from whichever soil moisture input is being used.
  subroutine calc_soil_moisture(this)
    class(DO3SE_State_t), intent(inout) :: this

    select case (this%SMD_conf%source)
    case ("disabled")
      ! Nothing to do
    case ("input SWP")
      this%V%SMD = soil_moisture_from_SWP(this%SMD_conf, this%V%SMD%SWP)
    case ("input SWC")
      this%V%SMD = soil_moisture_from_SWC(this%SMD_conf, this%V%SMD%Sn)
    case ("P-M")
      if (this%start_of_day()) then
        ! Apply change in soil moisture from previous day's P-M summary.
        this%V%SMD = soil_moisture_from_SWC(this%SMD_conf, this%V%SMD%Sn + this%V%PM%Sn_diff)
      end if
    case default
      UNKNOWN_STRING(this%SMD_conf%source)
    end select
  end subroutine calc_soil_moisture

  !> Calculate inputs to soil moisture, if necessary.  Should probably be called
  !! as the very last thing in the model loop.  Currently only implements the
  !! P-M water vapour stuff.
  subroutine calc_soil_moisture_changes(this)
    class(DO3SE_State_t), intent(inout) :: this

    select case (this%SMD_conf%source)
    case ("P-M")
      !
      ! Do hourly Penman-Monteith
      !

      ! Reset values at beginning of day
      if (this%start_of_day()) then
        call penman_monteith_reset(this%V%PM)
      end if

      ! Adapt multi-layer O3 resistance model to single-layer H2O resistance
      call init_ResistanceModel(this%V%rmodel_H2O, 1)
      this%V%rmodel_H2O%Ra = this%V%rmodel_O3%Ra
      this%V%rmodel_H2O%Rb = Rb(this%V%met%ustar, DIFF_H2O)
      ! Combine multi-layer into single-layer
      this%V%rmodel_H2O%Rinc(1) = 1.0/sum(1.0/this%V%rmodel_O3%Rinc(1:this%V%rmodel_O3%nL))
      this%V%rmodel_H2O%Rext(1) = 1.0/sum(1.0/this%V%rmodel_O3%Rext(1:this%V%rmodel_O3%nL))
      this%V%rmodel_H2O%Rsto(1) = DRATIO * 1.0/sum(1.0/this%V%rmodel_O3%Rsto(1:this%V%rmodel_O3%nL))
      this%V%rmodel_H2O%Rgs = this%V%rmodel_O3%Rgs
      this%V%rmodel_H2O%Rsur(1:this%nL) = Rsur(this%nL, &
                                               this%V%rmodel_H2O%Rb, &
                                               this%V%rmodel_H2O%Rsto(1:1), &
                                               this%V%rmodel_H2O%Rext(1:1), &
                                               [sum(this%MLMC(:,:)%LAI)], &
                                               [sum(this%MLMC(:,:)%SAI)])
      this%V%rmodel_O3%Rtotal(1:this%nL) = Rtotal(this%nL, &
                                                  this%V%rmodel_O3%Rsur(1:1), &
                                                  this%V%rmodel_O3%Rinc(1:1), &
                                                  this%V%rmodel_O3%Rgs)
      ! Is soil evaporation blocked?
      ! TODO: this assumes the first land cover is the only one that matters
      associate (gc => this%LCs(1)%gsto)
        select case (gc%f_SW_method)
        case ("fSWP exp", "fLWP exp")
          this%V%Es_blocked = this%V%SMD%SWP < inverse_f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, 1.0)
        case ("fSWP linear")
          this%V%Es_blocked = this%V%SMD%SWP < gc%SWP_max
        case ("fPAW")
          this%V%Es_blocked = this%V%SMD%ASW < inverse_f_PAW(this%SMD_conf%ASW_FC, gc%fmin, 1.0)
        case default
          this%V%Es_blocked = .true.
        end select
      end associate

      ! Run hourly accumulations of precipitation and evapotranspiration
      this%V%PM%precip_acc = this%V%PM%precip_acc + this%V%met%precip/1000
      call penman_monteith_hourly(this%V%met%Rn*1000000, this%V%met%P*1000, this%V%met%Ts_C, &
                                  this%V%met%esat*1000, this%V%met%eact*1000, this%V%met%VPD*1000, this%V%rmodel_H2O, &
                                  sum(this%MLMC(:,:)%LAI), logical(this%V%Es_blocked), this%V%PM)

      ! Calculate daily balance
      if (this%end_of_day()) then
        call penman_monteith_daily(this%V%PM, sum(this%MLMC(:,:)%LAI), this%SMD_conf%root, &
                                   this%SMD_conf%run_off_fraction, this%V%SMD%ASW, this%V%SMD%SMD)
      end if
    end select
  end subroutine calc_soil_moisture_changes

  !> Run everything that should be done in an hourly timestep.  Hourly nputs
  !! should be read before calling this, and hourly outputs written after.
  subroutine run_hour(this)
    class(DO3SE_State_t), intent(inout) :: this

    integer :: iL, iLC

    ! Ensure we have canopy height, LAI and SAI values
    call this%calc_phenology()

    ! Large-scale met data
    call this%calc_met_data()
    ! Canopy-dependent and per-layer met data
    call this%calc_micromet_data()

    ! Calculate soil moisture data (from input or using P-M)
    call this%calc_soil_moisture()

    ! Calculate gsto
    do iL = 1, this%nL
      do iLC = 1, this%nLC
        call this%calc_gsto_parameters(this%LCs(iLC), this%V, this%ML(iL), this%MLMC(iL,iLC))
        call this%calc_gsto(this%location, this%LCs(iLC), this%V, this%ML(iL), this%MLMC(iL,iLC))
      end do
    end do

    ! Canopy-level deposition (Vd, calculate O3 concentrations at layers)
    call this%calc_ozone_deposition()
    ! Leaf-level flux and ozone dose (FSt, POD, OT40, etc.)
    call this%calc_ozone_dose()


    ! TODO: work out how to combine the MLMC O3 flux values down to whole-canopy,
    !       whole-species and whole-layer values.

    ! Calculate soil moisture inputs from this hour for the next hour
    call this%calc_soil_moisture_changes()
  end subroutine run_hour

  !> Read non-land-cover configuration.
  logical function read_DO3SE_Config(config_file, location, LCC, met, SMD)
    character(len=*), intent(in) :: config_file
    type(Location_t), intent(out) :: location
    type(LandCoverConfig_t), target, intent(out) :: LCC
    type(MetConfig_t), intent(out) :: met
    type(SMDConfig_t), intent(out) :: SMD

    integer :: config_unit
    type(LandCover_t), pointer :: LC

    namelist /DO3SE_Config/ location, LCC, LC, met, SMD

    location = Location_t()
    LCC = LandCoverConfig_t()
    LC => LCC%LCs(1)
    met = MetConfig_t()
    SMD = SMDConfig_t()

    open (newunit=config_unit, file=config_file, status="old", action="read", position="rewind")
    read (unit=config_unit, nml=DO3SE_Config, end=100)
    read_DO3SE_Config = .true.
    close (config_unit)
    return
100 read_DO3SE_Config = .false.
    close (config_unit)
    return
  end function read_DO3SE_Config

end module DO3SE_ml
