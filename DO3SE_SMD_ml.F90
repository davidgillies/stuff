module DO3SE_SMD_ml

  use DO3SE_ModelConstants_ml, only: UNDEF
  use DO3SE_PhysicalConstants_ml, only: T0, DRATIO
  use DO3SE_Resistance_ml, only: ResistanceModel_t
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"
#include "interop_types.h"

  implicit none
  private

  !> SMD soil texture parameters.
  TYPE, public :: Soil_t
    REAL :: b = UNDEF         !< Texture dependent soil conductivity parameter
    REAL :: FC = UNDEF        !< Field capacity (m3 m-3)
    REAL :: SWP_AE = UNDEF    !< Water potential at air entry (MPa)
    REAL :: Ksat = UNDEF      !< Saturated soil conductance (s-2 MPa-1 mm-1)
  end type Soil_t

  ! Commonly used soil textures:
  !> Sandy loam
  type(Soil_t), parameter, public :: SOIL_SANDY_LOAM = &
    Soil_t(b = 3.31, FC = 0.16, SWP_AE = -0.00091, Ksat = 0.0009576)
  !> Silt loam
  type(Soil_t), parameter, public :: SOIL_SILT_LOAM = &
    Soil_t(b = 4.38, FC = 0.26, SWP_AE = -0.00158, Ksat = 0.0002178)
  !> Loam
  type(Soil_t), parameter, public :: SOIL_LOAM = &
    Soil_t(b = 6.58, FC = 0.29, SWP_AE = -0.00188, Ksat = 0.0002286)
  !> Clay loam (Ksat estimated)
  type(Soil_t), parameter, public :: SOIL_CLAY_LOAM = &
    Soil_t(b = 7.00, FC = 0.37, SWP_AE = -0.00588, Ksat = 0.00016)

  TYPE, public :: SMDConfig_t
    !> Soil texture:
    !!    - "sandy loam"
    !!    - "silt loam"
    !!    - "loam"
    !!    - "clay loam"
    !!    - "" or "custom":   Soil texture parameters set individually
    CHARACTER(len=16) :: soil_texture = "loam"
    !> Soil texture parameters
    type(Soil_t) :: soil = Soil_t()

    REAL :: root = 1.2        !< Root depth (m)
    ! TODO: seed PWP from SWP_min if available?
    REAL :: PWP = -4.0        !< "Permanent wilting point", minimum level for SWP (MPa)
    REAL :: ASW_FC = UNDEF    !< Available soil water at field capacity (m)
                              !! (calculated by check_SMDConfig)

    !> Soil water data source:
    !!    - "disabled":   No soil water data
    !!    - "input SWP":  Input soil water potential (MPa)
    !!    - "input SWC":  Input soil water volumetric content (m3 m-3)
    !!    - "P-M":        Use Penman-Monteith method to track soil water content
    CHARACTER(len=16) :: source = "disabled"
    ! Penman-Monteith method parameters
    REAL :: initial_SWC = UNDEF     !< Initial soil water content, defaults to soil%FC
    REAL :: run_off_fraction = 0.0  !< Irrigation inefficiency lost to run-off (fraction)

    !> (WEAP) "Management Allowable Depletion" (m), maximum SMD before irrigation is requested.
    REAL :: MAD = 0.0
  end type SMDConfig_t

  TYPE, public :: SMDData_t
    REAL :: Sn = UNDEF        !< Soil water content (m3 m-3)
    REAL :: SWP = UNDEF       !< Soil water potential (MPa)
    REAL :: ASW = UNDEF       !< Available soil water, above PWP (m)
    REAL :: SMD = UNDEF       !< Soil moisture deficit, below FC (m)
  end type SMDData_t

  TYPE, public :: PM_State_t
    ! Latest evapotranspiration values
    REAL :: Ei = 0.0          !< Evaporation from canopy (m)
    REAL :: Et = 0.0          !< Plant transpiration (m)
    REAL :: Es = 0.0          !< Soil surface evaporation (m)
    REAL :: Eat = 0.0         !< Evapotranspiration (m)

    ! Accumulated values (processed and cleared daily)
    REAL :: Ei_acc = 0.0      !< Accumulated evaporation from canopy (m)
    REAL :: Et_acc = 0.0      !< Accumulated plant transpiration (m)
    REAL :: Es_acc = 0.0      !< Accumulated soil surface evaporation (m)
    REAL :: Eat_acc = 0.0     !< Accumulated evapotranspiration (m)
    REAL :: precip_acc = 0.0  !< Accumulated precipitation (m)

    ! Soil water content tracking
    REAL :: Sn_diff = 0.0     !< Latest change in Sn (m3 m-3)

    ! Water destination tracking (updated daily)
    REAL :: input = 0.0           !< Input from rainfall + irrigation (m)
    REAL :: run_off = 0.0         !< Loss to run-off (m)
    REAL :: effective = 0.0       !< Effective irrigation (m)
    REAL :: intercepted_evaporated = 0.0  !< Loss to evaporation of intercepted (m)
    REAL :: evapotranspiration = 0.0      !< Loss to evapotranspiration (m)
    REAL :: percolated = 0.0      !< Loss to deep percolation (m)

    REAL :: run_off_acc = 0.0     !< Accumulated run-off (m)
    REAL :: percolated_acc = 0.0  !< Accumulated deep percolation (m)
  end type PM_State_t

  public :: check_SMDConfig
  public :: soil_moisture_from_SWC
  public :: soil_moisture_from_SWP
  public :: penman_monteith_hourly
  public :: penman_monteith_daily
  public :: penman_monteith_reset

contains

  !> Check that a Soil_t is valid.  The optional *name* argument can be used to
  !! select one of the preset values, passing name="" or name="custom" has the
  !! same effect as omitting the argument.
  subroutine check_Soil(soil, name)
    type(Soil_t), intent(inout) :: soil
    character(len=*), intent(in), optional :: name

    if (present(name)) then
      select case (name)
      case ("sandy loam")
        soil = SOIL_SANDY_LOAM
      case ("silt loam")
        soil = SOIL_SILT_LOAM
      case ("loam")
        soil = SOIL_LOAM
      case ("clay loam")
        soil = SOIL_CLAY_LOAM
      case ("", "custom")
        ! Parameters should have been set manually
      case default
        ERROR("unknown soil texture: "//trim(name))
      end select
    end if

    ! Check that parameters were set by the name or manually
    ASSERT_DEFINED(soil%b)
    ASSERT_DEFINED(soil%FC)
    ASSERT_DEFINED(soil%SWP_AE)
    ASSERT_DEFINED(soil%Ksat)
  end subroutine check_Soil

  !> Check that a SMDConfig_t is valid.
  !!
  !!    - Checks that the Soil_t is valid
  !!    - Applies named soil texture if necessary
  !!    - Calculates various soil-dependent constants
  subroutine check_SMDConfig(config)
    type(SMDConfig_t), intent(inout) :: config

    type(SMDData_t) :: tmp_data

    ! Check/set soil texture parameters
    call check_Soil(config%soil, config%soil_texture)

    ! Calculate available soil water at field capacity
    tmp_data = soil_moisture_from_SWC(config, config%soil%FC)
    config%ASW_FC = tmp_data%ASW

    ! Set up soil water data source parameters
    if (config%source == "P-M") then
      if (.not. is_def(config%initial_SWC)) then
        config%initial_SWC = config%soil%FC
      end if
    end if
  end subroutine check_SMDConfig

  !> Use soil water release curve to convert from soil water content (m3 m-3) to
  !! soil water potential (MPa).
  pure real function SWC_to_SWP(SWP_AE, b, SWC)
    real, intent(in) :: SWP_AE    !< Water potential at air entry (MPa)
    real, intent(in) :: b         !< Texture dependent soil conductivity parameter
    real, intent(in) :: SWC       !< Soil water content (m3 m-3)

    real, parameter :: SWC_sat = 0.4 ! Saturated soil water content for soil water release curve

    SWC_to_SWP = SWP_AE * ((SWC_sat / SWC)**b)
  end function SWC_to_SWP

  !> Use soil water release curve to convert from soil water potential (MPa) to
  !! soil water content (m3 m-3).
  pure real function SWP_to_SWC(SWP_AE, b, SWP)
    real, intent(in) :: SWP_AE    !< Water potential at air entry (MPa)
    real, intent(in) :: b         !< Texture dependent soil conductivity parameter
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    real, parameter :: SWC_sat = 0.4 ! Saturated soil water content for soil water release curve

    SWP_to_SWC = 1.0 / (((SWP/SWP_AE)**(1.0/b)) / SWC_sat)
  end function SWP_to_SWC

  !> Fill soil moisture data from a soil water content value.
  pure function soil_moisture_from_SWC(config, Sn) result(data)
    type(SMDConfig_t), intent(in) :: config
    real, intent(in) :: Sn      !< Soil water content (m3 m-3)

    type(SMDData_t) :: data

    real :: PWP_vol

    associate (soil => config%soil)
      ! Convert PWP to volumetric content to use as a minimum soil water content
      PWP_vol = SWP_to_SWC(soil%SWP_AE, soil%b, config%PWP)

      ! Constrain soil water content to be between field capacity and PWP
      data%Sn = max(PWP_vol, min(soil%FC, Sn))

      ! Calculate soil water potential (SWP)
      data%SWP = SWC_to_SWP(soil%SWP_AE, soil%b, data%Sn)

      ! Calculate available soil water (ASW)
      data%ASW = (data%Sn - PWP_vol) * config%root

      ! Calculate soil moisture deficit (SMD)
      data%SMD = (soil%FC - data%Sn) * config%root
    end associate
  end function soil_moisture_from_SWC

  !> Fill soil moisture data from a soil water potential value.
  pure function soil_moisture_from_SWP(config, SWP) result(data)
    type(SMDConfig_t), intent(in) :: config
    real, intent(in) :: SWP     !< Soil water potential (MPa)

    type(SMDData_t) :: data

    associate (soil => config%soil)
      data = soil_moisture_from_SWC(config, SWP_to_SWC(soil%SWP_AE, soil%b, SWP))
    end associate
  end function soil_moisture_from_SWP


  !> Hourly Penman-Monteith calculations for evaporation and transpiration.  
  !! Fills in values and updates accumulators in *state* for Ei, Et, Es and Eat.
  !!
  !! TODO: should we still be using Rsoil, instead of rm\%Rgs?
  subroutine penman_monteith_hourly(Rn, P, Ts_C, esat, eact, VPD, rm, LAI, Es_blocked, state)
    real, intent(in) :: Rn              !< Net radiation (J)
    real, intent(in) :: P               !< Atmospheric pressure (Pa)
    real, intent(in) :: Ts_C            !< Air temperature (degrees C)
    real, intent(in) :: esat            !< Saturated vapour pressure (Pa)
    real, intent(in) :: eact            !< Actual vapour pressure (Pa)
    real, intent(in) :: VPD             !< Vapour pressure deficit (Pa)
    type(ResistanceModel_t), intent(in) :: rm
    real, intent(in) :: LAI             !< Leaf area index (m2 m-2)
    logical, intent(in) :: Es_blocked   !< Is soil evaporation blocked?

    type(PM_State_t), intent(inout) :: state

    real :: Tvir, delta, lambda, psychro, Pair, Cair, G

    real :: Et_1, Et_2, Ei_3, Et_3
    real :: t, Es_Rn, Es_G, Es_1, Es_2, Es_3
    real :: SW_a, SW_s, SW_c, C_canopy, C_soil

    ! This model (probably) makes some one-layer assumptions, so don't allow 
    ! multi-layer resistance model.
    ASSERT(rm%nL == 1)


    associate (Ra => rm%Ra, &
               Rb_H2O => rm%Rb, &
               Rinc => rm%Rinc(1), &
               Rsto_H2O => rm%Rsto(1), &
               Rsoil => rm%Rgs, &
               Ei => state%Ei, &
               Et => state%Et, &
               Es => state%Es, &
               Eat => state%Eat)

      Tvir = (Ts_c+T0)/(1-(0.378*(eact/P)))
      delta= ((4098*esat)/((Ts_C+237.3)**2))
      lambda = (2501000-(2361*Ts_C))
      psychro = 1628.6 * (P/lambda)
      Pair = (0.003486*(P/Tvir))
      Cair = (0.622*((lambda*psychro)/P))

      G = 0.1 * Rn

      Et_1 = (delta * (Rn - G)) / lambda
      Et_2 = 3600 * Pair * Cair * VPD / Rb_H2O / lambda

      Ei_3 = delta + psychro
      Ei = (Et_1 + Et_2) / Ei_3 / 1000

      Et_3 = delta + psychro * (1 + Rsto_H2O / Rb_H2O)
      Et = (Et_1 + Et_2) / Et_3 / 1000

      if (Es_blocked) then
        Es = 0
      else
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lambda
        Es_2 = ((3600 * Pair * Cair * VPD) - (delta * Rinc * ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc + Rb_H2O) / lambda
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es = (Es_1 + Es_2) / Es_3 / 1000
      end if

      ! Calculate Eat from Et and Es (after Shuttleworth and Wallace, 1985)
      SW_a = (delta + psychro) * Rb_H2O
      SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
      SW_c = psychro * Rsto_H2O  ! Boundary layer resistance = 0
      C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
      C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
      if (Es <= 0) then
        Eat = Et
      else
        Eat = (C_canopy * Et) + (C_soil * Es)
      end if

    end associate

    ! Accumulate values
    state%Ei_acc = state%Ei_acc + state%Ei
    state%Et_acc = state%Et_acc + state%Et
    state%Es_acc = state%Es_acc + state%Es
    state%Eat_acc = state%Eat_acc + state%Eat
  end subroutine penman_monteith_hourly

  !> Daily soil water content update from accumulated Penman-Monteith values.  
  !! Sets state\%Sn_diff (change in soil water content) and clears daily
  !! accumulators.  state\%Sn_diff isn't constrained by field capacity here, so
  !! that consideration must be handled elsewhere.
  !!
  !! Also processes run-off and deep percolation according to SMD configuration.
  !! Both "most recent" and "accumulated" values for these are stored in state.
  pure subroutine penman_monteith_daily(state, LAI, root, run_off_fraction, ASW, SMD)
    type(PM_State_t), intent(inout) :: state
    real, intent(in) :: LAI     !< Leaf area index (m2 m-2)
    real, intent(in) :: root    !< Root depth (m)
    real, intent(in) :: run_off_fraction    !< Amount of precipitation that is lost as run-off
    real, intent(in) :: ASW     !< Current available soil water, above PWP (m)
    real, intent(in) :: SMD     !< Current soil moisture deficit, from FC (m)

    real :: max_ET, delta_SM

    ! Start with full amount of precipitation
    state%input = state%precip_acc

    ! Estimate loss to run-off
    state%run_off = run_off_fraction * state%input
    state%run_off_acc = state%run_off_acc + state%run_off
    ! Calculate "effective irrigation"
    state%effective = state%input - state%run_off

    ! Estimate loss of intercepted precipitation to evaporation.  Intercepted
    ! precipitation is estimated as 0.0001*LAI, which is therefore a limit on
    ! how much can be evaporated.
    state%intercepted_evaporated = min(state%effective, 0.0001 * LAI, state%Ei_acc)

    ! Can't lose water below PWP, so constrain evapotranspiration to ASW by
    ! restricting evapotranspiration.
    max_ET = ASW + state%effective - state%intercepted_evaporated
    state%evapotranspiration = min(max_ET, state%Eat_acc)

    ! Total balance = input - run_off - evaporated - evapotranspiration
    delta_SM = state%effective - state%intercepted_evaporated - state%evapotranspiration
    ! Converted to volumetric change using root depth.
    state%Sn_diff = delta_SM / root

    ! Amount that will go to deep percolation = remainder if water balance
    ! refills soil water, i.e. if it is greater than SMD.
    state%percolated = max(0.0, delta_SM - SMD)
    state%percolated_acc = state%percolated_acc + state%percolated
  end subroutine penman_monteith_daily

  !> Reset daily accumulated Penman-Monteith values.
  pure subroutine penman_monteith_reset(state)
    type(PM_State_t), intent(inout) :: state

    ! Clear accumulated variables
    state%Ei_acc = 0.0
    state%Et_acc = 0.0
    state%Es_acc = 0.0
    state%Eat_acc = 0.0
    state%precip_acc = 0.0
  end subroutine penman_monteith_reset

end module DO3SE_SMD_ml
