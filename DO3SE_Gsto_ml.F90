module DO3SE_Gsto_ml

  use DO3SE_ConfigTypes_ml
  use DO3SE_Util_ml
  use DO3SE_PhysicalConstants_ml, only: PAR_Wm2_to_photons
  use DO3SE_ModelConstants_ml, only: UNDEF, ASW_MIN, ASW_MAX
#include "interop_types.h"

  implicit none
  private

  !> Multiplicative stomatal conductance parameters.
  TYPE, public :: GstoParams_t
    REAL :: fmin = UNDEF        !< Minimum stomatal conductance (fraction)
    REAL :: gmax = UNDEF        !< Maximum stomatal conductance (mmol O3 m-2 PLA s-1)
    REAL :: gmorph = 1.0        !< Sun/shade morphology factor (fraction)

    REAL :: f_phen = 1.0        !< Phenology-related effect on gsto (fraction)
    REAL :: leaf_f_phen = 1.0   !< Phenology-related effect on leaf gsto (fraction)
    REAL :: f_light = 1.0       !< Irradiance effect on gsto (fraction)
    REAL :: leaf_f_light = 1.0  !< Irradiance effect on leaf gsto (fraction)
    REAL :: f_temp = 1.0        !< Temperature effect on gsto (fraction)
    REAL :: f_VPD = 1.0         !< VPD effect on gsto (fraction)
    REAL :: f_SW = 1.0          !< Soil water effect on gsto (fraction)
    REAL :: f_O3 = 1.0          !< O3 effect on gsto (fraction)
  end type GstoParams_t

  public :: f_phen_simple_PLF
  public :: f_phen_complex_PLF
  public :: leaf_f_phen_PLF
  public :: f_light
  public :: leaf_f_light
  public :: f_temp
  public :: f_temp_square_high
  public :: f_VPD_linear
  public :: inverse_f_VPD_linear
  public :: f_VPD_log
  public :: inverse_f_VPD_log
  public :: f_SWP_exp
  public :: inverse_f_SWP_exp
  public :: f_SWP_linear
  public :: f_PAW
  public :: inverse_f_PAW
  public :: gsto_leaf
  public :: gsto_mean
  public :: apply_VPD_crit

contains

  !> Phenology effect on stomatal conductance, according to a simple piecewise
  !! linear function.
  !!
  !! To handle all situations, including winter growing seasons, everything is
  !! re-indexed to be relative to SGS=0.
  !!
  !! ~~~
  !!
  !!                 f_phen_c
  !!              ______________
  !!             /              \\
  !!            /                \\
  !! f_phen_a  /                  \\
  !!          |<1>              <4>\  f_phen_e
  !!          |                     |
  !!   _______|                     |________
  !!         SGS                   EGS
  !! ~~~
  real function f_phen_simple_PLF(gc, SGS, EGS, dd)
    type(GstoConfig_t), intent(in) :: gc  !< gsto parameters
    integer, intent(in) :: SGS            !< Start of growing season (day of year)
    integer, intent(in) :: EGS            !< End of growing season (day of year)
    integer, intent(in) :: dd             !< Day of year

    real, dimension(2, 6) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      SGS, 0.0, &
      SGS, gc%f_phen_a, &
      (SGS + gc%f_phen_1), gc%f_phen_c, &
      (EGS - gc%f_phen_4), gc%f_phen_c, &
      EGS, gc%f_phen_e, &
      EGS, 0.0 /), shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    call reindex(func(1,:), real(SGS), 365.0)
    call assert(all(func(1,1:5) <= func(1,2:6)), "f_phen_simple_PLF: points not in order")
    dd_adj = dd - SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    f_phen_simple_PLF = PLF_value(func, real(dd_adj))
  end function f_phen_simple_PLF

  !> Phenology effect on stomatal conductance, according to a more complicated
  !! piecewise linear function.
  !!
  !! To handle all situations, including winter growing seasons, everything is
  !! re-indexed to be relative to SGS=0.
  !!
  !! ~~~
  !!                 f_phen_b                  f_phen_d
  !!                __________                __________
  !!               /         :\              /:         \\
  !!              /          : \  f_phen_c  / :          \\
  !!             /           :  \__________/  :           \\
  !!            /            :<2>          <3>:            \\
  !! f_phen_a  /             :                :             \\
  !!          | <1>          :                :          <4> \  f_phen_e
  !!          |              :                :               |
  !!   _______|              :                :               |________
  !!         SGS            limA            limB             EGS
  !! ~~~
  real function f_phen_complex_PLF(gc, SGS, EGS, dd)
    type(GstoConfig_t), intent(in) :: gc  !< gsto parameters
    integer, intent(in) :: SGS            !< Start of growing season (day of year)
    integer, intent(in) :: EGS            !< End of growing season (day of year)
    integer, intent(in) :: dd             !< Day of year

    real, dimension(2, 10) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      SGS, 0.0, &
      SGS, gc%f_phen_a, &
      (SGS + gc%f_phen_1), gc%f_phen_b, &
      gc%f_phen_limA, gc%f_phen_b, &
      (gc%f_phen_limA + gc%f_phen_2), gc%f_phen_c, &
      (gc%f_phen_limB - gc%f_phen_3), gc%f_phen_c, &
      gc%f_phen_limB, gc%f_phen_d, &
      (EGS - gc%f_phen_4), gc%f_phen_d, &
      EGS, gc%f_phen_e, &
      EGS, 0.0 /), shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    call reindex(func(1,:), real(SGS), 365.0)
    call assert(all(func(1,1:9) <= func(1,2:10)), "f_phen_complex_PLF: points not in order")
    dd_adj = dd - SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    f_phen_complex_PLF = PLF_value(func, real(dd_adj))
  end function f_phen_complex_PLF


  !> Phenology effect on leaf stomatal conductance, according to a simple
  !! piecewise linear function.
  !!
  !! To handle all situations, including winter growing seasons, everything is
  !! re-indexed to be relative to Astart=0.
  !!
  !! ~~~
  !!                    b
  !!              ______________
  !!             /              \\
  !!            /                \\
  !!        a  /                  \\
  !!          |<1>              <2>\  c
  !!          |                     |
  !!   _______|                     |________
  !!       Astart                  Aend
  !! ~~~
  real function leaf_f_phen_PLF(gc, Astart, Aend, dd)
    type(GstoConfig_t), intent(in) :: gc  !< gsto parameters
    integer, intent(in) :: Astart         !< Start of accumulation period (day of year)
    integer, intent(in) :: Aend           !< End of accumulation period (day of year)
    integer, intent(in) :: dd             !< Day of year

    real, dimension(2, 6) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      Astart, 0.0, &
      Astart, gc%leaf_f_phen_a, &
      (Astart + gc%leaf_f_phen_1), gc%leaf_f_phen_b, &
      (Aend - gc%leaf_f_phen_2), gc%leaf_f_phen_b, &
      Aend, gc%leaf_f_phen_c, &
      Aend, 0.0 /), shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before Astart to the end of the year
    call reindex(func(1,:), real(Astart), 365.0)
    call assert(all(func(1,1:5) <= func(1,2:6)), "leaf_f_phen_PLF: points not in order")
    dd_adj = dd - Astart
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    leaf_f_phen_PLF = PLF_value(func, real(dd_adj))
  end function leaf_f_phen_PLF


  !> Temperature effect on stomatal conductance.
  !!
  !! Describes a sort of bell-shaped curve with minima at T_min and T_max and a
  !! maximum at T_opt.
  pure real function f_temp(Ts_C, T_min, T_opt, T_max, fmin)
    real, intent(in) :: Ts_C    !< Air temperature (degrees C)
    real, intent(in) :: T_min   !< Minimum temperature (degrees C)
    real, intent(in) :: T_opt   !< Optimum temperature (degrees C)
    real, intent(in) :: T_max   !< Maximum temperature (degrees C)
    real, intent(in) :: fmin    !< Minimum f_temp

    real :: bt

    bt = (T_max - T_opt) / (T_opt - T_min)
    f_temp = ((Ts_C - T_min) / (T_opt - T_min)) * ((T_max - Ts_C) / (T_max - T_opt))**bt
    f_temp = max(fmin, min(1.0, f_temp))
  end function f_temp

  !> Temperature effect on stomatal conductance without reduction after T_opt.
  !!
  !! Follows f_temp exactly before T_opt, but the high end of the function is
  !! square - it continues at 1.0 until dropping to 0.0 at T_max.
  pure real function f_temp_square_high(Ts_C, T_min, T_opt, T_max, fmin)
    real, intent(in) :: Ts_C    !< Air temperature (degrees C)
    real, intent(in) :: T_min   !< Minimum temperature (degrees C)
    real, intent(in) :: T_opt   !< Optimum temperature (degrees C)
    real, intent(in) :: T_max   !< Maximum temperature (degrees C)
    real, intent(in) :: fmin    !< Minimum f_temp

    if (Ts_C >= T_max) then
      f_temp_square_high = 0.0
    else if (Ts_C >= T_opt) then
      f_temp_square_high = 1.0
    else
      f_temp_square_high = f_temp(Ts_C, T_min, T_opt, T_max, fmin)
    end if
  end function f_temp_square_high


  !> Calculate solar irradiance effect on canopy stomatal conductance.
  !!
  !! Leaf stomatal conductances for sun and shade leaves are combined according
  !! LAIsunfrac to give a canopy mean.
  pure real function f_light(f_lightfac, PARsun, PARshade, LAIsunfrac)
    real, intent(in) :: f_lightfac        !< Single leaf f_light coefficient
    real, intent(in) :: PARsun            !< PAR received by sunlit leaves (W m-2)
    real, intent(in) :: PARshade          !< PAR received by shaded leaves (W m-2)
    real, intent(in) :: LAIsunfrac        !< Fraction of canopy component that is sunlit

    real :: Flightsun, Flightshade

    ! TODO: does this need albedo?
    Flightsun = leaf_f_light(f_lightfac, PARsun)
    Flightshade = leaf_f_light(f_lightfac, PARshade)

    f_light = LAIsunfrac * Flightsun + (1.0 - LAIsunfrac) * Flightshade
  end function f_light


  !> Calculate solar irradiance effect on leaf stomatal conductance.
  pure real function leaf_f_light(f_lightfac, PAR)
    real, intent(in) :: f_lightfac        !< Single leaf f_light coefficient
    real, intent(in) :: PAR               !< Photosynthetically active radiation (W m-2)

    leaf_f_light = 1.0 - exp(-f_lightfac * (PAR * PAR_Wm2_to_photons))
  end function leaf_f_light


  !> Calculate the VPD effect on stomatal conductance from a linear function
  !! from (VPD_max, 1.0) to (VPD_min, fmin) - VPD_max should be smaller than
  !! VPD_min.
  pure real function f_VPD_linear(VPD, VPD_max, VPD_min, fmin)
    real, intent(in) :: VPD       !< Vapour pressure deficit (kPa)
    real, intent(in) :: VPD_max   !< VPD for maximum gsto (kPa)
    real, intent(in) :: VPD_min   !< VPD for minimum gsto (kPa)
    real, intent(in) :: fmin      !< Minimum f_VPD

    real :: slope

    slope = (1.0 - fmin) / (VPD_min - VPD_max)
    f_VPD_linear = max(fmin, min(1.0, (fmin + (VPD_min - VPD)*slope)))
  end function f_VPD_linear

  !> Calculate the VPD for a given f_VPD based on the linear VPD model.
  pure real function inverse_f_VPD_linear(f_VPD, VPD_max, VPD_min, fmin)
    real, intent(in) :: f_VPD
    real, intent(in) :: VPD_max   !< VPD for maximum gsto (kPa)
    real, intent(in) :: VPD_min   !< VPD for minimum gsto (kPa)
    real, intent(in) :: fmin      !< Minimum f_VPD

    real :: slope

    slope = (1.0 - fmin) / (VPD_min - VPD_max)
    inverse_f_VPD_linear = VPD_min - ((f_VPD - fmin) / slope)
  end function inverse_f_VPD_linear

  !> Calculate the VPD effect on stomatal conductance based on a simple
  !! logarithmic relationship.
  pure real function f_VPD_log(VPD, fmin)
    real, intent(in) :: VPD       !< Vapour pressure deficit (kPa)
    real, intent(in) :: fmin      !< Minimum f_VPD

    f_VPD_log = max(fmin, min(1.0, (1.0 - 0.6*log(VPD))))
  end function f_VPD_log

  !> Calculate the VPD for a given f_VPD based on the simple logarithmic VPD model.
  pure real function inverse_f_VPD_log(f_VPD, fmin)
    real, intent(in) :: f_VPD
    real, intent(in) :: fmin      !< Minimum f_VPD

    inverse_f_VPD_log = exp((1 - f_VPD) / 0.6)
  end function inverse_f_VPD_log


  !> Calculate the soil water effect on stomatal conductance using an
  !! exponential relationship.
  !!
  !! \f[
  !!    f_{SW} = \min(1, \max(f_{min}, a \times (-SWP)^b))
  !! \f]
  pure real function f_SWP_exp(a, b, fmin, SWP)
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: fmin      !< Minimum f_SWP
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    f_SWP_exp = min(1.0, max(fmin, a * (-SWP)**b))
  end function f_SWP_exp

  !> Calculate the SWP (MPa) for a given f_SWP based on the exponential
  !! relationship.
  !!
  !! Useful for determining SWP_max when using exponential f_SWP instead of
  !! linear f_SWP.
  pure real function inverse_f_SWP_exp(a, b, f_SWP)
    real, intent(in) :: a
    real, intent(in) :: b
    real, intent(in) :: f_SWP     !< Target f_SWP

    inverse_f_SWP_exp = -(f_SWP/a)**(1.0/b)
  end function inverse_f_SWP_exp

  !> Calculate the soil water effect on stomatal conductance with a linear
  !! function between (SWP_min, fmin) and (SWP_max, 1.0).
  pure real function f_SWP_linear(SWP_min, SWP_max, fmin, SWP)
    real, intent(in) :: SWP_min   !< SWP for minimum gsto (MPa)
    real, intent(in) :: SWP_max   !< SWP for maximum gsto (MPa)
    real, intent(in) :: fmin      !< Minimum f_SWP
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    f_SWP_linear = fmin + (1-fmin) * ((SWP_min - SWP)/(SWP_min - SWP_max))
    f_SWP_linear = min(1.0, max(fmin, f_SWP_linear))
  end function f_SWP_linear

  !> Calculate the soil water effect on stomatal conductance with a linear
  !! relationship based on available soil water (ASW).
  pure real function f_PAW(ASW_FC, fmin, ASW)
    real, intent(in) :: ASW_FC    !< Available soil water at field capacity (m)
    real, intent(in) :: fmin      !< Minimum f_PAW
    real, intent(in) :: ASW       !< Available soil water (m)

    f_PAW = fmin + (1-fmin) * ((100 * (ASW/ASW_FC)) - ASW_MIN)/(ASW_MAX - ASW_MIN)
    f_PAW = min(1.0, max(fmin, f_PAW))
  end function f_PAW

  !> Calculate the ASW (m) for a given f_PAW based on the f_PAW relationship.
  pure real function inverse_f_PAW(ASW_FC, fmin, f_PAW)
    real, intent(in) :: ASW_FC    !< Available soil water at field capacity (m)
    real, intent(in) :: fmin      !< Minimum f_PAW
    real, intent(in) :: f_PAW     !< Target f_PAW

    inverse_f_PAW = (ASW_MIN + ((f_PAW-fmin)/(1-fmin)) * (ASW_MAX-ASW_MIN)) * (ASW_FC/100)
  end function inverse_f_PAW


  !> Mean stomatal conductance of top leaf (mmol O3 m-2 PLA s-1).
  pure real function gsto_leaf(gp)
    type(GstoParams_t), intent(in) :: gp    !< Multiplicative gsto parameters

    gsto_leaf = gp%gmax * min(gp%leaf_f_phen, gp%f_O3) * gp%leaf_f_light * &
                max(gp%fmin, gp%f_temp * gp%f_VPD * gp%f_SW)
  end function gsto_leaf

  !> Mean canopy stomatal conductance (mmol O3 m-2 PLA s-1).
  !!
  !! TODO: can SMD canopy gsto just use this with f_SW=1.0 ?
  pure real function gsto_mean(gp)
    type(GstoParams_t), intent(in) :: gp    !< Multiplicative gsto parameters

    gsto_mean = gp%gmax * gp%gmorph * gp%f_phen * gp%f_light * &
                max(gp%fmin, gp%f_temp * gp%f_VPD * gp%f_SW)
  end function gsto_mean


  !> Stop stomatal conductance from increasing any more if the daily accumulated
  !! VPD threshold has been exceeded.
  pure real function apply_VPD_crit(VPD_crit, VPD_dd, old_gsto, new_gsto)
    real, intent(in) :: VPD_crit    !< Accumulated VPD threshold (kPa)
    real, intent(in) :: VPD_dd      !< Accumulated VPD for the day (kPa)
    real, intent(in) :: old_gsto    !< Previous stomatal conductance
    real, intent(in) :: new_gsto    !< New stomatal conductance

    if (VPD_dd >= VPD_crit) then
      apply_VPD_crit = min(old_gsto, new_gsto)
    else
      apply_VPD_crit = new_gsto
    end if
  end function apply_VPD_crit

end module DO3SE_Gsto_ml
