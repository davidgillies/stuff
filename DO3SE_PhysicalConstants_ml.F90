module DO3SE_PhysicalConstants_ml

  implicit none
  public

  real, parameter :: PI = 3.141592653589793238
  real, parameter :: DEG2RAD = 0.017453292519943295

  !> Atmospheric pressure at sea level (kPa)
  real, parameter :: seaP = 101.325

  !> 0 degrees Celsius in Kelvin
  real, parameter :: T0 = 273.15
  !> Stefan-Boltzmann constant (W m-2 K-4)
  real, parameter :: SBC = 5.670373e-8

  !> Universal gas constant (J K-1 mol-1)
  ! TODO: update this value to 8.3144621
  real, parameter :: R = 8.314472

  !> Approximate fraction of global radiation in PAR waveband
  real, parameter :: PARfrac = 0.45
  !> PAR conversion from W m-2 to umol photons m-2 s-1
  real, parameter :: PAR_Wm2_to_photons = 4.57
  !> Net radiation conversion from MJ m-2 h-1 to W m-2
  real, parameter :: Rn_MJ_to_W = 277.8

  !> Molecular diffusivity of O3 in air (m2 s-1)
  real, parameter :: DIFF_O3 = 0.000015
  !> Molecular diffusivity of H2O (m2 s-1)
  real, parameter :: DIFF_H2O = 0.000025
  !> Ratio between molecular diffusivity of O3 and H2O
  real, parameter :: DRATIO = 0.663

contains
end module DO3SE_PhysicalConstants_ml
