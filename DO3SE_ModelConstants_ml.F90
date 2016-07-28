module DO3SE_ModelConstants_ml
#include "interop_types.h"

  implicit none
  public

  REAL, parameter :: UNDEF = -999.0
  INTEGER, parameter :: IUNDEF = -999

  ! TODO: see how much we can remove these
  INTEGER, parameter :: MAX_LC = 3      !< Maximum number of land covers (used in some static allocations)
  INTEGER, parameter :: MAX_LAYERS = 5  !< Maximum number of layers (used in some static allocations)

  real, parameter :: DT = 60*60   !< Number of seconds in a timestep

  !> Canopy displacement (fraction of canopy height)
  real, parameter :: CANOPY_D = 0.7
  !> Canopy roughness length (fraction of canopy height)
  real, parameter :: CANOPY_Z0 = 0.1

  !> ASW (available soil water) for minimum gsto (percent of ASW at field capacity)
  real, parameter :: ASW_MIN = 0.0
  !> ASW (available soil water) for maximum gsto (percent of ASW at field capacity)
  real, parameter :: ASW_MAX = 50.0

  real, parameter :: LEAF_G_HEAT = 0.135
  real, parameter :: LEAF_G_H2O = 0.147
  real, parameter :: LEAF_G_CO2 = 0.110
  real, parameter :: LEAF_G_O3 = 0.105

contains
end module DO3SE_ModelConstants_ml
