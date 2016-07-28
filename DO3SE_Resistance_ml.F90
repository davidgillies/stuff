module DO3SE_Resistance_ml

  use DO3SE_ModelConstants_ml, only: UNDEF, MAX_LAYERS
  use DO3SE_Util_ml, only: alloc_init
#include "interop_types.h"

  implicit none
  private

  !> Resistance model.  All LAI-related resistances are bulk.
  TYPE, public :: ResistanceModel_t
    INTEGER :: nL = 0     !< Number of layers used
    REAL :: Ra_c = UNDEF  !< Aerodynamic resistance (s m-1) between 50m and
                          !! inside the canopy.
    REAL :: Ra = UNDEF    !< Aerodynamic resistance (s m-1) between 50m and the
                          !! reference height for this model (e.g. canopy
                          !! height or measurement height for O3).
    REAL :: Rb = UNDEF    !< Quasi-laminar boundary layer resistance (s m-1)
    REAL, dimension(MAX_LAYERS) :: &
      Rinc = UNDEF, &     !< In-canopy aerodynamic resistance (s m-1)
      Rext = UNDEF, &     !< External plant cuticle resistance (s m-1)
      Rsto = UNDEF        !< Stomatal resistance (s m-1)
    REAL :: Rgs = UNDEF   !< Ground surface resistance (s m-1)
    REAL, dimension(MAX_LAYERS) :: &
      Rsur = UNDEF, &     !< Combined surface resistance (s m-1)
      Rtotal = UNDEF      !< Total resistance for each layer downwards (s m-1)
  end type

  !> Leaf-level resistance model.  All LAI-related resistances are mean.
  TYPE, public :: LeafResistanceModel_t
    REAL :: Rb = UNDEF    !< Leaf boundary layer resistances (s m-1)
    REAL :: Rext = UNDEF  !< Leaf external plant cuticle resistance (s m-1)
    REAL :: Rsto = UNDEF  !< Leaf stomatal resistance (s m-1)
  end type

  public :: init_ResistanceModel
  public :: Ra_simple
  public :: Rb
  public :: leaf_gb
  public :: leaf_rb
  public :: Rinc
  public :: Rinc_prototype
  public :: Rext
  public :: Rsto
  public :: Rsur
  public :: Rtotal
  public :: deposition_velocity
  public :: stomatal_flux_rate

contains

  !> Initialise a ResistanceModel_t with the specified number of layers.
  pure subroutine init_ResistanceModel(rmodel, nL)
    type(ResistanceModel_t), intent(inout) :: rmodel
    integer, intent(in) :: nL

    rmodel%nL = nL
    rmodel%Ra_c = UNDEF
    rmodel%Ra = UNDEF
    rmodel%Rb = UNDEF
    rmodel%Rinc = UNDEF
    rmodel%Rext = UNDEF
    rmodel%Rsto = UNDEF
    rmodel%Rgs = UNDEF
    rmodel%Rsur = UNDEF
  end subroutine init_ResistanceModel

  !> Calculate aerodynamic resistance (Ra, s m-1) between two heights using a
  !! simple, neutral stability model.
  !!
  !! Must satisfy \f$z_2 \leq z_1\f$, \f$z_2 \gt d\f$ and \f$z_1 \gt d\f$.
  pure real function Ra_simple(ustar, z1, z2, d)
    real, intent(in) :: ustar   !< Friction velocity (m/s)
    real, intent(in) :: z1      !< Lower height (m)
    real, intent(in) :: z2      !< Upper height (m)
    real, intent(in) :: d       !< Zero displacement height (m)

    real, parameter :: K = 0.41 ! von Karman's constant

    Ra_simple = (1.0 / (ustar * K)) * log((z2 - d) / (z1 - d))
  end function Ra_simple

  !> Calculate quasi-laminar boundary layer resistance (Rb, s m-1) based on a
  !! given friction velocity and diffusivity.
  pure real function Rb(ustar, diff)
    real, intent(in) :: ustar   !< Friction velocity (m s-1)
    real, intent(in) :: diff    !< Molecular diffusivity in air (m2 s-1)

    real, parameter :: PR = 0.72    ! Prandtl number
    real, parameter :: K = 0.41     ! von Karman's constant
    real, parameter :: V = 0.000015 ! Kinematic viscosity of air at 20 C (m2 s-1)

    Rb = (2.0 / (K * ustar)) * (((V/diff)/PR)**(2.0/3.0))
  end function Rb

  !> Calculate leaf-level quasi-laminar boundary layer conductance
  !! (gb, mol m-2 s-1), for a particular kind of quantity specified by the base
  !! conductance of a single leaf surface, G.
  elemental real function leaf_gb(G, Lm, u)
    real, intent(in) :: G     !< Leaf surface conductance (mol m-2 s-1)
    real, intent(in) :: Lm    !< Cross-wind leaf dimension (m)
    real, intent(in) :: u     !< Wind speed (m s-1)

    ! G * 2 : from single surface to PLA (both sides of leaf)
    leaf_gb = (G * 2) * sqrt(u / Lm)
  end function leaf_gb

  !> Calculate leaf-level quasi-laminar boundary layer resistance (rb, s m-1)
  !! from conductance (gb, mol m-2 s-1).
  elemental real function leaf_rb(gb)
    real, intent(in) :: gb    !< Leaf boundary layer conductance (mol m-2 s-1)

    ! gb / 41 : 'mol m-2 s-1' to 'm s-1'
    ! 1 / gb  : 'm s-1' to 's m-1'
    leaf_rb = 41 / gb
  end function leaf_rb

  !> Estimate in-canopy aerodynamic resistance (Rinc, s m-1).
  !!
  !! This is the older single-layer DO3SE method.
  !!
  !! TODO: to use in a multilayer model, what does h represent?  Height above 
  !!       ground, or thickness of layer?
  pure real function Rinc(SAI, h, ustar)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)
    real, intent(in) :: h     !< Vegetation height (m)
    real, intent(in) :: ustar !< Friction velocity (m s-1)

    real, parameter :: Rinc_b = 14    ! Rinc coefficient

    Rinc = Rinc_b * SAI * h/ustar
  end function Rinc

  !> Estimate in-canopy aerodynamic resistance (Rinc, s m-1).
  !!
  !! This is the experimental method developed for the Keenley grassland
  !! multilayer model.
  !! 
  !! TODO: decide if we should keep this method
  pure real function Rinc_prototype(SAI, ustar)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)
    real, intent(in) :: ustar !< Friction velocity (m s-1)

    real, parameter :: Rinc_b = 14    ! Rinc coefficient

    Rinc_prototype = Rinc_b * SAI * Rinc_b/ustar
  end function Rinc_prototype

  !> Estimate external plant cuticle resistance (Rext, s m-1).
  pure real function Rext(SAI)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)

    real, parameter :: Rext_base = 2500

    Rext = Rext_base / SAI
  end function Rext

  !> Convert stomatal conductance to stomatal resistance (Rsto, s m-1).
  !!
  !! The maximum stomatal resistance is capped to prevent infinite values
  !! when the conductance is 0.
  elemental real function Rsto(Gsto)
    real, intent(in) :: Gsto    !< Stomatal conductance (mmol m-2 s-1)

    real, parameter :: MAX_RSTO = 100000

    ! (gsto in m s-1) = 41000 * (gsto in mmol m-2 s-1)
    ! (rsto in s m-1) = 1 / (gsto in m s-1)
    Rsto = min(MAX_RSTO, 41000.0 / Gsto)
  end function Rsto

  !> Calculate per-layer surface resistance - combined Rb, Rsto and Rext.
  ! TODO: per-layer Rb
  pure function Rsur(nL, Rb, Rsto, Rext, LAI, SAI)
    integer, intent(in) :: nL
    real, intent(in) :: Rb
    real, dimension(nL), intent(in) :: Rsto
    real, dimension(nL), intent(in) :: Rext
    real, dimension(nL), intent(in) :: LAI
    real, dimension(nL), intent(in) :: SAI

    real, dimension(nL) :: Rsur
    integer :: i

    do i = 1, nL
      ! TODO: let infinities happen and propagate through this? 1/Inf = 0?
      if (LAI(i) > 0) then
        ! LAI (and SAI) > 0, include Rsto and Rext components
        Rsur(i) = Rb + 1/(1/Rsto(i) + 1/Rext(i))
      else if (SAI(i) > 0) then
        ! Only SAI, omit the Rsto component
        Rsur(i) = Rb + Rext(i)
      else
        ! No foliage, very high resistance!
        ! TODO: find a justification for this, probably based on Rsto
        ! TODO: have an "R_INF" constant?
        Rsur(i) = 1000000
      end if
    end do
  end function Rsur

  !> Calculate multi-layer Rtotal - the total resistance for each layer and
  !! everything below that layer.
  pure function Rtotal(nL, Rsur, Rinc, Rgs)
    integer, intent(in) :: nL
    real, dimension(nL), intent(in) :: Rsur
    real, dimension(nL), intent(in) :: Rinc
    real, intent(in) :: Rgs

    real, dimension(nL) :: Rtotal

    real, dimension(nL+1) :: tmp
    integer :: i

    tmp(nL+1) = Rgs
    do i = nL, 1, -1
      tmp(i) = 1/(1/Rsur(i) + 1/(Rinc(i) + tmp(i+1)))
    end do
    Rtotal = tmp(1:nL)
  end function Rtotal

  !> Calculate deposition velocity (\f$V_d\f$) from a canopy resistance model.
  pure real function deposition_velocity(rmodel)
    type(ResistanceModel_t), intent(in) :: rmodel

    deposition_velocity = 1.0 / (rmodel%Ra_c + rmodel%Rtotal(1))
  end function deposition_velocity

  !> Calculate the rate of stomatal flux for a leaf resistance model.
  elemental real function stomatal_flux_rate(leaf_rmodel)
    type(LeafResistanceModel_t), intent(in) :: leaf_rmodel

    real :: leaf_r

    leaf_r = 1.0 / ((1.0 / leaf_rmodel%Rsto) + (1.0 / leaf_rmodel%Rext))
    stomatal_flux_rate = (1.0/leaf_rmodel%Rsto) * (leaf_r / (leaf_rmodel%Rb + leaf_r))
  end function stomatal_flux_rate

end module DO3SE_Resistance_ml
