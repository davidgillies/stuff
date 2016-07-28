module DO3SE_Met_ml

  use DO3SE_PhysicalConstants_ml, only: DEG2RAD, seaP, T0, PI, PARfrac, &
                                        PAR_Wm2_to_photons, DIFF_O3, &
                                        SBC
  use DO3SE_ModelConstants_ml, only: UNDEF, CANOPY_D, CANOPY_Z0
  use DO3SE_ConfigTypes_ml
  use DO3SE_Resistance_ml
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"
#include "interop_types.h"

  implicit none
  private

  !> "Global" meteorological data, i.e. external to and independent of the
  !! canopy.
  TYPE, public :: MetData_t
    REAL :: Ts_C = UNDEF      !< Surface air temperature (degrees C)
    REAL :: P = UNDEF         !< Atmospheric pressure (kPa)
    REAL :: precip = UNDEF    !< Precipitations (mm)
    REAL :: ustar = UNDEF     !< Friction velocity over target canopy (m s-1)
    REAL :: u = UNDEF         !< Measured windspeed (m s-1)
    REAL :: u_50 = UNDEF      !< Decoupled windspeed at 50m (m s-1)
    REAL :: O3 = UNDEF        !< Measured O3 concentration (ppb)
    REAL :: O3_50 = UNDEF     !< Decoupled O3 concentration at 50m (ppb)

    REAL :: sinB = UNDEF      !< sin() of solar elevation angle
    REAL :: Rn = UNDEF        !< Net radiation (MJ m-2 h-1)
    REAL :: R = UNDEF         !< Global radiation (W m-2)
    REAL :: PAR = UNDEF       !< Photosynthetically active radiation (W m-2)
    REAL :: PPFD = UNDEF      !< Photosynthetic photon flux density (umol m-2 s-1)
    REAL :: Idrctt = UNDEF    !< Direct PAR irradiance (W m-2)
    REAL :: Idfuse = UNDEF    !< Diffuse PAR irradiance (W m-2)

    REAL :: VPD = UNDEF       !< Vapour pressure deficit (kPa)
    REAL :: esat = UNDEF      !< Saturated vapour pressure (kPa)
    REAL :: eact = UNDEF      !< Actual vapour pressure (kPa)
    REAL :: RH = UNDEF        !< Relative humidity (fraction)

    REAL :: CO2 = UNDEF       !< CO2 concentration (ppm)
  end type MetData_t

  !> "Local" meteorological data, i.e. dependent on the canopy.
  !!
  !! Represents the value at the top of the canopy - should be used as an array
  !! for multiple layers, representing the value at the top of each layer.
  TYPE, public :: MicroMetData_t
    REAL :: PARsun = UNDEF    !< PAR received by sunlit leaves (W m-2)
    REAL :: PARshade = UNDEF  !< PAR received by shaded leaves (W m-2)

    REAL :: u = UNDEF         !< Windspeed (m s-1)
    REAL :: O3 = UNDEF        !< O3 concentration (ppb)

    REAL :: Tleaf_C = UNDEF   !< Leaf temperature (degrees C)
  end type MicroMetData_t

  public :: met_humidity
  public :: saturated_vapour_pressure
  public :: met_radiation
  public :: PAR_sun_shade
  public :: met_windspeed
  public :: multi_layer_windspeed
  public :: sunlit_LAI
  public :: MLMC_sunlit_LAI
  public :: met_O3
  public :: multi_layer_O3
  public :: O3_ppb_to_nmol_factor
  public :: met_CO2
  public :: leaf_net_radiation
  public :: leaf_temp_EB
  public :: leaf_temp_de_Boeck

contains

  !> Define VPD, relative humidity (RH), saturated vapour pressure and actual
  !! vapour pressure from either VPD or RH being set.
  pure subroutine met_humidity(met)
    type(metdata_t), intent(inout) :: met

    met%esat = saturated_vapour_pressure(met%Ts_C)
    if (is_undef(met%VPD) .and. is_undef(met%RH)) then
      ! TODO: error if neither VPD or RH supplied
    else if (is_undef(met%RH)) then
      ! Calculate relative humidity from VPD
      met%eact = met%esat - met%VPD
      met%RH = met%eact / met%esat
    else if (is_undef(met%VPD)) then
      ! Calculate VPD from relative humidity
      met%eact = met%esat * met%RH
      met%VPD = met%esat - met%eact
    end if
  end subroutine met_humidity


  !> Calculate saturated vapour pressure (kPa).
  pure real function saturated_vapour_pressure(Ts_C)
    real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)

    saturated_vapour_pressure = 0.611 * exp(17.27 * Ts_C / (Ts_C + 237.3))
  end function saturated_vapour_pressure


  pure subroutine met_radiation(met, loc, dd, hr)
    type(MetData_t), intent(inout) :: met
    type(Location_t), intent(in) :: loc
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)

    ! Solar elevation
    met%sinB = solar_elevation(loc%lat, loc%lon, dd, hr)

    ! Irradiance
    if (is_undef(met%PAR)) then
      ! Calculate PAR from some other available source
      if (is_def(met%Idrctt) .and. is_def(met%Idfuse)) then
        ! PAR is sum of direct and diffuse PAR
        met%PAR = met%Idrctt + met%Idfuse
      else if (is_def(met%PPFD)) then
        ! PAR from PPFD
        met%PAR = met%PPFD / PAR_Wm2_to_photons
      else if (is_def(met%R)) then
        ! Estimate PAR from global radiation
        met%PAR = met%R * PARfrac
      else
        ! TODO: error if no source of PAR is available
      end if
    end if
    if (is_undef(met%PPFD)) then
      ! If PPFD is missing, convert from PAR
      met%PPFD = met%PAR * PAR_Wm2_to_photons
    end if
    if (is_undef(met%Idrctt) .or. is_undef(met%Idfuse)) then
      ! If direct + diffuse are missing, estimate from PAR
      call PAR_direct_diffuse(met%PAR, met%sinB, met%P, met%Idrctt, met%Idfuse)
    end if
    if (is_undef(met%R)) then
      ! If global radiation is absent, estimate from PAR
      met%R = met%PAR / PARfrac
    end if

    ! Net radiation
    if (is_undef(met%Rn)) then
      met%Rn = net_radiation(loc%lat, loc%lon, loc%elev, loc%albedo, &
                             dd, hr, met%sinB, met%R, met%Ts_C, met%eact)
    end if
  end subroutine met_radiation


  pure function solar_elevation(lat, lon, dd, hr) result(sinB)
    real, intent(in) :: lat     !< Latitude (degrees North)
    real, intent(in) :: lon     !< Longitude (degrees East)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)

    real :: sinB      ! Output: sin() of solar elevation angle

    real :: t0_, dec, h

    t0_ = solar_noon(lon, dd)
    dec = solar_declination(dd)

    ! Hour-angle of the sun
    h = (15 * (hr - t0_)) * DEG2RAD

    ! sin() of solar elevation angle
    sinB = sin(lat * DEG2RAD)*sin(dec) + cos(lat * DEG2RAD)*cos(dec)*cos(h)
    ! TODO: should this line be removed? what effect does it have?  Does any
    !       use of sinB happen when sinB < 0?
    sinB = max(0.0, sinB)
  end function solar_elevation


  !> Calculate solar noon (hour of day).
  pure real function solar_noon(lon, dd)
    real, intent(in) :: lon     !< Longitude (degrees East)
    integer, intent(in) :: dd   !< Day of year

    real :: f, e, lonm, LC

    ! Solar noon correction for day of year
    f = (279.575 + (0.9856 * dd)) * DEG2RAD
    e = (-104.7*sin(f) + 596.2*sin(2*f) + 4.3*sin(3*f) - 12.7*sin(4*f) &
        - 429.3*cos(f) - 2.0*cos(2*f) + 19.3*cos(3*f)) / 3600
    ! Calculate the longitudinal meridian
    lonm = nint(lon / 15.0) * 15.0
    ! Solar noon, with day of year and longitudinal correction
    LC = (lon - lonm) / 15
    solar_noon = 12 - LC - e
  end function solar_noon


  !> Calculate solar declination (radians).
  pure real function solar_declination(dd)
    integer, intent(in) :: dd   !< Day of year

    solar_declination = (-23.4 * cos((360 * ((dd + 10) / 365.0))*DEG2RAD))*DEG2RAD
  end function solar_declination


  !> Estimate diffuse and direct PAR components.
  pure subroutine PAR_direct_diffuse(PAR, sinB, P, Idrctt, Idfuse)
    real, intent(in) :: PAR       !< Photosynthetically active radiation (W m-2)
    real, intent(in) :: sinB      !< sin() of solar elevation angle
    real, intent(in) :: P         !< Atmospheric pressure (kPa)

    real, intent(out) :: Idrctt   !< Direct PAR irradiance (W m-2)
    real, intent(out) :: Idfuse   !< Diffuse PAR irradiance (W m-2)

    real :: m, pPARdir, pPARdif, pPARtotal, ST, fPARdir, fPARdif

    if (sinB > 0.0) then
      m = 1.0 / sinB

      ! Potential direct PAR
      pPARdir = 600 * exp(-0.185 * (P / seaP) * m) * sinB
      ! Potential diffuse PAR
      pPARdif = 0.4 * (600 - pPARdir) * sinB
      ! Potential total PAR
      pPARtotal = pPARdir + pPARdif

      ! Sky transmissivity
      ST = max(0.21, min(0.9, PAR / pPARtotal))

      ! Direct and diffuse fractions
      fPARdir = (pPARdir / pPARtotal) * (1.0 - ((0.9 - ST) / 0.7)**(2.0/3.0))
      fPARdif = 1 - fPARdir

      ! Apply calculated direct and diffuse fractions to PARtotal
      Idrctt = fPARdir * PAR
      Idfuse = fPARdif * PAR
    else
      Idrctt = 0.0
      Idfuse = 0.0
    end if
  end subroutine PAR_direct_diffuse


  !> Estimate net radiation (MJ m-2 h-1).
  pure real function net_radiation(lat, lon, elev, albedo, dd, hr, sinB, R, Ts_C, eact)
    real, intent(in) :: lat     !< Latitude (degrees North)
    real, intent(in) :: lon     !< Longitude (degrees East)
    real, intent(in) :: elev    !< Elevation (m above sea level)
    real, intent(in) :: albedo  !< Surface albedo (fraction)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)
    real, intent(in) :: sinB    !< sin() of solar elevation angle
    real, intent(in) :: R       !< Global radiation (W m-2)
    real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)
    real, intent(in) :: eact    !< Actual vapour pressure (kPa)

    real, parameter :: Gsc = 0.082          ! Solar constant (MJ/m^2/min)
    real, parameter :: SBC = 4.903e-9 / 24  ! Stephan Boltzman constant

    real :: lat_rad, R_MJ, t0_, h, h1, h2, dr, dec, Re, pR, Rnl, Rns

    if (sinB <= 0) then
      net_radiation = 0.0
    else
      ! Latitude in radians
      lat_rad = lat * DEG2RAD

      ! Convert global radiation W m-2 to MJ m-2 s-1
      R_MJ = R * 0.0036

      ! Hour-angle of the sun
      t0_ = solar_noon(lon, dd)
      h = (15 * (hr - t0_)) * DEG2RAD
      h1 = h - (PI/24)
      h2 = h + (PI/24)

      dr = 1 + (0.033 * cos(((2 * PI) / 365) * dd))
      dec = solar_declination(dd)
      ! External radiation (with fix to stop div by zero)
      ! TODO: fix this to be less hackish
      Re = max(0.00000000001, &
               ((12 * 60) / PI) * Gsc * dr * ((h2 - h1) * sin(lat_rad) * sin(dec) &
               + cos(lat_rad) * cos(dec) * (sin(h2) - sin(h1))))
      ! TODO: what was this for?
      !Re = max(0.0, ((12*60)/PI)*Gsc*dr*sinB)

      ! Calculate net longwave radiation
      pR = (0.75 + (2e-5 * elev)) * Re

      Rnl = max(0.0, (SBC*((Ts_C + T0)**4)) * (0.34-(0.14*sqrt(eact))) &
                     * ((1.35*(min(1.0, R_MJ/pR)))-0.35))
      Rns = (1 - albedo) * R_MJ

      net_radiation = max(0.0, Rns - Rnl)
    end if
  end function net_radiation


  !> Estimate PAR received by sun and shade leaves within the canopy.
  pure subroutine PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI, PARsun, PARshade)
    real, intent(in) :: Idrctt        !< Direct PAR irradiance (W m-2)
    real, intent(in) :: Idfuse        !< Diffuse PAR irradiance (W m-2)
    real, intent(in) :: sinB          !< sin() of solar elevation angle
    real, intent(in) :: cosA          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)
    real, intent(in) :: LAI           !< Leaf area index (m2 m-2)

    real, intent(out) :: PARsun       !< PAR received by sunlit leaves (W m-2)
    real, intent(out) :: PARshade     !< PAR received by shaded leaves (W m-2)

    if (sinB > 0.0) then
      ! PAR flux densities evaluated using method of Norman (1982, p.79):
      ! "conceptually, 0.07 represents a scattering coefficient"
      PARshade = Idfuse * exp(-0.5 * LAI**0.8) + &
        0.07 * Idrctt * (1.1 - (0.1 * LAI)) * exp(-sinB)
      PARsun = Idrctt * 0.8 * (cosA / sinB) + PARshade
    else
      PARshade = 0.0
      PARsun = 0.0
    end if
  end subroutine PAR_sun_shade


  subroutine met_windspeed(met, umet, loc, h)
    type(MetData_t), intent(inout) :: met         !< Met data
    type(MicroMetData_t), intent(inout) :: umet   !< Top layer micromet data
    type(Location_t), intent(in) :: loc
    real, intent(in) :: h                         !< Canopy height

    real, parameter :: MIN_WINDSPEED = 0.01 ! Minimum windspeed value (m s-1)
    real, parameter :: MIN_USTAR = 0.0001   ! Minimum friction velocity value (m s-1)

    real :: h_u, z_u
    real :: u_d, u_z0, d, z0
    real :: ustar_ref

    if (loc%OTC) then
      h_u = h
      z_u = h
    else
      if (is_def(loc%h_u)) then
        h_u = loc%h_u
      else
        h_u = h
      end if
      z_u = loc%z_u
    end if

    u_d = h_u * CANOPY_D
    u_z0 = h_u * CANOPY_Z0
    d = h * CANOPY_D
    z0 = h * CANOPY_Z0

    ustar_ref = ustar_from_velocity(max(MIN_WINDSPEED, met%u), (z_u - u_d), u_z0)
    met%u_50 = max(MIN_WINDSPEED, velocity_from_ustar(ustar_ref, (50 - u_d), u_z0))
    met%ustar = max(MIN_USTAR, ustar_from_velocity(met%u_50, (50 - d), z0))
    umet%u = max(MIN_WINDSPEED, velocity_from_ustar(met%ustar, (h - d), z0))
  end subroutine met_windspeed


  subroutine multi_layer_windspeed(h, w, SAI, u_h, z, u_z)
    real, intent(in) :: h                 !< Canopy height (m)
    real, intent(in) :: w                 !< Leaf width (m)
    real, intent(in) :: SAI               !< Stand area index (m2 m-2)
    real, intent(in) :: u_h               !< Wind speed at canopy height (m s-1)
    real, dimension(:), intent(in) :: z   !< Heights within the canopy (m)

    real, dimension(size(z)), intent(out) :: u_z  !< Wind speed at each z (m s-1)

    real :: lm, a

    ! TODO: different definitions for lm
    lm = sqrt((4 * w * h) / (PI * SAI))

    a = sqrt((0.2 * SAI * h) / lm)

    u_z = u_h * exp(a * (z/h - 1))
  end subroutine multi_layer_windspeed


  !> Estimate windspeed velocity from friction velocity.
  pure elemental function velocity_from_ustar(ustar, z, z0) result (u)
    real, intent(in) :: ustar   ! Friction velocity (m/s)
    real, intent(in) :: z       ! Height above boundary, e.g. z - d (m)
    real, intent(in) :: z0      ! Roughness length, height at which u=0 (m)
    real :: u                   ! Output: velocity (m/s)

    real, parameter :: K = 0.41     ! von Karman's constant

    u = (ustar / K) * log(z / z0)
  end function velocity_from_ustar


  !> Estimate friction velocity from windspeed velocity.
  pure elemental function ustar_from_velocity(u, z, z0) result (ustar)
    real, intent(in) :: u   ! Velocity at height above boundary (m/s)
    real, intent(in) :: z   ! Height above boundary, e.g. z - d (m)
    real, intent(in) :: z0  ! Roughness length, height at which u=0 (m)
    real :: ustar           ! Output: friction velocity, ustar (m/s)

    real, parameter :: K = 0.41 ! von Karman's constant

    ustar = (u * K) / log(z / z0)
  end function ustar_from_velocity


  !> Calculate the sunlit LAI from total LAI and solar elevation angle.
  !!
  !! TODO: (2 * sinB) should be sinB/cosA
  pure real function sunlit_LAI(LAI, sinB)
    real, intent(in) :: LAI     !< Leaf area index (m^2/m^2)
    real, intent(in) :: sinB    !< sin() of solar elevation angle

    if (LAI > 0.0 .and. sinB > 0.0) then
      sunlit_LAI = ((1 - exp(-0.5 * LAI / sinB)) * (2 * sinB))
    else
      sunlit_LAI = 0.0
    end if
  end function sunlit_LAI

  !> Multi-layer multi-component sunlit LAI model.  LAI and LAIsunfrac must be
  !! the same dimension(nL,nLC).
  pure function MLMC_sunlit_LAI(LAI, sinB) result (LAIsunfrac)
    real, dimension(:,:), intent(in) :: LAI  !< Leaf area index (m^2/m^2)
    real, intent(in) :: sinB                 !< sin() of solar elevation angle

    real, dimension(size(LAI,1),size(LAI,2)) :: LAIsunfrac   !< Output: fraction of LAI that is sunlit

    real, dimension(0:size(LAI,1)) :: sunLAI_acc
    real :: sunLAI_layer
    real :: LAI_layer
    integer :: iL

    sunLAI_acc(0) = 0.0
    do iL = 1, size(LAI, 1)
      ! How much of "canopy so far" is sunlit?
      sunLAI_acc(iL) = sunlit_LAI(sum(LAI(1:iL,:)), sinB)
      ! How much of that is in this layer?
      sunLAI_layer = sunLAI_acc(iL) - sunLAI_acc(iL - 1)
      ! Fraction of LAI which is sunlit
      LAI_layer = sum(LAI(iL,:))
      if (LAI_layer > 0.0) then
        LAIsunfrac(iL,:) = sunLAI_layer / LAI_layer
      else
        LAIsunfrac(iL,:) = 0.0
      end if
    end do
  end function MLMC_sunlit_LAI


  !> Define an ozone concentration at the canopy
  subroutine met_O3(config, met, umet, loc, rmodel, h)
    type(MetConfig_t), intent(in) :: config         !< Met configuration
    type(MetData_t), intent(inout) :: met           !< Met data
    type(MicroMetData_t), intent(inout) :: umet     !< Top layer micromet data
    type(Location_t), intent(in) :: loc
    type(ResistanceModel_t), intent(in) :: rmodel   !< Resistance model for O3 over target canopy
    real, intent(in) :: h                           !< Canopy height

    real :: h_O3
    real :: O3_d, O3_z0
    real :: ustar_ref
    type(ResistanceModel_t) :: rmodel_ref

    select case (config%O3_method)
    case ("constant")
      ASSERT_DEFINED(config%O3_constant)
      met%O3 = config%O3_constant
    case ("offset")
      ASSERT_DEFINED(config%O3_constant)
      ASSERT_DEFINED(met%O3)
      met%O3 = max(0.0, met%O3 + config%O3_constant)
    case ("input")
      ASSERT_DEFINED(met%O3)
    case default
      UNKNOWN_STRING(config%O3_method)
    end select

    if (loc%OTC) then
      met%O3_50 = met%O3
      umet%O3 = met%O3
    else
      if (is_def(loc%h_O3)) then
        h_O3 = loc%h_O3
      else
        h_O3 = h
      end if

      O3_d = h_O3 * CANOPY_D
      O3_z0 = h_O3 * CANOPY_Z0

      ! TODO: this is unnecessary if h_O3 = h?
      ustar_ref = ustar_from_velocity(met%u_50, 50.0 - O3_d, O3_z0)
      rmodel_ref = rmodel
      rmodel_ref%Ra_c = Ra_simple(ustar_ref, O3_d + O3_z0, 50.0, O3_d)
      rmodel_ref%Ra = Ra_simple(ustar_ref, loc%z_O3, 50.0, O3_d)
      rmodel_ref%Rb = Rb(ustar_ref, DIFF_O3)

      met%O3_50 = O3_transfer_up(rmodel_ref, met%O3)
      umet%O3 = O3_transfer_down(rmodel, met%O3_50)
    end if
  end subroutine met_O3

  !> Calculate O3 concentration for all layers.  Requires that the value for the
  !! top layer (umet(1)%O3) is already known.
  subroutine multi_layer_O3(rmodel, O3_in, O3_out)
    type(ResistanceModel_t), intent(in) :: rmodel               !< Resistance model for O3 over target canopy
    real, intent(in) :: O3_in
    real, dimension(rmodel%nL), intent(out) :: O3_out

    real, dimension(rmodel%nL+1) :: bigR, smallR, C, ipiv_
    real, dimension(rmodel%nL+1,rmodel%nL+1) :: X

    integer :: i, j, info

    associate (nL => rmodel%nL, &
               Ra => rmodel%Ra, &
               Rb => rmodel%Rb, &
               Rinc => rmodel%Rinc(1:rmodel%nL), &
               Rsur => rmodel%Rsur(1:rmodel%nL), &
               Rgs => rmodel%Rgs)

      bigR(1) = Ra
      bigR(2:nL+1) = Rinc

      ! TODO: per-layer Rb
      smallR(1:nL) = Rsur(:)
      smallR(nL+1) = Rgs

      X = 0
      ! Iterate over columns
      do j = 1, nL + 1
        ! Fill in "R" values in upper triangle
        X(1:j,j) = bigR(1:j)
        ! Fill in "r" values on diagonal
        X(j,j) = X(j,j) + smallR(j)
        ! Fill in "r" values below diagonal
        if (j <= nL) then
          X(j+1,j) = -smallR(j)
        end if
      end do

      C = 0
      C(1) = O3_in

!#ifdef HAVE_LAPACK
      call SGESV(nL + 1, 1, X, nL + 1, ipiv_, C, nL + 1, info)
!#else
!      ERROR("Not built with LAPACK, multi_layer_O3 broken")
!#endif
      ASSERT(info == 0)
      C = smallR * C
      O3_out(:) = C(1:nL)

    end associate
  end subroutine multi_layer_O3


  !> Scale O3 concentration up from the resistance model's reference height
  !! to 50m.
  pure real function O3_transfer_up(rmodel, O3)
    type(ResistanceModel_t), intent(in) :: rmodel
    real, intent(in) :: O3

    real :: Vd

    Vd = deposition_velocity(rmodel)
    O3_transfer_up = O3 / (1.0 - (rmodel%Ra * Vd))
  end function O3_transfer_up


  !> Scale O3 concentration down from 50m to the resistance model's reference
  !! height.
  pure real function O3_transfer_down(rmodel, O3_50)
    type(ResistanceModel_t), intent(in) :: rmodel
    real, intent(in) :: O3_50

    real :: Vd

    Vd = deposition_velocity(rmodel)
    O3_transfer_down = O3_50 * (1.0 - (rmodel%Ra * Vd))
  end function O3_transfer_down


  !> Calculate the conversion factor between parts per billion and nmol m-3 for
  !! O3 at a given temperature and pressure.
  pure real function O3_ppb_to_nmol_factor(Ts_c, P)
    real, intent(in) :: Ts_C      !< Air temperature (degrees C)
    real, intent(in) :: P         !< Atmospheric pressure (kPa)

    real, parameter :: M_O3 = 48.0      ! Molecular weight of O3 (g)

    real :: Vn

    ! Specific molar volume of an ideal gas at this temperature and pressure
    Vn = 8.314510 * ((Ts_C + T0) / P)
    ! Conversion to nmol m-3 (1 microgram O3 = 20.833 nmol m-3)
    O3_ppb_to_nmol_factor = (1.0/Vn) * M_O3 * 20.833
  end function O3_ppb_to_nmol_factor


  !> Define a CO2 concentration, either using a supplied value or a configured
  !! constant value.
  subroutine met_CO2(config, met)
    type(MetConfig_t), intent(in) :: config
    type(MetData_t), intent(inout) :: met

    select case (config%CO2_method)
    case ("constant")
      met%CO2 = config%CO2_constant
    case ("input")
      ASSERT_DEFINED(met%CO2)
    case default
      UNKNOWN_STRING(config%CO2_method)
    end select
  end subroutine met_CO2


  !> Leaf net radiation calculation based on de Boeck (2012), which in turn is
  !! based on Campbell & Norman (1998).
  pure function leaf_net_radiation(R, e_a, T_air, T_leaf, albedo, cloud_cover) result(Rn)
    real, intent(in) :: R             !< Global radiation (W m-2)
    real, intent(in) :: e_a           !< Vapour pressure (kPa)
    real, intent(in) :: T_air         !< Air temperature (degrees C)
    real, intent(in) :: T_leaf        !< Leaf temperature (degrees C)
    real, intent(in) :: albedo        !< Soil albedo (fraction)
    real, intent(in) :: cloud_cover   !< Cloud cover (fraction)

    real :: Rn      !< Output: leaf net (short- and long-wave) radiation (W m-2)

    ! Leaf short-wave absorptivity (de Boeck (2012); from Campbell & Norman (1998), p153, in text)
    real, parameter :: alpha_s_leaf = 0.5
    ! Leaf long-wave absorptivity (de Boeck, 2012)
    real, parameter :: alpha_l_leaf = 0.97
    ! Soil long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_soil = 0.945
    ! Leaf long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_leaf = 0.97

    real :: T_soil, eps_ac, eps_l_sky, R_s_in, R_l_in, R_l_out

    ! TODO: should this be something else?
    T_soil = T_air

    ! Clear sky emissivity (p163, eq. 10.10)
    eps_ac = 1.72 * (e_a / (T_air + T0))**(1.0/7.0)
    ! Sky long-wave emmisivity (p164, eq. 10.12)
    eps_l_sky = (1 - 0.84 * cloud_cover) * eps_ac + 0.84 * cloud_cover

    ! Net radiation calculations from de Boeck (2012):
    ! Absorbed short-wave (direct on top of leaf, reflected from soil on underside)
    R_s_in = (0.5 * R * alpha_s_leaf) + (0.5 * R * albedo * alpha_s_leaf)
    ! Absorbed long-wave (emitted from soil on underside, emitted from sky on top)
    R_l_in = (0.5 * alpha_l_leaf * eps_l_soil * SBC * (T_soil + T0)**4) + &
             (0.5 * alpha_l_leaf * eps_l_sky * SBC * (T_air + T0)**4)
    ! Emitted long-wave
    R_l_out = eps_l_leaf * SBC * (T_leaf + T0)**4
    Rn = R_s_in + R_l_in - R_l_out
  end function leaf_net_radiation


  !> Leaf temperature from Environmental Biophysics (Campbell & Norman, 1998).
  !! TODO: make 1.4 "enhancement factor" optional
  pure function leaf_temp_EB(d, T_air, P, e_a, R_ni, u, g_vs) result(T_leaf)
    real, intent(in) :: d       !< Leaf characteristic dimension (m)
    real, intent(in) :: T_air   !< Air temperature (degrees C)
    real, intent(in) :: P       !< Air pressure (kPa)
    real, intent(in) :: e_a     !< Vapour pressure (kPa)
    real, intent(in) :: R_ni    !< Leaf net radiation (W m-2)
    real, intent(in) :: u       !< Wind speed (m s-1)
    real, intent(in) :: g_vs    !< Stomatal conductance to water vapour (mol m-2 s-1)

    real :: T_leaf    !< Output: leaf temperature (degrees C)

    ! Specific heat capacity of dry air (J mol-1 C-1; p279, table A1)
    real, parameter :: c_p = 29.3
    ! Latent heat of vapourisation of water (J mol-1; p37, in text)
    real, parameter :: lam = 44e3
    ! Emissivity of natural surface (fraction; p163, below eq. 10.9)
    real, parameter :: eps_s = 0.97

    real :: T_ak, g_r, g_Ha, g_Hr, g_va, g_v, a, b, c, e_s, VPD, delta, s, gam, gamstar

    T_ak = T_air + T0
    ! Radiative conductance to heat transfer (mol m-2 s-1; p188, eq. 12.7)
    g_r = 4 * eps_s * SBC * T_ak**3 / c_p
    ! Boundary layer (forced convection) conductance to heat transfer (mol m-2 s-1; p101, eq. 7.30)
    ! 1.4 enhancement factor for outdoor "turbulent wind" (p108, section 7.13)
    g_Ha = 1.4 * 0.135 * sqrt(u / d)
    ! Convective-radiative conductance to heat transfer (mol m-2 s-1; p188, eq. 12.9)
    g_Hr = g_Ha + g_r

    ! Boundary layer conductance to water vapour (mol m-2 s-1)
    ! 1.4 enhancement factor for outdoor "turbulent wind" (p108, section 7.13)
    g_va = 1.4 * 0.147 * sqrt(u / d)
    ! Whole-leaf (stomatal and boundary) conductance to water vapour (mol m-2 s-1; p224, eq. 14.2)
    ! (hypostomatus, adaxial stomatal conductance = 0)
    ! TODO: hypostomatus boolean switch
    g_v = leaf_conductance(g_vs, 0.0, g_va)

    ! Saturation vapour pressure coefficients (p41, below eq. 3.8)
    a = 0.611   ! (kPa)
    b = 17.502
    c = 240.97  ! (degrees C)
    ! Saturated vapour pressure (kPa; p41, eq. 3.8)
    e_s = a * exp(b * T_air / (T_air + c))
    ! Vapour deficit (kPa)
    VPD = e_s - e_a
    ! Slope of saturation vapour pressure function (p41, eq. 3.9)
    delta = b * c * e_s / (c + T_air)**2
    ! Slope of saturation mole fraction (p41, eq. 3.10)
    s = delta / P

    ! Thermodynamic psychrometer constant (C-1; p44, below eq. 3.16)
    gam = c_p / lam
    ! Apparent psychrometer constant (C-1; p225, below eq. 14.6)
    gamstar = gam * g_Hr / g_v
    ! Leaf temperature (degrees C; p225, eq. 14.6)
    T_leaf = T_air + (gamstar / (s + gamstar)) * (R_ni / (g_Hr * c_p) - VPD / (P * gamstar))
  contains
    pure real function leaf_conductance(gs_ab, gs_ad, ga)
      real, intent(in) :: gs_ab   !< Abaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: gs_ad   !< Adaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: ga      !< Boundary layer conductance (mol m-2 s-1)

      leaf_conductance = (0.5 * gs_ab * ga / (gs_ab + ga)) + (0.5 * gs_ad * ga / (gs_ad + ga))
    end function leaf_conductance
  end function leaf_temp_EB

  pure function leaf_temp_de_Boeck(R, e_a, T_air, initial_T_leaf, P, u, g_vs, hypostomatous, d, albedo, cloud_cover, &
                                   balance_threshold, adjustment_factor, max_iterations) result(T_leaf)
    real, intent(in) :: R                   !< Global radiation (W m-2)
    real, intent(in) :: e_a                 !< Vapour pressure (Pa)
    real, intent(in) :: T_air               !< Air temperature (degrees C)
    real, intent(in) :: initial_T_leaf      !< T_leaf starting value (degrees C)
    real, intent(in) :: P                   !< Air pressure (Pa)
    real, intent(in) :: u                   !< Wind speed (m s-1)
    real, intent(in) :: g_vs                !< Stomatal conductance to water vapour (mol m-2 s-1)
    logical, intent(in) :: hypostomatous    !< Are leaves hypostomatous?
    real, intent(in) :: d                   !< Leaf characteristic dimension (m)
    real, intent(in) :: albedo              !< Surface albedo (fraction)
    real, intent(in) :: cloud_cover         !< Cloud cover (fraction)
    real, intent(in) :: balance_threshold   !< Threshold within which to accept leaf energy balance
    real, intent(in) :: adjustment_factor   !< Multiplier for energy balance based T_leaf adjustments
    integer, intent(in) :: max_iterations   !< Maximum iteration count

    real :: T_leaf    !< Output: Leaf temperature (degrees C)

    ! Stefan-Boltzmann constant, W m-2 K-4 (Campbell & Norman (1998), p281, table A5)
    real, parameter :: SBC = 5.670373e-8
    ! Specific heat capacity of dry air, J mol-1 C-1 (Campbell & Norman (1998), p279, table A1)
    real, parameter :: c_p = 29.3

    ! Leaf short-wave absorptivity (de Boeck (2012); from Campbell & Norman (1998), p153, in text)
    real, parameter :: alpha_s_leaf = 0.5
    ! Leaf long-wave absorptivity (de Boeck, 2012)
    real, parameter :: alpha_l_leaf = 0.97
    ! Soil long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_soil = 0.945
    ! Leaf long-wave emissivity (de Boeck, 2012)
    real, parameter :: eps_l_leaf = 0.97

    real :: rho_s_soil, VPL, eps_ac, eps_l_sky, lam, T_soil, T_leaf_adjust, g_Ha, g_vb, g_v, e_s_Tleaf
    real :: R_s_in, R_l_in, R_l_out, H, lam_E, energy_balance
    integer :: i

    ! Soil short-wave reflectivity
    rho_s_soil = albedo
    ! Water vapour path length (de Boeck 2012)
    VPL = 46.5 * ((0.01 * e_a) / (T_air + T0))
    ! Clear sky emissivity (de Boeck 2012)
    eps_ac = 1 - (1 + VPL) * exp(-(1.2 + 3.0 * VPL)**0.5)
    ! Sky long-wave emmisivity (de Boeck (2012); from Campbell & Norman (1998), p164, eq. 10.12)
    eps_l_sky = (1 - 0.84 * cloud_cover) * eps_ac + 0.84 * cloud_cover
    ! Latent heat of vapourisation, J mol-1 (de Boeck (2012))
    lam = -42.575 * T_air + 44994

    ! Assume soil temperature is equal to air temperature
    ! TODO: use the banded estimate from de Boeck (2012)?
    T_soil = T_air

    ! Starting point
    T_leaf = initial_T_leaf
    T_leaf_adjust = 0.0

    do i = 1, max_iterations
      ! Apply adjustment (from previous iteration)
      T_leaf = T_leaf + T_leaf_adjust

      ! Boundary layer conductances (de Boeck (2012))
      ! XXX: Assume forced convection:
      g_Ha = 1.4 * 0.135 * sqrt(u / d)
      g_vb = 1.4 * 0.147 * sqrt(u / d)

      ! Total conductivity to water vapour
      if (hypostomatous) then
        g_v = combined_leaf_conductance(g_vs, 0.0, g_vb)
      else
        g_v = combined_leaf_conductance(g_vs, g_vs, g_vb)
      end if

      ! Saturated vapour pressure for leaf temperature, Pa
      e_s_Tleaf = saturated_vapour_pressure(T_leaf) * 1000  ! Converted from kPa to Pa

      ! Absorbed short-wave (direct on top of leaf, reflected from soil on underside) (de Boeck (2012))
      R_s_in = (0.5 * R * alpha_s_leaf) + (0.5 * R * rho_s_soil * alpha_s_leaf)

      ! Absorbed long-wave (emitted from soil on underside, emitted from sky on top) (de Boeck (2012))
      R_l_in = (0.5 * alpha_l_leaf * eps_l_soil * SBC * (T_soil + T0)**4) + &
               (0.5 * alpha_l_leaf * eps_l_sky * SBC * (T_air + T0)**4)

      ! Emitted long-wave (de Boeck (2012))
      R_l_out = eps_l_leaf * SBC * (T_leaf + T0)**4

      ! Sensible heat flux (de Boeck (2012))
      H = g_Ha * c_p * (T_leaf - T_air)

      ! Latent heat flux (de Boeck (2012))
      lam_E = lam * g_v * (e_s_Tleaf - e_a) / P

      ! Energy balance equation (de Boeck (2012))
      energy_balance = R_s_in + R_l_in - R_l_out - H - lam_E

      if (abs(energy_balance) <= balance_threshold) then
        exit
      else
        T_leaf_adjust = energy_balance * adjustment_factor
      end if
    end do
  contains
    pure function combined_leaf_conductance(gs_ab, gs_ad, ga) result(g)
      real, intent(in) :: gs_ab   !< Abaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: gs_ad   !< Adaxial stomatal conductance (mol m-2 s-1)
      real, intent(in) :: ga      !< Boundary layer conductance (mol m-2 s-1)
      real :: g     !< Output: combined leaf conductance (mol m-2 s-1)

      g = (0.5 * gs_ab * ga / (gs_ab + ga)) + (0.5 * gs_ad * ga / (gs_ad + ga))
    end function combined_leaf_conductance
  end function leaf_temp_de_Boeck

end module DO3SE_Met_ml
