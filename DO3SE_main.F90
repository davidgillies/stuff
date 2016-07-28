module DO3SE_main
contains
    subroutine run_do3se(T3D, PL, RINC, REXT, RB, RSTO, RSUR, RA, RA_C, RGS)

      use DO3SE_ml, only: DO3SE_State_t, read_DO3SE_Config, &

                          Location_t, LandCoverConfig_t, &

                          MetConfig_t, SMDConfig_t

      use DO3SE_Util_ml, only: assert, assert_not, inverted_assert, DO3SE_assert
      INTEGER LON,LAT,NIV
      PARAMETER (LON=128, LAT=64, NIV=31)
!     Model grid
      REAL, intent(in) :: T3D (LON,LAT,NIV) !!  Temperature in K

      REAL, intent(in) :: PL  (LON,LAT,NIV) !!  Pressuer in Pa

      REAL, intent(out) :: RINC (LON,LAT,NIV) !!  In-canopy aerodynamic resistance (s m-1)
      REAL, intent(out) :: REXT (LON,LAT,NIV) !!  External Plant Cuticle Resistance (s m-1)
      REAL, intent(out) :: RB   (LON,LAT,NIV) !!  Quai-laminar boundary layer resistance (s m-1)
      REAL, intent(out) :: RSTO (LON,LAT,NIV) !!  Stomatal Resistance (s m-1)
      REAL, intent(out) :: RSUR (LON,LAT,NIV) !!  Combined Surface Resistance (s m-1)
      REAL, intent(out) :: RA   (LON,LAT,NIV) !<  Aerodynamic Resistance (s m-1) between 50m
                              !!  reference height for this model.
      REAL, intent(out) :: RA_C (LON,LAT,NIV) !<  Aerodynamic resistance (s m-1) between 50m
                              !!  and inside the canopy.
      REAL, intent(out) :: RGS  (LON,LAT,NIV) !!  Ground Surface Resistance.





      INTEGER I,K,L


!     Number of tracers

      INTEGER NTRA

      PARAMETER (NTRA=100)


!     Dummy tomcat arrays



      ! DO3SE config (for now assuming 1 land cover type)

      type(Location_t) :: location  ! Used as initial values, overwrite per grid square

      type(LandCoverConfig_t), target :: LC_conf

      type(MetConfig_t), target :: met_conf

      type(SMDConfig_t), target :: SMD_conf

      logical :: DO3SE_success


      ! DO3SE model state (1 per lat/lon, only bottom layer)

      type(DO3SE_State_t), dimension(LON,LAT) :: DO3SE



      ! Wire up DO3SE assertion functions to implementations (used in config checking)

      assert => DO3SE_assert

      assert_not => inverted_assert

      ! Load configuration from file

      DO3SE_success = read_DO3SE_Config("do3se.nml", location, LC_conf, met_conf, SMD_conf)

      call assert(DO3SE_success, "Failed to read DO3SE config")

      do K=1,LAT

      do I=1,LON

      associate (this => DO3SE(I,K))

        ! Make a copy of location config to update per grid square

        allocate(this%location)

        this%location = location

        ! TODO: update lat/lon/elev of location

        this%location%lat = 50.0

        this%location%lon = 10.0

        this%location%elev = 100.0

        ! All states point at the same configuration otherwise (for now?)

        this%LC_conf => LC_conf

        this%met_conf => met_conf

        this%SMD_conf => SMD_conf

        ! Initialise the state, allocating DO3SE's variables

        call this%init()

      end associate

      end do

      end do


!     Set up some dummy values for TOMCAT arrays




!     Call DO3SE routines..


      ! Call DO3SE for the bottom layer of every LAT/LON

      ! TODO: where do timesteps happen?

      do K=1,LAT

      do I=1,LON

      associate (this => DO3SE(I,K))

        this%V%dd = 1               ! Day of year

        this%V%hr = 0               ! Hour of day (0-23)

        this%V%met%Ts_C = T3D(I,K,1) - 273.15   ! Surface air temperature (degrees C)

        this%V%met%P = 0.001 * PL(I,K,1)    ! Atmospheric pressure (kPa)

        this%V%met%precip = 0.1     ! Precipitation (mm)

        this%V%met%u = 0.8          ! Wind speed (m s-1)

        this%V%met%RH = 0.5         ! Relative humidity (fraction)

        this%V%met%O3 = 40.0        ! O3 concentration (ppb)

        this%V%met%CO2 = 392.0      ! CO2 concentration (ppm)

        this%V%met%Idrctt = 123.0   ! Direct PAR irradiance (W m-2)

        this%V%met%Idfuse = 123.0   ! Diffuse PAR irradiance (W m-2)

        call this%run_hour()

        RINC(I,K,1) = this%V%rmodel_O3%Rinc(1)
        REXT(I,K,1) = this%V%rmodel_O3%Rext(1)
        RB(I,K,1) = this%V%rmodel_O3%Rb
        RSTO(I,K,1) = this%V%rmodel_O3%Rsto(1)
        RSUR(I,K,1) = this%V%rmodel_O3%Rsur(1)
        RA(I,K,1) = this%V%rmodel_O3%Ra
        RA_C(I,K,1) = this%V%rmodel_O3%Ra_c
        RGS(I,K,1) = this%V%rmodel_O3%Rgs

        include "write.inc"

        ! TODO: save DO3SE values somewhere?

      end associate

      end do

      end do


!     Check values passed back


      END


end module
