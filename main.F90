program call_do3se
    use do3se_main
    implicit none
      INTEGER LON,LAT,NIV
      PARAMETER (LON=128, LAT=64, NIV=31)
      REAL, dimension(LON,LAT,NIV) :: T3D !!  Temperature in K

      REAL, dimension(LON,LAT,NIV) :: PL   !!  Pressuer in Pa
      REAL, dimension(LON,LAT,NIV) :: RINC  !!  In-canopy aerodynamic resistance (s m-1)
      REAL, dimension(LON,LAT,NIV) :: REXT  !!  External Plant Cuticle Resistance (s m-1)
      REAL, dimension(LON,LAT,NIV) :: RB    !!  Quai-laminar boundary layer resistance (s m-1)
      REAL, dimension(LON,LAT,NIV) :: RSTO  !!  Stomatal Resistance (s m-1)
      REAL, dimension(LON,LAT,NIV) :: RSUR  !!  Combined Surface Resistance (s m-1)
      REAL, dimension(LON,LAT,NIV) :: RA    !<  Aerodynamic Resistance (s m-1) between 50m
                              !!  reference height for this model.
      REAL, dimension(LON,LAT,NIV) :: RA_C  !<  Aerodynamic resistance (s m-1) between 50m
                              !!  and inside the canopy.
      REAL, dimension(LON,LAT,NIV) :: RGS   !!  Ground Surface Resistance.
      INTEGER I,K,L
      DO L=1,NIV

      DO K=1,LAT

      DO I=1,LON

        T3D(I,K,L)=200.0

      ENDDO

      ENDDO

      ENDDO
    call run_do3se(T3D, PL, RINC, REXT, RB, RSTO, RSUR, RA, RA_C, RGS)

end program call_do3se
