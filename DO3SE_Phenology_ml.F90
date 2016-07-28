module DO3SE_Phenology_ml

  use DO3SE_ModelConstants_ml, only: UNDEF
  use DO3SE_ConfigTypes_ml, only: Season_t
  use DO3SE_Util_ml

  implicit none
  private

  public :: LAI_day_PLF
  public :: SAI_wheat

contains

  !> Estimate LAI based on a piecewise linear function of the day of year.
  !!
  !! To handle all situations, including winter growing seasons, everything is
  !! re-indexed to be relative to SGS=0.
  !!
  !! Interpolation wraps around between SGS and EGS to handle LAI seasons that
  !! aren't a simple "bump".
  real function LAI_day_PLF(s, dd)
    type(Season_t), intent(in) :: s   !< Season configuration
    integer, intent(in) :: dd         !< Day of year

    real, dimension(2, 5) :: func
    integer :: dd_adj

    ! Build function
    func = reshape((/ real :: &
      s%SGS, s%LAI_a, &
      (s%SGS + s%LAI_1), s%LAI_b, &
      (s%EGS - s%LAI_2), s%LAI_c, &
      s%EGS, s%LAI_d, &
      (s%SGS + 365), s%LAI_a /),  shape(func))
    ! Re-index everything to SGS = 0, wrapping dates before SGS to the end of the year
    call reindex(func(1,:), real(s%SGS), 365.0)
    call assert(all(func(1,1:4) <= func(1,2:5)), "LAI_day_PLF: points not in order")
    dd_adj = dd - s%SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if
    ! Lookup value in PLF
    LAI_day_PLF = PLF_value(func, real(dd_adj))
  end function LAI_day_PLF


  !> Wheat stand area index based on the growing season.
  real function SAI_wheat(s, dd, LAI)
    type(Season_t), intent(in) :: s   !< Season configuration
    integer, intent(in) :: dd         !< Day of year
    real, intent(in) :: LAI           !< Leaf area index (m2 m-2)

    real, dimension(2, 4) :: func
    integer :: i, dd_adj

    SAI_wheat = UNDEF

    func = reshape((/ real :: &
      s%SGS, LAI, &
      (s%SGS + s%LAI_1), LAI + ((5.0/3.5) - 1) * LAI, &
      s%EGS + 1, LAI + 1.5, &
      s%SGS + 365, LAI /), shape(func))
    call reindex(func(1,:), real(s%SGS), 365.0)
    dd_adj = dd - s%SGS
    if (dd_adj < 0) then
      dd_adj = dd_adj + 365
    end if

    do i = 1, 4
      if (dd_adj < func(1,i)) then
        SAI_wheat = func(2,i)
        exit
      end if
    end do

    call assert(is_def(SAI_wheat), "SAI_wheat: no result")
  end function SAI_wheat

end module DO3SE_Phenology_ml
