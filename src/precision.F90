module precision
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! set precision and use constants for gtce_cfd code
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! rdf = read real default (would be nice if it was shorter)
  !
  integer, parameter :: rdf = selected_real_kind(p=6) ! precision = 6 same as real*4

  ! idf = integer default
  !
  integer, parameter :: idf = selected_int_kind(r=9)  ! range = 9; same as integer*4

  ! define constants
  ! 
  real (kind = rdf), parameter :: zero = 0.0_rdf
  real (kind = rdf), parameter :: one  = 1.0_rdf
  real (kind = rdf), parameter :: two  = 2.0_rdf
  real (kind = rdf), parameter :: three= 3.0_rdf
  real (kind = rdf), parameter :: four = 4.0_rdf
  real (kind = rdf), parameter :: five = 5.0_rdf
  real (kind = rdf), parameter :: six  = 6.0_rdf
  real (kind = rdf), parameter :: seven= 7.0_rdf
  real (kind = rdf), parameter :: eight= 8.0_rdf
  real (kind = rdf), parameter :: nine = 9.0_rdf

  real (kind = rdf), parameter :: ten     = 10.0_rdf
  real (kind = rdf), parameter :: sixteen = 16.0_rdf

  real (kind = rdf), parameter :: one_fourth  = one / four
  real (kind = rdf), parameter :: pt1         = one / ten
  real (kind = rdf), parameter :: pt25        = one / four
  real (kind = rdf), parameter :: pt5         = one / two
  real (kind = rdf), parameter :: one_half    = one / two
  real (kind = rdf), parameter :: one_third   = one / three
  real (kind = rdf), parameter :: one_pt_five = three / two
  real (kind = rdf), parameter :: onept5      = three / two
  real (kind = rdf), parameter :: four_third  = four/ three

  real (kind = rdf), parameter :: pi = 3.141592653589793_rdf

end module precision


