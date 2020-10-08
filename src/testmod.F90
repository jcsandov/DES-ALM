

module rs_test

  implicit none 

  public :: rs_driver
  private

  integer :: rs_im
contains

  subroutine rs_driver(im, arr)

    integer, intent(in) :: im
    real(4), dimension(1:im), intent(inout) :: arr

    rs_im = im

    print *, 'im =', im

    call rs_set()

  end subroutine rs_driver

  subroutine rs_set ()

    real(4), dimension(1:rs_im) :: arr

    arr = 21
    print *, 'arr =', arr

  end subroutine rs_set

end module rs_test

program driver

  use rs_test

  real (kind = 4), dimension(1:3) :: var

  var = 1

  print *, 'before var =', var

  call rs_driver (3, var)
    
  print *, 'after var =', var

end program driver
