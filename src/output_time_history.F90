
subroutine output_time_history
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global
  use global_param
  implicit none

  integer :: l, m, i, j, k

  integer :: p
  integer :: num_points

  integer, parameter :: ione = 1

  ! 20 is hard-wired
  integer, dimension(3,100) :: isave  ! 1 = i, j = 2, k = 3

!   ! des & duct, interpolation coefficients for D1 & D2 cross-sections
!   real (kind = rdf), parameter :: cd1d = 0.02641_rdf
!   real (kind = rdf), parameter :: cd1u = 0.05663_rdf
!   real (kind = rdf), parameter :: cd2d = 0.07162_rdf
!   real (kind = rdf), parameter :: cd2u = 0.04246_rdf

!   real (kind = rdf), dimension(4,jmg(ns),kmg(ns)) :: qd1, qd2
!   real (kind = rdf), dimension(jmg(ns),kmg(ns)) :: xnud1, xnud2

  ! file names
  character (len = 50, kind = 1) :: filename

!   ! convert notation to new code
!   num_points = monitor_num_points
!   do p = 1, num_points
!      isave(:,p) = monitor_point_ijk(p,:)
!   end do

!   ! store the whole flowfield for the for computing metrics
!   if ((ntime / icnw) * icnw == ntime) then
!      filename = ''
!      write(filename,'(a16,i6.6)') 'usolufiles/usolu',ntime
!      open  (1000, file = filename, form = 'unformatted')
!      write (1000) ((q(m,l), l = ls(1), le(1)), m = 1, 4)
!      if (turbulence) &
!      write (1000) (xnut(l), l = ls(1), le(1))
!      close (1000)
!   end if

!   ! store the flowfield on D1 & D2 cross-section for comparison
!   if (unsteady .and. abut) then
!      k = kmg(ns)
!      filename = ''
!      write(filename,'(a15,i6.6)') 'sym_plane/usolu',ntime
!      open  (unit = 11, file = filename, form = 'unformatted')
!      write (unit = 11) (((q(m,ln(i,j,k,ns)), i=1,img(ns)), j=1,jmg(ns)), m=1,4)
!      if ( turbulence ) &
!      write (unit = 11) ((xnut(ln(i,j,k,ns)), i=1,img(ns)), j=1,jmg(ns))
!      close (unit = 11)
!   end if

!   if (duct .and. des) then

!      do k = 1, kmg(ns)
!      do j = 1, jmg(ns)

!         l = ln(94,j,k,ns)
!         qd1(1:4,j,k) = (q(1:4,l) * cd1d + q(1:4,l-icg(ns)) * cd1u) &
!                      / (cd1d + cd1u)
!         xnud1(j,k)   = ( xnut(l) * cd1d +  xnut(l-icg(ns)) * cd1u) &
!                      / (cd1d + cd1u)

!         l = ln(133,j,k,ns)
!         qd2(1:4,j,k) = (q(1:4,l) * cd2d + q(1:4,l-icg(ns)) * cd2u) &
!                      / (cd2d + cd2u)
!         xnud2(j,k)   = ( xnut(l) * cd2d +  xnut(l-icg(ns)) * cd2u) &
!                      / (cd2d + cd2u)

!      end do
!      end do

!      filename = ''
!      write(filename,'(a14,i6.6)') 'station_1/sta1',ntime
!      open (unit = 21, file = filename, form = 'unformatted')
!      write(unit = 21) (((qd1(m,j,k), j = 1, jmg(ns)), k = 1, kmg(ns)), m = 1, 4)
!      write(unit = 21) ((xnud1(j,k), j = 1, jmg(ns)), k = 1, kmg(ns))
!      close(unit = 21)

!      filename = ''
!      write(filename,'(a14,i6.6)') 'station_2/sta2',ntime
!      open (unit = 22, file = filename, form = 'unformatted')
!      write(unit = 22) (((qd2(m,j,k), j = 1, jmg(ns)), k = 1, kmg(ns)), m = 1, 4)
!      write(unit = 22) ((xnud2(j,k), j = 1, jmg(ns)), k = 1, kmg(ns))
!      close(unit = 22)

!   end if

  ! performance monitor
  open  (unit = 1100, file = 'performance', position = 'append')
  write (unit = 1100, fmt = '(2(1x, i5.5),2(1X,G15.7))')  &
                   ntime, itc, time_step_time, total_time
  close (unit = 1100)


  ! history monitoring of solution may need to be performed in
  ! a distributed fashion by saving the stuff to files in process-sized
  ! chunks ???

!   ! time history at select points read from the input file
!   open  (unit = 1200, file = 'history', position = 'append')
!   do p = 1, num_points
!      l = ln(isave(1,p), isave(2,p), isave(3,p), 1)
!      write(unit = 1200, fmt = '(i2.2,6(1X,F15.8))') p, time, q(1,l), q(2,l), &
!                                             q(3,l), q(4,l), xnut(l)
!   end do
!   close (unit = 1200)

end subroutine output_time_history



