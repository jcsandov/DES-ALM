module rfg

  use precision
  !use global
  !use global_param
  !use global_mpi
  !use global_turb
  !use global_app
  !use global_osg


  implicit none

  ! public routines
  public :: initialize_rfg, rfg_vel

  ! everything else (subroutines, function, variables, parameters)
  ! are private by default
  private

  ! type parameters
  integer, parameter :: k4b = selected_int_kind(9)

  ! external variables

  ! for random number generator
  ! (record if you need to generate same sample)
  integer (kind = k4b), save :: seed 

  ! total number of wavenumbers to sample
  integer, save, private :: N

  real (kind = rdf), allocatable, save, private :: p_rfg(:,:)
  real (kind = rdf), allocatable, save, private :: q_rfg(:,:)
  real (kind = rdf), allocatable, save, private :: k_rfg(:,:)
  real (kind = rdf), allocatable, save, private :: omega_rfg(:)

  logical, save :: initialized = .false.

  contains
    
    subroutine initialize_rfg (wavenumbers, iseed)

      integer (kind = k4b) :: iseed

      real (kind = rdf), allocatable :: zet_rfg(:,:)
      real (kind = rdf), allocatable :: csi_rfg(:,:)

      integer :: wavenumbers
      integer :: i 

      ! N determine by calling routine
      N = wavenumbers
      
      ! seed converted
      seed = iseed

      ! NR random number generator requires
      ! negative integer for seed
      if (iseed >= 0) then
         write (*,*) 'NR random number generator requires '
         write (*,*) 'a negative integer for seed value'
      end if
      
      ! allocate arrays
      allocate ( p_rfg(3,N), &
                 q_rfg(3,N), &
                 k_rfg(3,N), &
               omega_rfg(N), &
               zet_rfg(3,N), &
               csi_rfg(3,N) )

      ! obtain pseudo-random numbers
      call get_random_numbers (zet_rfg(1,:), zero, one)
      call get_random_numbers (zet_rfg(2,:), zero, one)
      call get_random_numbers (zet_rfg(3,:), zero, one)
      call get_random_numbers (csi_rfg(1,:), zero, one)
      call get_random_numbers (csi_rfg(2,:), zero, one)
      call get_random_numbers (csi_rfg(3,:), zero, one)
      call get_random_numbers (k_rfg(1,:), zero, pt5)
      call get_random_numbers (k_rfg(2,:), zero, pt5)
      call get_random_numbers (k_rfg(3,:), zero, pt5)
      call get_random_numbers (omega_rfg(:), zero, one)

      ! perform cross products
      do i = 1, N
         p_rfg(:,i) = cross_product(zet_rfg(:,i), k_rfg(:,i))
         q_rfg(:,i) = cross_product(csi_rfg(:,i), k_rfg(:,i))
      end do

      ! zet and csi are no longer needed
      deallocate (zet_rfg, csi_rfg)      

      initialized = .true.

    end subroutine initialize_rfg
    
    subroutine get_random_numbers(vec, mean, stdev)

      ! routine needed because my random number generator
      ! can't be called with a vector

      integer :: i

      real (kind = rdf), dimension(N) :: vec
      real (kind = rdf) :: mean
      real (kind = rdf) :: stdev

      do i = 1, N
         vec(i) = box_muller(seed, mean, stdev)
      end do

    end subroutine get_random_numbers
    
    function cross_product (a, b) result (d)

      real (kind = rdf), dimension(3) :: d
      real (kind = rdf), intent(in),  dimension(3) :: a, b

      d(1) = a(2) * b(3) - a(3) * b(2)
      d(2) = a(3) * b(1) - a(1) * b(3)
      d(3) = a(1) * b(2) - a(2) * b(1)

    end function cross_product
    
    function box_muller(idum, mean, stdev)
    !======================================================================
    !
    ! function : box_muller
    ! ========
    !
    ! returns a normally distributed random deviate with mean
    ! and standard deviation (stdev) using the box_muller
    ! transformation
    !
    ! from Numerical Receipes in Fortran (2ed) p 280
    !
    !======================================================================    

      real (kind = rdf) :: box_muller

      integer (kind = k4b), intent(inout) :: idum
      integer, save :: iset = 0

      real (kind = rdf), intent(in) :: mean
      real (kind = rdf), intent(in) :: stdev

      real (kind = rdf), save :: gset
      real (kind = rdf) :: fac
      real (kind = rdf) :: rsq
      real (kind = rdf) :: v1
      real (kind = rdf) :: v2


!!$      print *, 'idum =', idum

      if (iset == 0) then
1        v1 = two * ran(idum) - one
         v2 = two * ran(idum) - one
         rsq = v1**2 + v2**2
         if (rsq >= one .or. rsq == zero) goto 1
         fac = sqrt(-two * log(rsq) / rsq)
         gset = v1 * fac
         box_muller = v2 * fac
         iset = 1
      else
         box_muller = gset
         iset = 0
      end if

      box_muller = mean + stdev * box_muller

    end function box_muller
    

    function ran(idum)
    !======================================================================
    ! function : ran
    ! ========
    !
    ! source   : p. 1142, Chap. B7, Numerical Receipes for Fortran 90
    ! ======
    !            http://www.ulib.org/webRoot/Books/Numerical_Recipes \
    !                               /bookf90pdf/chap7f9.pdf
    !
    !            also see http://www.nr.com
    !
    ! Description
    ! ===========
    !
    ! "Minimal" random number generator of Park and Miller combined with a
    ! Marsaglia shift sequence. Returns a uniform random deviate between 
    ! 0.0 and 1.0 (exclusive of the endpoint values). This fully portable,
    ! scalar generator has the "traditional" (not Fortran 90) calling 
    ! sequence with a random deviate as the returned function value: call 
    ! with idum a negative integer to initialize; thereafter, do not alter 
    ! idum expect to reinitialize. The period of this generator is
    ! about 3.1 x 10^18.
    ! 
    !======================================================================

      real (kind = rdf) :: ran

      integer (kind = k4b), intent(inout) :: idum

      integer (kind = k4b), parameter :: ia = 16807
      integer (kind = k4b), parameter :: im = 2147483647
      integer (kind = k4b), parameter :: iq = 127773
      integer (kind = k4b), parameter :: ir = 2836

      real (kind = rdf), save :: am

      integer (kind = k4b), save :: ix = -1
      integer (kind = k4b), save :: iy = -1
      integer (kind = k4b), save :: k

      if (idum <= 0 .or. iy < 0) then   ! initialize
         am = nearest(1.0, -1.0) / im
         iy = ior(ieor(888889999, abs(idum)), 1)
         ix = ieor(777755555, abs(idum))
         idum = abs(idum) + 1           ! set idum positive
      end if

      ix = ieor(ix, ishft(ix,  13))    ! Marsaglia shift sequence with
      ix = ieor(ix, ishft(ix, -17))    ! period 2^32 - 1
      ix = ieor(ix, ishft(ix,   5))

      k = iy / iq                       ! Park-Miller sequence by
      iy = ia * (iy - k * iq) - ir * k  ! Schrage's method, period 2^31 - 1

      if (iy < 0) iy = iy + im

      ran = am * ior(iand(im, ieor(ix,iy)), 1)   ! combine the two generators
                                                 ! with masking to ensure
                                                 ! nonzero value
    end function ran

    function get_iso_velocity (t, x, c, upsilon) result (v)
      
      ! we can scale x() & t before they enter; however, we can't
      ! scale k() until we call routine because scale_c may
      ! change in space and we must preserve k() for all
      ! subsequently generated velocities

      real (kind = rdf), dimension(3) :: v

      real (kind = rdf), dimension(3), intent(in) :: x
      real (kind = rdf), dimension(3), intent(in) :: c
      real (kind = rdf), intent(in) :: upsilon
      real (kind = rdf), intent(in) :: t

      ! local variables
      real (kind = rdf), dimension(3) :: sum1
      real (kind = rdf), dimension(3,N) :: k_scale
      real (kind = rdf) :: term1
      real (kind = rdf) :: const
      
      integer :: i

      ! confirm that we've initialized spectrum generation
      if (.not. initialized) then
         print *, ' Problem ... '
         print *, ''
         print *, 'initialize spectrum via initialize_spectrum'
         print *, 'before calling get_spectrum'
         stop
      end if

      ! scale wave numbers
      do i = 1, N
         k_scale(:,i) = k_rfg(:,i) * upsilon / c(:)
      end do

      ! perform summation using f90 array notation
      sum1(:) = zero
      do i = 1, N
         term1 = dot_product(k_scale(:,i), x(:)) + omega_rfg(i) * t
         sum1(:) = sum1(:) + p_rfg(:,i) * cos(term1) + q_rfg(:,i) * sin(term1)
      end do
      
      const = sqrt(two / real(N, kind = rdf))
      v(:) = const * sum1(:)
      
    end function get_iso_velocity
    
    subroutine diagonalize (r, m, c)
      ! uses lapack
      real (kind = rdf), intent(in), dimension(3,3) :: r
      real (kind = rdf), intent(out), dimension(3) :: c
      real (kind = rdf), intent(out), dimension(3,3) :: m

      real (kind = rdf), dimension(3) :: b
      real (kind = rdf), dimension(200) :: work

      integer :: ok

      ! call correct routine based on single
      ! or double precision
      if (rdf == 4) then
         call ssyev ('V', 'U', 3, r, 3, b, work, 102, ok)
      else if (rdf == 8) then
         call dsyev ('V', 'U', 3, r, 3, b, work, 102, ok)
      else
         print *, 'precision problem with lapack'
      end if

      ! error check
      if (ok /= 0) print *, 'Lapack returned bad result'

      m(:,:) = r(:,:)
      c(:) = sqrt(b(:))

    end subroutine diagonalize
    
    subroutine anti_diagonalize (m, c, upsilon, v)
      
      real (kind = rdf), dimension(3,3), intent(in) :: m
      real (kind = rdf), dimension(3), intent(in) :: c
      real (kind = rdf), intent(in) :: upsilon
      real (kind = rdf), dimension(3), intent(inout) :: v      
      
      ! see eqs. 8 & 9 of Smirnov et al. 2001
!!$      v(:) = c(:) * upsilon(:) * v(:)
      v(:) = c(:) * v(:)
      v(:) = matmul(m(:,:), v(:))

    end subroutine anti_diagonalize
    
    function rfg_vel (t, x, tau, eta, r) result (v)
      ! an external function to generate random velocity
      ! at a given point in space and time
      
      ! variables for external program
      ! 
      ! t    := time
      ! x    := location in space, vector, x(3)
      ! tau  := turbulence time scale
      ! eta  := turbulence length scale, vector, eta(3)
      ! r    := velocity correlation tensor, array, r(3,3)
      !         but only upper part used (because symmetric)

      real (kind = rdf), intent(in) :: t, tau, eta
      real (kind = rdf), intent(in), dimension(3) :: x
      real (kind = rdf), intent(in), dimension(3,3) :: r

      real (kind = rdf), dimension(3) :: v

      ! local variables
      !
      ! c       := fluctuations in diagonalized coordinates, vector, c(3)
      ! upsilon := turbulence velocity scale, vector upsilon(3)
      !         :=

      real (kind = rdf), dimension(3) :: c
      real (kind = rdf) :: upsilon

      real (kind = rdf), dimension(3,3) :: m
      real (kind = rdf) :: t1
      real (kind = rdf), dimension(3) :: x1, v1
      
      call diagonalize (r, m, c)

      ! scale variables
      t1 = t / tau
      x1(:) = x(:) / eta
      upsilon = eta / tau
      
      v1(:) = get_iso_velocity(t1, x1, c, upsilon)

      call anti_diagonalize (m, c, upsilon, v1)

      v(:) = v1(:)

    end function rfg_vel
    

  end module rfg

