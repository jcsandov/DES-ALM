module global_param

  use precision
  implicit none

  ! Reynolds number
  !
  real (kind = rdf) :: ren

  ! C_des constant & working variable limit in des
  !
  real (kind = rdf) :: cdes, cvin

  ! k-epslion & k-omega limits
  ! 
  real (kind = rdf) :: ckin, cein

  ! (real) time step
  !
  real (kind = rdf) :: delti
  real (kind = rdf) :: time0
  real (kind = rdf) :: time
  integer :: nt1, nt2

  ! Runge-Kutta parameters
  !
  integer :: irk
  real (kind = rdf) :: alfa(4)

  ! parameter for unsteady term
  !
  real (kind = rdf) :: e_source

  ! local time-stepping
  !
  real (kind = rdf), dimension(:), allocatable :: cfl1, vnn1, cfl2, vnn2

  ! artificial compressibility
  !
  real (kind = rdf) :: beta

  ! pressure dissipation (continuity eq)
  !
  real (kind = rdf) :: pdiss_coef

  ! matrix dissipation coefficients
  !
  real (kind = rdf), dimension(4) :: eps

  ! residual smoothing coefficients
  !
  real (kind = rdf), dimension(:,:,:), allocatable :: ep

  ! iteration information
  !
  integer :: it, itmax, itk

  ! minimum number of pseudo-iterations per time step
  !
  integer :: it_min

  ! monitoring of the program
  !
  integer :: monitor_num_points ! read in input program

  ! hardwired variable
  !
  integer, dimension(10,4) :: monitor_point_ijk  

  integer :: checkpoint ! = 10 ! should be in input file

  ! convergence tolerances
  ! 
  real (kind = rdf) :: er_min, eo_min

  ! particle program
  !
  integer :: partc_skip

  ! unsteady ?
  !
  logical :: unsteady = .true.

  ! solution method in time
  !
  logical :: rkm = .false.  ! multistage Runge-Kutta method
  logical :: daf = .true.   ! diagonalized version of 
                            ! approximate factorization

  ! switch for jet, cavity , abutment or particles
  !
  logical :: cavity = .false.
  logical :: jet    = .false.
  logical :: abut   = .false.
  logical :: duct   = .false.
  logical :: draft  = .false.

  ! temp ; convective terms
  ! 
  logical :: quick = .true.

  ! for initial guess
  logical :: upwind = .false.

  ! more global info!
  !
  integer :: itc, itr
  integer :: icn
  integer :: icnw
  integer :: ntime

  real (kind = rdf) :: time_step_time
  real (kind = rdf) :: total_time

      
  ! Turbulence models
  ! -----------------
  logical :: turbulence = .true. ! turbulence switch
                                 ! --> used in init and other places
  logical :: les = .false.
  logical :: des = .true.
  logical :: s_a = .false.
  logical :: k_w = .false.       ! k-omega turbulence model
  logical :: k_e = .false.       ! k-epsilon turbulence model

  logical :: nlinc = .false.      ! Non-linear closure of turbulence model
  logical :: craft = .false.      ! Craft-Launder k-w nonlinear turbulence model

  ! consider roughness? (Aupoix & Spalart 2003, IJHFF)
  logical :: rough_wall = .false.
  real (kind = rdf), dimension(:,:), allocatable :: hsw

  ! Random forcing generation
  ! -------------------------
  logical :: rfgt = .true.


end module global_param





