module rfg_inlet
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! module  : rfg_inlet
  !
  ! purpose : provide random, turbuence based forcing for inlet
  !           of yaras duct flow
  !
  ! date    : 12 March 2002
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Force the inlet of the duct using the random flow generator
  ! module, which was written by S Casey Jones in March 2002 based
  ! on the following paper:
  !
  !    A Smirnov, S. Shi, & I. Celik
  !    Random flow generator technique for large eddy simulations
  !      and particle dynamics modeling
  !    ASME J. Fluids Engineering
  !    v. 123, June 2001, pp. 359-371
  !    
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use precision
  use rfg       
  use global
  use global_mpi
  use global_param

  implicit none

  public :: rfg_init, rfg_mpi_inlet_forcing

  private
  ! rfg module provides
  !
  ! subroutine : rfg_init
  ! function   : rfg_vel

  ! inlet flow field properties
  real (kind = rdf), allocatable, dimension(:,:,:), save :: rfg_um
  real (kind = rdf), allocatable, dimension(:,:,:), save :: rfg_uij
  real (kind = rdf), allocatable, dimension(:,:), save   :: rfg_t_ls
  real (kind = rdf), allocatable, dimension(:,:), save   :: rfg_t_ts

contains

  subroutine rfg_init()
    integer :: seed
    integer :: samples

    ! y, z in binary input file, not needed
    real (kind = rdf), allocatable, dimension(:,:) :: rfg_y, rfg_z

    ! allocate variables
    allocate (rfg_y(jmg(1,1),kmg(1,1)), &
              rfg_z(jmg(1,1),kmg(1,1)))

    allocate (rfg_um(3,jmg(1,1),kmg(1,1)),  &
              rfg_uij(6,jmg(1,1),kmg(1,1)), &
              rfg_t_ts(jmg(1,1),kmg(1,1)),  &
              rfg_t_ls(jmg(1,1),kmg(1,1)))

    if ( myid == root ) then

       ! get flow field info (read arrays in directly)
       open  (unit = 40, file = 'inlet.bin.dat', form = 'unformatted')
       read  (unit = 40) rfg_y !coordenadas y?
       read  (unit = 40) rfg_z !coordenadas z?
       read  (unit = 40) rfg_um !velocidades medias?
       read  (unit = 40) rfg_t_ts !escala de tiempo?
       read  (unit = 40) rfg_t_ls !escalas de longitud?
       read  (unit = 40) rfg_uij !correlaciones?
       close (unit = 40)

       ! get seed and number of random numbers to
       ! use to sample the isotropic/homogeneous
       ! turbulence spectrum
       open  (unit = 45, file = 'rfg_inlet_seed.dat', form = 'formatted')
       read  (unit = 45, fmt = *) samples !Numero de series (N)
       read  (unit = 45, fmt = *) seed !¿?
       close (unit = 45)

    end if

    ! broadcast to all nodes
    ! 
    call mpi_bcast (rfg_um, 3*jmg(1,1)*kmg(1,1), mpi_real, root,&
         & mpi_comm_world, ierr)
    call mpi_bcast (rfg_t_ts, jmg(1,1)*kmg(1,1), mpi_real, root,&
         & mpi_comm_world, ierr)
    call mpi_bcast (rfg_t_ls, jmg(1,1)*kmg(1,1), mpi_real, root,&
         & mpi_comm_world, ierr)
    call mpi_bcast (rfg_uij, 6*jmg(1,1)*kmg(1,1), mpi_real, root,&
         & mpi_comm_world, ierr)
    call mpi_bcast (samples, 1, mpi_integer, root, mpi_comm_world,&
         & ierr)
    call mpi_bcast (seed, 1, mpi_integer, root, mpi_comm_world, ierr)

    ! initialize rfg
    !
    call initialize_rfg (samples, seed) 

    deallocate (rfg_y, rfg_z)

  end subroutine rfg_init

  subroutine rfg_mpi_inlet_forcing ()

    ! local variables
    real (kind = rdf), dimension(3) :: x_rfg    
    real (kind = rdf), dimension(3) :: v_rfg
    real (kind = rdf), dimension(3,3) :: corr
    real (kind = rdf) :: ts_rfg, ls_rfg

    integer :: i, j, k, l

    integer :: j_mysta, j_myend
    integer :: k_mysta, k_myend

    ! must handle the boundaries correctly
    j_mysta = gi_ja(1)
    j_myend = gi_jb(1)

    k_mysta = gi_ka(1)
    k_myend = gi_kb(1)

    if (myleft == mpi_proc_null)  j_mysta = gi_ja(1) + 1
    if (myright == mpi_proc_null) j_myend = gi_jb(1) - 1

    if (mydown == mpi_proc_null)  k_mysta = gi_ka(1) + 1
    if (myup    == mpi_proc_null) k_myend = gi_kb(1) - 1

    if (gi_ia(1).eq.1) then !solo forzar la entrada si el procesador contiene el inlet
       do k = k_mysta, k_myend
       do j = j_mysta, j_myend

          l = gi_2_le_idx(1,j,k,1)

          x_rfg(1) = x(l)
          x_rfg(2) = y(l)
          x_rfg(3) = z(l)
          
          ! turbulence time scale
          !
          ts_rfg = rfg_t_ts(j,k)

          ! turbulence length scale
          !
          ls_rfg = rfg_t_ls(j,k)

          ! turbulence reynolds stress tensor
          !
          corr(1,1) = rfg_uij(1,j,k)
          corr(1,2) = rfg_uij(2,j,k)
          corr(1,3) = rfg_uij(3,j,k)
          corr(2,2) = rfg_uij(4,j,k)
          corr(2,3) = rfg_uij(5,j,k)
          corr(3,3) = rfg_uij(6,j,k)

          v_rfg(:) = rfg_vel(time, x_rfg(:), ts_rfg, ls_rfg, corr(:,:))
          q(2:4,l) = rfg_um(:,j,k) + v_rfg(:)
        
        ! if ((j.eq.181) .and. (k.eq.61)) then
        !   open (unit = 1100, file = 'fluctcenter.dat', position = 'append')
        !   write (unit =1100, fmt = '(3(1x,G20.7))') v_rfg(:)
        !   close(unit = 1100)
        ! end if

       end do
       end do
     end if

  end subroutine rfg_mpi_inlet_forcing

end module rfg_inlet

