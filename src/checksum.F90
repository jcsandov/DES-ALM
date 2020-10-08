!
!
! Module : Checksum
!
!        for use with the gtce_cfd code
!
! a mpi aware module that provides bitwise-accurate checksums for
! vectors of various ranks. At present, this module is limited to
! single precision floating points numbers but it could easily be
! extended to double precision numbers or to integers
!
! S. Casey Jones <casey.jones@ce.gatech.edu>
!
! written : 18 Oct 2002
!
! future ideas
! ------------
!
! o use generic interface feature of f90 modules to simplify the
!   interface to the checksum routine, that is, regardless of the
!   rank of the var argument just have one call to the routine
!
module checksum
  use precision
  use global
  use global_mpi

  implicit none 

  private                       ! default to prevent namespace pollution

  public :: checksum_1d_par, &
            checksum_2d_par, &
            checksum_3d_par, &
            checksum_4d_par, &
            checksum_5d_par

  public :: cksum_4d_par, cksum_3d_par, cksum_2d_par, cksum_5d_par
  public :: init_cksum

  ! integer , parameter :: idf = selected_int_kind (R = 9) ! R = range

  real (kind = rdf) :: rsum
  real (kind = rdf) :: ssum

  integer (kind = idf) :: iv
  integer (kind = idf) :: iw
  integer (kind = idf) :: ix
  integer (kind = idf) :: iy
  integer (kind = idf) :: iz

  integer (kind = idf) :: isum
  integer (kind = idf) :: cksum
  
  integer, save :: ighost             ! i-direction # ghost points
  integer, save :: jghost             ! j-direction # ghost points
  integer, save :: kghost             ! k-direction # ghost points

  logical, save :: initialized = .false.

contains
  
  subroutine init_cksum (igp, jgp, kgp)

    integer, intent (in) :: igp
    integer, intent (in) :: jgp
    integer, intent (in) :: kgp

    ighost = igp
    jghost = jgp
    kghost = kgp

    !print *, 'init-cksum - 1', igp, jgp, kgp, ighost, jghost, kghost    

    initialized = .true.

  end subroutine init_cksum

  ! --

  subroutine print_result (comment)

    character (len = *)  :: comment
    character (len = 25) :: out

    ! '(a9,2x,a25,z,2x,g)', 'checksum:',comment, cksum, ssum

!    write (*, fmt = "('checksum:',1x,a,t25,z,t35,g20.10)" ) trim(comment), cksum, ssum

  end subroutine print_result

  ! --
  
  subroutine cksum_5d_par (comment, var)
    
    implicit none
    character (len = *) :: comment

    real (kind = rdf), dimension (:,:,:,:,:), intent(in)  :: var
    real (kind = rdf), dimension (:,:,:,:,:), allocatable :: vartmp

    integer :: mb
    integer :: ib
    integer :: jb
    integer :: kb
    integer :: nb
    
    ! check if init_checksum has been called to define
    ! number of ghost points in each direction
    ! 
    if (.not. initialized) then
       print *, 'cksum_par_4d not initialized'
       print *, 'call init_cksum(igp,jgp,kgp) before'
       print *, 'cksum_par_4d'
       return
    end if
    
    ! *** note ***
    !
    ! F90 standard default: dummy arrays default to lbound 1 unless
    ! otherwise specified; hence we could get away without ma,ia,ja,ka
    !

    ! assign limits 
    ! 
    nb = ubound(var,1)
    mb = ubound(var,2)
    ib = ubound(var,3)
    jb = ubound(var,4)
    kb = ubound(var,5)


    ! remove ghost cells
    !
    allocate (vartmp(1:nb,1:mb,1:ib-2*ighost,1:jb-2*jghost,1:kb-2*kghost))

    vartmp = var(1:nb,1:mb,1+ighost:ib-ighost,1+jghost:jb-jghost,1+kghost:kb-kghost)

    ! compute checksum
    ! 
    call checksum_5d_par (comment,  vartmp)

    ! free mem
    !
    deallocate (vartmp)

  end subroutine cksum_5d_par


  subroutine cksum_4d_par (comment, var)
    
    implicit none
    character (len = *) :: comment

    real (kind = rdf), dimension (:,:,:,:), intent(in)  :: var
    real (kind = rdf), dimension (:,:,:,:), allocatable :: vartmp

    integer :: mb
    integer :: ib
    integer :: jb
    integer :: kb
    
    ! check if init_checksum has been called to define
    ! number of ghost points in each direction
    ! 
    if (.not. initialized) then
       print *, 'cksum_par_4d not initialized'
       print *, 'call init_cksum(igp,jgp,kgp) before'
       print *, 'cksum_par_4d'
       return
    end if
    
    ! *** note ***
    !
    ! F90 standard default: dummy arrays default to lbound 1 unless
    ! otherwise specified; hence we could get away without ma,ia,ja,ka
    !

    ! assign limits 
    ! 
    mb = ubound(var,1)
    ib = ubound(var,2)
    jb = ubound(var,3)
    kb = ubound(var,4)

    ! remove ghost cells
    !
    allocate (vartmp(1:mb,1:ib-2*ighost,1:jb-2*jghost,1:kb-2*kghost))

    vartmp = var(1:mb,1+ighost:ib-ighost,1+jghost:jb-jghost,1+kghost:kb-kghost)

    ! compute checksum
    ! 
    call checksum_4d_par (comment,  vartmp)

    ! free mem
    !
    deallocate (vartmp)

  end subroutine cksum_4d_par

  ! --

  subroutine cksum_3d_par (comment, var)
    
    implicit none
    character (len = *) :: comment

    real (kind = rdf), dimension (:,:,:), intent(in)  :: var
    real (kind = rdf), dimension (:,:,:), allocatable :: vartmp

    integer :: ib
    integer :: jb
    integer :: kb
    
    ! check if init_checksum has been called to define
    ! number of ghost points in each direction
    ! 
    if (.not. initialized) then
       print *, 'cksum_par_3d not initialized'
       print *, 'call init_cksum(igp,jgp,kgp) before'
       print *, 'cksum_par_3d'
       return
    end if
    
    ! *** note ***
    !
    ! F90 standard default: dummy arrays default to lbound 1 unless
    ! otherwise specified; hence we could get away without ma,ia,ja,ka
    !

    ! assign limits 
    ! 
    ib = ubound(var,1)
    jb = ubound(var,2)
    kb = ubound(var,3)

    ! remove ghost cells
    !
    allocate (vartmp(1:ib-2*ighost,1:jb-2*jghost,1:kb-2*kghost))

    vartmp = var(1+ighost:ib-ighost,1+jghost:jb-jghost,1+kghost:kb-kghost)

    ! compute checksum
    ! 
    call checksum_3d_par (comment,  vartmp)

    ! free mem
    !
    deallocate (vartmp)

  end subroutine cksum_3d_par

! --

  subroutine cksum_2d_par (comment, var)
    
    implicit none
    character (len = *) :: comment

    real (kind = rdf), dimension (:,:), intent(in)  :: var
    real (kind = rdf), dimension (:,:), allocatable :: vartmp

    integer :: jb
    integer :: kb
    
    ! check if init_checksum has been called to define
    ! number of ghost points in each direction
    ! 
    if (.not. initialized) then
       print *, 'cksum_par_3d not initialized'
       print *, 'call init_cksum(igp,jgp,kgp) before'
       print *, 'cksum_par_3d'
       return
    end if
    
    ! *** note ***
    !
    ! F90 standard default: dummy arrays default to lbound 1 unless
    ! otherwise specified; hence we could get away without ma,ia,ja,ka
    !

    ! assign limits 
    ! 
    jb = ubound(var,1)
    kb = ubound(var,2)

    ! remove ghost cells
    !
    allocate (vartmp(1:jb-2*jghost,1:kb-2*kghost))

    vartmp = var(1+jghost:jb-jghost,1+kghost:kb-kghost)

    ! compute checksum
    ! 
    call checksum_2d_par (comment,  vartmp)

    ! free mem
    !
    deallocate (vartmp)

  end subroutine cksum_2d_par

!--

  subroutine checksum_1d_par (comment, var)
    implicit none
    character (len = *) :: comment
    real (kind = rdf), dimension (:), intent(in)  :: var
    integer (kind = idf), dimension (:), allocatable :: mold

    cksum = 0
    rsum = zero
    ssum = zero

    ix = size(var,1)

    allocate (mold(ix))
    
    isum = sum(transfer(var(1:ix), mold, size(mold)))

    call mpi_reduce (isum, cksum, 1, mpi_integer4, MPI_SUM&
         &, root, mpi_comm_world, ierr)

    ! used to show value of the checksum that is not bitwise accurate
    ! 
    rsum = sum(var(1:ix))
    
    call mpi_reduce(rsum, ssum, 1, mpi_real4, mpi_sum, root,&
         & mpi_comm_world, ierr)

    if (myid == root) call print_result (comment)

    deallocate(mold)

  end subroutine checksum_1d_par

  ! --

  subroutine checksum_2d_par (comment, var)
    implicit none
    character (len = *) :: comment
    real (kind = rdf), dimension (:,:), intent(in)  :: var
    integer (kind = idf), dimension (:), allocatable :: mold


    cksum = 0
    rsum = zero
    ssum = zero

    ix = size(var,1)
    iy = size(var,2)

    allocate (mold(ix*iy))
    
    isum = sum(transfer(var(1:ix,1:iy), mold, size(mold)))

    call mpi_reduce (isum, cksum, 1, mpi_integer4, MPI_SUM&
         &, root, mpi_comm_world, ierr)

    ! used to show value of the checksum that is not bitwise accurate
    ! 
    rsum = sum(var(1:ix,1:iy))
    
    call mpi_reduce(rsum, ssum, 1, mpi_real4, mpi_sum, root,&
         & mpi_comm_world, ierr)

    if (myid == root) call print_result (comment)

    deallocate(mold)

  end subroutine checksum_2d_par

  ! --

  subroutine checksum_3d_par (comment, var)
    implicit none
    character (len = *) :: comment
    real (kind = rdf), dimension (:,:,:), intent(in)  :: var
    integer (kind = idf), dimension (:), allocatable :: mold

    cksum = 0
    rsum = zero
    ssum = zero

    ix = size(var,1)
    iy = size(var,2)
    iz = size(var,3)

    allocate (mold(ix*iy*iz))
    
    isum = sum(transfer(var(1:ix,1:iy,1:iz), mold, size(mold)))

    call mpi_reduce (isum, cksum, 1, mpi_integer4, MPI_SUM&
         &, root, mpi_comm_world, ierr)

    ! used to show value of the checksum that is not bitwise accurate
    ! 
    rsum = sum(var(1:ix,1:iy,1:iz))
    
    call mpi_reduce(rsum, ssum, 1, mpi_real4, mpi_sum, root,&
         & mpi_comm_world, ierr)

    if (myid == root) call print_result(comment)

    deallocate(mold)

  end subroutine checksum_3d_par

  ! --

  subroutine checksum_4d_par (comment, varr)
    implicit none
    character (len = *) :: comment
    real (kind = rdf), dimension (:,:,:,:), intent(in)  :: varr
    integer (kind = idf), dimension (:), allocatable :: mold

    cksum = 0
    rsum = zero
    ssum = zero

    ix = size(varr,1)
    iy = size(varr,2)
    iz = size(varr,3)
    iw = size(varr,4)

    allocate (mold(ix*iy*iz*iw))
    
    isum = sum(transfer(varr(1:ix,1:iy,1:iz,1:iw), mold, size(mold)))

    call mpi_reduce (isum, cksum, 1, mpi_integer4, MPI_SUM&
         &, root, mpi_comm_world, ierr)

    ! used to show value of the checksum that is not bitwise accurate
    ! 
    rsum = sum(varr(1:ix,1:iy,1:iz,1:iw))
    
    call mpi_reduce(rsum, ssum, 1, mpi_real4, mpi_sum, root,&
         & mpi_comm_world, ierr)

    if (myid == root) call print_result (comment)

    deallocate(mold)

  end subroutine checksum_4d_par

  subroutine checksum_5d_par (comment, varr)
    implicit none
    character (len = *) :: comment
    real (kind = rdf), dimension (:,:,:,:,:), intent(in)  :: varr
    integer (kind = idf), dimension (:), allocatable :: mold

    cksum = 0
    rsum = zero
    ssum = zero

    ix = size(varr,1)
    iy = size(varr,2)
    iz = size(varr,3)
    iw = size(varr,4)
    iv = size(varr,5)

    allocate (mold(ix*iy*iz*iw*iv))
    
    isum = sum(transfer(varr(1:ix,1:iy,1:iz,1:iw,1:iv), mold, size(mold)))

    call mpi_reduce (isum, cksum, 1, mpi_integer4, MPI_SUM&
         &, root, mpi_comm_world, ierr)

    ! used to show value of the checksum that is not bitwise accurate
    ! 
    rsum = sum(varr(1:ix,1:iy,1:iz,1:iw,1:iv))
    
    call mpi_reduce(rsum, ssum, 1, mpi_real4, mpi_sum, root,&
         & mpi_comm_world, ierr)

    if (myid == root) call print_result (comment)

    deallocate(mold)

  end subroutine checksum_5d_par


end module checksum

