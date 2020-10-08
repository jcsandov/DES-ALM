
module global_mpi

  use precision

  include 'mpif.h'

  !use mpi   !--->won't work in pgi and ScaliMPI
  

  integer, parameter :: root = 0

  integer :: myid
  integer :: nproc
  integer :: ierr
  integer :: status(MPI_STATUS_SIZE)
  integer :: errcode

  ! 3-D domain decomposition
  !integer :: comm3d
  integer, dimension(:), allocatable :: comm3d

  ! neighoring processes
  ! 
  integer :: myleft
  integer :: myright
  integer :: myup
  integer :: mydown
  integer :: myfront
  integer :: myback

  ! overset grid
  !
  integer, dimension(:), allocatable :: comm_nz
  integer :: myzone
  integer, dimension(:),   allocatable :: nproc_nz

  ! integer variables for domain decomosition
  integer, dimension(:,:), allocatable :: dims
  integer, dimension(:,:), allocatable :: coords

  ! logical variables for domain decomosition
  logical, dimension(:,:), allocatable :: periods
  logical, dimension(:,:), allocatable :: reorder



contains

  subroutine init_mpi

    call mpi_init (ierr)
    call mpi_comm_rank (mpi_comm_world, myid,  ierr)
    call mpi_comm_size (mpi_comm_world, nproc, ierr)

  end subroutine init_mpi

  ! --

  subroutine test_mpi

    print *, 'myid is ', myid, 'of', nproc, 'processors'

  end subroutine test_mpi
  
  ! --

  subroutine mpe_decomp1d (n, numprocs, myid, s, e)

    ! This file contains a routine for producing a decomposition of a
    ! 1-d array when given a number of processors.  It may be used in
    ! "direct" product decomposition.  The values returned assume
    ! a "global" domain in [1:n]

    ! From the book:
    ! 
    ! "Using MPI", Gropp, Lusk, & Skjellum
    ! 2nd edition, MIT Press, 1999

    ! global variables
    !
    integer, intent(in) :: n
    integer, intent(in) :: numprocs
    integer, intent(in) :: myid
    integer, intent(out) :: s
    integer, intent(out) :: e

    ! local variables
    ! 
    integer :: nlocal
    integer :: deficit

    nlocal = n / numprocs
    s = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s = s + min(myid,deficit)
    if (myid < deficit) then
       nlocal = nlocal + 1
    end if
    
    e = s + nlocal - 1
    if (e > n .or. myid == numprocs-1) e = n

  end subroutine mpe_decomp1d

  ! --

  subroutine length_1d(n1, n2, nproc, l)

    ! global variables
    !
    integer, intent(in) :: n1
    integer, intent(in) :: n2
    integer, intent(in) :: nproc

    integer, intent(out), dimension(0:nproc-1) :: l

    ! local variables
    !
    integer :: iwork1
    integer :: iwork2
    integer :: n

    iwork1 = (n2 - n1 + 1) / nproc
    iwork2 = mod(n2 - n1 + 1, nproc)
    do n = 0, nproc-1
       l(n) = iwork1
       if (n < iwork2) l(n) = l(n) + 1
    end do
       
  end subroutine length_1d
  
  ! --

  subroutine finish_mpi

    call mpi_finalize (ierr)

  end subroutine finish_mpi
  
  ! --

  subroutine checksum_1p (comment, var)

    ! compute checksum for a 3d array on a single processor

    implicit none

    character (len = *) :: comment

    real (kind = rdf), dimension (:,:,:), intent(in) :: var
    real (kind = rdf), dimension (:,:,:), allocatable :: vartmp

    real (kind = rdf) :: rsum, ssum

    integer (kind = 4), dimension(:), allocatable :: mold

    integer (kind = 4) :: isum
    integer (kind = 4) :: cksum

    integer :: ix, iy, iz

    cksum = 0
    ssum = zero

    ix = size(var,1)
    iy = size(var,2)
    iz = size(var,3)

    allocate (vartmp(1:ix,1:iy,1:iz), mold(ix*iy*iz))

    cksum = sum(transfer(var(1:ix,1:iy,1:iz),mold,size(mold)))
    ssum  = sum(var(1:ix,1:iy,1:iz))

    !print '(a,z,2x,g)', 'Checksum: '//comment, cksum, ssum

    deallocate(vartmp, mold)

  end subroutine checksum_1p

end module global_mpi
