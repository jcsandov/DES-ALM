
  subroutine mg_bintp_mom
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Simple or Mass flux based interpolation methods
  ! Send and receive donor and donee informations
  ! J. Paik, Sep, 2004

  ! - - - - - - - - - - - - - - - - - - - - - - - - 

  implicit none

  integer :: ll
  integer :: ni, np, nz
  integer :: ihost, jhost, khost, nzhost
  integer :: nhost1, nhost2, nhost3, nhost4, nhost5, nhost6, nhost7, nhost8

  real (kind = rdf) :: rei
  real (kind = rdf) :: xco, yco, zco
  real (kind = rdf) :: v000, v100, v010, v001, v101, v011, v110, v111
  real (kind = rdf) :: c000, c100, c010, c001, c101, c011, c110, c111

  real (kind = rdf), dimension(:), allocatable :: sbuf
  real (kind = rdf), dimension(:), allocatable :: rbuf

  integer, dimension(0:nproc-1) :: sendcounts
  integer, dimension(0:nproc-1) :: recvcounts
  integer, dimension(0:nproc-1) :: sdispls
  integer, dimension(0:nproc-1) :: rdispls

  integer :: nosend
  integer :: norecv

  integer :: nme

  rei = one / ren

  nme = 4

  ! counts
  !
  nosend = 0
  norecv = 0
  do np = 0, nproc - 1
     sendcounts(np) = nme * donor(np)
     recvcounts(np) = nme * donee(np)
     nosend = nosend + sendcounts(np)
     norecv = norecv + recvcounts(np)
  end do

  ! displacements (start at zero like myid == root == 0)
  !
  sdispls(0) = 0
  rdispls(0) = 0
  do np = 1, nproc - 1
     sdispls(np) = sdispls(np-1) + sendcounts(np-1)
     rdispls(np) = rdispls(np-1) + recvcounts(np-1)
  end do

  ! send and receive information
  !
  allocate (sbuf(1:nosend))
  allocate (rbuf(1:norecv))

  sbuf = zero
  rbuf = zero
  
  ll = 0
  do np = 0, nproc - 1

     do ni = 1, donor(np)

        ihost = hosts(1,ni,np)
        jhost = hosts(2,ni,np)
        khost = hosts(3,ni,np)
        nzhost= hosts(4,ni,np)

        if (nzhost /= myzone) print*, 'something wrong host in mg_bintp_mom'

        xco = coef(1,ni,np)
        yco = coef(2,ni,np)
        zco = coef(3,ni,np)

        !nhost1 = gi_2_le_idx(ihost  ,jhost  ,khost  ,ns)
        !nhost2 = gi_2_le_idx(ihost+1,jhost  ,khost  ,ns)
        !nhost3 = gi_2_le_idx(ihost  ,jhost+1,khost  ,ns)
        !nhost4 = gi_2_le_idx(ihost  ,jhost  ,khost+1,ns)
        !nhost5 = gi_2_le_idx(ihost+1,jhost  ,khost+1,ns)
        !nhost6 = gi_2_le_idx(ihost  ,jhost+1,khost+1,ns)
        !nhost7 = gi_2_le_idx(ihost+1,jhost+1,khost  ,ns)
        !nhost8 = gi_2_le_idx(ihost+1,jhost+1,khost+1,ns)
        nhost1 = le_idx(ihost  ,jhost  ,khost  ,ns)
        nhost2 = le_idx(ihost+1,jhost  ,khost  ,ns)
        nhost3 = le_idx(ihost  ,jhost+1,khost  ,ns)
        nhost4 = le_idx(ihost  ,jhost  ,khost+1,ns)
        nhost5 = le_idx(ihost+1,jhost  ,khost+1,ns)
        nhost6 = le_idx(ihost  ,jhost+1,khost+1,ns)
        nhost7 = le_idx(ihost+1,jhost+1,khost  ,ns)
        nhost8 = le_idx(ihost+1,jhost+1,khost+1,ns)

        c000 = (one-xco)*(one-yco)*(one-zco)
        c100 =      xco *(one-yco)*(one-zco)
        c010 = (one-xco)*     yco *(one-zco)
        c001 = (one-xco)*(one-yco)*     zco
        c101 =      xco *(one-yco)*     zco
        c011 = (one-xco)*     yco *     zco
        c110 =      xco *     yco *(one-zco)
        c111 =      xco *     yco *     zco

        do m = 1, 4

           v000 = q(m,nhost1)
           v100 = q(m,nhost2)
           v010 = q(m,nhost3)
           v001 = q(m,nhost4)
           v101 = q(m,nhost5)
           v011 = q(m,nhost6)
           v110 = q(m,nhost7)
           v111 = q(m,nhost8)

           ll = ll + 1
           sbuf(ll) = v000*c000 + v100*c100 + v010*c010 + v001*c001 + &
                      v101*c101 + v011*c011 + v110*c110 + v111*c111
        end do

     end do

  end do

  if (ll /= nosend) print*, 'll /= nosend:', myid, myzone, ll, nosend

  ! Sends a distinct message from each task to every task.
  ! Messages can have different sizes and displacements.
  !
  call mpi_alltoallv (sbuf, sendcounts, sdispls, mpi_real, &
                      rbuf, recvcounts, rdispls, mpi_real, mpi_comm_world, ierr)

  ! receive interface information and interpolate data
  !
  ll = 0
  do np = 0, nproc - 1

     do ni = 1, donee(np)

        i = nodes(1,ni,np)
        j = nodes(2,ni,np)
        k = nodes(3,ni,np)
        nz= nodes(4,ni,np)

       !l = gi_2_le_idx(i,j,k,ns)
        l = le_idx(i,j,k,ns)

        if (nz /= myzone) print*, 'something wrong nz in mg_bintp_mom'

        do m = 1, 4

           ll  = ll + 1
           q(m,l) = rbuf(ll)

        end do

     end do

  end do

  ! exchnge variables
  !call mg_exchng3_2d (ns, q)

  if (ll /= norecv) print*, 'll /= norecv:', myid, myzone, ll, norecv

  deallocate ( sbuf, rbuf )
 
end subroutine mg_bintp_mom


