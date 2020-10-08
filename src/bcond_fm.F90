
subroutine bcond_fm (il, iu, jl, ju, kl, ku, igp, jgp, kgp, q, gi_ia, gi_ja, gi_ka)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! general version !Kim's curved duct (full 3d geometry)
  ! 
  ! fine-mesh
  ! boundary conditions
  ! 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  use global_app
  use global_mpi
  use global_turb

  implicit none
  
  ! extents
  ! 
  integer :: il
  integer :: jl
  integer :: kl
  integer :: iu
  integer :: ju
  integer :: ku

  ! ghost points
  ! 
  integer :: igp
  integer :: jgp
  integer :: kgp

  ! solution vector
  ! 
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku) :: q

  ! fixed pressure
  ! 
  real (kind = rdf) :: pfix

  integer :: i, j, k, b
  
  integer :: i_mysta
  integer :: j_mysta
  integer :: k_mysta

  integer :: i_myend
  integer :: j_myend
  integer :: k_myend

  integer :: n_t, counter, number_disk_nodes
  real (kind = rdf) :: radius

  integer ::  gi_ia, gi_ja, gi_ka, gi_jb, gi_kb, i_t

  ! make definition of i_mysta etc., consistent with
  ! all other uses, e.g., solver_daf & rhs_kw_bcond
  !
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! processes on the domain boundaries
  ! 
  !if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  !if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  !if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  !if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  !if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  !if (myup    == mpi_proc_null) k_myend = ku - kgp - 1


  ! boundary in csi-direction (i = 1)
  !
  if (myback == mpi_proc_null) then
     b = 1
     i = i_mysta

     csi1: select case(btype(b,myzone))
     case(0) ! -> interface

     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i+1,j,k) + &
                          sb(1:4,b) * q(1:4,i+2,j,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1,i,j,k) = (one + rat(j,k,b)) * q(1,i+1,j,k) - &
                               rat(j,k,b)  * q(1,i+2,j,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i+1,j,k) + &
           !             sb(1,b) * q(1,i+2,j,k)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        !do k = k_mysta, k_myend
        !do j = j_mysta, j_myend
        !   q(1:4,i,j,k) = q(1:4,i_myend,j,k)
        !end do
        !end do
     end select csi1
  end if
     
  ! boundary in csi-direction (i = imax)
  !
  if (myfront == mpi_proc_null) then
     b = 2
     i = i_myend

     csi2: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i-1,j,k) + &
                          sb(1:4,b) * q(1:4,i-2,j,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1,i,j,k) = (one + rat(j,k,b)) * q(1,i-1,j,k) - &
                               rat(j,k,b)  * q(1,i-2,j,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i-1,j,k) + &
           !             sb(1,b) * q(1,i-2,j,k)
        end do
        end do
     case(5) ! exit (MOC)
		!j = j_mysta
        !do k = k_mysta, k_myend
		!   q( 1 ,i,j,k) = q( 1 ,i-1,j,k)
        !   q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
        !                  sb(2:4,b) * q(2:4,i-2,j,k)
		!end do

		!j = j_myend
        !do k = k_mysta, k_myend
		!   q( 1 ,i,j,k) = q( 1 ,i-1,j,k)
        !   q(2:4,i,j,k) = sa(2:4,b) * q(2:4,i-1,j,k) + &
        !                  sb(2:4,b) * q(2:4,i-2,j,k)
		!end do

     case(6) ! periodic condition
        !do k = k_mysta, k_myend
        !do j = j_mysta, j_myend
        !   q(1:4,i,j,k) = q(1:4,i_mysta,j,k)
        !end do
        !end do

case(8) ! free surface
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           q(1:4,i,j,k) = q(1:4,i-1,j,k) !sa(1:4,b)*q(1:4,i-1,j,k) + sb(1:4,b)*q(1:4,i-2,j,k)

        end do
        end do


     end select csi2
  end if
     
  ! boundary in eta-direction (j = 1)
  ! 
  if (myleft == mpi_proc_null) then
     b = 3
     j = j_mysta

     eta1: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j+1,k) + &
                          sb(1:4,b) * q(1:4,i,j+2,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,k,b)) * q(1,i,j+1,k) - &
                               rat(i,k,b)  * q(1,i,j+2,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j+1,k) + &
           !             sb(1,b) * q(1,i,j+2,k)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j_myend,k)
        end do
        end do
     end select eta1
  end if

  ! boundary in eta-direction (j = jmax)
  ! 
  if (myright == mpi_proc_null) then
     b = 4
     j = j_myend

     eta2: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j-1,k) + &
                          sb(1:4,b) * q(1:4,i,j-2,k)
        end do
        end do
     case(4) ! inflow
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,k,b)) * q(1,i,j-1,k) - &
                               rat(i,k,b)  * q(1,i,j-2,k)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j-1,k) + &
           !             sb(1,b) * q(1,i,j-2,k)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do k = k_mysta, k_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j_mysta,k)
        end do
        end do
     end select eta2
  end if

  ! boundary in zet-direction (k = 1)
  ! 
  if (mydown == mpi_proc_null) then
     b = 5
     k = k_mysta

     zet1: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k+1) + &
                          sb(1:4,b) * q(1:4,i,j,k+2)
        end do
        end do
     case(4) ! inflow
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,j,b)) * q(1,i,j,k+1) - &
                               rat(i,j,b)  * q(1,i,j,k+2)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j,k+1) + &
           !             sb(1,b) * q(1,i,j,k+2)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j,k_myend)
        end do
        end do
     end select zet1
  end if

  ! boundary in zet-direction (k = kmax)
  ! 
  if (myup == mpi_proc_null) then
     b = 6
     k = k_myend

     zet2: select case(btype(b,myzone))
     case(0) ! -> interface
     case(1:3) ! wall, symmetric plane & freestream
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = sa(1:4,b) * q(1:4,i,j,k-1) + &
                          sb(1:4,b) * q(1:4,i,j,k-2)
        end do
        end do
     case(4) ! inflow
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1,i,j,k) = (one + rat(i,j,b)) * q(1,i,j,k-1) - &
                               rat(i,j,b)  * q(1,i,j,k-2)
           !q(1,i,j,k) = sa(1,b) * q(1,i,j,k-1) + &
           !             sb(1,b) * q(1,i,j,k-2)
        end do
        end do
     case(5) ! exit (MOC)
     case(6) ! periodic condition
        do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           q(1:4,i,j,k) = q(1:4,i,j,k_mysta)
        end do
        end do
     end select zet2
  end if

  ! probe default pressure value
  !
  if (myid == pfix_proc) pfix = q(1,local_ifix,local_jfix,local_kfix)

  ! distribute fixed pressure to all processes
  ! 
  call mpi_bcast(pfix, 1, mpi_real, pfix_proc, mpi_comm_world, ierr)

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     q(1,i,j,k) = q(1,i,j,k) - pfix
  end do
  end do
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Zero velocity at hub
  ! WARNING: every x-plane must be contained in the same proc.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (turbines) then
    do n_t=1,n_tur
      number_disk_nodes = disk_node_coordinates(1,1,n_t)

      do counter = 1, number_disk_nodes
        i        = disk_node_coordinates(counter+1,1,n_t)
        j        = disk_node_coordinates(counter+1,2,n_t)
        k        = disk_node_coordinates(counter+1,3,n_t)
        radius   = disk_dat(j,k,1,n_t)*length_scale

        if ( ((i-gi_ia+1).ge.i_mysta.and.(i-gi_ia+1).le.i_myend).and.((j-gi_ja+1).ge.j_mysta.and.(j-gi_ja+1).le.j_myend) &
          .and.((k-gi_ka+1).ge.k_mysta.and.(k-gi_ka+1).le.k_myend).and.(radius.le.hub_rad(n_t)) ) then
          !print *, myid, 'proc'
          ! if (myid==1) print *, 'i, j, k originales: ',i,j,k
          j = j - gi_ja + 1
          k = k - gi_ka + 1
          
          ! if (myid==1) print *, 'gi_ia, gi_ja, gi_ka', gi_ia, gi_ja, gi_ka 
          ! if (myid==1) print *, 'i, j, k cambiados: ',i,j,k


          !call mpi_barrier(mpi_comm_world,ierr)
          do i_t = itur1(n_t),itur2(n_t)
            i = i_t - gi_ia + 1
            q(2,i,j,k) = 0.0  
          end do
!
        end if
      end do
      !print *, 'saliendo de todos los nodos en la turbina'
    end do
  end if


end subroutine bcond_fm

