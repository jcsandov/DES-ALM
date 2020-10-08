
subroutine rhs_sa_bcond ()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! boundary conditions for the modified eddy viscosity in des
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! local parameter for S-A
  real (kind = rdf) :: chi, chic, fv1
  real (kind = rdf) :: xnut_org

  ! local dummy variables
  real (kind = rdf) :: rei

  integer :: b

  integer :: ista, iend, &
             jsta, jend, &
             ksta, kend

  rei  = one / ren

  ista = i_mysta
  iend = i_myend
  jsta = j_mysta
  jend = j_myend
  ksta = k_mysta
  kend = k_myend

  ! include points on solid wall
  if (myback  == mpi_proc_null) ista = i_mysta - 1
  if (myleft  == mpi_proc_null) jsta = j_mysta - 1
  if (mydown  == mpi_proc_null) ksta = k_mysta - 1

  if (myfront == mpi_proc_null) iend = i_myend + 1
  if (myright == mpi_proc_null) jend = j_myend + 1
  if (myup    == mpi_proc_null) kend = k_myend + 1

  ! boundary in csi-direction (i = 1)
  !
  if ( myback == mpi_proc_null ) then
      b = 1
      i = i_mysta - 1
     csi1: select case(btype(b,myzone))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = sa(1,b) * q(5,i+1,j,k) + &
                        sb(1,b) * q(5,i+2,j,k)
        end do
        end do
     case(4) ! inlet
     case(5)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = (one + rat(j,k,b)) * q(5,i+1,j,k) - &
                               rat(j,k,b)  * q(5,i+2,j,k)
        end do
        end do
     case(6)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = q(5,i_myend,j,k)
        end do
        end do
     end select csi1
  end if

  ! boundary in csi-direction (i = imax)
  ! 
  if ( myfront == mpi_proc_null ) then
     b = 2
     i = i_myend + 1
     csi2: select case(btype(b,myzone))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = sa(1,b) * q(5,i-1,j,k) + &
                        sb(1,b) * q(5,i-2,j,k)
        end do
        end do
     case(4)
     case(5)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = (one + rat(j,k,b))* q(5,i-1,j,k) - &
                               rat(j,k,b) * q(5,i-2,j,k)
        end do
        end do
     case(6)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = q(5,i_mysta,j,k)
        end do
        end do
case(8) ! free surface
        do k = ksta, kend
        do j = jsta, jend
         q(5,i,j,k) = q(5,i-1,j,k) !sa(1:4,b)*q(1:4,i-1,j,k) + sb(1:4,b)*q(1:4,i-2,j,k)

        end do
        end do
 

    end select csi2
  end if

  ! boundary in eta-direction (j = 1)
  !
  if ( myleft == mpi_proc_null ) then
      b = 3
      j = j_mysta - 1
     eta1: select case(btype(b,myzone))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j+1,k) + &
                        sb(1,b) * q(5,i,j+2,k)
        end do
        end do
     case(4)
     case(5)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,k,b))* q(5,i,j+1,k) - &
                               rat(i,k,b) * q(5,i,j+2,k)
        end do
        end do
     case(6)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j_myend,k)
        end do
        end do
     end select eta1
  end if

  ! boundary in eta-direction (j = jmax)
  !
  if ( myright == mpi_proc_null ) then
      b = 4
      j = j_myend + 1
     eta2: select case(btype(b,myzone))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j-1,k) + &
                        sb(1,b) * q(5,i,j-2,k)
        end do
        end do
     case(4)
     case(5)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,k,b))* q(5,i,j-1,k) - &
                               rat(i,k,b) * q(5,i,j-2,k)
        end do
        end do
     case(6)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j_mysta,k)
        end do
        end do
     end select eta2
  end if

  ! boundary in zet-direction (k = 1)
  !
  if ( mydown == mpi_proc_null ) then
      b = 5
      k = k_mysta - 1
     zet1: select case(btype(b,myzone))
     case(0)
     case(1) ! solid
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j,k+1) + &
                        sb(1,b) * q(5,i,j,k+2)
        end do
        end do
     case(4)
     case(5)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,j,b))* q(5,i,j,k+1) - &
                               rat(i,j,b) * q(5,i,j,k+2)
        end do
        end do
     case(6)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j,k_myend)
        end do
        end do
     end select zet1
  end if

  ! boundary in zet-direction (k = kmax)
  !
  if ( myup == mpi_proc_null ) then
      b = 6
      k = k_myend + 1
     zet2: select case(btype(b,myzone))
     case(0)
     case(1) ! solid
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j,k-1) + &
                        sb(1,b) * q(5,i,j,k-2)
        end do
        end do
     case(4)
     case(5)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,j,b))* q(5,i,j,k-1) - &
                               rat(i,j,b) * q(5,i,j,k-2)
        end do
        end do
     case(6)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j,k_mysta)
        end do
        end do
     end select zet2
  end if

  !==================================================
  !
  ! calculate eddy viscosity : xnut = xnut_bar * f_v1
  ! exclude ghost points & update gp after returing
  ! to des eddy
  ! 
  !==================================================
  !
  do k = ksta, kend
  do j = jsta, jend
  do i = ista, iend

     q(5,i,j,k) = max(1.0E-07, q(5,i,j,k))
     chi = q(5,i,j,k) / rei
     chic = chi*chi*chi
     fv1  = chic / (chic + cv1c)
     xnut(i,j,k) = q(5,i,j,k) * fv1

  end do
  end do
  end do


end subroutine rhs_sa_bcond


