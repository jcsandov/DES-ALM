
subroutine brhs_convec_quick (decide_upwind)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Orthogonal, curvlinear, cartesian coordinates

  ! Calculate the convective terms on the exit plane with the
  ! quick scheme in a non-conservative form.

  ! input
  !     decide_upwind
  !     csi(3,ijk)
  !     ucn_j(3,ijk)
  !     aj(ijk)
  !     q(2-4,ijk)
  !    rh(2-4,ijk)

  ! output
  !     rh(2-4,ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! global choice
  integer :: decide_upwind

  ! for 2nd order upwind coef = 1/2
  ! for 3rd order upwind coef = 1/6
  ! for QUICK scheme coef = 1/8 
  real (kind = rdf), parameter :: coef = 0.125_rdf ! 1/8

  ! local arrays
  real (kind = rdf), dimension(:,:), allocatable   :: up
  real (kind = rdf), dimension(:,:), allocatable   :: um
  real (kind = rdf), dimension(:,:,:), allocatable :: rh_dum

  ! local dummy variables
  real (kind = rdf) :: ucon

  integer :: &
       jsta, jend, &
       ksta, kend

  jsta = j_mysta
  ksta = k_mysta

  jend = j_myend
  kend = k_myend

  if (myleft == mpi_proc_null)  jsta = j_mysta + 1
  if (mydown == mpi_proc_null)  ksta = k_mysta + 1

  if (myright == mpi_proc_null) jend = j_myend - 1
  if (myup    == mpi_proc_null) kend = k_myend - 1


  allocate (up(jl:ju,kl:ku), &
            um(jl:ju,kl:ku), &
    rh_dum(1:3,jl:ju,kl:ku))

  ! eta direction

  do k = kl, ku
  do j = jl, ju
     ucon = one_half * ucn_j(2,j,k)
     up(j,k) = ucon + abs(ucon)
     um(j,k) = ucon - abs(ucon)
  end do
  end do

  ! calculate convective term using either a 
  ! a first-order upwind scheme
  ! or a second-order quick scheme

  if (decide_upwind == 1) then
     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
        rh(2:4,j,k) = rh(2:4,j,k) + de &
                    * (up(j,k) * (q(2:4,j  ,k) - q(2:4,j-1,k)) + &
                       um(j,k) * (q(2:4,j+1,k) - q(2:4,j  ,k))  )
     end do
     end do

  else ! if (decide_upwind == 2) then ! leave out if there is only one choice
     do k = k_mysta, k_myend
     do j = jsta, jend
        rh_dum(1:3,j,k) = coef * de &
                           * (up(j,k) * (three * q(2:4,j+1,k)  &
                                      +  three * q(2:4,j,k)    &
                                      -  seven * q(2:4,j-1,k)  &
                                      +          q(2:4,j-2,k)) &

                           +  um(j,k) * (-three * q(2:4,j-1,k) &
                                      -   three * q(2:4,j,k)   &
                                      +   seven * q(2:4,j+1,k) &
                                      -           q(2:4,j+2,k)))
     end do
     end do

     ! Nodes near boundary
     if ( myleft == mpi_proc_null ) then
        do k = k_mysta, k_myend
           j = j_mysta
           rh_dum(1:3,j,k) = de &
                    * (up(j,k) * (q(2:4,j  ,k) - q(2:4,j-1,k)) + &
                       um(j,k) * (q(2:4,j+1,k) - q(2:4,j  ,k))  )
        end do
     end if

     if ( myright == mpi_proc_null ) then
        do k = k_mysta, k_myend
           j = j_myend
           rh_dum(1:3,j,k) = de &
                    * (up(j,k) * (q(2:4,j  ,k) - q(2:4,j-1,k)) + &
                       um(j,k) * (q(2:4,j+1,k) - q(2:4,j  ,k))  )
        end do
     end if

     ! update
     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
        rh(2,j,k) = rh(2,j,k) + rh_dum(1,j,k)
        rh(3,j,k) = rh(3,j,k) + rh_dum(2,j,k)
        rh(4,j,k) = rh(4,j,k) + rh_dum(3,j,k)
     end do
     end do

  end if
  
  ! zeta direction
  !
  do k = kl, ku
  do j = jl, ju
     ucon = one_half * ucn_j(3,j,k)
     up(j,k) = ucon + abs(ucon)
     um(j,k) = ucon - abs(ucon)
  end do
  end do

  ! calculate convective term using either a 
  ! a first-order upwind scheme or a second-order quick scheme

  if (decide_upwind == 1) then
     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
        rh(2:4,j,k) = rh(2:4,j,k) + dz &
                        * (up(j,k) * (q(2:4,j,k  ) - q(2:4,j,k-1)) + &
                           um(j,k) * (q(2:4,j,k+1) - q(2:4,j,k  ))  )
     end do
     end do

  else ! if (decide_upwind == 2) then ! leave out if there is only one choice

     do k = ksta, kend
     do j = j_mysta, j_myend
        rh_dum(1:3,j,k) = coef * dz &
                           * (up(j,k) * (three * q(2:4,j,k+1)  &
                                      +  three * q(2:4,j,k)    &
                                      -  seven * q(2:4,j,k-1)  &
                                      +          q(2:4,j,k-2)) &

                           +  um(j,k) * (-three * q(2:4,j,k-1) &
                                      -   three * q(2:4,j,k)   &
                                      +   seven * q(2:4,j,k+1) &
                                      -           q(2:4,j,k+2)) )
     end do
     end do

     ! Nodes near boundary
     !
     if ( mydown == mpi_proc_null ) then
        do j = j_mysta, j_myend
           k = k_mysta
           rh_dum(1:3,j,k) = dz &
                        * (up(j,k) * (q(2:4,j,k  ) - q(2:4,j,k-1)) + &
                           um(j,k) * (q(2:4,j,k+1) - q(2:4,j,k  ))  )
        end do
     end if

     if ( myup == mpi_proc_null ) then
        do j = j_mysta, j_myend
           k = k_myend
           rh_dum(1:3,j,k) = dz &
                        * (up(j,k) * (q(2:4,j,k  ) - q(2:4,j,k-1)) + &
                           um(j,k) * (q(2:4,j,k+1) - q(2:4,j,k  ))  )
        end do
     end if

     ! update
     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
        rh(2,j,k) = rh(2,j,k) + rh_dum(1,j,k)
        rh(3,j,k) = rh(3,j,k) + rh_dum(2,j,k)
        rh(4,j,k) = rh(4,j,k) + rh_dum(3,j,k)
     end do
     end do

  end if
  
  ! zero dissipation for the momentum equations,
  ! a kludge but keep it for now---could eliminate this
  ! variable
  diss(2:4,:,:) = zero

  !
  deallocate (up, um, rh_dum)


end subroutine brhs_convec_quick

