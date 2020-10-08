
subroutine rhs_convec_quick (decide_upwind)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Orthogonal, curvlinear, cartesian coordinates

  ! Calculate the dissipation part of the quick scheme
  ! or other first- or second-order upwind for the momentum
  ! equations

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
  !
  integer :: decide_upwind

  ! for 2nd order upwind coef = 1/2
  ! for 3rd order upwind coef = 1/6
  ! for QUICK scheme coef = 1/8 
  !
  real (kind = rdf), parameter :: coef = 0.125_rdf ! 1/8

  ! local arrays
  !
  real (kind = rdf), dimension(:,:,:), allocatable :: up
  real (kind = rdf), dimension(:,:,:), allocatable :: um

  ! local dummy variables
  !
  real (kind = rdf) :: ucon
  real (kind = rdf) :: dc2, de2, dz2

  integer :: &
       ista, iend, &
       jsta, jend, &
       ksta, kend

  ! local integers
  !
  integer :: i, j, k

  allocate (up(il:iu,jl:ju,kl:ku), &
            um(il:iu,jl:ju,kl:ku))

  ! csi direction

  ! calculate the positive and negative contravariant velocities
  ! 
  ! Note -----
  !
  ! (U/J)+ & (U/J)- not just U+ or U-

  !============================================================
  !
  ! MPI issues
  !
  ! out-of-array bounds access ---> need
  ! switch for process who own the boundary
  !
  ! 
  !============================================================

  ! use ista, iend to keep from having
  ! a array-bound access error when process
  ! owns the boundary

  ista=i_mysta; jsta=j_mysta; ksta=k_mysta
  iend=i_myend; jend=j_myend; kend=k_myend
  
  if (myback == mpi_proc_null)  ista=i_mysta+1
  if (myleft == mpi_proc_null)  jsta=j_mysta+1
  if (mydown == mpi_proc_null)  ksta=k_mysta+1

  if (myfront == mpi_proc_null) iend=i_myend-1
  if (myright == mpi_proc_null) jend=j_myend-1
  if (myup    == mpi_proc_null) kend=k_myend-1

  do k=kl,ku
  do j=jl,ju
  do i=il,iu
     ucon=one_half*ucn_j(1,i,j,k)
     up(i,j,k)=ucon+abs(ucon)
     um(i,j,k)=ucon-abs(ucon)
  end do
  end do
  end do

  ! calculate convective term using either a 
  ! a first-order upwind scheme
  ! or a second-order quick scheme

  if (decide_upwind == 1) then

     do k=k_mysta,k_myend
     do j=j_mysta,j_myend
     do i=i_mysta,i_myend

        fv(1:3,i,j,k) = dc &
                        * (up(i,j,k)*(q(2:4,i  ,j,k)-q(2:4,i-1,j,k)) + &
                           um(i,j,k)*(q(2:4,i+1,j,k)-q(2:4,i  ,j,k))  )
     end do
     end do
     end do

  else ! if (decide_upwind == 2) then ! leave out if there is only one choice
 
     do k=k_mysta,k_myend
     do j=j_mysta,j_myend
     do i=ista,iend

        fv(1:3,i,j,k) = coef * dc &
                         * (up(i,j,k) * (three * q(2:4,i+1,j,k)  &
                                      +  three * q(2:4,i,j,k)    &
                                      -  seven * q(2:4,i-1,j,k)  &
                                      +          q(2:4,i-2,j,k)) &

                         +  um(i,j,k) * (-three * q(2:4,i-1,j,k) &
                                      -   three * q(2:4,i,j,k)   &
                                      +   seven * q(2:4,i+1,j,k) &
                                      -           q(2:4,i+2,j,k)))
     end do
     end do
     end do


     ! Near-boundary nodes
     ! 

     if (myback == mpi_proc_null) then
        do k=k_mysta,k_myend
        do j=j_mysta,j_myend
           i=i_mysta
           fv(1:3,i,j,k) = dc &
                        * (up(i,j,k)*(q(2:4,i  ,j,k)-q(2:4,i-1,j,k)) + &
                           um(i,j,k)*(q(2:4,i+1,j,k)-q(2:4,i  ,j,k))  )

        end do
        end do
     end if
     
     if (myfront == mpi_proc_null) then
        do k=k_mysta,k_myend
        do j=j_mysta,j_myend
           i=i_myend
           fv(1:3,i,j,k) = dc &
                        * (up(i,j,k)*(q(2:4,i  ,j,k)-q(2:4,i-1,j,k)) + &
                           um(i,j,k)*(q(2:4,i+1,j,k)-q(2:4,i  ,j,k))  )
        end do
        end do
     end if

     ! blanking interface
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        i = li_blk_ia(n,nb) - 1
        if ( blktype(1,nb,myzone) == 0 .and. i > i_mysta ) then
           do k = li_blk_ka(n,nb), li_blk_kb(n,nb) 
           do j = li_blk_ja(n,nb), li_blk_jb(n,nb) 
              fv(1:3,i,j,k) = dc &
                           * (up(i,j,k) * (q(2:4,i  ,j,k) - q(2:4,i-1,j,k)) + &
                              um(i,j,k) * (q(2:4,i+1,j,k) - q(2:4,i  ,j,k)))
           end do
           end do
        end if
        i = li_blk_ib(n,nb) + 1
        if ( blktype(2,nb,myzone) == 0 .and. i < i_myend ) then
           do k = li_blk_ka(n,nb), li_blk_kb(n,nb) 
           do j = li_blk_ja(n,nb), li_blk_jb(n,nb) 
              fv(1:3,i,j,k) = dc &
                           * (up(i,j,k) * (q(2:4,i  ,j,k) - q(2:4,i-1,j,k)) + &
                              um(i,j,k) * (q(2:4,i+1,j,k) - q(2:4,i  ,j,k)))
           end do
           end do
        end if
     end do
     end if

  end if


  do k=k_mysta,k_myend
  do j=j_mysta,j_myend
  do i=i_mysta,i_myend
     rh(2,i,j,k)=rh(2,i,j,k)+fv(1,i,j,k)
     rh(3,i,j,k)=rh(3,i,j,k)+fv(2,i,j,k)
     rh(4,i,j,k)=rh(4,i,j,k)+fv(3,i,j,k)
  end do
  end do
  end do


  ! eta direction

  do k=kl,ku
  do j=jl,ju
  do i=il,iu
     ucon=one_half*ucn_j(2,i,j,k)
     up(i,j,k)=ucon+abs(ucon)
     um(i,j,k)=ucon-abs(ucon)
  end do
  end do
  end do

  ! calculate convective term using either a first-order upwind scheme
  ! or a second-order quick scheme

  if (decide_upwind == 1) then

      do k=k_mysta,k_myend
      do j=j_mysta,j_myend
      do i=i_mysta,i_myend

         fv(1:3,i,j,k) = de &
                      * (up(i,j,k)*(q(2:4,i,j  ,k)-q(2:4,i,j-1,k)) + &
                         um(i,j,k)*(q(2:4,i,j+1,k)-q(2:4,i,j  ,k))  )
     end do
     end do
     end do

  else ! if (decide_upwind == 2) then ! leave out if there is only one choice

     do k=k_mysta, k_myend
     do j=jsta, jend
     do i=i_mysta, i_myend
        fv(1:3,i,j,k) = coef * de &
                         * (up(i,j,k) * ( three*q(2:4,i,j+1,k)  &
                                      +   three*q(2:4,i,j,k)    &
                                      -   seven*q(2:4,i,j-1,k)  &
                                      +         q(2:4,i,j-2,k)) &

                         +  um(i,j,k) * (-three*q(2:4,i,j-1,k) &
                                      -   three*q(2:4,i,j,k)   &
                                      +   seven*q(2:4,i,j+1,k) &
                                      -         q(2:4,i,j+2,k)))
     end do
     end do
     end do

     ! Nodes near boundary

     if (myleft == mpi_proc_null) then
        do k=k_mysta,k_myend
        do i=i_mysta,i_myend
           j=j_mysta
           fv(1:3,i,j,k) = de &
                        * (up(i,j,k)*(q(2:4,i,j  ,k)-q(2:4,i,j-1,k)) + &
                           um(i,j,k)*(q(2:4,i,j+1,k)-q(2:4,i,j  ,k)))

        end do
        end do
     end if

     if (myright == mpi_proc_null) then
        do k=k_mysta,k_myend
        do i=i_mysta,i_myend
           j=j_myend
           fv(1:3,i,j,k) = de &
                        * (up(i,j,k)*(q(2:4,i,j  ,k)-q(2:4,i,j-1,k)) + &
                           um(i,j,k)*(q(2:4,i,j+1,k)-q(2:4,i,j  ,k)))
        end do
        end do
     end if

     ! blanking interface
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        j = li_blk_ja(n,nb) - 1
        if ( blktype(3,nb,myzone) == 0 .and. j > j_mysta ) then
           do k=li_blk_ka(n,nb),li_blk_kb(n,nb) 
           do i=li_blk_ia(n,nb),li_blk_ib(n,nb) 
              fv(1:3,i,j,k) = de &
                            * (up(i,j,k)*(q(2:4,i,j  ,k)-q(2:4,i,j-1,k)) + &
                               um(i,j,k)*(q(2:4,i,j+1,k)-q(2:4,i,j  ,k)))
           end do
           end do
        end if
        j = li_blk_jb(n,nb) + 1
        if ( blktype(4,nb,myzone) == 0 .and. j < j_myend ) then
           do k = li_blk_ka(n,nb), li_blk_kb(n,nb) 
           do i = li_blk_ia(n,nb), li_blk_ib(n,nb) 
              fv(1:3,i,j,k) = de &
                            * (up(i,j,k) * (q(2:4,i,j  ,k) - q(2:4,i,j-1,k)) + &
                               um(i,j,k) * (q(2:4,i,j+1,k) - q(2:4,i,j  ,k)))
           end do
           end do
        end if
     end do
     end if

  end if

     do k=k_mysta,k_myend
     do j=j_mysta,j_myend
     do i=i_mysta,i_myend
        rh(2,i,j,k)=rh(2,i,j,k)+fv(1,i,j,k)
        rh(3,i,j,k)=rh(3,i,j,k)+fv(2,i,j,k)
        rh(4,i,j,k)=rh(4,i,j,k)+fv(3,i,j,k)
     end do
     end do
     end do
  
  ! zeta direction

  do k=kl,ku
  do j=jl,ju
  do i=il,iu
     ucon=one_half*ucn_j(3,i,j,k)
     up(i,j,k)=ucon+abs(ucon)
     um(i,j,k)=ucon-abs(ucon)
  end do
  end do
  end do

  ! calculate convective term using either a 
  ! a first-order upwind scheme
  ! or a second-order quick scheme

  if (decide_upwind == 1) then

     do k=k_mysta,k_myend
     do j=j_mysta,j_myend
     do i=i_mysta,i_myend
        fv(1:3,i,j,k) = dz &
                     * (up(i,j,k)*(q(2:4,i,j,k  )-q(2:4,i,j,k-1)) + &
                        um(i,j,k)*(q(2:4,i,j,k+1)-q(2:4,i,j,k  ))  )
     end do
     end do
     end do

  else ! if (decide_upwind == 2) then ! leave out if there is only one choice

     do k=ksta, kend
     do j=j_mysta,j_myend
     do i=i_mysta,i_myend
        fv(1:3,i,j,k) = coef * dz &
                         * (up(i,j,k) * ( three*q(2:4,i,j,k+1)  &
                                      +   three*q(2:4,i,j,k)    &
                                      -   seven*q(2:4,i,j,k-1)  &
                                      +         q(2:4,i,j,k-2)) &

                         +  um(i,j,k) * (-three*q(2:4,i,j,k-1) &
                                      -   three*q(2:4,i,j,k)   &
                                      +   seven*q(2:4,i,j,k+1) &
                                      -         q(2:4,i,j,k+2)) )
     end do
     end do
     end do

     ! Nodes near boundary

     if (mydown == mpi_proc_null) then

        do j=j_mysta,j_myend
        do i=i_mysta,i_myend
           k=k_mysta
           fv(1:3,i,j,k) = dz &
                        * (up(i,j,k)*(q(2:4,i,j,k  )-q(2:4,i,j,k-1)) + &
                           um(i,j,k)*(q(2:4,i,j,k+1)-q(2:4,i,j,k  ))  )

        end do
        end do
     end if

     if (myup == mpi_proc_null) then
        do j=j_mysta,j_myend
        do i=i_mysta,i_myend
           k=k_myend
           fv(1:3,i,j,k) = dz &
                        * (up(i,j,k)*(q(2:4,i,j,k  )-q(2:4,i,j,k-1)) + &
                           um(i,j,k)*(q(2:4,i,j,k+1)-q(2:4,i,j,k  ))  )
        end do
        end do
     end if
  
     ! blanking interface
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        k = li_blk_ka(n,nb) - 1
        if ( blktype(5,nb,myzone) == 0 .and. k > k_mysta ) then
           do j = li_blk_ja(n,nb), li_blk_jb(n,nb) 
           do i = li_blk_ia(n,nb), li_blk_ib(n,nb) 
              fv(1:3,i,j,k) = dz &
                            * (up(i,j,k)*(q(2:4,i,j,k  )-q(2:4,i,j,k-1)) + &
                               um(i,j,k)*(q(2:4,i,j,k+1)-q(2:4,i,j,k  )))
           end do
           end do
        end if
        k = li_blk_kb(n,nb) + 1
        if ( blktype(6,nb,myzone) == 0 .and. k < k_myend ) then
           do j = li_blk_ja(n,nb), li_blk_jb(n,nb) 
           do i = li_blk_ia(n,nb), li_blk_ib(n,nb) 
              fv(1:3,i,j,k) = dz &
                            * (up(i,j,k) * (q(2:4,i,j,k  ) - q(2:4,i,j,k-1)) + &
                               um(i,j,k) * (q(2:4,i,j,k+1) - q(2:4,i,j,k  )))
           end do
           end do
        end if
     end do
     end if

  end if

     do k=k_mysta,k_myend
     do j=j_mysta,j_myend
     do i=i_mysta,i_myend
        rh(2,i,j,k)=rh(2,i,j,k)+fv(1,i,j,k)
        rh(3,i,j,k)=rh(3,i,j,k)+fv(2,i,j,k)
        rh(4,i,j,k)=rh(4,i,j,k)+fv(3,i,j,k)
     end do
     end do
     end do


  ! zero dissipation for the momentum equations, 
  ! (implicit in non-conservative quick scheme)
  ! a kludge but keep it for now---could eliminate this
  ! variable

  diss(2:4,:,:,:) = zero

  deallocate (up, um)


end subroutine rhs_convec_quick



