
subroutine rhs_kw_bcond ()

  use global_app

  implicit none

  
  ! The definitions of i_mysta etc., match solver_daf & kw_eddy and
  ! not the definitions in bcond_fm 
  !
  ! *** This is a terrible practice ---> inconsistency ****
  !
  !  The two routines should be consistent ---> change bcond_fm
  
  ! kw parameters
  !
  real (kind = rdf), parameter :: sixty = 60.00_rdf
  real (kind = rdf), parameter :: pt075 =  0.075_rdf

  ! local dummy variables
  !
  real (kind = rdf) :: yd
  real (kind = rdf) :: srenp

  integer :: ista, iend, &
             jsta, jend, &
             ksta, kend

  ! dummy
  srenp = sixty / ren / pt075

  ista = i_mysta
  iend = i_myend
  jsta = j_mysta
  jend = j_myend
  ksta = k_mysta
  kend = k_myend

  ! exclude points on solid wall
  if (myback  == mpi_proc_null .and. btype(1) /= 1) ista = i_mysta - 1
  if (myleft  == mpi_proc_null .and. btype(3) /= 1) jsta = j_mysta - 1
  if (mydown  == mpi_proc_null .and. btype(5) /= 1) ksta = k_mysta - 1

  if (myfront == mpi_proc_null .and. btype(2) /= 1) iend = i_myend + 1
  if (myright == mpi_proc_null .and. btype(4) /= 1) jend = j_myend + 1
  if (myup    == mpi_proc_null .and. btype(6) /= 1) kend = k_myend + 1


  !==================================================
  !
  ! k & omega
  !
  !==================================================

  ! boundary in csi-direction (i = 1)
  !
  if ( myback == mpi_proc_null ) then
      b = 1
      i = i_mysta - 1
     csi1: select case(btype(b))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do j = jsta, jend
           yd = wd(i+1,j,k)
           q(6,i,j,k) = srenp/yd/yd
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = sa(1,b) * q(5,i+1,j,k) + &
                        sb(1,b) * q(5,i+2,j,k)
           q(6,i,j,k) = sa(1,b) * q(6,i+1,j,k) + &
                        sb(1,b) * q(6,i+2,j,k)
           ! bound k and omega
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(4) ! inlet
     case(5) ! outlet
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = (one + rat(j,k,b)) * q(5,i+1,j,k) - &
                               rat(j,k,b)  * q(5,i+2,j,k)
           q(6,i,j,k) = (one + rat(j,k,b)) * q(6,i+1,j,k) - &
                               rat(j,k,b)  * q(6,i+2,j,k)
           ! bound k and omega
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(6)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = q(5,i_myend,j,k)
           q(6,i,j,k) = q(6,i_myend,j,k)
           ! bound k and omega
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     end select csi1
  end if

  ! boundary in csi-direction (i = imax)
  !
  if ( myfront == mpi_proc_null ) then
     b = 2
     i = i_myend + 1
     csi2: select case(btype(b))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do j = jsta, jend
           yd = wd(i-1,j,k)
           q(6,i,j,k) = srenp/yd/yd
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = sa(1,b) * q(5,i-1,j,k) + &
                        sb(1,b) * q(5,i-2,j,k)
           q(6,i,j,k) = sa(1,b) * q(6,i-1,j,k) + &
                        sb(1,b) * q(6,i-2,j,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(4)
     case(5)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = (one + rat(j,k,b))* q(5,i-1,j,k) - &
                               rat(j,k,b) * q(5,i-2,j,k)
           q(6,i,j,k) = (one + rat(j,k,b))* q(6,i-1,j,k) - &
                               rat(j,k,b) * q(6,i-2,j,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(6)
        do k = ksta, kend
        do j = jsta, jend
           q(5,i,j,k) = q(5,i_mysta,j,k)
           q(6,i,j,k) = q(6,i_mysta,j,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     end select csi2
  end if

  ! boundary in eta-direction (j = 1)
  !
  if ( myleft == mpi_proc_null ) then
      b = 3
      j = j_mysta - 1
     eta1: select case(btype(b))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do i = ista, iend
           yd = wd(i,j+1,k)
           q(6,i,j,k) = srenp/yd/yd
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j+1,k) + &
                        sb(1,b) * q(5,i,j+2,k)
           q(6,i,j,k) = sa(1,b) * q(6,i,j+1,k) + &
                        sb(1,b) * q(6,i,j+2,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(4)
     case(5)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,k,b))* q(5,i,j+1,k) - &
                               rat(i,k,b) * q(5,i,j+2,k)
           q(6,i,j,k) = (one + rat(i,k,b))* q(6,i,j+1,k) - &
                               rat(i,k,b) * q(6,i,j+2,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(6)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j_myend,k)
           q(6,i,j,k) = q(6,i,j_myend,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     end select eta1
  end if


  ! boundary in eta-direction (j = jmax)
  !
  if ( myright == mpi_proc_null ) then
      b = 4
      j = j_myend + 1
     eta2: select case(btype(b))
     case(0)
     case(1) ! solid
        do k = ksta, kend
        do i = ista, iend
           yd = wd(i,j-1,k)
           q(6,i,j,k) = srenp/yd/yd
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j-1,k) + &
                        sb(1,b) * q(5,i,j-2,k)
           q(6,i,j,k) = sa(1,b) * q(6,i,j-1,k) + &
                        sb(1,b) * q(6,i,j-2,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(4)
     case(5)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,k,b))* q(5,i,j-1,k) - &
                               rat(i,k,b) * q(5,i,j-2,k)
           q(6,i,j,k) = (one + rat(i,k,b))* q(6,i,j-1,k) - &
                               rat(i,k,b) * q(6,i,j-2,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(6)
        do k = ksta, kend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j_mysta,k)
           q(6,i,j,k) = q(6,i,j_mysta,k)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     end select eta2
  end if


  ! boundary in zet-direction (k = 1)
  !
  if ( mydown == mpi_proc_null ) then
      b = 5
      k = k_mysta - 1
     zet1: select case(btype(b))
     case(0)
     case(1) ! solid
        do j = jsta, jend
        do i = ista, iend
           yd = wd(i,j,k+1)
           q(6,i,j,k) = srenp/yd/yd
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j,k+1) + &
                        sb(1,b) * q(5,i,j,k+2)
           q(6,i,j,k) = sa(1,b) * q(6,i,j,k+1) + &
                        sb(1,b) * q(6,i,j,k+2)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(4)
     case(5)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,j,b))* q(5,i,j,k+1) - &
                               rat(i,j,b) * q(5,i,j,k+2)
           q(6,i,j,k) = (one + rat(i,j,b))* q(6,i,j,k+1) - &
                               rat(i,j,b) * q(6,i,j,k+2)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(6)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j,k_myend)
           q(6,i,j,k) = q(6,i,j,k_myend)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     end select zet1
  end if


  ! boundary in zet-direction (k = kmax)
  !
  if ( myup == mpi_proc_null ) then
      b = 6
      k = k_myend + 1
     zet2: select case(btype(b))
     case(0)
     case(1) ! solid
        do j = jsta, jend
        do i = ista, iend
           yd = wd(i,j,k-1)
           q(6,i,j,k) = srenp/yd/yd
           q(5,i,j,k) = zero
        end do
        end do
     case(2:3)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = sa(1,b) * q(5,i,j,k-1) + &
                        sb(1,b) * q(5,i,j,k-2)
           q(6,i,j,k) = sa(1,b) * q(6,i,j,k-1) + &
                        sb(1,b) * q(6,i,j,k-2)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(4)
     case(5)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = (one + rat(i,j,b))* q(5,i,j,k-1) - &
                               rat(i,j,b) * q(5,i,j,k-2)
           q(6,i,j,k) = (one + rat(i,j,b))* q(6,i,j,k-1) - &
                               rat(i,j,b) * q(6,i,j,k-2)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     case(6)
        do j = jsta, jend
        do i = ista, iend
           q(5,i,j,k) = q(5,i,j,k_mysta)
           q(6,i,j,k) = q(6,i,j,k_mysta)
           q(5,i,j,k) = max(ckin, q(5,i,j,k))
           q(6,i,j,k) = max(cein, q(6,i,j,k))
        end do
        end do
     end select zet2
  end if

  !==================================================
  !
  ! calculate eddy viscosity (include ghost points)
  !
  !==================================================
  !
  call rhs_kw_coef ()

  if ( nlinc ) then
     ! calcuate xnut using definition by Craft et al. (1995)
     ! in kw_eddy subroutine
  else
     do k = kl, ku
     do j = jl, ju
     do i = il, iu
        xnut(i,j,k) = As(i,j,k) * q(5,i,j,k) / q(6,i,j,k)
     end do
     end do
     end do
  end if


  !==================================================
  !
  ! average omega at the corners (c)
  !
  !==================================================

  ! c = 1 --> i = 1; j = 1
  ! 
  if ( myback == mpi_proc_null .and. myleft == mpi_proc_null .and. &
       btype(1) == 1 .and btype(3) == 1) then
     i = i_mysta - 1
     j = j_mysta - 1
     do k = k_mysta - 1, k_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i+1,j,k) + q(6,i,j+1,k))
     end do
  end if

  ! c = 2 --> i = 1; j = jmax
  ! 
  if ( myback == mpi_proc_null .and. myright == mpi_proc_null .and. &
       btype(1) == 1 .and. btype(4) == 1) then
     i = i_mysta - 1
     j = j_myend + 1
     do k = k_mysta - 1, k_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i+1,j,k) + q(6,i,j-1,k))
     end do
  end if
  
  ! c = 3 --> i = 1; k = 1
  ! 
  if ( myback == mpi_proc_null .and. mydown == mpi_proc_null .and. &
       btype(1) == 1 .and. btype(5) == 1) then
     i = i_mysta - 1
     k = k_mysta - 1
     do j = j_mysta - 1, j_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i+1,j,k) + q(6,i,j,k+1))
     end do
  end if

  ! c = 4 --> i = 1; k = kmax
  ! 
  if ( myback == mpi_proc_null .and. myup == mpi_proc_null .and. &
       btype(1) == 1 .and. btype(6) == 1) then
     i = i_mysta - 1
     k = k_myend + 1
     do j = j_mysta - 1, j_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i+1,j,k) + q(6,i,j,k-1))
     end do
  end if
  
  ! c = 5 --> i = imax; j = 1
  ! 
  if ( myfront == mpi_proc_null .and. myleft == mpi_proc_null .and. &
       btype(2) == 1 .and. btype(3) == 1) then
     i = i_myend + 1
     j = j_mysta - 1
     do k = k_mysta - 1, k_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i-1,j,k) + q(6,i,j+1,k))
     end do
  end if

  ! c = 6 --> i = imax; j = jmax
  ! 
  if ( myfront == mpi_proc_null .and. myright == mpi_proc_null .and. &
       btype(2) == 1 .and. btype(4) == 1) then
     i = i_myend + 1
     j = j_myend + 1
     do k = k_mysta - 1, k_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i-1,j,k) + q(6,i,j-1,k))
     end do
  end if
  
  ! c = 7 --> i = imax; k = 1
  ! 
  if ( myfront == mpi_proc_null .and. mydown == mpi_proc_null .and. &
       btype(2) == 1 .and. btype(5) == 1) then
     i = i_myend + 1
     k = k_mysta - 1
     do j = j_mysta - 1, j_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i-1,j,k) + q(6,i,j,k+1))
     end do
  end if

  ! c = 8 --> i = imax; k = kmax
  ! 
  if ( myfront == mpi_proc_null .and. myup == mpi_proc_null .and. &
       btype(2) == 1 .and. btype(6) == 1) then
     i = i_myend + 1
     k = k_myend + 1
     do j = j_mysta - 1, j_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i-1,j,k) + q(6,i,j,k-1))
     end do
  end if

  ! c = 9 --> j = 1; k = 1
  ! 
  if ( myleft == mpi_proc_null .and. mydown == mpi_proc_null .and. &
       btype(3) == 1 .and. btype(5) == 1) then
     j = j_mysta - 1
     k = k_mysta - 1
     do i = i_mysta - 1, i_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i,j+1,k) + q(6,i,j,k+1))
     end do
  end if

  ! c = 10 --> j = 1; k = kmax
  ! 
  if ( myleft == mpi_proc_null .and. myup == mpi_proc_null .and. &
       btype(3) == 1 .and. btype(6) == 1) then
     j = j_mysta - 1
     k = k_myend + 1
     do i = i_mysta - 1, i_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i,j+1,k) + q(6,i,j,k-1))
     end do
  end if
  
  ! c = 11 --> j = jm; k = 1
  ! 
  if ( myright == mpi_proc_null .and. mydown == mpi_proc_null .and. &
       btype(4) == 1 .and. btype(5) == 1) then
     j = j_myend + 1
     k = k_mysta - 1
     do i = i_mysta - 1, i_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i,j-1,k) + q(6,i,j,k+1))
     end do
  end if
  
  ! c = 12 --> j = jmax; k = kmax
  !
  if ( myright == mpi_proc_null .and. myup == mpi_proc_null .and. &
       btype(4) == 1 .and. btype(6) == 1) then
     j = j_myend + 1
     k = k_myend + 1
     do i = i_mysta - 1, i_myend + 1
        q(6,i,j,k) = pt5 * (q(6,i,j-1,k) + q(6,i,j,k-1))
     end do
  end if

end subroutine rhs_kw_bcond


