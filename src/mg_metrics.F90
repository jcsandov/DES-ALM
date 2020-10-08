
subroutine mg_metrics

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates
  !
  ! for mp program

  ! Calculate metrics of the geometric transformation for
  ! for the fine grid and each of the coarse grids

  ! input
  !     x(nmxg)
  !     y(nmxg)
  !     z(nmxg)

  ! output
  !     csi(1:3,nmxg) csi_x, csi_y and csi_z
  !     eta(1:3,nmxg) eta_x, eta_y and eta_z
  !     zet(1:3,nmxg) zet_x, zet_y and zet_z
  !      aj(nmxg)   Jacobian
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global
  use global_param
  use global_app
  use global_mpi

  implicit none

  ! local temporary arrays
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: xc
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: xe
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: xz
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: yc
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: ye
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: yz
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: zc
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: ze
  !real (kind = rdf), dimension (le_idx_a(1):le_idx_b(ng)) :: zz

  ! local variables
  real (kind = rdf) :: gj, dum, dum1
  real (kind = rdf) :: Sc, Se, Sz, SL

  integer :: n, i, j, k, l, p, m
  integer :: ic, ie, iz

  integer :: ka
  integer :: kb
  integer :: ja
  integer :: jb
  integer :: ia
  integer :: ib


  do n = 1, ng

     ic = 1
     ie = le_ix(n)
     iz = le_ix(n) * le_jx(n)

     ! csi direction
     !
     dum1 = dc(n)
     dum  = one_half * dc(n)

     ka = li_ka(n)
     kb = li_kb(n)
     ja = li_ja(n)
     jb = li_jb(n)
     ia = li_ia(n)
     ib = li_ib(n)

     if (myback  == mpi_proc_null) ia = li_ia(n) + 1
     if (myfront == mpi_proc_null) ib = li_ib(n) - 1

     do k = ka, kb
     do j = ja, jb
     do i = ia, ib
        l = le_idx(i,j,k,n)
        xc(l) = dum * (x(l+ic) - x(l-ic))
        yc(l) = dum * (y(l+ic) - y(l-ic))
        zc(l) = dum * (z(l+ic) - z(l-ic))
     end do
     end do
     end do

     ! i boundaries
     !
     if (myback == mpi_proc_null) then
        do k = ka, kb
        do j = ja, jb
           l = le_idx(li_ia(n),j,k,n)
           xc(l) = dum1 * (x(l+ic) - x(l))
           yc(l) = dum1 * (y(l+ic) - y(l))
           zc(l) = dum1 * (z(l+ic) - z(l))
        end do
        end do
     end if

     if (myfront == mpi_proc_null) then
        do k = ka, kb
        do j = ja, jb
           l = le_idx(li_ib(n),j,k,n)
           xc(l) = dum1 * (x(l) - x(l-ic))
           yc(l) = dum1 * (y(l) - y(l-ic))
           zc(l) = dum1 * (z(l) - z(l-ic))
        end do
        end do
     end if

     ! eta direction (not 2d ready)
     !
     dum1 = de(n)
     dum = one_half * de(n)

     ka = li_ka(n)
     kb = li_kb(n)
     ja = li_ja(n)
     jb = li_jb(n)
     ia = li_ia(n)
     ib = li_ib(n)

     if (myleft  == mpi_proc_null) ja = li_ja(n) + 1
     if (myright == mpi_proc_null) jb = li_jb(n) - 1

     do k = ka, kb
     do j = ja, jb
     do i = ia, ib
        l = le_idx(i,j,k,n)
        xe(l) = dum * (x(l+ie) - x(l-ie))
        ye(l) = dum * (y(l+ie) - y(l-ie))
        ze(l) = dum * (z(l+ie) - z(l-ie))
     end do
     end do
     end do

     ! j boundaries
     !
     if (myleft  == mpi_proc_null) then
        do k = ka, kb
        do i = ia, ib
           l = le_idx(i,li_ja(n),k,n)
           xe(l) = dum1 * (x(l+ie) - x(l))
           ye(l) = dum1 * (y(l+ie) - y(l))
           ze(l) = dum1 * (z(l+ie) - z(l))
        end do
        end do
     end if

     if (myright == mpi_proc_null) then
        do k = ka, kb
        do i = ia, ib
           l = le_idx(i,li_jb(n),k,n)
           xe(l) = dum1 * (x(l) - x(l-ie))
           ye(l) = dum1 * (y(l) - y(l-ie))
           ze(l) = dum1 * (z(l) - z(l-ie))
        end do
        end do
     end if

     ! zet direction
     !
     dum1 = dz(n)
     dum = one_half * dz(n)

     ka = li_ka(n)
     kb = li_kb(n)
     ja = li_ja(n)
     jb = li_jb(n)
     ia = li_ia(n)
     ib = li_ib(n)

     if (mydown == mpi_proc_null) ka = li_ka(n) + 1
     if (myup   == mpi_proc_null) kb = li_kb(n) - 1

     do k = ka, kb
     do j = ja, jb
     do i = ia, ib
        l = le_idx(i,j,k,n)
        xz(l) = dum * (x(l+iz) - x(l-iz))
        yz(l) = dum * (y(l+iz) - y(l-iz))
        zz(l) = dum * (z(l+iz) - z(l-iz))
     end do
     end do
     end do

     ! k boundaries
     !
     if (mydown == mpi_proc_null) then
        do j = ja, jb
        do i = ia, ib
           l = le_idx(i,j,li_ka(n),n)
           xz(l) = dum1 * (x(l+iz) - x(l))
           yz(l) = dum1 * (y(l+iz) - y(l))
           zz(l) = dum1 * (z(l+iz) - z(l))
        end do
        end do
     end if
     
     if (myup == mpi_proc_null) then
        do j = ja, jb
        do i = ia, ib
           l = le_idx(i,j,li_kb(n),n)
           xz(l) = dum1 * (x(l) - x(l-iz))
           yz(l) = dum1 * (y(l) - y(l-iz))
           zz(l) = dum1 * (z(l) - z(l-iz))
        end do
        end do
     end if
     
     ! calculate Jacobian and the inverse transormation
     ! on the interior nodes only; we will communicate the
     ! metrics & jacobian for the ghost nodes
     !
     ka = li_ka(n)
     kb = li_kb(n)
     ja = li_ja(n)
     jb = li_jb(n)
     ia = li_ia(n)
     ib = li_ib(n)

     do k = ka, kb
     do j = ja, jb
     do i = ia, ib
        l = le_idx(i,j,k,n)

        gj=xc(l)*(ye(l)*zz(l)-yz(l)*ze(l))+ &
           xe(l)*(yz(l)*zc(l)-yc(l)*zz(l))+ &
           xz(l)*(yc(l)*ze(l)-ye(l)*zc(l))

        if (gj == 0.0) print*, 'zero jacobian at ',n,i,j,k,gj

        aj(l) = one / gj

        csi(1,l)=aj(l)*(ye(l)*zz(l)-yz(l)*ze(l))
        csi(2,l)=aj(l)*(xz(l)*ze(l)-xe(l)*zz(l))
        csi(3,l)=aj(l)*(xe(l)*yz(l)-xz(l)*ye(l))

        eta(1,l)=aj(l)*(yz(l)*zc(l)-yc(l)*zz(l))
        eta(2,l)=aj(l)*(xc(l)*zz(l)-xz(l)*zc(l))
        eta(3,l)=aj(l)*(xz(l)*yc(l)-xc(l)*yz(l))

        zet(1,l)=aj(l)*(yc(l)*ze(l)-ye(l)*zc(l))
        zet(2,l)=aj(l)*(xe(l)*zc(l)-xc(l)*ze(l))
        zet(3,l)=aj(l)*(xc(l)*ye(l)-xe(l)*yc(l))
     end do
     end do
     end do
  end do

  n = 1

  ka = li_ka(ns)
  kb = li_kb(ns)
  ja = li_ja(ns)
  jb = li_jb(ns)
  ia = li_ia(ns)
  ib = li_ib(ns)

  ! calculate computational time step for des equaions
  if ( turbulence ) then

     do k = ka, kb
     do j = ja, jb
     do i = ia, ib
        l = le_idx(i,j,k,ns)

        Sc = sqrt(xc(l)*xc(l) + yc(l)*yc(l) + zc(l)*zc(l))
        Se = sqrt(xe(l)*xe(l) + ye(l)*ye(l) + ze(l)*ze(l))
        Sz = sqrt(xz(l)*xz(l) + yz(l)*yz(l) + zz(l)*zz(l))

        SL = min(Sc, Se, Sz)

        dtev(l) = SL * min(cfl2(myzone), vnn2(myzone)*ren*SL)

        ! Arnone et al. (1995)
        if ( unsteady ) dtev(l) = min(dtev(l),delti*two/three)

     end do
     end do
     end do

     ! blanking area
     if (nblk /= 0) then
     do nb = 1, nblk
        do k = li_blk_ka(ns,nb), li_blk_kb(ns,nb)
        do j = li_blk_ja(ns,nb), li_blk_jb(ns,nb)
        do i = li_blk_ia(ns,nb), li_blk_ib(ns,nb)
           l = le_idx(i,j,k,ns)

           dtev(l) = zero
           wd(l) = zero
        end do
        end do
        end do
     end do
     end if

  end if


  !==================================================
  !
  ! distance read in grid file not calculated because
  ! of limitations with mpi
  !
  !==================================================
  !
  ! Calculate variables determined only by the
  ! metrics of the grid
  !
  !if (daf) then

     ! calculate computational time step for daf only
     !
     !call met_daf_dtau

     ! calculate model matrices of Jacobian
     !
     !call met_model

  !end if

!contains

  !include 'met_model.F90'

  !include 'met_daf_dtau.F90'
  !include 'met_distance.F90'
  !include 'met_model_matrices.F90'

end subroutine mg_metrics

