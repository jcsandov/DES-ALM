
subroutine met_daf_dtau
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculates dtau(ijk) at each node (local, pseudo time stepping)
  ! eigenvalues of Jacobian matrices

  ! input
  !     ren
  !     cfl1
  !     vnn1
  !     csi(3,ijk)
  !     aj(ijk)
  !     ucn(3,ijk)

  ! output
  !     dtau(ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  real (kind = rdf) :: e1, e2, e3
  real (kind = rdf) :: g1, g2, g3
  real (kind = rdf) :: dti1, dti2, dti3
  real (kind = rdf) :: dtv1, dtv2, dtv3
  real (kind = rdf) :: dcsq, desq, dzsq

  real (kind = rdf) :: cfl, vnn
  real (kind = rdf) :: dtemp

  integer :: ka
  integer :: kb
  integer :: ja
  integer :: jb
  integer :: ia
  integer :: ib

  integer :: n
  integer :: l
  integer :: i
  integer :: j
  integer :: k

  cfl = cfl1
  vnn = vnn1

  ! we have made the inverse of the grid spacing
  ! universal; therefore, dsp = dc^(-1) in old code
  ! this caused a bug in the translations

  do n = ns, ng
     ! Arnone et al. (1995)
     dtemp = delti / (one_pt_five * two ** (ng - 1))

     dcsq = dc(n) * dc(n)
     desq = de(n) * de(n)
     dzsq = dz(n) * dz(n)

     ! Note, that each direction is treated independently

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

        g1  = csi(1,l) * csi(1,l) + &
              csi(2,l) * csi(2,l) + &
              csi(3,l) * csi(3,l)
        g2  = eta(1,l) * eta(1,l) + &
              eta(2,l) * eta(2,l) + &
              eta(3,l) * eta(3,l)
        g3  = zet(1,l) * zet(1,l) + &
              zet(2,l) * zet(2,l) + &
              zet(3,l) * zet(3,l)

        e1 = sqrt(g1)
        e2 = sqrt(g2)
        e3 = sqrt(g3)

        dti1 = cfl / e1 / dc(n)
        dti2 = cfl / e2 / de(n)
        dti3 = cfl / e3 / dz(n)

        dtv1 = vnn * ren / dcsq / g1
        dtv2 = vnn * ren / desq / g2
        dtv3 = vnn * ren / dzsq / g3

        dtau(l) = min(dti1, dti2, dti3, dtv1, dtv2, dtv3)
        !if (unsteady) dtau(l) = min(dtau(l), dtemp)

     end do
     end do
     end do

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
        do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
        do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
           l = le_idx(i,j,k,n)
           dtau(l) = zero
        end do
        end do
        end do
     end do
     end if

  end do

end subroutine met_daf_dtau


