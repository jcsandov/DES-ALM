
subroutine rhs_daf_dtau
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
  real (kind = rdf) :: rei, ret

  cfl=cfl1(myzone)
  vnn=vnn1(myzone)
  rei=one / ren

  ! we have made the inverse of the grid spacing
  ! universal; therefore, dsp = dc^(-1) in old code
  ! this caused a bug in the translations

  ! Arnone et al. (1995)
  dtemp = delti / (one_pt_five * two ** (3 - 1))

  dcsq=dc*dc
  desq=de*de
  dzsq=dz*dz

  ! Note, that each direction is treated independently

  do k=kl,ku
  do j=jl,ju
  do i=il,iu
     g1=csi(1,i,j,k)*csi(1,i,j,k) + &
        csi(2,i,j,k)*csi(2,i,j,k) + &
        csi(3,i,j,k)*csi(3,i,j,k)
     g2=eta(1,i,j,k)*eta(1,i,j,k) + &
        eta(2,i,j,k)*eta(2,i,j,k) + &
        eta(3,i,j,k)*eta(3,i,j,k)
     g3=zet(1,i,j,k)*zet(1,i,j,k) + &
        zet(2,i,j,k)*zet(2,i,j,k) + &
        zet(3,i,j,k)*zet(3,i,j,k)

     e1=sqrt(g1)
     e2=sqrt(g2)
     e3=sqrt(g3)

     dti1=cfl/e1/dc
     dti2=cfl/e2/de
     dti3=cfl/e3/dz

     ret=rei+xnut(i,j,k)

     dtv1=vnn/ret/dcsq/g1
     dtv2=vnn/ret/desq/g2
     dtv3=vnn/ret/dzsq/g3

     dtau(i,j,k)=min(dti1,dti2,dti3,dtv1,dtv2,dtv3)
     !if (unsteady) dtau(i,j,k) = min(dtau(i,j,k), dtemp)

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
        dtau(i,j,k) = zero
     end do
     end do
     end do
  end do
  end if


end subroutine rhs_daf_dtau


