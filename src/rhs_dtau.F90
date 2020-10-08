
subroutine rhs_dtau
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
  use global_param
  implicit none
  
  real (kind = rdf) :: e1, e2, e3
  real (kind = rdf) :: g1, g2, g3
  real (kind = rdf) :: abu, abv, abw
  real (kind = rdf) :: dti1, dti2, dti3
  real (kind = rdf) :: dtv1, dtv2, dtv3
  real (kind = rdf) :: dcsq, desq, dzsq

  real (kind = rdf) :: cfl, vnn
  real (kind = rdf) :: dtemp

  ! we have made the inverse of the grid spacing
  ! universal; therefore, dsp = dc^(-1) in old code
  ! this caused a bug in the translations

  dcsq = dc * dc
  desq = de * de
  dzsq = dz * dz

  cfl = cfl1
  vnn = vnn1

  ! Arnone et al (1995)
  dtemp = delti / (one_pt_five * two ** (3 - 1))

  ! Note, that each direction is treated independently
  do k = 1, km
  do j = 1, jm
  do i = 1, im
     g1  = csi(1,i,j,k) * csi(1,i,j,k) + csi(2,i,j,k) * csi(2,i,j,k) &
         + csi(3,i,j,k) * csi(3,i,j,k)
     g2  = eta(1,i,j,k) * eta(1,i,j,k) + eta(2,i,j,k) * eta(2,i,j,k) &
         + eta(3,i,j,k) * eta(3,i,j,k)
     g3  = zet(1,i,j,k) * zet(1,i,j,k) + zet(2,i,j,k) * zet(2,i,j,k) &
         + zet(3,i,j,k) * zet(3,i,j,k)

     abu = abs(ucn_j(1,i,j,k) * aj(i,j,k))
     abv = abs(ucn_j(2,i,j,k) * aj(i,j,k))
     abw = abs(ucn_j(3,i,j,k) * aj(i,j,k))
     
     e1 = abu + sqrt(abu * abu + g1)
     e2 = abv + sqrt(abv * abv + g2)
     e3 = abw + sqrt(abw * abw + g3)

     dti1 = cfl / e1 / dc
     dti2 = cfl / e2 / de
     dti3 = cfl / e3 / dz
     dtv1 = vnn * ren / dcsq / g1
     dtv2 = vnn * ren / desq / g2
     dtv3 = vnn * ren / dzsq / g3

     dtau(i,j,k) = min(dti1, dti2, dti3, dtv1, dtv2, dtv3)
!!$     dtau(i,j,k) = min(dti1, dti2, dti3)
     if (unsteady) dtau(i,j,k) = min(dtau(i,j,k), dtemp)

  end do
  end do
  end do

end subroutine rhs_dtau


