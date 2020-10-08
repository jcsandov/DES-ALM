subroutine rhs_modal_matrices
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate matrics of the geometric transformation for
  ! for the fine grid and each of the coarse grids

  ! input
  !     csi(1:3,nmxg) csi_x, csi_y and csi_z
  !     eta(1:3,nmxg) eta_x, eta_y and eta_z
  !     zet(1:3,nmxg) zet_x, zet_y and zet_z

  ! output
  !     mai(4,4,nmxg)
  !     mc (4,4,nmxg)
  !     n1i(4,4,nmxg)
  !     n2i(4,4,nmxg)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local dummy reals
  real(kind=rdf) :: g1, g2, g3
  real(kind=rdf) :: sg1, sg2, sg3
  real(kind=rdf) :: alpha, sigma, theta

  real(kind=rdf) :: csi1, csi2, csi3
  real(kind=rdf) :: eta1, eta2, eta3
  real(kind=rdf) :: zet1, zet2, zet3

  real(kind=rdf) :: csi1s, csi2s, csi3s
  real(kind=rdf) :: eta1s, eta2s, eta3s
  real(kind=rdf) :: zet1s, zet2s, zet3s

  real(kind=rdf) :: scsi, seta, szet
  real(kind=rdf) :: dcsi, deta, dzet
  real(kind=rdf) :: rc1, rc2, rc3
  real(kind=rdf) :: re1, re2, re3
  real(kind=rdf) :: rz1, rz2, rz3

  real(kind=rdf), dimension(4,4) :: ma
  real(kind=rdf), dimension(4,4) :: mb
  real(kind=rdf), dimension(4,4) :: mbi
  real(kind=rdf), dimension(4,4) :: mci
  real(kind=rdf), dimension(4,4) :: rp

  real(kind=rdf) :: ssgd, sgd2

  ! model matrices of Jacobian matrices
  !
  do k = kl, ku
  do j = jl, ju
  do i = il, iu

     csi1=csi(1,i,j,k)
     csi2=csi(2,i,j,k)
     csi3=csi(3,i,j,k)
     eta1=eta(1,i,j,k)
     eta2=eta(2,i,j,k)
     eta3=eta(3,i,j,k)
     zet1=zet(1,i,j,k)
     zet2=zet(2,i,j,k)
     zet3=zet(3,i,j,k)

     csi1s=csi1*csi1
     csi2s=csi2*csi2
     csi3s=csi3*csi3
     eta1s=eta1*eta1
     eta2s=eta2*eta2
     eta3s=eta3*eta3
     zet1s=zet1*zet1
     zet2s=zet2*zet2
     zet3s=zet3*zet3
 
     g1=csi1s+csi2s+csi3s
     g2=eta1s+eta2s+eta3s
     g3=zet1s+zet2s+zet3s

     sg1=sqrt(g1)
     sg2=sqrt(g2)
     sg3=sqrt(g3)

     alpha=one+e_source*onept5*dtau(i,j,k)/delti
     sigma=sqrt(beta*alpha)
     theta=sqrt(beta/alpha)

     ! spectral radius
     spr(1,i,j,k)=theta*sg1
     spr(2,i,j,k)=theta*sg2
     spr(3,i,j,k)=theta*sg3

     ! base
     scsi=sqrt(two*g1+(csi2-csi3)**two)
     seta=sqrt(two*g2+(eta2-eta3)**two)
     szet=sqrt(two*g3+(zet2-zet3)**two)

     dcsi=g1+two*csi2*csi3+csi1*csi1
     deta=g2+two*eta2*eta3+eta1*eta1
     dzet=g3+two*zet2*zet3+zet1*zet1

     rc1=csi2 +csi3
     rc2=csi2 -csi3
     rc3=csi2s-csi3s

     re1=eta2 +eta3
     re2=eta2 -eta3
     re3=eta2s-eta3s

     rz1=zet2 +zet3
     rz2=zet2 -zet3
     rz3=zet2s-zet3s

     ! Ma
     !
     ma(1,1)=zero
     ma(1,2)=zero
     ma(1,3)=sigma
     ma(1,4)=sigma

     ma(2,1)=-(csi2 + csi3) / scsi
     ma(2,2)= csi1*(csi3-csi2) / (scsi * sg1)
     ma(2,3)= csi1 / sg1
     ma(2,4)=-csi1 / sg1

     ma(3,1)= csi1 / scsi
     ma(3,2)= (g1 + csi2*(csi3-csi2)) / (scsi * sg1)
     ma(3,3)= csi2 / sg1
     ma(3,4)=-csi2 / sg1

     ma(4,1)= csi1 / scsi
     ma(4,2)=-(g1 - csi3*(csi3-csi2)) / (scsi * sg1)
     ma(4,3)= csi3 / sg1
     ma(4,4)=-csi3 / sg1

     ! Mai (= Ma^-1) matrix
     !

     mai(1,1,i,j,k)= zero
     mai(1,2,i,j,k)=-scsi*(csi2 + csi3) / dcsi
     mai(1,3,i,j,k)= scsi*csi1 / dcsi
     mai(1,4,i,j,k)= mai(1,3,i,j,k)

     ssgd=scsi / sg1 / dcsi
     sgd2=one / sg1 / dcsi / two

     mai(2,1,i,j,k)= zero
     mai(2,2,i,j,k)= csi1*(csi3 - csi2) * ssgd
     mai(2,3,i,j,k)= (csi2*csi3 + csi1s + csi3s) * ssgd
     mai(2,4,i,j,k)=-(csi2*csi3 + csi1s + csi2s) * ssgd

     mai(3,1,i,j,k)= pt5 / sigma
     mai(3,2,i,j,k)= csi1 * (two*g1 + two*csi2*csi3 - csi2s - csi3s) * sgd2
     mai(3,3,i,j,k)= (g1*rc1 + csi1s*rc2 + csi3*rc3) * sgd2
     mai(3,4,i,j,k)= (g1*rc1 - csi1s*rc2 - csi2*rc3) * sgd2

     mai(4,1,i,j,k)= pt5 / sigma
     mai(4,2,i,j,k)= -mai(3,2,i,j,k)
     mai(4,3,i,j,k)= -mai(3,3,i,j,k)
     mai(4,4,i,j,k)= -mai(3,4,i,j,k)

     ! Mb
     !
     mb(1,1)=zero
     mb(1,2)=zero
     mb(1,3)=sigma
     mb(1,4)=sigma

     mb(2,1)=-(eta2 + eta3) / seta
     mb(2,2)= eta1*(eta3-eta2) / (seta * sg2)
     mb(2,3)= eta1 / sg2
     mb(2,4)=-eta1 / sg2

     mb(3,1)= eta1 / seta
     mb(3,2)= (g2 + eta2*(eta3-eta2)) / (seta * sg2)
     mb(3,3)= eta2 / sg2
     mb(3,4)=-eta2 / sg2

     mb(4,1)= eta1 / seta
     mb(4,2)=-(g2 - eta3*(eta3-eta2)) / (seta * sg2)
     mb(4,3)= eta3 / sg2
     mb(4,4)=-eta3 / sg2

     ! Mbi (= Mb^-1) matrix
     !

     mbi(1,1)= zero
     mbi(1,2)=-seta*(eta2 + eta3) / deta
     mbi(1,3)= seta* eta1 / deta
     mbi(1,4)= mbi(1,3)

     ssgd=seta / sg2 / deta
     sgd2=one / sg2 / deta / two

     mbi(2,1)= zero
     mbi(2,2)= eta1 * (eta3 - eta2) * ssgd
     mbi(2,3)= (eta2*eta3 + eta1s + eta3s) * ssgd
     mbi(2,4)=-(eta2*eta3 + eta1s + eta2s) * ssgd

     mbi(3,1)= pt5 / sigma
     mbi(3,2)= eta1 * (two*g2 + two*eta2*eta3 - eta2s - eta3s) * sgd2
     mbi(3,3)= (g2*re1 + eta1s*re2 + eta3*re3) * sgd2
     mbi(3,4)= (g2*re1 - eta1s*re2 - eta2*re3) * sgd2

     mbi(4,1)= pt5 / sigma
     mbi(4,2)= -mbi(3,2)
     mbi(4,3)= -mbi(3,3)
     mbi(4,4)= -mbi(3,4)

     ! N1^-1 (= Mb^-1*Ma)
     !
     rp(1:4,1:4) = matmul(mbi,ma)
     n1i(1:4,1:4,i,j,k) = rp(1:4,1:4)

     ! Mc
     !
     mc(1,1,i,j,k)= zero
     mc(1,2,i,j,k)= zero
     mc(1,3,i,j,k)= sigma
     mc(1,4,i,j,k)= sigma

     mc(2,1,i,j,k)=-(zet2 + zet3) / szet
     mc(2,2,i,j,k)= zet1*(zet3-zet2) / (szet * sg3)
     mc(2,3,i,j,k)= zet1 / sg3
     mc(2,4,i,j,k)=-zet1 / sg3

     mc(3,1,i,j,k)= zet1 / szet
     mc(3,2,i,j,k)= (g3 + zet2*(zet3-zet2)) / (szet * sg3)
     mc(3,3,i,j,k)= zet2 / sg3
     mc(3,4,i,j,k)=-zet2 / sg3

     mc(4,1,i,j,k)= zet1 / szet
     mc(4,2,i,j,k)=-(g3 - zet3*(zet3-zet2)) / (szet * sg3)
     mc(4,3,i,j,k)= zet3 / sg3
     mc(4,4,i,j,k)=-zet3 / sg3

     ! Mci (= Mc^-1) matrix
     !

     mci(1,1)= zero
     mci(1,2)=-szet*(zet2 + zet3) / dzet
     mci(1,3)= szet* zet1 / dzet
     mci(1,4)= mci(1,3)

     ssgd=szet / sg3 / dzet
     sgd2=one / sg3 / dzet / two

     mci(2,1)= zero
     mci(2,2)= zet1*(zet3 - zet2) * ssgd
     mci(2,3)= (zet2*zet3 + zet1s + zet3s) * ssgd
     mci(2,4)=-(zet2*zet3 + zet1s + zet2s) * ssgd

     mci(3,1)= pt5 / sigma
     mci(3,2)= zet1 * (two*g3 + two*zet2*zet3 - zet2s - zet3s) * sgd2
     mci(3,3)= (g3*rz1 + zet1s*rz2 + zet3*rz3) * sgd2
     mci(3,4)= (g3*rz1 - zet1s*rz2 - zet2*rz3) * sgd2

     mci(4,1)= pt5 / sigma
     mci(4,2)= -mci(3,2)
     mci(4,3)= -mci(3,3)
     mci(4,4)= -mci(3,4)

     ! N2i = Mc^-1 * Mb
     !
     rp(1:4,1:4) = matmul(mci, mb)
     n2i(1:4,1:4,i,j,k) = rp(1:4,1:4)

   end do
   end do
   end do

end subroutine rhs_modal_matrices

