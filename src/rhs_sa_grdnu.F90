subroutine rhs_sa_grdnu ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, nonorthogonal curvilinear coodinates

  ! Calculate gradient of eddy viscosity in S-A model
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! local dummy variables
  real (kind = rdf) :: dc2, de2, dz2
  real (kind = rdf) :: dnu1, dnu2, dnu3
  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku) :: q5j
  real (kind = rdf) :: aj2

  ! dum --> dc2 or de2 or dz2
  dc2=pt5*dc
  de2=pt5*de
  dz2=pt5*dz 

  do k=kl,ku
  do j=jl,ju
  do i=il,iu
     q5j(i,j,k)=q(5,i,j,k)/aj(i,j,k)
  end do
  end do
  end do


  ! interior nodes only
  do k=k_mysta,k_myend
  do j=j_mysta,j_myend
  do i=i_mysta,i_myend

     ! csi direction
             
     dnu1=dc2*(q5j(i+1,j,k)*csi(1,i+1,j,k)-q5j(i-1,j,k)*csi(1,i-1,j,k))
     dnu2=dc2*(q5j(i+1,j,k)*csi(2,i+1,j,k)-q5j(i-1,j,k)*csi(2,i-1,j,k))
     dnu3=dc2*(q5j(i+1,j,k)*csi(3,i+1,j,k)-q5j(i-1,j,k)*csi(3,i-1,j,k))

     ! eta direction

     dnu1=dnu1+de2*(q5j(i,j+1,k)*eta(1,i,j+1,k)-q5j(i,j-1,k)*eta(1,i,j-1,k))
     dnu2=dnu2+de2*(q5j(i,j+1,k)*eta(2,i,j+1,k)-q5j(i,j-1,k)*eta(2,i,j-1,k))
     dnu3=dnu3+de2*(q5j(i,j+1,k)*eta(3,i,j+1,k)-q5j(i,j-1,k)*eta(3,i,j-1,k))

     ! zet direction

     dnu1=dnu1+dz2*(q5j(i,j,k+1)*zet(1,i,j,k+1)-q5j(i,j,k-1)*zet(1,i,j,k-1))
     dnu2=dnu2+dz2*(q5j(i,j,k+1)*zet(2,i,j,k+1)-q5j(i,j,k-1)*zet(2,i,j,k-1))
     dnu3=dnu3+dz2*(q5j(i,j,k+1)*zet(3,i,j,k+1)-q5j(i,j,k-1)*zet(3,i,j,k-1))

     aj2=aj(i,j,k)*aj(i,j,k)

     dnu(i,j,k) = aj2*(dnu1*dnu1 + dnu2*dnu2 + dnu3*dnu3)

  end do
  end do
  end do


end subroutine rhs_sa_grdnu


