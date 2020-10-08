
subroutine rhs_sa_prod ()
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, nonorthogonal curvilinear coordinates

  ! Calculate rhs terms of S-A (1994), 
  ! which represent the production of eddy viscosity,
  !                             turbulence diffusion,
  !             and near-wall turbulence dustruction.

  ! terms of fv2 & fv3 are based on Squires et al.
  ! (AIAA Paper 02-1021) except trim terms (ft1 & ft2)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! cut of limit for production term
  real (kind = rdf), parameter :: vtc_col = 0.00001_rdf
  real (kind = rdf), parameter :: one_sixth = one / six

  real (kind = rdf) :: rei
  real (kind = rdf) :: chi, chic
  real (kind = rdf) :: fv1, fv2, fv3
  real (kind = rdf) :: vtcm, cawd, vtc_org, vtc
  real (kind = rdf) :: rr, rr6 
  real (kind = rdf) :: grr_org, grr, grr6
  real (kind = rdf) :: vwd

  rei = one / ren

  do k=k_mysta,k_myend
  do j=j_mysta,j_myend
  do i=i_mysta,i_myend
     
     if ( wd(i,j,k) == zero ) then

     psa(i,j,k) = zero

     else
     !
     ! production term
     !
     chi  = q(5,i,j,k) / rei
     chic = chi * chi * chi
     fv1  = chic / (chic + cv1c)
     fv2  = (one + chi / cv2) ** (-three)
     fv3  = ( (one + chi * fv1) * (one - fv2) ) / chi

     ! org: fv2  = one - chi / (one + chi * fv1)
     ! ft2 = ct3 * exp(-ct4 * chi * chi)

     ! original version by Spalart-Allmaras (1994)
     vtcm = sqrt(two*wo(i,j,k))

     !  modified version by Dacles-Mariani et al. (1995)
     !  vtcw=sqrt(two*wo(i,j,k))
     !  vtcs=sqrt(two*so(i,j,k))
     !  asmw=vtcs-vtcw
     !  vtcm=vtcw + two * min(zero,asmw)

     cawd    = capa * capa * wd(i,j,k) * wd(i,j,k)
     vtc_org = fv3 * vtcm + q(5,i,j,k) * fv2 / cawd
     ! org: vtc_org = vtcm + q(5,i,j,k) * fv2 / cawd

     ! production term, S, is limited to be positive, with a cut-off limit
     vtc     = max(vtc_org, vtc_col)

     psa(i,j,k) = cb1 * vtc * q(5,i,j,k)
     ! psa(i,j,k) = cb1 * (one - ft2) * vtc * q(5,i,j,k)

     !
     ! destruction term
     !

     rr  = q(5,i,j,k) / (vtc * cawd)
     rr6 = rr * rr * rr * rr * rr * rr

     grr_org   = rr + cw2 * (rr6 - rr)
     ! the closure function, g, cannot be negative
     grr  = max(zero, grr_org)
     grr6 = grr * grr * grr * grr * grr * grr

     fw(i,j,k) = grr * ((one + cw3h) / (grr6 + cw3h))**one_sixth

     vwd = q(5,i,j,k) / wd(i,j,k)

     psa(i,j,k) = psa(i,j,k) - cw1 * fw(i,j,k) * vwd * vwd

     !
     ! turbulence diffustion term
     !

     psa(i,j,k) = psa(i,j,k) + cb2 * dnu(i,j,k) / sig

     end if

  end do
  end do
  end do


end subroutine rhs_sa_prod


