
subroutine nonreflect_ibc (inout, jl, ju, kl, ku, &
                           igp, jgp, kgp, &
                            dtau, &
                            csi, &
                            ucn_j, &
                            aj, &
                            q, &
                            rh)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Non-reflecting implicit boundary condition based on MOC

  ! Characteristics-based boundary conditions for all boundaries
  ! of the computational domain. Requires the calculation of the
  ! right-hand side for the planes and edges of the computational
  ! domain (but not the corners).

  ! LODI (Locally one-dimensional inviscid) boundary conditions

  ! Nonorthogonal curvilinear coordinates

  ! scj, 28 Jun 2002 added non-orthogonal terms
  ! scj, 30 May 2002 all boundaries included
  ! jp,  25 Jan 2004 modified for implicit scheme

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global_param
  use global_mpi

  implicit none

  ! extents
  ! -------
  integer :: il, jl, kl, iu, ju, ku

  ! ghost points
  ! ------------
  integer, intent(in) :: igp, jgp, kgp

  integer, intent(in) :: inout

  ! main variables
  ! --------------
  real (kind = rdf), dimension(jl:ju,kl:ku) ,intent(in) :: dtau, aj, ucn_j
  real (kind = rdf), dimension(1:3,jl:ju,kl:ku), intent(in) :: csi
  real (kind = rdf), dimension(1:4,1:2,jl:ju,kl:ku), intent(in) :: q
  real (kind = rdf), dimension(1:4,jl:ju,kl:ku), intent(inout)  :: rh

  integer :: i_mysta, j_mysta, k_mysta
  integer :: i_myend, j_myend, k_myend
  real (kind = rdf) :: dc, de, dz
  integer :: i, j, k

  ! local constants
  real (kind = rdf) :: lam1, lam2, lam3, lam4
  real (kind = rdf) :: L1, L2, L3, L4
  real (kind = rdf) :: u, v, w
  real (kind = rdf) :: cj
  real (kind = rdf) :: c1, c2, c3
  real (kind = rdf) :: ucon, gjj, sg

  real (kind = rdf) :: rh1, rh2, rh3, rh4

  real (kind = rdf), dimension(1:4) :: l_gen
  real (kind = rdf), dimension(1:4) :: dq

  ! process boundaries
  ! ------------------
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! processes on the domain boundaries
  ! ----------------------------------
  if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  if (myup    == mpi_proc_null) k_myend = ku - kgp - 1


  ! i = constant sides
  !
  if (inout == 1 .or. inout == 2) then
        i=inout
     do k=k_mysta,k_myend
     do j=j_mysta,j_myend
          !call i_boundary 

          u=q(2,i,j,k)
          v=q(3,i,j,k)
          w=q(4,i,j,k)
          ucon = ucn_j(j,k) * aj(j,k)
          c1=csi(1,j,k)
          c2=csi(2,j,k)
          c3=csi(3,j,k)
          gjj=c1*c1+c2*c2+c3*c3 
          sg =sqrt(gjj)
          cj =sqrt(ucon * ucon + gjj)
          if (i == 1) then 
             dq(1)=(q(1,i+1,j,k)-q(1,i,j,k))
             dq(2)=(q(2,i+1,j,k)-q(2,i,j,k))
             dq(3)=(q(3,i+1,j,k)-q(3,i,j,k))
             dq(4)=(q(4,i+1,j,k)-q(4,i,j,k))
             call calc_L_general
             call calc_L_min
          else
             dq(1)=(q(1,i,j,k)-q(1,i-1,j,k))
             dq(2)=(q(2,i,j,k)-q(2,i-1,j,k))
             dq(3)=(q(3,i,j,k)-q(3,i-1,j,k))
             dq(4)=(q(4,i,j,k)-q(4,i-1,j,k))
             call calc_L_general
             call calc_L_max
          end if

         call calc_rh

     end do
     end do
  end if

contains

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine calc_L_general

    ! Routine is general for all three directions

    real (kind = rdf), dimension(4,4) :: si  ! s inverse
    real (kind = rdf) :: d1, d2, d3, d4, d5, d6
    real (kind = rdf) :: kj

    lam1 = ucon - cj  ! lam 1 < 0 by defn. always
    lam2 = ucon + cj  ! lam 2 > 0 by defn.
    lam3 = ucon
    lam4 = ucon

    kj = c1 + c2 + c3

    d1 = (c2 - c1 + lam1 * (v - u)) / kj
    d2 = (c2 - c1 + lam2 * (v - u)) / kj
    d3 = (c3 - c1 + lam1 * (w - u)) / kj
    d4 = (c3 - c1 + lam2 * (w - u)) / kj

    d5 = (c2 - c1 + lam3 * (v - u)) / kj * two
    d6 = (c3 - c1 + lam4 * (w - u)) / kj * two
  
    ! inverse of the modal matrix
    si(1,1) = -lam2
    si(2,1) = -lam1
    si(3,1) =  lam2 * d1 + lam1 * d2
    si(4,1) =  lam2 * d3 + lam1 * d4

    si(1,2) = c1
    si(2,2) = c1
    si(3,2) = -two * cj * cj / kj - c1 * d5
    si(4,2) = -two * cj * cj / kj - c1 * d6

    si(1,3) =  c2
    si(2,3) =  c2
    si(3,3) =  two * cj * cj / kj - c2 * d5
    si(4,3) = -c2 * d6

    si(1,4) =  c3
    si(2,4) =  c3
    si(3,4) = -c3 * d5
    si(4,4) =  two * cj * cj / kj - c3 * d6

    si = si * sg / (two * cj * cj)      ! an element-by-element
                                        ! multiplication
    l_gen = matmul(si, dq)

    l_gen(1) = lam1 * l_gen(1)
    l_gen(2) = lam2 * l_gen(2)
    l_gen(3) = lam3 * l_gen(3)
    l_gen(4) = lam4 * l_gen(4)

  end subroutine calc_L_general

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine calc_L_min
    
    if (ucon < zero) then
       L1 = l_gen(1)
       L2 = zero
       L3 = l_gen(3)
       L4 = l_gen(4)
    else               ! ucon >= zero
       L1 = l_gen(1)
       L2 = zero
       L3 = zero
       L4 = zero
    end if

  end subroutine calc_L_min
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine calc_L_max

    if (ucon < zero) then
       L1 = zero
       L2 = l_gen(2)
       L3 = zero
       L4 = zero
    else
       L1 = zero
       L2 = l_gen(2)
       L3 = l_gen(3)
       L4 = l_gen(4)
    end if

  end subroutine calc_L_max

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine calc_rh

    ! most variables defined in calling routine

    ! local variables
    real (kind = rdf) :: rh1, rh2, rh3, rh4

    rh1=(-cj        *L1 +  cj        *L2                              ) / sg
    rh2=((c1+u*lam1)*L1 + (c1+u*lam2)*L2 -  c2      *L3 -  c3      *L4) / sg
    rh3=((c2+v*lam1)*L1 + (c2+v*lam2)*L2 + (c1 + c3)*L3 -  c3      *L4) / sg
    rh4=((c3+w*lam1)*L1 + (c3+w*lam2)*L2 -  c2      *L3 + (c1 + c2)*L4) / sg

    rh(1,j,k)=rh(1,j,k)+rh1/aj(j,k)
    rh(2,j,k)=rh(2,j,k)+rh2/aj(j,k)
    rh(3,j,k)=rh(3,j,k)+rh3/aj(j,k)
    rh(4,j,k)=rh(4,j,k)+rh4/aj(j,k)

  end subroutine calc_rh

end subroutine nonreflect_ibc



