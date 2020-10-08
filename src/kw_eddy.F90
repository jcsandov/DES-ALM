
subroutine kw_eddy ( il, iu, jl, ju, kl, ku, & !  1
                     igp, jgp, kgp,          & !  2
                     dc, de, dz,             & !  3
                     dtev,                   & !  4
                     q,                      & !  5
                     qn,                     & !  6
                     qnm1,                   & !  7
                     csi,                    & !  8
                     eta,                    & !  9
                     zet,                    & ! 10
                     aj,                     & ! 11
                     wd,                     & ! 12
                     xnut )                  & ! 13

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, nonorthogonal curvilinear coordinates

  ! Calculate the eddy viscosity according to k-w model

  ! only works if dc = de = dz = 1.0 ; that is, on the fine grid only

  ! input
  ! q(5,ijk)
  ! ucn_j(3,ijk)
  ! aj(i,j,k)
  ! csi(3,ijk), eta(3,ijk), zet(3,ijk)
  ! wd(i,j,k)
  ! itk
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  use global_param
  use global_app
  use global_mpi
  use checksum
  
  implicit none

  integer :: il, &
             jl, &
             kl, &
             iu, &
             ju, &
             ku

  integer :: igp, &
             jgp, &
             kgp

  integer :: i,  j,  k, iterk

  integer :: i_mysta, &
             j_mysta, &
             k_mysta, &
             i_myend, &
             j_myend, &
             k_myend

  ! constants in low-Reynolds number k-w model (Wilcox 1994)
  !
  real (kind = rdf), parameter :: cmu = one
  real (kind = rdf), parameter :: ck  = 0.09_rdf
  real (kind = rdf), parameter :: cw1 = five / nine
  real (kind = rdf), parameter :: cw2 = 0.075_rdf
  real (kind = rdf), parameter :: sk  = two
  real (kind = rdf), parameter :: sw  = two
  
  ! local arrays
  !
  real (kind = rdf), dimension(6,il:iu,jl:ju,kl:ku) :: q
  real (kind = rdf), dimension(3,il:iu,jl:ju,kl:ku) :: csi, eta, zet, ucn_j
  real (kind = rdf), dimension(2,il:iu,jl:ju,kl:ku) :: qn, qnm1, qold, visc, rh, pkw
  real (kind = rdf),   dimension(il:iu,jl:ju,kl:ku) :: xnut, aj, wd, dtev
  real (kind = rdf),   dimension(il:iu,jl:ju,kl:ku) :: ptke, dk, dw, wo, so
  real (kind = rdf), dimension(2,il:iu,jl:ju,kl:ku) :: fv

  real (kind = rdf),   dimension(il:iu,jl:ju,kl:ku) :: As, Ac, Bs

  ! local reals
  !
  real (kind = rdf) :: dc, de, dz
  real (kind = rdf) :: qe5, qe6, dum

  integer :: nm

#ifdef DEBUG
  call init_cksum (igp, jgp, kgp)
#endif

  ! set up like solver_daf.f90
  !
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! interior nodes only; processes on the domain boundaries
  ! 
  if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  if (myup    == mpi_proc_null) k_myend = ku - kgp - 1


  ! save solution vector for iteration (from delta form)  
  ! 
  qold (1:2,:,:,:) = q(5:6,:,:,:)

  ! calculate contravariant velocities (that is, U/J)
  call rhs_contra_j ()

  ! calculate the production of turbulence kinetic energy
  ! also, calculate mean strain rate and mean rate of rotation

  call rhs_kw_ptke ()

#ifdef DEBUG
  call cksum_4d_par ('kw.qold', qold)
  call cksum_3d_par ('kw.dtev', dtev)
  call cksum_3d_par ('kw.xnut', xnut)
  call cksum_3d_par ('kw.ptke', ptke)
  call cksum_3d_par ('kw.wd', wd)
  call cksum_3d_par ('kw.so', so)
  call cksum_3d_par ('kw.wo', wo)

  !call cksum_4d_par ('kw.qold', qold)
#endif
  

  ! update k-omega coefficients As, Ac, Bs
  !
  ! (included it here because I didn't want to include
  ! it in the init routine---I don't want these to be
  ! global variables if we can avoid it --- we can do this
  ! because the boundary condition routine is also not
  ! global for kw model)
  ! 
  call rhs_kw_coef ()

  do iterk = 1, itk     

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_coef'
     call cksum_3d_par ('kw.As', As)
     call cksum_3d_par ('kw.Ac', Ac)
     call cksum_3d_par ('kw.Bs', Bs)     
#endif

     ! update tke and omega on ghost points
     ! 
     if ( iterk > 1 ) call rhs_exchng3_4d (q(5:6,:,:,:))

     ! calculate production-destruction term
     !
     call rhs_kw_prod ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_prod'
     call cksum_4d_par ('kw.pkw', pkw)
#endif

     ! nonlinear kw model by Craft et al., 1995
     !
     if (nlinc .and. craft) call rhs_kw_eddy_craft

     ! update eddy viscosity on ghost cells
     !
     call rhs_exchng3_3d (xnut)

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.eddy_craft'
     call cksum_3d_par ('kw.xnut', xnut)
#endif

     ! calculate viscouse terms
     !
     call rhs_kw_viscous ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_viscous'
     call cksum_4d_par ('kw.visc', visc)
#endif

     ! calculate convective flux 
     !
     call rhs_kw_flux ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_flux'
     call cksum_4d_par ('kw.flux.rh', rh)
#endif

     ! add convective, viscous, production, and unsteady terms to rhs
     !
     call rhs_kw_unst_visc ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_unst_visc'
     call cksum_4d_par ('kw.unst.rh', rh)
#endif

     ! implicit treatment of positive contribution of source term
     !
     call rhs_kw_implicit_source ()

#ifdef DEBUG
  call cksum_3d_par ('kw.imp.dk', dk)
  call cksum_3d_par ('kw.imp.dw', dw)
  call cksum_4d_par ('kw.imp.rh', rh)  
#endif

     ! adi_solution
     !
     !call rhs_kw_adi_solver ()

     ! maximum grid size in each direction
     !
     !nm = max(im, jm, km)

     nm = max(i_myend - i_mysta + 3, &
              j_myend - j_mysta + 3, &
              k_myend - k_mysta + 3  )

     call rhs_kw_adi_solver_ijk ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_adi_solver_ij'
     call cksum_4d_par ('kw.adi.ijk.rh', rh)
#endif

     ! update
     !
     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
     do i = i_mysta, i_myend

        qe5 = qold(1,i,j,k) + rh(1,i,j,k)
        qe6 = qold(2,i,j,k) + rh(2,i,j,k)

        dum = one + onept5 * e_source * dtev(i,j,k) / delti

        q(5,i,j,k) = q(5,i,j,k) + (-q(5,i,j,k) + qe5) / dum
        q(6,i,j,k) = q(6,i,j,k) + (-q(6,i,j,k) + qe6) / dum

     end do
     end do
     end do

     ! boundary condition
     !
#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.k & omega n+1'
     call cksum_4d_par ('kw.new.q', q)
#endif

     ! boundary conditions
     ! 
     call rhs_kw_bcond ()


#ifdef DEBUG
     !print *, 'myid =', myid, 'kw_eddy.kw_bcond', ' itk = ', iterk
     call cksum_4d_par ('kw.bcond', q)
#endif

     ! save solution for iterations
     !

     qold (1:2,:,:,:) = q(5:6,:,:,:)

!==================================================
!
! MPI need exchange for q(5:6,:,:,:) at least
!
!==================================================


  end do

contains

  include 'rhs_contra_j.F90'
  include 'rhs_kw_ptke.F90'
  include 'rhs_kw_prod.F90'
  include 'rhs_kw_viscous.F90'
  include 'rhs_kw_flux.F90'
  include 'rhs_kw_unst_visc.F90'
  include 'rhs_kw_implicit_source.F90'
  include 'rhs_kw_eddy_craft.F90'
  include 'rhs_kw_adi_solver_ijk.F90'
  include 'rhs_kw_bcond.F90'
  include 'rhs_kw_coef.F90'
  include 'rhs_exchng3_4d.F90'
  include 'rhs_exchng3_3d.F90'

end subroutine kw_eddy

!include 'rhs_kw_adi_solver.F90'
