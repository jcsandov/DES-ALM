
subroutine des_eddy (il, iu, jl, ju, kl, ku, & !  1
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
                     xnut )                    ! 13


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


  ! constants in Spalart-Allmaras one eq. eddy viscosity model
  ! cw1 = cb1/(kappa*kappa) + (1 + cb2)/sig
  real (kind = rdf), parameter :: cb1  = 0.1355_rdf
  real (kind = rdf), parameter :: cb2  = 0.622_rdf
  real (kind = rdf), parameter :: sig  = 2.0_rdf / 3.0_rdf
  real (kind = rdf), parameter :: capa = 0.41_rdf
  real (kind = rdf), parameter :: cw1  = 3.23906782_rdf
  real (kind = rdf), parameter :: cw2  = 0.3_rdf
  real (kind = rdf), parameter :: cw3  = 2.0_rdf
  real (kind = rdf), parameter :: cw3h = 64.0_rdf
  real (kind = rdf), parameter :: cv1  = 7.1_rdf
  real (kind = rdf), parameter :: cv1c = 357.911_rdf
  real (kind = rdf), parameter :: cv2  = 5.0_rdf
  
  ! local arrays
  !
  real (kind = rdf), dimension(1:5,il:iu,jl:ju,kl:ku) :: q
  real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) :: csi, eta, zet
  real (kind = rdf), dimension(    il:iu,jl:ju,kl:ku) :: qn, qnm1, rh
  real (kind = rdf),     dimension(il:iu,jl:ju,kl:ku) :: xnut, aj, wd, dtev
  ! eliminated fv(1:2,il:iu ... from kw_eddy

  ! allocatable variables
  real (kind = rdf), dimension(:,:,:,:), allocatable :: ucn_j
  real (kind = rdf), dimension(:,:,:),   allocatable :: qold, visc
  real (kind = rdf), dimension(:,:,:),   allocatable :: dnu, psa, dk, wo, fw
  real (kind = rdf), dimension(:,:,:),   allocatable :: fv

  ! local reals
  !
  real (kind = rdf) :: dc, de, dz
  real (kind = rdf) :: qe5, dum

#ifdef DEBUG
  call init_cksum (igp, jgp, kgp)
#endif

  ! allocate varibles 
  allocate (ucn_j(1:3,il:iu,jl:ju,kl:ku), &
                 qold(il:iu,jl:ju,kl:ku), &
                 visc(il:iu,jl:ju,kl:ku), &
                  dnu(il:iu,jl:ju,kl:ku), &
                  psa(il:iu,jl:ju,kl:ku), &
                   dk(il:iu,jl:ju,kl:ku), &
                   wo(il:iu,jl:ju,kl:ku), &
                   fw(il:iu,jl:ju,kl:ku), &
                   fv(il:iu,jl:ju,kl:ku) )

  ucn_j=zero; qold=zero; visc=zero; dnu=zero; psa=zero
  dk=zero; wo=zero; fw=zero; fv=zero

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
  qold (:,:,:) = q(5,:,:,:)

  ! calculate contravariant velocities (that is, U/J)
  call rhs_contra_j ()

  ! calculate the production of turbulence kinetic energy
  ! also, calculate mean strain rate and mean rate of rotation

  call rhs_sa_source ()

#ifdef DEBUG
  call cksum_3d_par ('sa.qold', qold)
  call cksum_3d_par ('sa.dtev', dtev)
  call cksum_3d_par ('sa.xnut', xnut)
  call cksum_3d_par ('sa.wd', wd)
  call cksum_3d_par ('sa.wo', wo)

  !$call cksum_3d_par ('sa.ptke', ptke)
  !$call cksum_3d_par ('sa.so', so)
  !$call cksum_3d_par ('sa.qold', qold)
#endif
  
  do iterk = 1, itk     

     ! update tke and omega on ghost points
     ! 
     if ( iterk > 1 ) call rhs_exchng3_4d (q(5:5,:,:,:))

     ! calculate production-destruction term
     !
     call rhs_sa_grdnu ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'sa_eddy.sa_grdnu'
     call cksum_3d_par ('sa.dnu', dnu)
#endif

     ! calculate production-destruction term
     !
     call rhs_sa_prod ()

#ifdef DEBUG
     !print *, 'myid =', myid, 'sa_eddy.sa_prod'
     call cksum_3d_par ('sa.psa', psa)
#endif

#ifdef DEBUG
!      if (iterk == 1) then
!   ! dump rh
!   ! 
!   open (unit = 11, file = 'sa.par.psa', form = 'formatted')
!   do k = k_mysta - 1, k_myend + 1
!   do j = j_mysta - 1, j_myend + 1
!   do i = i_mysta - 1, i_myend + 1        
!      write (unit = 11, fmt = '(3(i3.3,1x),6(e15.7,1x))') &
!           i,j,k,psa(i,j,k)
!   end do
!   end do
!   end do
!   close (unit = 11)
!   end if
#endif

     ! calculate viscouse terms
     !
     call rhs_sa_viscous ()

#ifdef DEBUG
     call cksum_3d_par ('sa.visc', visc)
#endif

     ! calculate convective flux 
     !
     call rhs_sa_flux ()

#ifdef DEBUG
     call cksum_3d_par ('sa.flux.rh', rh)
#endif

     ! add convective, viscous, production, and unsteady terms to rhs
     !
     call rhs_sa_unst_visc ()

#ifdef DEBUG
     call cksum_3d_par ('sa.unst.rh', rh)
#endif

     ! implicit treatment of positive contribution of source term
     !
     call rhs_sa_implicit_source ()

#ifdef DEBUG
  !call cksum_3d_par ('sa.imp.dk', dk)
  call cksum_3d_par ('sa.imp.rh', rh)  
#endif


     if (myback  == mpi_proc_null) rh(i_mysta-1,:,:) = zero
     if (myfront == mpi_proc_null) rh(i_myend+1,:,:) = zero
     if (myleft  == mpi_proc_null) rh(:,j_mysta-1,:) = zero
     if (myright == mpi_proc_null) rh(:,j_myend+1,:) = zero
     if (mydown  == mpi_proc_null) rh(:,:,k_mysta-1) = zero
     if (myup    == mpi_proc_null) rh(:,:,k_myend+1) = zero

     if (nblk /= 0) then
     do nb = 1, nblk
        rh(li_blk_ia(1,nb):li_blk_ib(1,nb), &
           li_blk_ja(1,nb):li_blk_jb(1,nb), &
           li_blk_ka(1,nb):li_blk_kb(1,nb)) = zero
     end do
     end if

     ! adi_solution
     !
     call rhs_sa_adi_solver ()

#ifdef DEBUG
     call cksum_3d_par ('sa.adi.ijk.rh', rh)
#endif

!      do k = 2, km - 1
!      do j = 2, jm - 1
!      do i = 2, im - 1

     ! update
     !
     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
     do i = i_mysta, i_myend

        qe5 = qold(i,j,k) + rh(i,j,k)

        dum = one + onept5 * e_source * dtev(i,j,k) / delti

        q(5,i,j,k) = q(5,i,j,k) + (-q(5,i,j,k) + qe5) / dum

     end do
     end do
     end do

#ifdef DEBUG
     call cksum_4d_par ('sa.new.q', q)
#endif

     ! boundary conditions
     ! 
     call rhs_sa_bcond ()

     ! update eddy viscosity on ghost cells
     !
     call rhs_exchng3_3d (xnut)

#ifdef DEBUG
     call cksum_3d_par ('sa.xnut', xnut)
     call cksum_4d_par ('sa.bcond', q)
#endif

     ! save solution for iterations
     !
     qold = q(5,:,:,:)

  end do

  call rhs_exchng3_4d (q(5:5,:,:,:))

  ! deallocate variables
  deallocate (ucn_j, qold, visc, dnu, psa, dk, wo, fw, fv)

contains

  include 'rhs_contra_j.F90'
  include 'rhs_sa_source.F90'
  include 'rhs_sa_prod.F90'
  include 'rhs_sa_grdnu.F90'
  include 'rhs_sa_viscous.F90'
  include 'rhs_sa_flux.F90'
  include 'rhs_sa_unst_visc.F90'
  include 'rhs_sa_implicit_source.F90'
  include 'rhs_sa_adi_solver.F90'
  include 'rhs_sa_bcond.F90'
  include 'rhs_exchng3_4d.F90'
  include 'rhs_exchng3_3d.F90'

end subroutine des_eddy

!include 'rhs_kw_adi_solver.F90'
