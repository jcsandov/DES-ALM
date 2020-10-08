
subroutine solver_daf (me, decide_grid_level, decide_recalc_rh, & ! 1
     decide_calc_pk, il, iu, jl, ju, kl, ku, & ! 2
     igp, jgp, kgp, &           ! 3
     dc, de, dz, &              ! 4
     q, &                       ! 5
     qn, &                      ! 6
     qnm1, &                    ! 7
     csi, &                     ! 8
     eta, &                     ! 9
     zet, &                     ! 10
     aj, &                      ! 11
     xnut, &                    ! 12
     pk, &                      ! 13
     rh, &                      ! 14
     gi_ia, gi_ja, gi_ka, gi_jb, gi_kb)   		! 20
 !   sec_blade, blade_dat) !remove ) from disk_dat and ad , 


!     uij )                             ! 31

  ! Solve the system of equations using approximate factorization
  ! implicit solver with local time stepping and residual smoothing

  ! input
  !     decide_grid_level
  !     decide_recalc_rh
  !     q(4,ijk)
  !     qn(4,ijk)
  !     qnm1(4,ijk)
  !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
  !     aj(ijk)
  !     xnut(ijk)
  !     dc, de, dz
  !     dc2, de2, dz2
  !     dcsq, desq, dzsq

  ! output
  !     xnut(ijk)
  !     rh(4,ijk)
  !     
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  use global_param
  use global_app
  use global_mpi
  use checksum
  use global_turb
  

  implicit none

  integer :: n_t
  integer :: me

  integer :: il                 ! i lower bound
  integer :: jl                 ! j lower bound
  integer :: kl                 ! k lower bound
  integer :: iu                 ! i upper bound
  integer :: ju                 ! j upper bound
  integer :: ku                 ! k upper bound

  integer :: i_mysta            ! first interior i node
  integer :: j_mysta            ! first interior j node
  integer :: k_mysta            ! first interior k node
  integer :: i_myend            ! last interior i node
  integer :: j_myend            ! last interior j node
  integer :: k_myend            ! last interior k node

  integer :: igp                ! i-direction ghost points
  integer :: jgp                ! j-direction ghost points
  integer :: kgp                ! k-direction ghost points

  real (kind = rdf) :: dc       ! csi-direction grid spacing
  real (kind = rdf) :: de       ! eta-direction grid spacing
  real (kind = rdf) :: dz       ! zeta-direction grid spacing

  ! arrays for rhs & solver routines
  !
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: q
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(in)    :: qn, qnm1
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(out)   :: rh
  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku), intent(inout) :: pk
  real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku), intent(in) :: csi,eta,zet
  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku)    , intent(in) :: aj

  ! turbulence stuff
  ! 
  real (kind = rdf), dimension(il:iu,jl:ju,kl:ku), intent(in) :: xnut

  ! Reynolds stresses
  !
  !real (kind = rdf), dimension(1:6,il:iu,jl:ju,kl:ku) :: uij
  !real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) :: rstr

  ! dummy variables
  real (kind = rdf), dimension(:,:,:,:), allocatable   :: qold_af, &
                                                             diss, &
                                                             visc, &
                                                            ucn_j, &
                                                              spr, &
                                                               fv

  real (kind = rdf), dimension(:,:,:,:,:), allocatable ::  mai, &
                                                           n1i, &
                                                           n2i, &
                                                            mc

  real (kind = rdf), dimension(:,:,:), allocatable     ::  dtau

  ! control from previous program
  ! 
  integer :: decide_grid_level ! actually just contains current grid level
  integer :: decide_calc_pk    ! calculate fine-grid forcing function
  integer :: decide_recalc_rh

  integer :: i, j, k, idend, m
  integer :: istage

  ! start using this again;
  ! decide prefix is just confusing
  ! 
  integer :: n
  
  integer :: inout, ibs, ibe, ibse

  ! Needed to call the turbines forces !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
   integer ::  gi_ia, gi_ja, gi_ka, gi_jb, gi_kb
   !real (kind = rdf), dimension(gi_ja:gi_jb,gi_ka:gi_kb,1:4), intent(in)    :: disk_dat

!   integer :: ntur, n_t, sec_blade     
!   integer ::  gi_ia, gi_ja, gi_ka, gi_jb, gi_kb   
!   integer dimension(ntur) ::  itur1, itur2 
!   real (kind = rdf), dimension(gi_ja:gi_jb,gi_ka:gi_kb,1:3,ntur), intent(in)    :: disk_dat
!   real (kind = rdf), dimension(1:5,1:sec_blade), intent(in)    :: disk_dat
         
  n = decide_grid_level

  ! allocate variables: dummy, flux variables
  allocate (qold_af(1:4,il:iu,jl:ju,kl:ku), &
               diss(1:4,il:iu,jl:ju,kl:ku), &
               visc(1:3,il:iu,jl:ju,kl:ku), &
              ucn_j(1:3,il:iu,jl:ju,kl:ku), &
                 fv(1:3,il:iu,jl:ju,kl:ku), &
                   dtau(il:iu,jl:ju,kl:ku) )
  qold_af=zero; visc=zero; diss=zero; ucn_j=zero; dtau=zero
  rh=zero

  ! daf, arrays for model matrices
  allocate ( mai(1:4,1:4,il:iu,jl:ju,kl:ku), &
             n1i(1:4,1:4,il:iu,jl:ju,kl:ku), &
             n2i(1:4,1:4,il:iu,jl:ju,kl:ku), &
              mc(1:4,1:4,il:iu,jl:ju,kl:ku), &
                 spr(1:3,il:iu,jl:ju,kl:ku) )
  mai=zero; n1i=zero; n2i=zero; mc =zero; spr=zero

#ifdef DEBUG
  call init_cksum (igp, jgp, kgp)
#endif



  ! set first and last ** interior ** grid nodes for this process
  ! 
  ! *** note ***
  !
  ! using the 'my' notation because this varies depending on
  ! whether or not the process owns any nodes on the boundary
  ! of the computational domain
  !
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! processes on the domain boundaries
  ! 
  if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

  ! save solution for iterations; include ghost points because we
  ! don't know if the boundary is included (and needed or not)
  ! because of the characteristics-based boundary conditions or the
  ! fact that we might be on an interior process
  ! 
  qold_af = q

  ! calculate contravariant velocities (that is, U/J)
  ! 
  call rhs_contra_j ()

#ifdef DEBUG
  call cksum_4d_par ('daf.ucn', ucn_j)
#endif

  ! calculate local time step
  !
  call rhs_daf_dtau
#ifdef DEBUG
  call cksum_3d_par ('daf.dtau', dtau)
#endif

  ! calculate expensive terms only during first stage
  ! viscous terms
  !
  call rhs_viscous ()

#ifdef DEBUG
  call cksum_4d_par ('daf.visc', visc)
#endif


  ! dissipation terms
  !
!  if (quick) then
     call rhs_diss_p ()

#ifdef DEBUG
     call cksum_4d_par ('daf.diss', diss)
#endif

     call rhs_flux_sans_convec ()

#ifdef DEBUG
     call cksum_4d_par ('daf.flux', rh)
#endif

     if ( n == 1 ) then
        call rhs_convec_quick (2)
     else
        call rhs_convec_quick (1)
     end if

#ifdef DEBUG
     call cksum_4d_par ('daf.quick', rh)
#endif

!  else
!     call rhs_diss_matrix ()
!     call rhs_flux ()
!#ifdef DEBUG
!     call cksum_4d_par ('daf.matrix', rh)
!#endif
!
!  end if


  ! compute non-linear closure terms
  ! 
!   if ( nlinc ) then
!      if ( craft .and. n == 1 ) then
!         call rhs_ns_nlin_craft ()
        
! #ifdef DEBUG
!         call cksum_4d_par ('daf.uij', uij)
! #endif

!      end if
!      call rhs_ns_rstress ()

! #ifdef DEBUG
!      call cksum_4d_par ('daf.rstr', rstr)
! #endif
  
!   end if
  
  ! add viscous, dissipation, and unsteady terms to rhs
  ! and reynolds stresses on interior nodes only
  !
  !call rhs_ns_unst_visc_diss ()
  call rhs_unst_visc_diss ()

#ifdef DEBUG
     call cksum_4d_par ('daf.unst', rh)
#endif

  ! add forcing term on coarse grids
  !
  if ( n > 1 ) then
     if (decide_calc_pk == 1) then
        call rhs_pk (1)
     else
        call rhs_pk (0)
     end if
  end if

#ifdef DEBUG
     call cksum_4d_par ('daf.pk', rh)
#endif

  ! BOUNDARIES
  if (myback  == mpi_proc_null) rh(:,i_mysta-1,:,:) = zero
  if (myfront == mpi_proc_null) rh(:,i_myend+1,:,:) = zero
  if (myleft  == mpi_proc_null) rh(:,:,j_mysta-1,:) = zero
  if (myright == mpi_proc_null) rh(:,:,j_myend+1,:) = zero
  if (mydown  == mpi_proc_null) rh(:,:,:,k_mysta-1) = zero
  if (myup    == mpi_proc_null) rh(:,:,:,k_myend+1) = zero

  ! BLANKING
  if (nblk /= 0) then
  do nb = 1, nblk
     rh(:,li_blk_ia(n,nb):li_blk_ib(n,nb), &
          li_blk_ja(n,nb):li_blk_jb(n,nb), &
          li_blk_ka(n,nb):li_blk_kb(n,nb)) = zero
  end do
  end if

  !     idend = i_myend + 1
  ! right-hand-side on outlet boundary
  !
  if ( n == 1 ) then
  do i = 1, 2
     ibs = 0
     ibe = 0
     ibse= 0
     inout=0
     if (i == 1 .and. myback  == mpi_proc_null .and. btype(i,myzone) == 5) then
        ibs =i_mysta-1
        ibe =i_mysta
        ibse=ibs
        inout=i
     end if
     if (i == 2 .and. myfront == mpi_proc_null .and. btype(i,myzone) == 5) then
        ibs =i_myend
        ibe =i_myend+1
        ibse=ibe
        inout=i
     end if
     if (inout /= 0) then
        call mg_brhs  (jl, ju, kl, ku, & !  1
                   dc, de, dz, igp, jgp, kgp, & !  2
                   eta(1:3,ibse,jl:ju,kl:ku), & !  4
                   zet(1:3,ibse,jl:ju,kl:ku), & !  5
                 ucn_j(1:3,ibse,jl:ju,kl:ku), & !  6
                        aj(ibse,jl:ju,kl:ku), & !  7
                      dtau(ibse,jl:ju,kl:ku), & !  8
                      xnut(ibse,jl:ju,kl:ku), & !  9
                     q(1:4,ibse,jl:ju,kl:ku), & ! 10
                    qn(1:4,ibse,jl:ju,kl:ku), & ! 11
                  qnm1(1:4,ibse,jl:ju,kl:ku), & ! 12
                    rh(1:4,ibse,jl:ju,kl:ku)  ) ! 13

        call nonreflect_ibc (inout, jl, ju, kl, ku,  &
                                     igp, jgp, kgp, &
                            dtau(ibse,jl:ju,kl:ku), &
                         csi(1:3,ibse,jl:ju,kl:ku), &
                         ucn_j(1,ibse,jl:ju,kl:ku), &
                              aj(ibse,jl:ju,kl:ku), &
                        q(1:4,ibs:ibe,jl:ju,kl:ku), &
                          rh(1:4,ibse,jl:ju,kl:ku) )
     end if

  end do
  end if

#ifdef DEBUG
     call cksum_4d_par ('daf.exit', rh)
#endif

  ! include boundary grid plane for characteristics-based boundary conditions
  !
  idend = i_myend
  if ( myfront == mpi_proc_null .and. n == 1 .and. btype(2,myzone) == 5) &
       idend = i_myend + 1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Turbine force calling 
if (turbines) then
  call rhs_lines()
  call rhs_exchng3_4d(rh(2:4,:,:,:))
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! calculate model matrices
  !


  call rhs_modal_matrices

#ifdef DEBUG
     call cksum_5d_par ('daf.mai', mai)
     call cksum_5d_par ('daf.n1i', n1i)
     call cksum_5d_par ('daf.n2i', n2i)
     call cksum_5d_par ('daf.mc', mc)
     call cksum_4d_par ('daf.spr', spr)
#endif

  ! implicit diagonal solver
  !
  call rhs_diag_solver ()

#ifdef DEBUG
     call cksum_4d_par ('daf._ijk', rh)
#endif

  ! update interior nodes
  !
  do k=k_mysta,k_myend
  do j=j_mysta,j_myend
  do i=i_mysta,idend ! i_myend
     q(1,i,j,k)=qold_af(1,i,j,k)+rh(1,i,j,k)
     q(2,i,j,k)=qold_af(2,i,j,k)+rh(2,i,j,k)
     q(3,i,j,k)=qold_af(3,i,j,k)+rh(3,i,j,k)
     q(4,i,j,k)=qold_af(4,i,j,k)+rh(4,i,j,k)
  end do
  end do
  end do


 if (n == 1) call bcond_fm (il, iu, jl, ju, kl, ku, igp, jgp, kgp, q, gi_ia, gi_ja, gi_ka)


#ifdef DEBUG
     call cksum_4d_par ('daf.new q', q)
#endif



 ! recompute rh for forcing computation in n level routines
 !
 if (decide_recalc_rh == 1) then

     ! update q on ghost points
     ! 
     call rhs_exchng3_4d (q)      ! update ghost points: q
     call rhs_contra_j
     call rhs_viscous

!     if (quick) then
        call rhs_flux_sans_convec ()
        if ( n > 1 ) then
           call rhs_convec_quick (1)
        else
           call rhs_convec_quick (2)
        end if
!     else
!        call rhs_flux ()
!     end if

     ! compute non-linear closure terms
     ! 
!      if ( nlinc ) then
!         if ( craft .and. n == 1 ) call rhs_ns_nlin_craft ()
!         call rhs_ns_rstress ()
!      end if

     !call rhs_ns_unst_visc_diss ()
     call rhs_unst_visc_diss ()


     if ( n > 1) call rhs_pk (0)

     if (myback  == mpi_proc_null) rh(:,i_mysta-1,:,:) = zero
     if (myfront == mpi_proc_null) rh(:,i_myend+1,:,:) = zero
     if (myleft  == mpi_proc_null) rh(:,:,j_mysta-1,:) = zero
     if (myright == mpi_proc_null) rh(:,:,j_myend+1,:) = zero
     if (mydown  == mpi_proc_null) rh(:,:,:,k_mysta-1) = zero
     if (myup    == mpi_proc_null) rh(:,:,:,k_myend+1) = zero

     if (nblk /= 0) then
     do nb = 1, nblk
        rh(:,li_blk_ia(n,nb):li_blk_ib(n,nb), &
             li_blk_ja(n,nb):li_blk_jb(n,nb), &
             li_blk_ka(n,nb):li_blk_kb(n,nb)) = zero
     end do
     end if


#ifdef DEBUG
     call cksum_4d_par ('daf.new rh', rh)
#endif

  end if

  ! deallocate dummy variables
  deallocate (qold_af, diss, visc, ucn_j, dtau, fv)
  deallocate (mai, n1i, n2i, mc, spr)

contains

  include 'rhs_contra_j.F90'
  include 'rhs_daf_dtau.F90'
  include 'rhs_viscous.F90'
  include 'rhs_diss_p.F90'
  include 'rhs_flux_sans_convec.F90'
  include 'rhs_convec_quick.F90'
  include 'rhs_pk.F90'
  include 'rhs_exchng3_4d.F90'
  include 'rhs_modal_matrices.F90'
  include 'rhs_diag_solver.F90'
  include 'rhs_unst_visc_diss.F90'
  include 'rhs_lines.F90' 
!  include 'rhs_flux.F90'
!  include 'rhs_diss_quick.F90'
!  include 'rhs_diss_matrix.F90'
!  include 'rhs_ns_nlin_craft.F90'
!  include 'rhs_ns_rstress.F90'
!  include 'rhs_ns_unst_visc_diss.F90'

end subroutine solver_daf


