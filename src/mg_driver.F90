
subroutine mg_driver
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! Multigrid driver (v-cycle only) using solver_rk (runge-kutta)
  ! as a smoother. The logic is intertwined between these two routines.
  ! And although, it might be nice I am not sure separate the logic
  ! so that the routines can be truely modular.
  !
  ! Improved version using solver_daf (daigonal approx. factorization)
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  use global
  use global_param
  use global_mpi
  use global_osg
  use global_app
  use wf_mpi
  use global_turb
  

  implicit none

  ! integers
  ! 
  integer :: decide_recalc_rh

  ! counters & placeholders
  ! 
  integer :: n
  integer :: l, m
  integer :: i, j, k
  integer :: la
  integer :: lb

  logical :: converged

  ! defaults for v-cycle multigrid
  ! 
  nstop = ng
  itc = 0
  itmax = itm(1)               ! remove connect to input file

  ! store 'n' and 'n-1' levels for new time step
  ! 
  la = le_idx_a(1)
  lb = le_idx_b(1)
  qnm1(:,la:lb) = qn(:,la:lb)
    qn(:,la:lb) =  q(:,la:lb)

  ! initialize q on coarse grids (mp should work)
  !
  do m = 1, me
     call mg_inject (qn(m,:))
     call mg_inject (qnm1(m,:))
  end do

n = 1
la = le_idx_a(1)
lb = le_idx_b(1)

call conservation(le_ia(n),le_ib(n),&
                le_ja(n),le_jb(n),&
                le_ka(n),le_kb(n),&
                igp(n), jgp(n), kgp(n), & 
                xc(la:lb), yc(la:lb), zc(la:lb),      &
                xe(la:lb), ye(la:lb), ze(la:lb),      &
                xz(la:lb), yz(la:lb), zz(la:lb),      &
                csi(1:3,la:lb),                       & ! 8
                eta(1:3,la:lb),                       & ! 9
                zet(1:3,la:lb),                      & ! 10
                aj(la:lb),      & ! 11
                q(1:4,la:lb) )



  ! pseudo-time iteration
  !
  inner_iteration : do it = 1, itmax
     converged = .false.

!      print *, ''
!      print *, 'proc ', myid, 'starting iteration ', it
!      print *, ''

     itc = itc + 1
     !if (myid == root) then
     !  print*,'.................'
     !  print*,'pseudo-time iteration Nº=', itc,'myid',myid
     !  print*,'.................'
     !end if

     ! save fine grid solution for
     ! various residual calculations
     !
     do l = le_idx_a(1), le_idx_b(1)
        qold_mg(:,l) = q(:,l)
     end do

     ! multigrid v-cycle
     do n = ns, nstop

    !    if (myid == root) then
    !      print*,'-> multigrid v-cycle Nº=',n
    !    end if

        la = le_idx_a(n)
        lb = le_idx_b(n)

!        if (myid == root) print '(/a,1x,i2/17("*"))', 'grid level =', n

        if (n  > 1) then
           ! exchange ghost points for q
           !
            call mg_exchng3_2d (n-1, rh)

           ! mg_inject; mg_calc_pk
           ! 
           call mg_nlevel (n)

           ! update ghost points for uij
           ! 
           if (turbulence .and. nlinc) then
              call mg_exchng3_2d(n, uij)
           end if
           
#ifdef DEBUG
           call mg_cksum ('mg.nlevel', n, q)
#endif


           call mg_zero_pk_boundary (n, pk)
!           call mg_cksum ('nl-pk',n, pk(:,:))
        end if

        do itr = 1, iter(n)

            !if (myid == root) then
            !  print*,'---> itr=',itr
            !end if
            
           if (n /= ng) then
              if (itr /= iter(n)) then
                 decide_recalc_rh = 0
              else
                 decide_recalc_rh = 1
              end if
           end if

           ! exchange ghost points for q
           !
           call mg_exchng3_2d (n, q)

           if (turbulence) then
              call mg_exchng3_1d (n, xnut)
           end if

#ifdef DEBUG
           call mg_cksum ('mg.exchg_q', n, q)
#endif

           if ( des .and. n == ns ) then

              call des_eddy (le_ia(n),le_ib(n),  & ! 1
                le_ja(n),le_jb(n),            & ! 1
                le_ka(n),le_kb(n),            & ! 1
                igp(n), jgp(n), kgp(n),       & ! 2
                dc(n), de(n), dz(n),          & ! 3
                dtev(la:lb),           & ! 4
                q(1:me,la:lb),         & ! 5
                qn(me,la:lb),          & ! 6
                qnm1(me,la:lb),        & ! 7
                csi(1:3,la:lb),        & ! 8
                eta(1:3,la:lb),        & ! 9
                zet(1:3,la:lb),        & ! 10
                aj(la:lb),             & ! 11
                wd(la:lb),             & ! 12
                xnut(la:lb) )            ! 13

              call mg_bintp_tur
#ifdef DEBUG
              call mg_cksum ('mg.osg.tur', n, q)
#endif

           end if  ! end of des_eddy

!print *, 'vmpi data dif numero nodos', gi_ia(n),gi_ib(n),gi_ja(n), gi_jb(n),gi_ka(n),gi_kb(n)

!if(myid==root) print *, 'antes de solver_daf en mg_driver'
           call solver_daf (me, n, decide_recalc_rh, itr, & ! 1
                le_ia(n),le_ib(n),le_ja(n),le_jb(n),le_ka(n),le_kb(n), & ! 2
                igp(n), jgp(n), kgp(n),        & ! 3
                dc(n), de(n), dz(n),           & ! 4
                q(1:4,la:lb),           & ! 5
                qn(1:4,la:lb),          & ! 6
                qnm1(1:4,la:lb),        & ! 7
                csi(1:3,la:lb),         & ! 8
                eta(1:3,la:lb),         & ! 9
                zet(1:3,la:lb),         & ! 10
                aj(la:lb),              & ! 11
                xnut(la:lb),            & ! 12
                pk(1:4,la:lb),          & ! 13
                rh(1:4,la:lb),           & ! 14
                gi_ia(1),gi_ja(1),gi_ka(1),gi_jb(1),gi_kb(1))
  ! Turbine case 2, with a turbine info file
	  	!disk_dat(gi_ja(1):gi_jb(1),gi_ka(1):gi_kb(1),1:3)
		!sec_blade, blade_dat(1:5,1:sec_blade)
!if(myid==root) print *, 'DESPUES de solver_daf en mg_driver'


#ifdef DEBUG
        call mg_cksum ('mg.daf-q', n, q)
        call mg_cksum ('mg.daf-rh', n, rh)
#endif

           !==============================================================
           ! FOR DARGAHI's CASE ONLY
           ! WALL FUNCTIONS
           if (n == ns .and. myzone == 1 .and. myleft  == mpi_proc_null) &
             call wf_eta (1)
           if (n == ns .and. myzone == 1 .and. myright == mpi_proc_null) &
             call wf_eta (2)

           !if (n == ns .and. myzone == 6 .and. myleft  == mpi_proc_null) &
           !   call wf_eta (1)
           !if (n == ns .and. myzone == 6 .and. myright == mpi_proc_null) &
           !   call wf_eta (2)

           !if (n == ns .and. myzone == 7 .and. myleft  == mpi_proc_null) &
           !   call wf_eta (1)
           !if (n == ns .and. myzone == 7 .and. myright == mpi_proc_null) &
           !   call wf_eta (2)

           !==============================================================


           ! boundary interpolation in overset grid
           !
           if (n == ns .and. nzone > 1) then

              call mg_bintp_mom
#ifdef DEBUG
              call mg_cksum ('mg.osg', n, q)
#endif
           end if

        end do ! itr

        call mg_zero_pk_boundary (n, pk)

     end do ! end n

     ! prolong from grid n+1 to grid n
     !
     do n = nstop-1, ns, -1

        ! prolong from coarse grid to fine grid
        ! 
        call mg_prolong (n)
        !
#ifdef DEBUG        
        call mg_cksum ('mg.prolong', n, q)
#endif

     end do

!      call mpi_barrier (mpi_comm_world, ierr)
!      call mpi_abort (mpi_comm_world, errcode, ierr)
!      stop

     ! finished V-cycle; apply boundary conditions
     ! 
     n  = 1                      
     la = le_idx_a(1)
     lb = le_idx_b(1)

!    if(myid==root) print *, 'antes de bcond_fm en mg_driver'

!     call bcond_fm (le_ia(n),le_ib(n),&
!                    le_ja(n),le_jb(n),&
!                    le_ka(n),le_kb(n),&
!                    igp(n), jgp(n), kgp(n), &
!                    q(1:4,la:lb))

          call bcond_fm (le_ia(n),le_ib(n),&
                    le_ja(n),le_jb(n),&
                    le_ka(n),le_kb(n),&
                    igp(n), jgp(n), kgp(n), &
                    q(1:4,la:lb), &
                    gi_ia(1),gi_ja(1),gi_ka(1))

!    if(myid==root) print *, 'DESPUES de bcond_fm en mg_driver'


#ifdef DEBUG        
     call mg_cksum ('mg.bc', n, q)
#endif

     ! boundary interpolation in overset grid
     ! exchange ghost points for q
     !
     call mg_exchng3_2d (n, q)

     !call mg_bintp_mom

#ifdef DEBUG
        call mg_cksum ('mg.osg', n, q)
#endif

     ! check for convergence during this time step
     ! 
     call mg_ck_conver ( converged )

     if ( converged ) exit inner_iteration

     !if ( it > 4 ) exit inner_iteration

!      call mpi_barrier (mpi_comm_world, ierr)
!      call mpi_abort (mpi_comm_world, errcode, ierr)

  end do inner_iteration
  
contains

  include 'mg_exchng3_2d.F90'
  include 'mg_exchng3_1d.F90'
  include 'mg_cksum.F90'
  include 'mg_prolong.F90'
  include 'mg_residual.F90'
  include 'mg_bintp_tur.F90'
  include 'mg_bintp_mom.F90'

  subroutine mg_ck_conver (converged)

    implicit none 

    !real (kind = rdf) :: erp
    real (kind = rdf) :: eru
    real (kind = rdf) :: erv
    real (kind = rdf) :: erw
    !real (kind = rdf) :: erx

    real (kind = rdf), save :: ermax
    real (kind = rdf), save :: eomax

    ! residual norm calcuations
    !
    real (kind = rdf), dimension(1:me), save :: init_resd
    real (kind = rdf), dimension(1:me) :: residual

    logical :: converged

    ! if (itc == 1 .or. ((itc / icn) * icn == itc)) then

    ! residual calculations
    !
    call mg_residual (1, residual)

    if (myid == root) then
       open  (unit = 8, file = 'conver', position = 'append', form = 'formatted')
       write (unit = 8, fmt = '(i5.5,1x,6(g13.6,1x))') itc, residual(:)
       close (unit = 8)
       !end if

       ! used mpi_reduce_scatter () in mg_resid; thus
       ! compute logical converged on all processes

       ! save initial residual
       !
       if (itc == 1) then

          init_resd(:) = residual(:)

!           init_resd(1) = erp
!           init_resd(2) = eru
!           init_resd(3) = erv
!           init_resd(4) = erw
!           init_resd(5) = erx
       end if

       ! check converged & go to next time step?
       !
       eomax = max(residual(2),residual(3),residual(4))

       eru = residual(2) - init_resd(2)
       erv = residual(3) - init_resd(3)
       erw = residual(4) - init_resd(4)

       ermax = max(eru, erv, erw)

       ! if (duct) ermax = max(eru, erv)
       !$ermax = erv

       converged = .false.
       if ((itc > it_min) .and. (ermax < er_min)) converged = .true.
       if ((itc > it_min) .and. (eomax < eo_min)) converged = .true.

    end if

    call mpi_bcast (converged, 1, mpi_logical, root, mpi_comm_world, ierr)
    !call mpi_bcast (converged, 1, mpi_logical, root, comm3d, ierr)

  end subroutine mg_ck_conver
  
  subroutine mg_zero_pk_boundary(n, var)

    integer , intent(in) :: n
    real (kind = rdf), dimension (:,:), intent(inout) :: var

    if (myback == mpi_proc_null)  then
       do k = le_ka(n), le_kb(n)
       do j = le_ja(n), le_jb(n)
          l = le_idx(li_ia(n),j,k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myfront == mpi_proc_null) then
       do k = le_ka(n), le_kb(n)
       do j = le_ja(n), le_jb(n)
          l = le_idx(li_ib(n),j,k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myleft == mpi_proc_null)  then
       do k = le_ka(n), le_kb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,li_ja(n),k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myright == mpi_proc_null) then
       do k = le_ka(n), le_kb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,li_jb(n),k,n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (mydown == mpi_proc_null)  then
       do j = le_ja(n), le_jb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,j,li_ka(n),n)
          var(:,l) = zero
       end do
       end do
    end if
    
    if (myup == mpi_proc_null) then
       do j = le_ja(n), le_jb(n)
       do i = le_ia(n), le_ib(n)
          l = le_idx(i,j,li_kb(n),n)
          var(:,l) = zero
       end do
       end do
    end if
    
  end subroutine mg_zero_pk_boundary


end subroutine mg_driver

!
! save for later use
!
!



!!$           if (rkm) then
!!$           call solver_rk (me, n, decide_recalc_rh, itr, & 
!!$                img(n), jmg(n), kmg(n), &
!!$                dc, de, dz,             &
!!$                dtau(ls(n):le(n)),      &
!!$                dtev(ls(n):le(n)),      &
!!$                q(1:me,ls(n):le(n)),     &
!!$                qn(1:me,ls(n):le(n)),    &
!!$                qnm1(1:me,ls(n):le(n)),  &
!!$                csi(1:3,ls(n):le(n)),   &
!!$                eta(1:3,ls(n):le(n)),   &
!!$                zet(1:3,ls(n):le(n)),   &
!!$                aj(ls(n):le(n)),        &
!!$                wd(ls(n):le(n)),        &
!!$                xnut(ls(n):le(n)),      &
!!$                pk(1:4,ls(n):le(n)),    &
!!$                rh(1:4,ls(n):le(n)) )
!!$           else if (daf) then
!            call solver_daf (me, n, decide_recalc_rh, itr, & 
!                 img(n), jmg(n), kmg(n), &
!                 dc, de, dz,             &
!                 dtau(ls(n):le(n)),      &
!                 dtev(ls(n):le(n)),      &
!                 q(1:me,ls(n):le(n)),     &
!                 qn(1:me,ls(n):le(n)),    &
!                 qnm1(1:me,ls(n):le(n)),  &
!                 csi(1:3,ls(n):le(n)),   &
!                 eta(1:3,ls(n):le(n)),   &
!                 zet(1:3,ls(n):le(n)),   &
!                 mai(1:4,1:4,ls(n):le(n)),  &
!                 n1i(1:4,1:4,ls(n):le(n)),  &
!                 n2i(1:4,1:4,ls(n):le(n)),  &
!                 mc(1:4,1:4,ls(n):le(n)),   &
!                 spr(1:3,ls(n):le(n)),   &
!                 aj(ls(n):le(n)),        &
!                 wd(ls(n):le(n)),        &
!                 xnut(ls(n):le(n)),      &
!                 pk(1:4,ls(n):le(n)),    &
!                 rh(1:4,ls(n):le(n)) )
! !!$           end if















