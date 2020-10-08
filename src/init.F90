
subroutine init ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  use global
  use global_param
  use global_app
  use global_mpi
  use global_osg
  use checksum
  use global_turb

  implicit none
  
  ! mpi-osg
  integer, dimension(:,:), allocatable :: id_ijkz

  real (kind = rdf) :: tmp

  integer :: icycle    ! left over from multiple
                       ! types of multigrid; eliminate?

  integer :: i, j, k, l, m, cont
  integer :: n, p
  integer :: nz
  integer :: ib
  integer :: n_t, n_t_aux

  ! mpi-osg
  allocate ( id_ijkz(1:7,0:nproc-1) )

  call init_input_file ()       ! get input info and distribute

! call init_mg_old ()           ! keep until mg_driver, etc fixed


  call init_decomp_3d ()        ! decide grid distribution among proc.

  call init_mg ()               ! allocate local memory for each proc.

  call init_blanking


  call init_grid ()             ! root read grid and distribute to proc.

if (myid==root) print *,'before init_lines ...'

!===================================================================================================
!Allocation of disk variables and initialization of the actuator disks.
! WARNING:
! 	- maximum 5000 points in every x-plane
!	- maximum 20 x-planes inside every disk
!===================================================================================================
do nz=1,nzone
  !allocate(disk_dat(1:jmg(1,nz),1:kmg(1,nz),1:7,1:n_tur))
  allocate(disk_dat(1:jmg(1,nz),1:kmg(1,nz),1:8,1:n_tur)) !JS,2018
  allocate(disk_dat_ker(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz),1:4,1:n_tur)) !JS,2019
  allocate(airfoil_before(1:jmg(1,nz),1:kmg(1,nz),1:n_tur))
  allocate(airfoil_after(1:jmg(1,nz),1:kmg(1,nz),1:n_tur))
  !allocate(disk_node_coordinates(1:5000,1:3,1:n_tur))
  allocate(disk_node_coordinates(1:50000,1:3,1:n_tur)) !Changed to save the position too. JS, 2018
  allocate(disk_node_coordinates_pos(1:50000,1:3,1:n_tur)) !Changed to save the position too. JS, 2018

  allocate(disk_node_coordinates_ker(1:100000,1:3,1:n_tur)) !Created to apply the forces with ker. JS, 2019
  allocate(disk_node_coordinates_pos_ker(1:100000,1:3,1:n_tur)) !Created to apply the forces with ker. JS, 2019

  allocate(C_T(1:n_tur))
  allocate(C_P(1:n_tur))
  allocate(kernel(1:20,1:n_tur))
  !================================================
  ! ALM variables allocation (J. Sandoval, 2018)  !
  ! WARNING:                                      !
  !   * maximum of 20000 dtheta discretizations   !
  !================================================
  allocate(dtheta(1:n_tur))
  allocate(nangles(1:n_tur))
  allocate(dtheta_interp(1:n_tur))
  allocate(blade_elemts(1:n_tur,1:maxval(nblades),1:maxval(nelems),1:3))
  allocate(blade_angles(1:n_tur,1:maxval(nblades)))
  allocate(node_indexes(1:n_tur,1:maxval(nelems),1:900000,1:3)) !*
  allocate(node_indexes_pos(1:n_tur,1:maxval(nelems),1:900000,1:27)) ! Now we will also write the positions. JS, 2018
!  allocate(kerALM(1:n_tur,1:maxval(nelems),1:20000,1:5000)) !*
!  Now this variable also incorporate the position of adjacent nodes to interpolate the velocity
!  in rhs_lines (Jorge Sandoval, 2019. ALM_v2_trilinear_interpolation)
  allocate(forces_lines(1:n_tur,1:maxval(nblades),1:maxval(nelems),1:6)) ! JS, 2018
end do

if (turbines) then
  call init_lines ()
end if

if (myid==root) print *, 'After init_lines ...'

  call mg_inject (x)  ! inject x to coarser grids
  call mg_inject (y)  ! inject y
  call mg_inject (z)  ! inject z

!print*, 'After mg_inject', myid

  do n = 2, ng                  ! exchange ghost points for grid
                                ! on coarse grid levels
     call exchng3_mg_1d (n, x)  
     call exchng3_mg_1d (n, y)
     call exchng3_mg_1d (n, z)
  end do

!print*, 'After exchng', myid


  call mg_metrics ()            ! calculate metrics and daf init
                                ! on interior nodes
!print*, 'After mg_metrics', myid

  call init_osg_interface       ! interface information in overset grid

!print*, 'After init_osg_interface', myid
  call init_outlet_bc ()


call mpi_barrier (mpi_comm_world, ierr)

  do n = 1, ng
     call exchng3_mg_2d(n, csi)
     call exchng3_mg_2d(n, eta)
     call exchng3_mg_2d(n, zet)

     call exchng3_mg_1d(n,dtev)
     call exchng3_mg_1d(n,  wd)
     call exchng3_mg_1d(n,  aj)
  end do
  
  call init_solu ()             ! root read solu and distribute to proc.

  call init_daf_bc ()           ! initialize bcs for implicit method

#ifdef DEBUG
  call ck_variables ()            ! temp routine
#endif

contains

  subroutine init_input_file ()
    
    integer :: myunit

    integer :: nproc_total

    integer :: blka, blkb, nzb

    character (len=256) :: filename

    ! allocate variables for boundary type and blanking area
    allocate (                  ep(2,3,nzone), & 
                                  cfl1(nzone), &
                                  vnn1(nzone), &
                                  cfl2(nzone), &
                                  vnn2(nzone), &
                               btype(6,nzone), &
                                bdir(6,nzone), &
                            nzblanking(nzone), &
                     blktype(6,nba_max,nzone), &
                  blanking(3,2*nba_max,nzone) )

    ! overset grid
    allocate ( nproc_nz(nzone), &
                 dims(3,nzone) )

    ! if (myid == root) then  ! ELIMINATED IN CORFU.ECE.GATECH.EDU

       ! parameter for unsteady term
       e_source = zero
       if ( unsteady ) e_source = one

       ! read input data
       open (1, file = 'ind3dmg.dat')

       read (1,*) ns, icycle
       read (1,*) (img(ns,nz), jmg(ns,nz), kmg(ns,nz), nz = 1, nzone)
       read (1,*) (icrs(nz), jcrs(nz), kcrs(nz), nz = 1, nzone)
       read (1,*) nt2, delti, it_min, er_min, eo_min
       read (1,*) (iter(n), n = 1, ng)
       read (1,*) (itm(n), n = 1, ng)
       read (1,*) ren, beta, icn, icnw
       read (1,*) (cfl1(nz), vnn1(nz), nz = 1, nzone)
       read (1,*) irk, (alfa(i), i = 1, irk)
       read (1,*) (eps(i), i = 1, 4)
       read (1,*) itk, cdes, ckin, cein
       read (1,*) (cfl2(nz), vnn2(nz), nz = 1, nzone)
       read (1,*) (((ep(j,k,nz), k = 1, 3), j = 1, 2), nz = 1, nzone)
       read (1,*) checkpoint
       read (1,*) pdiss_coef
       read (1,*) monitor_num_points
       read (1,*) ((monitor_point_ijk(j,k),k=1,4),j=1,monitor_num_points)
       read (1,*) ((btype(i,nz), i = 1, 6), nz = 1, nzone)
       read (1,*) (( bdir(i,nz), i = 1, 6), nz = 1, nzone)
       read (1,*) ifix, jfix, kfix, nzfix
       read (1,*) length_scale
       read (1,*) vel_infty!, num_blades
!       Turbine parameters
       read (1,*) turbines
       read (1,*) node_indexes_load
       read (1,*) n_tur
       allocate (itur1(n_tur),itur2(n_tur), jtur1(n_tur), jtur2(n_tur), ktur1(n_tur), ktur2(n_tur), tur_rad(n_tur))
       allocate (hub_rad(n_tur))
       allocate (TSR(n_tur))
       allocate (nblades(n_tur))
       allocate (nelems(n_tur))
       do n_t=1,n_tur
       read (1,*) itur1(n_t), itur2(n_t), jtur1(n_t), jtur2(n_t), ktur1(n_t), ktur2(n_t), tur_rad(n_t)
       read (1,*) hub_rad(n_t),TSR(n_t)
	     read (1,*) nblades(n_t),nelems(n_t) !modified by J. Sandoval, 18/05/2018
       print *, ''
       end do

      ! ! Memory allocation for blade_angles. Modified by J. Sandoval 22/05/2018
      ! allocate (blade_angles(n_tur,maxval(nblades)))

       ! read blanking data
       do nz = 1, nzone
          read (1,*) nzblanking(nz)
          if (nzblanking(nz) /= 0) then
             do nzb = 1, nzblanking(nz)
                blka = 1 + 2 * (nzb - 1)
                blkb = 2 + 2 * (nzb - 1)
                read (1,*) (blktype(n,nzb,nz), n = 1,6)
                read (1,*) ((blanking(i,j,nz), j = blka, blkb), i = 1, 3)
             end do
          else
             blanking(:,:,nz) = 0
          end if
       end do
	   ! read end of bed-load layer


       close (1)

       open (2, file = 'nproc.dat')
       read (2,*) (nproc_nz(nz), nz = 1, nzone)
       read (2,*) ((dims(i,nz), i = 1, 3), nz = 1, nzone)
       close (2)

       nproc_total = 0
       do nz = 1, nzone
          nproc_total = nproc_total + nproc_nz(nz)
       end do
       if (nproc /= nproc_total) &
          print *, 'nproc not match with nproc_total ',nproc ,':',nproc_total

       if ( des .or. s_a ) cvin = pt1 / ren

       ! Calculate img(n), jmg(n), kmg(n) instead of relying on input file
       ! 
       ! icrs, jcrs, = 0 ---> semi-coarsenzng
       ! kcrs        = 1 ---> full-coarsenzng
       !
       do nz=1,nzone
       do n =2,ng
          img(n,nz)=(img(n-1,nz)+icrs(nz))/2**icrs(nz)
          jmg(n,nz)=(jmg(n-1,nz)+jcrs(nz))/2**jcrs(nz)
          kmg(n,nz)=(kmg(n-1,nz)+kcrs(nz))/2**kcrs(nz)
       end do
       end do

       ! read file state
       !
       allocate(angle(n_tur))
       open (unit = 11, file = 'filestat')
       read (unit = 11, fmt = *) time0, nt1
       do n_t=1,n_tur
       read (unit = 11, fmt = *) angle(n_t)
       if (angle(n_t)>=2*pi) then
         angle(n_t)=angle(n_t)-2*pi
        end if
       end do
       close(unit = 11)

    ! end if 		  ! if (myid == root)

    ! broadcast variables from input file to all processes
    ! 
    ! call mpi_bcast (ns, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (icycle, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (me, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (img, ng*nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (jmg, ng*nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (kmg, ng*nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (nt2, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (delti, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (e_source, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (it_min, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (er_min, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (eo_min, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (iter, ng, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (itm, ng, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (ren, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (cfl1, nzone, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (vnn1, nzone, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (beta, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (icn, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (icnw, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (irk, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (alfa, irk, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (eps, 4, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (itk, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (cfl2, nzone, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (vnn2, nzone, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (cdes, 1, mpi_real, root, mpi_comm_world, ierr)    
    ! call mpi_bcast (ckin, 1, mpi_real, root, mpi_comm_world, ierr)    
    ! call mpi_bcast (cein, 1, mpi_real, root, mpi_comm_world, ierr)    
    ! call mpi_bcast (cvin, 1, mpi_real, root, mpi_comm_world, ierr)    
    ! call mpi_bcast (ep, 6*nzone, mpi_real, root, mpi_comm_world, ierr)
    !
    ! partc_skip not needed; should eliminate
    !
    ! call mpi_bcast (pdiss_coef, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (icrs, nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (jcrs, nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (kcrs, nzone, mpi_integer, root, mpi_comm_world, ierr)
    !
    ! monitor num_points needed only by root
    ! monitor_point_ijk() needed only by root
    ! 
    ! call mpi_bcast (btype, 6*nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (bdir , 6*nzone, mpi_integer, root, mpi_comm_world, ierr)

    ! call mpi_bcast (ifix, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (jfix, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (kfix, 1, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (nzfix, 1, mpi_integer, root, mpi_comm_world, ierr)

    ! call mpi_bcast (checkpoint, 1, mpi_integer, root, mpi_comm_world, ierr)

    ! call mpi_bcast (time0, 1, mpi_real, root, mpi_comm_world, ierr)
    ! call mpi_bcast (nt1, 1, mpi_integer, root, mpi_comm_world, ierr)

    ! nprocs in each direction
    ! call mpi_bcast (dims, 3*nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (nproc_nz, nzone, mpi_integer, root, mpi_comm_world, ierr)

    ! overset grid
    ! call mpi_bcast (nzblanking, nzone, mpi_integer, root, mpi_comm_world, ierr)
    ! call mpi_bcast (blanking, 6*nba_max* nzone, mpi_integer, root, &
    !                 mpi_comm_world, ierr)
    ! call mpi_bcast (blktype, 6*nba_max* nzone, mpi_integer, root, &
    !                mpi_comm_world, ierr)

    ! monitoring
    ! call mpi_bcast (monitor_num_points, 1, mpi_integer, root, &
    !                 mpi_comm_world, ierr)
    ! call mpi_bcast (monitor_point_ijk,4, mpi_integer, root, &
    !                 mpi_comm_world, ierr)

    ! call mpi_barrier(mpi_comm_world, ierr)
	
	! in case of simulation with sediment transport:

	! if (concentration) call mpi_bcast (bedlayer, 1, mpi_integer, root, mpi_comm_world, ierr)

#ifdef DEBUG

    ! check via output file
    !
    myunit = myid + 50

    write (filename, fmt = '(a,i2.2)') 'ck_input', myunit

    open (unit = myunit, file = trim(filename), form = 'formatted')

    write (unit = myunit, fmt = '(a15,3x,g)') 'ns', ns
    write (unit = myunit, fmt = '(a15,3x,g)') 'icycle', icycle
    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'img', img(:,nz)
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'jmg', jmg(:,nz)
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'kmg', kmg(:,nz)
    end do
    write (unit = myunit, fmt = '(a15,3x,g)') 'nt2', nt2
    write (unit = myunit, fmt = '(a15,3x,g)') 'delti', delti
    write (unit = myunit, fmt = '(a15,3x,g)') 'e_source', e_source
    write (unit = myunit, fmt = '(a15,3x,g)') 'it_min', it_min
    write (unit = myunit, fmt = '(a15,3x,g)') 'er_min', er_min
    write (unit = myunit, fmt = '(a15,3x,g)') 'eo_min', eo_min
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'iter', iter
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'itm', itm
    write (unit = myunit, fmt = '(a15,3x,g)') 'ren', ren
    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,g)') 'cfl1', cfl1(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'vnn1', vnn1(nz)
    end do
    write (unit = myunit, fmt = '(a15,3x,g)') 'beta', beta
    write (unit = myunit, fmt = '(a15,3x,g)') 'icn', icn
    write (unit = myunit, fmt = '(a15,3x,g)') 'icnw', icnw
    write (unit = myunit, fmt = '(a15,3x,g)') 'irk', irk
    write (unit = myunit, fmt = '(a15,3x,4(g,1x))') 'alfa', alfa
    write (unit = myunit, fmt = '(a15,3x,4(g,1x))') 'eps', eps
    write (unit = myunit, fmt = '(a15,3x,g)') 'itk', itk
    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,g)') 'cfl2', cfl2(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'vnn2', vnn2(nz)
    end do
    write (unit = myunit, fmt = '(a15,3x,g)') 'cdes', cdes
    write (unit = myunit, fmt = '(a15,3x,g)') 'ckin', ckin
    write (unit = myunit, fmt = '(a15,3x,g)') 'cein', cein
    write (unit = myunit, fmt = '(a15,3x,g)') 'cvin', cvin
    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'ep', ep(:,:,nz)
    end do
    write (unit = myunit, fmt = '(a15,3x,g)') 'pdiss_coef', pdiss_coef
    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,g)') 'icrs', icrs(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'jcrs', jcrs(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'kcrs', kcrs(nz)
    end do
    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'btype', btype(:,nz)
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'bdir', bdir(:,nz)
    end do
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'ifix', ifix
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'jfix', jfix
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'kfix', kfix

    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'dims', dims(:,nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'nproc_nz', nproc_nz(nz)
    end do

    do nz = 1, nzone
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'blanking', blanking(:,:,nz)
    end do

    close (unit = myunit)

#endif

  end subroutine init_input_file
  
  ! --

  subroutine init_daf_bc ()

    implicit none

    integer :: np

    integer :: myunit
    character (len=255) :: filename

    ! find process that own pfix grid point
    ! 
    ! limited to a 1d (k) decomposition
    ! 
    if (gi_ka(1) <= kfix .and. kfix <= gi_kb(1) .and. &
        gi_ja(1) <= jfix .and. jfix <= gi_jb(1) .and. &
        gi_ia(1) <= ifix .and. ifix <= gi_ib(1) .and. &
        myzone == nzfix ) then
        pfix_proc = myid
        local_kfix = kfix - gi_ka(1) + 1
        local_jfix = jfix - gi_ja(1) + 1
        local_ifix = ifix - gi_ia(1) + 1
    end if

#ifdef DEBUG
    myunit = myid + 60
    write (filename, fmt = '(a,i2.2)') 'pfix_input', myunit
    open (unit = myunit, file = trim(filename), form = 'formatted')
    write (unit = myunit, fmt = '(a15,3x,g)') 'pfix_proc', pfix_proc
    write (unit = myunit, fmt = '(a15,3x,g)') 'local_ifix', local_ifix
    write (unit = myunit, fmt = '(a15,3x,g)') 'local_jfix', local_jfix
    write (unit = myunit, fmt = '(a15,3x,g)') 'local_kfix', local_kfix
    close(unit = myunit)
#endif
    
    ! set boundary treatment coefficients, sa and sb
    ! btyep	= 0, interface
    !		= 1, solid wall
    !		= 2, symmetric plane
    !		= 3, free stream
    !		= 4, inflow
    !		= 5, outflow
    !		= 6, periodic

    allocate (sa(5,6), &
              sb(5,6))

    ! p    -- > first-order extrapolate
    ! wall -- > no-slip; no-flux

    ! default ---> first-order extrapolate
    sa = four_third
    sb = -one_third

    do ib = 1, 6
       if (btype(ib,myzone) == 1) then
          sa(2:4,ib) = zero
          sb(2:4,ib) = zero
       end if
       if (btype(ib,myzone) == 2) then
          if (bdir(ib,myzone) == 1) then
             sa(2,ib) = zero
             sb(2,ib) = zero
          else if (bdir(ib,myzone) == 2) then
             sa(3,ib) = zero
             sb(3,ib) = zero
          else if (bdir(ib,myzone) == 3) then
             sa(4,ib) = zero
             sb(4,ib) = zero
          end if
       end if
    end do

  end subroutine init_daf_bc

  ! --

  subroutine init_decomp_3d ()

    integer :: ii, jj, kk

    integer :: myunit

    integer :: np,  &
               npa, &
               npb

    integer :: base_grp, grp

    integer :: color, key

    integer :: myid_nz
    integer, dimension(:), allocatable :: myzone_start

    character (len=255) :: filename

    ! allocate place holders
    ! 
    allocate(li_ia(ng), li_ja(ng), li_ka(ng), &
         le_ia(ng), le_ja(ng), le_ka(ng), &
         gi_ia(ng), gi_ja(ng), gi_ka(ng), &
         ge_ia(ng), ge_ja(ng), ge_ka(ng), &
         li_ib(ng), li_jb(ng), li_kb(ng), &
         le_ib(ng), le_jb(ng), le_kb(ng), &
         gi_ib(ng), gi_jb(ng), gi_kb(ng), &
         ge_ib(ng), ge_jb(ng), ge_kb(ng), &
         li_ix(ng), li_jx(ng), li_kx(ng), &
         le_ix(ng), le_jx(ng), le_kx(ng), &
         gi_ix(ng), gi_jx(ng), gi_kx(ng), &
         ge_ix(ng), ge_jx(ng), ge_kx(ng), &
         igp(ng), jgp(ng), kgp(ng), &
         li_idx_a(ng), &
         le_idx_a(ng), &
         gi_idx_a(ng), &
         li_idx_b(ng), &
         le_idx_b(ng), &
         gi_idx_b(ng), &
         li_idx_mx(ng), &
         le_idx_mx(ng), &
         gi_idx_mx(ng) )

    ! Get a new communzcator for a decompositon of the domain.
    ! Let MPI find a "good" decomposition:
    !
    !dim(1) = 0 !need to be read in input date
    !dim(2) = 0 !if dims(i) is zero, use mpi_dims_create for i-direction only
    !dim(3) = 0

    allocate ( coords(3,nzone), &
              periods(3,nzone), &
              reorder(3,nzone) )

    allocate (comm_nz(nzone), &
               comm3d(nzone), &
         myzone_start(nzone) )

    ! need to be moved into init_input routine
    periods = .false.
    reorder = .true.

    ! establish the myzone to which this processor belongs
    !
    do np = 0, nproc - 1
       if (myid == np) then
          npa = 0
          npb = -1
          do nz = 1, nzone
             npb = npb + nproc_nz(nz)
             if ( np >= npa .and. np <= npb ) then
                myzone = nz
             end if
             npa = npb + 1
          end do
       end if
    end do

    myzone_start(1) = 0
    do nz = 2, nzone
       myzone_start(nz) = myzone_start(nz-1) + nproc_nz(nz-1)
    end do

    ! build nz communzcators
    !
    color = myzone
    key   = myid
    call mpi_comm_split (mpi_comm_world, color, key, comm_nz(myzone), ierr)

    ! mpi_dims_create works for dims(i) = 0 only
    !
    call mpi_dims_create (nproc_nz(myzone), 3, dims(1:3,myzone), ierr)
    call mpi_cart_create (comm_nz(myzone), 3, dims(1:3,myzone), &
                          periods(1:3,myzone), .true., comm3d(myzone), ierr)
                          
    ! get my position in this sub-communzcator:
    !
    call mpi_comm_rank (comm3d(myzone), myid_nz, ierr)
    myid = myid_nz + myzone_start(myzone)

    ! my neighbors are now +/- 1 with my rank.
    ! handle the case of the boundaries by using mpi_proc_null:
    !
    ! this routine determines the neighbors in a 3d decomposition of the domain
    ! this assumes that mpi_cart_create has already been called:
    !

    call mpi_cart_shift (comm3d(myzone), 0, 1, myback, myfront, ierr)
    call mpi_cart_shift (comm3d(myzone), 1, 1, myleft, myright, ierr)
    call mpi_cart_shift (comm3d(myzone), 2, 1, mydown, myup   , ierr)

    ! does it need?
    if (myback  /= mpi_proc_null) myback  = myback  + myzone_start(myzone)
    if (myfront /= mpi_proc_null) myfront = myfront + myzone_start(myzone)
    if (myleft  /= mpi_proc_null) myleft  = myleft  + myzone_start(myzone)
    if (myright /= mpi_proc_null) myright = myright + myzone_start(myzone)
    if (mydown  /= mpi_proc_null) mydown  = mydown  + myzone_start(myzone)
    if (myup    /= mpi_proc_null) myup    = myup    + myzone_start(myzone)

    ! decomposition:
    !
    ! 3-d domain decomposition using mpe_decomp1d

    ! retrieves cartesian topology information from a communzcator
    !
    call mpi_cart_get (comm3d(myzone), 3, dims(1:3,myzone), &
                       periods(1:3,myzone), coords(1:3,myzone), ierr)

    ! produce a decomposition of a 1-d array when given a number of processors
    if ( dims(1,myzone) /= 1) &
       call mpe_decomp1d(img(1,myzone), dims(1,myzone), coords(1,myzone), &
                         gi_ia(1), gi_ib(1))
    if ( dims(2,myzone) /= 1) &
       call mpe_decomp1d(jmg(1,myzone), dims(2,myzone), coords(2,myzone), &
                         gi_ja(1), gi_jb(1))
    if ( dims(3,myzone) /= 1) &
       call mpe_decomp1d(kmg(1,myzone), dims(3,myzone), coords(3,myzone), &
                         gi_ka(1), gi_kb(1))

    ! print dims & coords data for checking
    !
    myunit = myid! + 10 !70
    write (filename, fmt = '(a,i3.3)') 'decofiles/dims_coords', myunit
    open (unit = myunit, file = trim(filename), form = 'formatted')
    write(unit = myunit, fmt = *) dims(1:3,myzone)
    write(unit = myunit, fmt = *) coords(1:3,myzone)
    close(unit = myunit)

    ! sync
    call mpi_barrier (mpi_comm_world, ierr)

    ! on the coarse grids
    if (dims(1,myzone) == 1) then

       ! number ghost points
       igp(1:ng) = 0

       ! i direction fixed, global coordinates
       gi_ia(1:ng) = 1
       gi_ib(1:ng) = img(1:ng,myzone)

    else 

       ! number ghost points
       igp(1:ng) = 2

       ! set up i-direction decomposition on the coarse grids
       !
       ! if process owns grid node on finest grid then it owns
       ! grid node on all coarser grids
       !
       do n = 2, ng
          if (icrs(myzone) == 0) then
             gi_ia(n) = gi_ia(1)
             gi_ib(n) = gi_ib(1)
          else
             if (mod(gi_ia(n-1),2) == 0) then
                gi_ia(n) = (gi_ia(n-1) + 2) / 2
             else
                gi_ia(n) = (gi_ia(n-1) + 1) / 2
             end if
             if (mod(gi_ib(n-1),2) == 0) then
                gi_ib(n) = gi_ib(n-1) / 2
             else
                gi_ib(n) = (gi_ib(n-1) + 1) / 2
             end if
          end if
       end do

       ! check that decomposition is possible on all grids
       ! 
       if (minval(gi_ib(1:ng) - gi_ia(1:ng)) <= 0) then
          print *, 'Process ', myid, ', has only 1 layer. 2 ghost layers needed'
          print *, 'in i-direction'
          call mpi_abort(mpi_comm_world, errcode, ierr)
          stop
       end if
    end if

    if (dims(2,myzone) == 1) then

       ! number ghost points
       jgp(1:ng) = 0

       ! j direction fixed, global coordinates
       gi_ja(1:ng) = 1
       gi_jb(1:ng) = jmg(1:ng,myzone)

    else 

       ! number ghost points
       jgp(1:ng) = 2

       ! set up j-direction decomposition on the coarse grids
       !
       ! if process owns grid node on finest grid then it owns
       ! grid node on all coarser grids
       !
       do n = 2, ng
          if (jcrs(myzone) == 0) then
             gi_ja(n) = gi_ja(1)
             gi_jb(n) = gi_jb(1)
          else
             if (mod(gi_ja(n-1),2) == 0) then
                gi_ja(n) = (gi_ja(n-1) + 2) / 2
             else
                gi_ja(n) = (gi_ja(n-1) + 1) / 2
             end if
             if (mod(gi_jb(n-1),2) == 0) then
                gi_jb(n) = gi_jb(n-1) / 2
             else
                gi_jb(n) = (gi_jb(n-1) + 1) / 2
             end if
          end if
       end do

       ! check that decomposition is possible on all grids
       ! 
       if (minval(gi_jb(1:ng) - gi_ja(1:ng)) <= 0) then
          print *, 'Process ', myid, ', has only 1 layer. 2 ghost layers needed'
          print *, 'in j-direction'
          call mpi_abort(mpi_comm_world, errcode, ierr)
          stop
       end if
    end if

    if (dims(3,myzone) == 1) then

       ! number ghost points
       kgp(1:ng) = 0

       ! k direction fixed, global coordinates
       gi_ka(1:ng) = 1
       gi_kb(1:ng) = kmg(1:ng,myzone)

    else 

       ! number ghost points
       kgp(1:ng) = 2

       ! set up k-direction decomposition on the coarse grids
       !
       ! if process owns grid node on finest grid then it owns
       ! grid node on all coarser grids
       !
       do n = 2, ng
          if (kcrs(myzone) == 0) then
             gi_ka(n) = gi_ka(1)
             gi_kb(n) = gi_kb(1)
          else
             if (mod(gi_ka(n-1),2) == 0) then
                gi_ka(n) = (gi_ka(n-1) + 2) / 2
             else
                gi_ka(n) = (gi_ka(n-1) + 1) / 2
             end if
             if (mod(gi_kb(n-1),2) == 0) then
                gi_kb(n) = gi_kb(n-1) / 2
             else
                gi_kb(n) = (gi_kb(n-1) + 1) / 2
             end if
          end if
       end do

       ! check that decomposition is possible on all grids
       ! 
       if (minval(gi_kb(1:ng) - gi_ka(1:ng)) <= 0) then
          print *, 'Process ', myid, ', has only 1 layer. 2 ghost layers needed'
          print *, 'in k-direction'
          call mpi_abort(mpi_comm_world, errcode, ierr)
          stop
       end if
    end if

    ! mpi-osg only
    ! print ijk range in each processor
    !
    myunit = myid! + 10 !50
    write (filename, fmt = '(a,i3.3)') 'decofiles/gi.ijkz', myunit
    open (unit = myunit, file = trim(filename), form = 'formatted')
    write(unit = myunit, fmt = *) gi_ia(1), gi_ib(1), &
                                  gi_ja(1), gi_jb(1), &
                                  gi_ka(1), gi_kb(1), myzone

    close(unit = myunit)

    ! local processor, interior == exclude ghost points
    !
    li_ia(:) = 1 ; li_ib(:) = gi_ib(:) - gi_ia(:) + 1
    li_ja(:) = 1 ; li_jb(:) = gi_jb(:) - gi_ja(:) + 1
    li_ka(:) = 1 ; li_kb(:) = gi_kb(:) - gi_ka(:) + 1

    ! local processor, exterior == include ghost points
    !
    le_ia(:) = 1 - igp(:) ; le_ib(:) = li_ib(:) + igp(:)
    le_ja(:) = 1 - jgp(:) ; le_jb(:) = li_jb(:) + jgp(:)
    le_ka(:) = 1 - kgp(:) ; le_kb(:) = li_kb(:) + kgp(:)
    
    ! global numbering, exterior == include ghost points
    !
    ge_ia(:) = gi_ia(:) - igp(:)
    ge_ja(:) = gi_ja(:) - jgp(:)
    ge_ka(:) = gi_ka(:) - kgp(:)

    ge_ib(:) = gi_ib(:) + igp(:)
    ge_jb(:) = gi_jb(:) + jgp(:)
    ge_kb(:) = gi_kb(:) + kgp(:)

    ! number (x = max) of grid nodes for each
    ! way of numbering grids
    ! 
    li_ix(:) = li_ib(:) - li_ia(:) + 1
    li_jx(:) = li_jb(:) - li_ja(:) + 1
    li_kx(:) = li_kb(:) - li_ka(:) + 1

    le_ix(:) = le_ib(:) - le_ia(:) + 1
    le_jx(:) = le_jb(:) - le_ja(:) + 1
    le_kx(:) = le_kb(:) - le_ka(:) + 1
    
    ! recalculate gi_kx(1) ; change order if eliminate kmg(:)
    ! 
    gi_ix(:) = gi_ib(:) - gi_ia(:) + 1
    gi_jx(:) = gi_jb(:) - gi_ja(:) + 1
    gi_kx(:) = gi_kb(:) - gi_ka(:) + 1

    ge_ix(:) = ge_ib(:) - ge_ia(:) + 1
    ge_jx(:) = ge_jb(:) - ge_ja(:) + 1
    ge_kx(:) = ge_kb(:) - ge_ka(:) + 1

    ! make indexes (idx == index);
    ! 
    ! the lower bound of these indexes must be set explicity
    ! otherwise this they won't work --- and it is a very hard
    ! bug to find!
    !
    allocate( &
         li_idx(li_ia(1):li_ib(1),li_ja(1):li_jb(1),li_ka(1):li_kb(1),1:ng), &
         le_idx(le_ia(1):le_ib(1),le_ja(1):le_jb(1),le_ka(1):le_kb(1),1:ng), &
         gi_idx(minval(gi_ia(:)):maxval(gi_ib(:)),&
                minval(gi_ja(:)):maxval(gi_jb(:)),&
                minval(gi_ka(:)):maxval(gi_kb(:)),&
                1:ng), &
         gi_2_le_idx(minval(gi_ia(:)):maxval(gi_ib(:)),&
                     minval(gi_ja(:)):maxval(gi_jb(:)),&
                     minval(gi_ka(:)):maxval(gi_kb(:)),&
                     1:ng) ) !, &
         !ge_2_le_idx(minval(ge_ia(:)):maxval(ge_ib(:)),&
         !            minval(ge_ja(:)):maxval(ge_jb(:)),&
         !            minval(ge_ka(:)):maxval(ge_kb(:)),&
         !            1:ng) )

    ! index; local; interior nodes no ghost points
    ! 
    li_idx_a(1) = 1
    li_idx_b(1) = li_ix(1) * li_jx(1) * li_kx(1)
    li_idx_mx(1) = 0
    
    do n = 2, ng
       li_idx_a(n) = li_idx_b(n-1) + 1
       li_idx_b(n) = li_idx_a(n) + li_ix(n) * li_jx(n) * li_kx(n) - 1
       li_idx_mx(n) = li_idx_mx(n-1) + li_ix(n-1) * li_jx(n-1) * li_kx(n-1)
    end do
    
    l = 0
    do n = 1, ng
       do k = li_ka(n), li_kb(n)
          do j = li_ja(n), li_jb(n)
             do i = li_ia(n), li_ib(n)
                l = l + 1
                li_idx(i,j,k,n) = l
             end do
          end do
       end do
    end do
    
    ! index; local; exterior nodes too includes ghost points
    !
    ! this index must be numbered sequentially instead of
    ! calculating it using the formula because le_ka(:) doesn't
    ! start from 1; which is assumed by le_idx_a(1) = 1 and
    ! le_idx_b(1) = le_ix * le_jx * le_kx
    ! 
    le_idx_a(1) = 1
    le_idx_b(1) = le_ix(1) * le_jx(1) * le_kx(1)
    le_idx_mx(1) = 0
    
    do n = 2, ng
       le_idx_a(n) = le_idx_b(n-1) + 1
       le_idx_b(n) = le_idx_a(n) + le_ix(n) * le_jx(n) * le_kx(n) - 1
       le_idx_mx(n) = le_idx_mx(n-1) + le_ix(n-1) * le_jx(n-1) * le_kx(n-1)
    end do
    
    l = 0
    do n = 1, ng
       do k = le_ka(n), le_kb(n)
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1
                le_idx(i,j,k,n) = l
             end do
          end do
       end do
    end do

    ! index; global; interior nodes no ghost points;
    !
    ! this index must be numbered sequentially instead of
    ! calculating it using the formula because gi_ka(:) doesn't
    ! start from 1; which is assumed by gi_idx_a(1) = 1 and
    ! gi_idx_b(1) = gi_ix * gi_jx * gi_kx
    ! 
    gi_idx_a(1) = 1
    gi_idx_b(1) = gi_ix(1) * gi_jx(1) * gi_kx(1)
    gi_idx_mx(1) = 0
    
    do n = 2, ng
       gi_idx_a(n) = gi_idx_b(n-1) + 1
       gi_idx_b(n) = gi_idx_a(n) + gi_ix(n) * gi_jx(n) * gi_kx(n) - 1
       gi_idx_mx(n) = gi_idx_mx(n-1) + gi_ix(n-1) * gi_jx(n-1) * gi_kx(n-1)
    end do
    
    l = 0
    do n = 1, ng
       do k = gi_ka(n), gi_kb(n)
          do j = gi_ja(n), gi_jb(n)
             do i = gi_ia(n), gi_ib(n)
                l = l + 1
                gi_idx(i,j,k,n) = l
             end do
          end do
       end do
    end do

    ! conversion vectors from global numbering to
    ! local numbering

    do n = 1, ng
       do k = gi_ka(n), gi_kb(n)
          do j = gi_ja(n), gi_jb(n)
             do i = gi_ia(n), gi_ib(n)
                kk = k - gi_ka(n) + 1
                jj = j - gi_ja(n) + 1
                ii = i - gi_ia(n) + 1
                gi_2_le_idx(i,j,k,n) = le_idx(ii,jj,kk,n)
             end do
          end do
       end do
    end do

    !do n = 1, ng
    !   do k = ge_ka(n), ge_kb(n)
    !      do j = ge_ja(n), ge_jb(n)
    !         do i = ge_ia(n), ge_ib(n)
    !            kk = k - ge_ka(n) + 1
    !            jj = j - ge_ja(n) + 1
    !            ii = i - ge_ia(n) + 1
    !            ge_2_le_idx(i,j,k,n) = le_idx(ii,jj,kk,n)
    !         end do
    !      end do
    !   end do
    !end do



  ! check via output file
  !
#ifdef DEBUG
  myunit = myid + 50

  write (filename, fmt = '(a,i2.2)') 'nabor', myunit

  open (unit = myunit, file = trim(filename), form = 'formatted')

  write (unit = myunit, fmt = '(a15,3x,g)') 'myid',   myid
  write (unit = myunit, fmt = '(a15,3x,g)') 'myzone', myzone
  write (unit = myunit, fmt = '(a15,3x,g)') 'myback', myback
  write (unit = myunit, fmt = '(a15,3x,g)') 'myfront',myfront
  write (unit = myunit, fmt = '(a15,3x,g)') 'myleft', myleft
  write (unit = myunit, fmt = '(a15,3x,g)') 'myright',myright
  write (unit = myunit, fmt = '(a15,3x,g)') 'mydown', mydown
  write (unit = myunit, fmt = '(a15,3x,g)') 'myup',   myup

  close (unit = myunit)
#endif

  end subroutine init_decomp_3d


  ! --

 subroutine init_mg ()

    integer :: la
    integer :: lb

    ! allocate variables locally for each process
    ! need to include ghost points

    la = le_idx_a(1)
    lb = le_idx_b(ng)

!     print *, myid, la, lb

    allocate(x(la:lb),            &
             y(la:lb),            &
             z(la:lb),            &
             aj(la:lb),           &
             csi(1:3,la:lb),      &
             eta(1:3,la:lb),      &
             zet(1:3,la:lb),      &
             q(1:me,la:lb),       &
             qn(1:me,la:lb),      &
             qnm1(1:me,la:lb),    &
             qold_mg(1:me,la:lb), &
             rh(1:4,la:lb),       &
             pk(1:4,la:lb))
!             dtau(la:lb),         &

    allocate(xc(la:lb),           &
                            xe(la:lb),                       &
                            xz(la:lb),                       &
                            yc(la:lb),                       &
                            ye(la:lb),                       &
                            yz(la:lb),                       &
                            zc(la:lb),                       &
                            ze(la:lb),                       &
                            zz(la:lb))
           xc = zero; xe = zero; xz = zero
           yc = zero; ye = zero; yz = zero
           zc = zero; ze = zero; zz = zero


    !if (daf) allocate(mai(1:4,1:4,la:lb), &
    !                  n1i(1:4,1:4,la:lb), &
    !                  n2i(1:4,1:4,la:lb), &
    !                   mc(1:4,1:4,la:lb), &
    !                      spr(1:3,la:lb))
    
    if ( turbulence ) allocate (wd(la:lb), &
                              dtev(la:lb), &
                              xnut(la:lb))
    
    if ( nlinc ) allocate (uij(1:6,la:lb))



    ! initialize variables
    ! 
    x = zero; y = zero; z = zero
    csi = zero; eta = zero; zet = zero
    aj = zero

!     ret = one / ren

    q = zero; qn = zero; qnm1 = zero; qold_mg = zero

    rh = zero; pk = zero

    !if (daf) then
    !   mai = zero; n1i = zero; n2i = zero; mc = zero; spr = zero
    !end if

    if ( turbulence ) then
       wd = zero; dtev = zero; xnut = zero
    end if

    if ( nlinc ) uij = zero


       
    ! allocate vars needed for input files;
    !
    allocate ( convec(ng), &
                   dc(ng), &
                   de(ng), &
                   dz(ng)  )

    ! move to input file and init_input
    ! 
    convec(1) = 'nc_quick'
    convec(2:ng) = 'nc_first'

    do n = 1, ng
       dc(n) = one / two**real(((n - 1) * icrs(myzone)), kind = rdf)
       de(n) = one / two**real(((n - 1) * jcrs(myzone)), kind = rdf)
       dz(n) = one / two**real(((n - 1) * kcrs(myzone)), kind = rdf)
    end do


  end subroutine init_mg

!===================================================================================================
!init_disk:
!----------
! Initializates disk variables.
! ..................................................................................................
! cy (-)                = j-coordinate of disk's center.
! cz (-)                = k-coordinate of disk's center.
! radius (-)            = Distance from node to disk center in every x-plane.  
! dA (-)                = Node's frontal area element.
! posx (-)              = Node number of disk center in i-direction.
! counter (-)           = Dummy variable. Counts blade point being evaluated.
! counter2 (-)          = Dummy variable. Counts number of grid points inside each actuator disk.
! rad1 (m)              = Dummy variable for airfoil calculations.
! rad2 (m)              = Dummy variable for airfoil calculations.
! aux1 (-)              = Dummy variable for airfoil calculations.
! aux1 (-)              = Dummy variable for airfoil calculations.
! m (-)                 = Slope for chord and twist interpolations.
! xd (-)                = Grid in i-direction.
! yd (-)                = Grid in j-direction.
! zd (-)                = Grid in k-direction.
! blade_geom(n,1) (m)   = Radius of nth blade point in blade.txt.
! blade_geom(n,2) (deg) = Twist angle of nth blade point in blade.txt.
! blade_geom(n,3) (m)   = Chord length of nth blade point in blade.txt.
! blade_geom(n,4) (-)   = Relative thickness of nth blade point in blade.txt.
! airfoil_name(n)       = Airfoil of nth blade point in blade.txt.
! LX 					= Length of the turbine cylinder volume [adim]
! Xcent					= i-coordinate of the cylinder volume center [adim]
! SIGMA					= Standard deviation for regularization kernel [adim]
! ..................................................................................................
!===================================================================================================
  subroutine init_lines ()
  use global_turb

  real :: cy, cz, radius, dA, LX, SIGMA!, Xcent, SIGMA
  integer :: posx, ref_vel, counter, counter2, counter3, i
  real :: rad1, rad2, aux1, aux2, m  
  real (kind = rdf), dimension(:,:,:), allocatable :: xd
  real (kind = rdf), dimension(:,:,:), allocatable :: zd
  real (kind = rdf), dimension(:,:,:), allocatable :: yd
  real (kind = rdf), dimension(:,:), allocatable :: blade_geom
  character airfoil_name(30)*10


  !======================================
  ! ALM local variables for init_lines  !
  ! J. Sandoval, 2018                   !
  !======================================

  real :: theta, xelem, yelem, zelem, xcent, ycent, zcent, dis(6)
  integer :: ielem, jelem, kelem, nang, elem, bld, inside
  integer :: isearch, jsearch, ksearch, nodes_search!, nangles_aux
  real :: epsALM, dV, rdis, xker, yker, zker
  !integer :: counter_nodes, iker, jker, kker
  real::point(3),cell(8,3)

!  print *,'Reading grid...',myid
  open  (unit = 51, file = 'grid', form='unformatted')
  do nz = 1, nzone
    allocate (xd(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz)), &
              yd(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz)), &
              zd(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz)) )

    read  (unit = 51) (((xd(i,j,k),i=1,img(1,nz)), &
                                     j=1,jmg(1,nz)),k=1,kmg(1,nz))
    read  (unit = 51) (((yd(i,j,k),i=1,img(1,nz)), &
                                     j=1,jmg(1,nz)),k=1,kmg(1,nz))
    read  (unit = 51) (((zd(i,j,k),i=1,img(1,nz)), &
                                     j=1,jmg(1,nz)),k=1,kmg(1,nz))
  end do 
  close (unit = 51)

  nz = 1
  !----------------------------------------------------------
  ! Read blade geometry (blade.txt), including airfoil type.
  ! Note: Text file must contain no more than 30 lines.
  !----------------------------------------------------------
  allocate (blade_geom(30,4))
  call read_geom(blade_geom, airfoil_name)

  !----------------------------------------------------------
  ! Loop for initializing every actuator disk.
  !----------------------------------------------------------
  do n_t=1,n_tur

    counter2 = 0
    counter3 = 0

    !----------------------------------------------------------
    ! Geometry of disk center
    !----------------------------------------------------------
    ! Xcent as average position given by itu1(n_t) and itur2(n_t) (local variable)
    xcent = 0.5*(xd(itur1(n_t),jtur1(n_t),ktur1(n_t))+xd(itur2(n_t),jtur1(n_t),ktur1(n_t))) 
    LX = 0.5*(xd(itur2(n_t)+1,jtur1(n_t),ktur1(n_t))+xd(itur2(n_t),jtur1(n_t),ktur1(n_t)) - &
    			xd(itur1(n_t),jtur1(n_t),ktur1(n_t))-xd(itur1(n_t)-1,jtur1(n_t),ktur1(n_t)))
    SIGMA = 0.5*LX

    ! print*,'---------'
    ! print*,'Xcent:',Xcent,'LX',LX,'SIGMA',SIGMA
    ! print*,'---------'

    !----------------------------------------------------------
    !Regularization kernel in [1/m] 
    !----------------------------------------------------------
    do i=itur1(n_t),itur2(n_t)

    	kernel(i-itur1(n_t)+1,n_t) = exp(-0.5*((xd(i,jtur1(n_t),ktur1(n_t))-xcent)/SIGMA)**2)/sqrt(2*pi*(SIGMA*length_scale)**2)

      ! print*,'i:',i,'kernel',kernel(i-itur1(n_t),n_t)
    end do

    !----------------------------------------------------------
    !Since every x-plane has the same j,k distribution we can use itur1 for definzng
    !every other coordinate
    !----------------------------------------------------------
    !posx = itur1(n_t)
    posx = NINT((itur1(n_t)+itur2(n_t))/2.0)

    !print *,'itur1,itur2,posx',itur1(n_t),itur2(n_t),posx

    cy   = (yd(posx,jtur1(n_t),ktur1(n_t))+yd(posx,jtur2(n_t),ktur1(n_t)))/2
    cz   = (zd(posx,jtur1(n_t),ktur1(n_t))+zd(posx,jtur1(n_t),ktur2(n_t)))/2
      
    do j=1,jmg(1,nz)
      do k=1,kmg(1,nz)
        radius = sqrt((yd(posx,j,k)-cy)**2+(zd(posx,j,k)-cz)**2)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Indexation of nodes and positions for forces applying in the
        ! kernel formulation for ALM. JS, 2019. ALM_src_v3.
        ! The volume considered is a cylinder of 1.2*tur_rad(n_t).
        ! With this correction we seek to incoporate the forces
        ! generated at the tip in the radial direction

!        if (radius.le.(1.2*(tur_rad(n_t)/length_scale))) then
        if (radius.le.(1.3*(tur_rad(n_t)/length_scale))) then
            do i=itur1(n_t),itur2(n_t)
              
                !---------------------------------------------------------------------
                ! By definition first element of disk_node_coordinates is the number
                ! of nodes inside disk area.
                !---------------------------------------------------------------------
                counter3                                = counter3 + 1
                disk_node_coordinates_ker(1,1,n_t)      = counter3
                disk_node_coordinates_ker(counter3+1,1,n_t) = i
                disk_node_coordinates_ker(counter3+1,2,n_t) = j
                disk_node_coordinates_ker(counter3+1,3,n_t) = k
      
                disk_node_coordinates_pos_ker(1,1,n_t)          = counter3
                disk_node_coordinates_pos_ker(counter3+1,1,n_t) = xd(i,j,k) ! For ALM_v3. JS, 2019
                disk_node_coordinates_pos_ker(counter3+1,2,n_t) = yd(i,j,k) ! For ALM_v3. JS, 2019
                disk_node_coordinates_pos_ker(counter3+1,3,n_t) = zd(i,j,k) ! For ALM_v3. JS, 2019
      
                disk_dat_ker(i,j,k,1,n_t) = radius  
      
                !frontal area element
                disk_dat_ker(i,j,k,2,n_t) = ABS(yd(i,j-1,k)-yd(i,j+1,k))*ABS(zd(i,j,k-1)-zd(i,j,k+1))/4
                dA = ABS(yd(i,j-1,k)-yd(i,j+1,k))*ABS(zd(i,j,k-1)-zd(i,j,k+1))/4
                !---------------------------------------------------------------------
                ! Angle between node and disk center is positive clockwise following
                ! David Ingram's definition.
                !---------------------------------------------------------------------
                !if ( yd(i,j,k)-cy .ge. 0 .and. zd(i,j,k)-cz .ge. 0  ) then 
                !    disk_dat_ker(i,j,k,3,n_t)= asin((zd(i,j,k)-cz )/disk_dat_ker(i,j,k,1,n_t))
                !else if ( yd(i,j,k)-cy .le. 0 .and. zd(i,j,k)-cz .ge. 0  ) then   
                !    disk_dat_ker(i,j,k,3,n_t)= pi - asin((zd(i,j,k)-cz )/disk_dat_ker(i,j,k,1,n_t))
                !else if ( yd(i,j,k)-cy .le. 0 .and. zd(i,j,k)-cz .le. 0  ) then   
                !    disk_dat_ker(i,j,k,3,n_t)= pi + asin((cz-zd(i,j,k))/disk_dat_ker(i,j,k,1,n_t))
                !else 
                !   disk_dat_ker(i,j,k,3,n_t)= 2*pi - asin((cz-zd(i,j,k))/disk_dat_ker(i,j,k,1,n_t))
                !end if
                ! the 3rd position of this variable is not used (JS, 2019)
                disk_dat_ker(i,j,k,3,n_t)= 0.0 ! [-]. JS, 2018
                disk_dat_ker(i,j,k,4,n_t)= (dA*ABS(xd(i+1,j,k)-xd(i-1,j,k))/2) ! [-]. JS, 2018

            end do
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      	if (radius.le.(tur_rad(n_t)/length_scale)) then

          !---------------------------------------------------------------------
          ! By definition first element of disk_node_coordinates is the number
          ! of nodes inside disk area.
          !---------------------------------------------------------------------
          counter2                                = counter2 + 1
          disk_node_coordinates(1,1,n_t)          = counter2
          disk_node_coordinates(counter2+1,1,n_t) = posx
          disk_node_coordinates(counter2+1,2,n_t) = j
          disk_node_coordinates(counter2+1,3,n_t) = k

          disk_node_coordinates_pos(1,1,n_t)          = counter2
          disk_node_coordinates_pos(counter2+1,1,n_t) = xd(posx,j,k) ! For ALM. JS, 2018
          disk_node_coordinates_pos(counter2+1,2,n_t) = yd(posx,j,k) ! For ALM. JS, 2018
          disk_node_coordinates_pos(counter2+1,3,n_t) = zd(posx,j,k) ! For ALM. JS, 2018


          disk_dat(j,k,1,n_t) = radius	

          !frontal area element
          disk_dat(j,k,2,n_t) = ABS(yd(posx,j-1,k)-yd(posx,j+1,k))*ABS(zd(posx,j,k-1)-zd(posx,j,k+1))/4
      	  dA = ABS(yd(posx,j-1,k)-yd(posx,j+1,k))*ABS(zd(posx,j,k-1)-zd(posx,j,k+1))/4
          !---------------------------------------------------------------------
          ! Angle between node and disk center is positive clockwise following
          ! David Ingram's definition.
          !---------------------------------------------------------------------
          if ( yd(posx,j,k)-cy .ge. 0 .and. zd(posx,j,k)-cz .ge. 0  ) then 
              disk_dat(j,k,3,n_t)= asin((zd(posx,j,k)-cz )/disk_dat(j,k,1,n_t))
          else if ( yd(posx,j,k)-cy .le. 0 .and. zd(posx,j,k)-cz .ge. 0  ) then   
              disk_dat(j,k,3,n_t)= pi - asin((zd(posx,j,k)-cz )/disk_dat(j,k,1,n_t))
        	else if ( yd(posx,j,k)-cy .le. 0 .and. zd(posx,j,k)-cz .le. 0  ) then   
              disk_dat(j,k,3,n_t)= pi + asin((cz-zd(posx,j,k))/disk_dat(j,k,1,n_t))
        	else 
        	   disk_dat(j,k,3,n_t)= 2*pi - asin((cz-zd(posx,j,k))/disk_dat(j,k,1,n_t))
        	end if 
        
          disk_dat(j,k,4,n_t)= (dA*ABS(xd(posx+1,j,k)-xd(posx-1,j,k))/2) ! [-]. JS, 2018
        	!disk_dat(j,k,4,n_t) = 0.0

          !---------------------------------------------------------------------
          ! Calculations for chord and twist using blade geometry.
          ! Note: radius is scaled to meters for comparison with information in
          ! blade.txt, but data in disk_dat is stored non-dimensionally.
          !---------------------------------------------------------------------
          counter = 1
          rad1 = 10
          rad2 = 10
          radius = radius*length_scale

          !---------------------------------------------------------------------
          ! This loops finds in between which blade points is the node that's
          ! being evaluated. Stops when rad2 = -1 (end of blade.txt).
          ! For now, if radius is smaller thant first data point in blade.txt
          ! it assings the value of blade.txt's first row. Should include
          ! nacelle in the future (ask Cristian about it).
          !---------------------------------------------------------------------
          do while(rad2.gt.0)
            rad1 = blade_geom(counter,1)
            rad2 = blade_geom(counter+1,1)
            if (rad2.eq.-1) then
              print*,'Error calculating blade geometry.'
            else if (counter.eq.1 .and. radius.le.rad1) then
              !---------------------------------------------------------------------
              ! HANDLING POINT WITH RADIUS LESS THAN FIRST DATA IN BLADE.TXT
              !---------------------------------------------------------------------
              disk_dat(j,k,5,n_t)     = blade_geom(counter,3) !chord (me)
              !disk_dat(j,k,5,n_t)     = 0 !chord (me) !(IF IN NACELLE CHORD = 0)
              disk_dat(j,k,6,n_t)     = blade_geom(counter,2) !twist (deg)
              airfoil_before(j,k,n_t) = airfoil_name(counter)
              airfoil_after(j,k,n_t)  = airfoil_name(counter)
              rad2= -1
            else if ( (radius <= rad2) .and. (radius >= rad1) ) then

              !Interpolation for chord length [m]
              aux1                    = blade_geom(counter,3)
              aux2                    = blade_geom(counter+1,3)
              m                       = (aux2-aux1)/(rad2-rad1)
              disk_dat(j,k,5,n_t)     = m*(radius-rad1)+aux1
              
              !Interpolation for relative thickness [-] (J.Sandoval, 2018)

              aux1                    = blade_geom(counter,4)
              aux2                    = blade_geom(counter+1,4)
              m                       = (aux2-aux1)/(rad2-rad1)
              disk_dat(j,k,8,n_t)     = m*(radius-rad1)+aux1

              !Interpolation for twist angle
              aux1                    = blade_geom(counter,2)
              aux2                    = blade_geom(counter+1,2)
              m                       = (aux2-aux1)/(rad2-rad1)
              disk_dat(j,k,6,n_t)     = m*(radius-rad1)+aux1
              
              !Storing of the interpolation coefficient.
              disk_dat(j,k,7,n_t)     = (radius-rad1)/(rad2-rad1)
              
              !Storing name of airfoils before and after.
              airfoil_before(j,k,n_t) = airfoil_name(counter)
              airfoil_after(j,k,n_t)  = airfoil_name(counter+1)
	!	if (myid == root) then
	!		print*,'DEBUG INSIDE R1<R<ÂR2'
	!		print*,'a.bef:',airfoil_before(j,k,n_t),'|',airfoil_name(counter)
	!		print*,'a.aft:',airfoil_after(j,k,n_t),'|',airfoil_name(counter+1)
	!	end if
              rad2=-1
            end if
            counter = counter + 1
          end do
!	if (myid == root) then		
!		print*,'Point #:',counter2
!		print*,'Air bef:',airfoil_before(j,k,n_t)
!		print*,'Air aft:Â',airfoil_after(j,k,n_t)
!	end if


      	else  	
          !If node is not in disk surface disk_dat still stores the distance.
      		disk_dat(j,k,1,n_t) = radius
      	end if
      end do
    end do

   ! print*,'--------------------------------------------------------'
   ! print*,'Disk number ',n_t,' initialized with ', counter2, 'nodes inside it.'
  end do
 ! print*,'--------------------------------------------------------'

  if (myid == root) then
    print*, 'Total number of nodes for ker application: ', counter3
    print*, 'Total number of nodes inside the disk: ', counter2
  endif


!======================================
! ALM blades initialization           !
! J. Sandoval, 2018                   !
!======================================

  ! Initialization of blade elements positions
  do n_t=1,n_tur
      theta = angle(n_t)
      do bld=1,nblades(n_t)
        do elem = 1,nelems(n_t)
          dtheta(n_t) = TSR(n_t)*delti*length_scale/tur_rad(n_t)
          
          !posx = itur1(n_t)
          posx = NINT((itur1(n_t)+itur2(n_t))/2.0)

          xcent = 0.5*(xd(itur1(n_t),jtur1(n_t),ktur1(n_t))+xd(itur2(n_t),jtur1(n_t),ktur1(n_t)))
          ycent = (yd(posx,jtur1(n_t),ktur1(n_t))+yd(posx,jtur2(n_t),ktur1(n_t)))/2
          zcent = (zd(posx,jtur1(n_t),ktur1(n_t))+zd(posx,jtur1(n_t),ktur2(n_t)))/2

          blade_elemts(n_t,bld,elem,1)=xcent
          blade_elemts(n_t,bld,elem,2)=ycent + ((tur_rad(n_t)/length_scale)/nelems(n_t))*elem*cos(theta)
          blade_elemts(n_t,bld,elem,3)=zcent + ((tur_rad(n_t)/length_scale)/nelems(n_t))*elem*sin(theta)          
        end do ! do elem = 1,nelems(n_t)
        blade_angles(n_t,bld)=theta
        theta=theta+2*pi/nblades(n_t)
      end do ! do bld=1,nblades(n_t)
  end do ! do n_t=1,n_tur

  ! Node location for diferent blade position indexation

  do n_t=1,n_tur
    dtheta(n_t) = TSR(n_t)*delti*length_scale/tur_rad(n_t)
    nangles(n_t) = NINT(100*2*pi/dtheta(n_t))

    dtheta_interp(n_t) = 2*pi/nangles(n_t)


    if (myid == root) then
    	print *, 'nangles: ', nangles(n_t)
    	print *, 'dtheta: ', dtheta(n_t)
    	print *, 'dtheta_interp', dtheta_interp(n_t)
    endif


    !posx = itur1(n_t)
    posx = NINT((itur1(n_t)+itur2(n_t))/2.0)

    xcent = xd(posx,jtur1(n_t),ktur1(n_t))! 0.5*(xd(itur1(n_t),jtur1(n_t),ktur1(n_t))+xd(itur2(n_t),jtur1(n_t),ktur1(n_t)))
    ycent = (yd(posx,jtur1(n_t),ktur1(n_t))+yd(posx,jtur2(n_t),ktur1(n_t)))/2
    zcent = (zd(posx,jtur1(n_t),ktur1(n_t))+zd(posx,jtur1(n_t),ktur2(n_t)))/2

    theta = 0.0

    !nangles_aux=nangles(n_t)

    if(node_indexes_load) then

     if(myid==root) print *, 'Leyendo archivo node_indexes y node_indexes_pos...'
     open  (unit = 18, file = 'node_indexes', form='unformatted')
      do nz = 1, nzone
        
        !do n_t_aux=1,n_tur
          read  (unit = 18) (((node_indexes(n_t,elem,nang,cont),elem=1,nelems(n_t)), &
                                           nang=1,nangles(n_t)),cont=1,3)
        !end do
      end do 
    close (unit = 18)

    open  (unit = 21, file = 'node_indexes_pos', form='unformatted')
      do nz = 1, nzone
        
        !do n_t_aux=1,n_tur
          read  (unit = 21) (((node_indexes_pos(n_t,elem,nang,cont),elem=1,nelems(n_t)), &
                                           nang=1,nangles(n_t)),cont=1,27)
        !end do
      end do 
    close (unit = 21)

    if(myid==root) print *, 'Lectura archivo node_indexes y node_indexes_pos finalizada...'
    
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !if(myid==root) print *, 'Leyendo archivo kerALM'
    ! open  (unit = 19, file = 'node_indexes', form='unformatted')
    !  do nz = 1, nzone
    !    
    !    !do n_t_aux=1,n_tur
    !      read  (unit = 19) (((kerALM(n_t,elem,nang,cont),elem=1,nelems(n_t)), &
    !                                       nang=1,nangles(n_t)),cont=1,5000)
    !    !end do
    !  end do 
    !close (unit = 19)
    !if(myid==root) print *, 'Lectura archivo kerALM finalizada...'


    else

    do nang = 1, nangles(n_t)
      theta = theta + dtheta_interp(n_t)

      do elem = 1, nelems(n_t)
        xelem = xcent
        yelem = ycent + ((tur_rad(n_t)/length_scale)/nelems(n_t))*elem*cos(theta)
        zelem = zcent + ((tur_rad(n_t)/length_scale)/nelems(n_t))*elem*sin(theta)

      !if (myid == root) then
      !  print *, 'xelem: ', xelem
      !  print *, 'yelem: ', yelem
      !  print *, 'zelem: ', zelem
      !endif

        point(1)=xelem
        point(2)=yelem
        point(3)=zelem
        
        searchloop:do nodes_search = 1, disk_node_coordinates(1,1,n_t)
          isearch = disk_node_coordinates(nodes_search+1,1,n_t)
          jsearch = disk_node_coordinates(nodes_search+1,2,n_t)
          ksearch = disk_node_coordinates(nodes_search+1,3,n_t)

          !if (myid == root) then
          !	print *, 'nodos de busqueda: ', isearch, jsearch, ksearch
          !	print *, '_______________________________'
          !end if

          cell(1,1)=xd(isearch,jsearch,ksearch)
          cell(1,2)=yd(isearch,jsearch,ksearch)
          cell(1,3)=zd(isearch,jsearch,ksearch)
          cell(2,1)=xd(isearch+1,jsearch,ksearch)
          cell(2,2)=yd(isearch+1,jsearch,ksearch)
          cell(2,3)=zd(isearch+1,jsearch,ksearch)
          cell(3,1)=xd(isearch+1,jsearch+1,ksearch)
          cell(3,2)=yd(isearch+1,jsearch+1,ksearch)
          cell(3,3)=zd(isearch+1,jsearch+1,ksearch)
          cell(4,1)=xd(isearch,jsearch+1,ksearch)
          cell(4,2)=yd(isearch,jsearch+1,ksearch)
          cell(4,3)=zd(isearch,jsearch+1,ksearch)
          cell(5,1)=xd(isearch,jsearch,ksearch+1)
          cell(5,2)=yd(isearch,jsearch,ksearch+1)
          cell(5,3)=zd(isearch,jsearch,ksearch+1)
          cell(6,1)=xd(isearch+1,jsearch,ksearch+1)
          cell(6,2)=yd(isearch+1,jsearch,ksearch+1)
          cell(6,3)=zd(isearch+1,jsearch,ksearch+1)
          cell(7,1)=xd(isearch+1,jsearch+1,ksearch+1)
          cell(7,2)=yd(isearch+1,jsearch+1,ksearch+1)
          cell(7,3)=zd(isearch+1,jsearch+1,ksearch+1)
          cell(8,1)=xd(isearch,jsearch+1,ksearch+1)
          cell(8,2)=yd(isearch,jsearch+1,ksearch+1)
          cell(8,3)=zd(isearch,jsearch+1,ksearch+1)

          call search(point,cell,inside,dis)

          if(inside.eq.1) then

          	!if (myid == root) then
          	!	print *, 'inside = 1'
          	!	print *, isearch
          	!	print *, jsearch
          	!	print *, ksearch
        		!print *, '__________________________'
          	!endif

            node_indexes(n_t,elem,nang,1)=isearch
            node_indexes(n_t,elem,nang,2)=jsearch
            node_indexes(n_t,elem,nang,3)=ksearch

            node_indexes_pos(n_t,elem,nang,1)=xelem
            node_indexes_pos(n_t,elem,nang,2)=yelem
            node_indexes_pos(n_t,elem,nang,3)=zelem

            ! Indexation of adjacent nodes to interpolate the velocity
            ! in rhs_lines (Jorge Sandoval, 2019. ALM_v2_trilinear_interp)

            node_indexes_pos(n_t,elem,nang,4)  = xd(isearch,jsearch,ksearch)
            node_indexes_pos(n_t,elem,nang,5)  = yd(isearch,jsearch,ksearch)
            node_indexes_pos(n_t,elem,nang,6)  = zd(isearch,jsearch,ksearch)
            node_indexes_pos(n_t,elem,nang,7)  = xd(isearch+1,jsearch,ksearch)
            node_indexes_pos(n_t,elem,nang,8)  = yd(isearch+1,jsearch,ksearch)
            node_indexes_pos(n_t,elem,nang,9)  = zd(isearch+1,jsearch,ksearch)
            node_indexes_pos(n_t,elem,nang,10) = xd(isearch+1,jsearch+1,ksearch)
            node_indexes_pos(n_t,elem,nang,11) = yd(isearch+1,jsearch+1,ksearch)
            node_indexes_pos(n_t,elem,nang,12) = zd(isearch+1,jsearch+1,ksearch)
            node_indexes_pos(n_t,elem,nang,13) = xd(isearch,jsearch+1,ksearch)
            node_indexes_pos(n_t,elem,nang,14) = yd(isearch,jsearch+1,ksearch)
            node_indexes_pos(n_t,elem,nang,15) = zd(isearch,jsearch+1,ksearch)
            node_indexes_pos(n_t,elem,nang,16) = xd(isearch,jsearch,ksearch+1)
            node_indexes_pos(n_t,elem,nang,17) = yd(isearch,jsearch,ksearch+1)
            node_indexes_pos(n_t,elem,nang,18) = zd(isearch,jsearch,ksearch+1)
            node_indexes_pos(n_t,elem,nang,19) = xd(isearch+1,jsearch,ksearch+1)
            node_indexes_pos(n_t,elem,nang,20) = yd(isearch+1,jsearch,ksearch+1)
            node_indexes_pos(n_t,elem,nang,21) = zd(isearch+1,jsearch,ksearch+1)
            node_indexes_pos(n_t,elem,nang,22) = xd(isearch+1,jsearch+1,ksearch+1)
            node_indexes_pos(n_t,elem,nang,23) = yd(isearch+1,jsearch+1,ksearch+1)
            node_indexes_pos(n_t,elem,nang,24) = zd(isearch+1,jsearch+1,ksearch+1)
            node_indexes_pos(n_t,elem,nang,25) = xd(isearch,jsearch+1,ksearch+1)
            node_indexes_pos(n_t,elem,nang,26) = yd(isearch,jsearch+1,ksearch+1)
            node_indexes_pos(n_t,elem,nang,27) = zd(isearch,jsearch+1,ksearch+1)


            !print *, 'i: ', node_indexes(n_t,elem,nang,1)
            !print *, 'j: ', node_indexes(n_t,elem,nang,2)
            !print *, 'k: ', node_indexes(n_t,elem,nang,3)
            !print *, '________________________________________________'
            !print *, 'x: ', node_indexes_pos(n_t,elem,nang,1)
            !print *, 'y: ', node_indexes_pos(n_t,elem,nang,2)
            !print *, 'z: ', node_indexes_pos(n_t,elem,nang,3)
            !print *, '////////////////////////////////////////////////////////////'

          exit searchloop
          end if !if(inside.eq.1) then
        end do searchloop !searchloop:do nodes_search = 1, disk_node_coordinates(1,1,n_t)
      end do ! do elem = 1, nelems(n_t)
    end do ! do nang = 1, nangles(n_t)

    ! Writing of node indexes file after search algorithm

    if(myid==root) print *, 'Escribiendo archivo node_indexes y node_indexes_pos...'

    open  (unit = 11, file = 'node_indexes', form='unformatted')
      do nz = 1, nzone

       !do n_t=1,n_tur
          write  (unit = 11) (((node_indexes(n_t,elem,nang,cont),elem=1,nelems(n_t)), &
                                           nang=1,nangles(n_t)),cont=1,3)      
       !end do
      end do
    close(unit = 11)

    open  (unit = 17, file = 'node_indexes_pos', form='unformatted')
      do nz = 1, nzone

       !do n_t=1,n_tur
          write  (unit = 17) (((node_indexes_pos(n_t,elem,nang,cont),elem=1,nelems(n_t)), &
                                           nang=1,nangles(n_t)),cont=1,27)      
       !end do
      end do
    close(unit = 17)

    if(myid==root) print *, 'Archivo node_indexes y node_indexes_pos finalizado...'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !if(myid==root) print *, 'Escribiendo archivo kerALM...'

    !open  (unit = 12, file = 'kerALM', form='unformatted')
    !  do nz = 1, nzone

    !   !do n_t=1,n_tur
    !      write  (unit = 12) (((kerALM(n_t,elem,nang,cont),elem=1,nelems(n_t)), &
    !                                       nang=1,nangles(n_t)),cont=1,5000)      
    !   !end do
    !  end do
    !close(unit = 12)

    !if(myid==root) print *, 'Archivo kerALM finalizado...'



    end if !if(node_indexes_load) then. Boolean to load node_indexes file
  end do ! do n_t=1,n_tur

  deallocate (xd, yd, zd, blade_geom)

	!if (myid == root) then
  	!	print*,"Finished initializing actuator disks."
 	!	do n_t=1,n_tur
  	!		print*,'Disk Nº',n_t, ' has ', disk_node_coordinates(1,1,n_t), ' points.'
  	!	end do
  	!	print*,'--------------------------------------------------------'
!
  	!	print *, node_indexes
	!end if
  
end subroutine init_lines
 
!===================================================================================================
!read_geom:
!----------
! Reads blade.txt and stores it's data.
! ..................................................................................................
! counter (-)           = Dummy variable. Counts which data point is being evaluated.
! rad (m)               = Radius of data point being read.
! twist (deg)           = Twist angle of data point being read.
! chord (m)             = Chords length of data point being read.
! thick (-)             = Relative thickness of data point being read.
! name1 (-)             = Airfoil that defines data point being read.
! name2 (-)             = Airfoil that defines data point being read.
! blade_geom(n,1) (m)   = Radius of nth blade point in blade.txt.
! blade_geom(n,2) (deg) = Twist angle of nth blade point in blade.txt.
! blade_geom(n,3) (m)   = Chord length of nth blade point in blade.txt.
! blade_geom(n,4) (-)   = Relative thickness of nth blade point in blade.txt.
! airfoil_name(n)       = Airfoil of nth blade point in blade.txt.
! ..................................................................................................
!===================================================================================================
 subroutine read_geom(blade_geom, airfoil_name)

  integer :: counter
  real :: rad, twist, chord, thick
  real :: blade_geom(30,4)
  character(LEN=10) :: name1, name2
  character  airfoil_name(30)*10

  counter = 1
  rad = 10

  open  (unit = 61, FILE = 'blade.txt', STATUS='OLD')
  do while(rad.gt.0)
    read(61,*) rad, twist, chord, thick,name1, name2
    blade_geom(counter,1) = rad
    blade_geom(counter,2) = twist
    blade_geom(counter,3) = chord
    blade_geom(counter,4) = thick
    airfoil_name(counter) = name1
    counter = counter + 1
!	if (myid == root) then
!		print*,'DEBUG READGEOM'
!		print*,'rad:',rad,'twist',twist,'chord',chord
!		print*,'airfoil ',name1
!	end if
  end do

  close(unit=61)
  return

  end subroutine
 
  ! --

  subroutine init_blanking ()
  !
  ! specify blanking infomation in mpi+overset grid
  !

    integer :: blk_ka, blk_kb, &
               blk_ja, blk_jb, &
               blk_ia, blk_ib

    integer :: blka, blkb
    integer :: nzb

    logical :: i_in , &
               j_in, &
               k_in

    character (len = 255) :: filename
    integer :: myunit

    integer, dimension(:,:), allocatable :: gi_blk_ia, gi_blk_ib, &
                                            gi_blk_ja, gi_blk_jb, &
                                            gi_blk_ka, gi_blk_kb

    allocate (gi_blk_ia(ng,nba_max), &
              gi_blk_ib(ng,nba_max), &
              gi_blk_ja(ng,nba_max), &
              gi_blk_jb(ng,nba_max), &
              gi_blk_ka(ng,nba_max), &
              gi_blk_kb(ng,nba_max) )

    if (nzblanking(myzone) /= 0) then

       nblk = 0
       do nzb = 1, nzblanking(myzone)

          blka = 1 + 2 * (nzb - 1)
          blkb = 2 + 2 * (nzb - 1)

          i_in = .false.
          j_in = .false.
          k_in = .false.

          ! i inside ?
          if (gi_ia(1) <= blanking(1,blka,myzone) .and. &
                          blanking(1,blka,myzone) <= gi_ib(1) ) i_in = .true.
          if (gi_ia(1) <= blanking(1,blkb,myzone) .and. &
                          blanking(1,blkb,myzone) <= gi_ib(1) ) i_in = .true.
          if (gi_ia(1) >= blanking(1,blka,myzone) .and. &
                          blanking(1,blkb,myzone) >= gi_ib(1) ) i_in = .true.

          ! j inside ?
          if (gi_ja(1) <= blanking(2,blka,myzone) .and. &
                          blanking(2,blka,myzone) <= gi_jb(1) ) j_in = .true.
          if (gi_ja(1) <= blanking(2,blkb,myzone) .and. &
                          blanking(2,blkb,myzone) <= gi_jb(1) ) j_in = .true.
          if (gi_ja(1) >= blanking(2,blka,myzone) .and. &
                          blanking(2,blkb,myzone) >= gi_jb(1) ) j_in = .true.

          ! j inside ?
          if (gi_ka(1) <= blanking(3,blka,myzone) .and. &
                          blanking(3,blka,myzone) <= gi_kb(1) ) k_in = .true.
          if (gi_ka(1) <= blanking(3,blkb,myzone) .and. &
                          blanking(3,blkb,myzone) <= gi_kb(1) ) k_in = .true.
          if (gi_ka(1) >= blanking(3,blka,myzone) .and. &
                          blanking(3,blkb,myzone) >= gi_kb(1) ) k_in = .true.

          if (i_in .and. j_in .and. k_in) then

             blk_ia = max(gi_ia(1),blanking(1,blka,myzone))
             blk_ib = min(gi_ib(1),blanking(1,blkb,myzone))
             blk_ja = max(gi_ja(1),blanking(2,blka,myzone))
             blk_jb = min(gi_jb(1),blanking(2,blkb,myzone))
             blk_ka = max(gi_ka(1),blanking(3,blka,myzone))
             blk_kb = min(gi_kb(1),blanking(3,blkb,myzone))

             if( blk_ia > blk_ib .or. &
                 blk_ja > blk_jb .or. &
                 blk_ka > blk_kb ) print*, 'something wrong in blanking'

             nblk = nblk + 1

             gi_blk_ia(1,nblk) = blk_ia
             gi_blk_ib(1,nblk) = blk_ib
             gi_blk_ja(1,nblk) = blk_ja
             gi_blk_jb(1,nblk) = blk_jb
             gi_blk_ka(1,nblk) = blk_ka
             gi_blk_kb(1,nblk) = blk_kb

             ! on the coarse grids
             do n = 2, ng

                ! i-dir
                gi_blk_ia(n,nblk) = (gi_blk_ia(n-1,nblk) + icrs(myzone)) / &
                                    (2**icrs(myzone))
                if (blktype(2,nzb,myzone) /= 0) &
                gi_blk_ib(n,nblk) = (gi_blk_ib(n-1,nblk) + icrs(myzone)  ) / &
                                    (2**icrs(myzone))
                if (blktype(2,nzb,myzone) == 0) &
                gi_blk_ib(n,nblk) = (gi_blk_ib(n-1,nblk) + icrs(myzone)+1) / &
                                    (2**icrs(myzone))
                gi_blk_ia(n,nblk) = max(gi_blk_ia(n,nblk),1)
                gi_blk_ib(n,nblk) = min(gi_blk_ib(n,nblk),gi_ib(n))

                ! j-dir
                gi_blk_ja(n,nblk) = (gi_blk_ja(n-1,nblk) + jcrs(myzone)) / &
                                    (2**jcrs(myzone))
                if (blktype(4,nzb,myzone) /= 0) &
                gi_blk_jb(n,nblk) = (gi_blk_jb(n-1,nblk) + jcrs(myzone)  ) / &
                                    (2**jcrs(myzone))
                if (blktype(4,nzb,myzone) == 0) &
                gi_blk_jb(n,nblk) = (gi_blk_jb(n-1,nblk) + jcrs(myzone)+1) / &
                                    (2**jcrs(myzone))
                gi_blk_ja(n,nblk) = max(gi_blk_ja(n,nblk),1)
                gi_blk_jb(n,nblk) = min(gi_blk_jb(n,nblk),gi_jb(n))

                ! k-dir
                gi_blk_ka(n,nblk) = (gi_blk_ka(n-1,nblk) + kcrs(myzone)) / &
                                    (2**kcrs(myzone))
                if (blktype(6,nzb,myzone) /= 0) &
                gi_blk_kb(n,nblk) = (gi_blk_kb(n-1,nblk) + kcrs(myzone)  ) / &
                                    (2**kcrs(myzone))
                if (blktype(6,nzb,myzone) == 0) &
                gi_blk_kb(n,nblk) = (gi_blk_kb(n-1,nblk) + kcrs(myzone)+1) / &
                                    (2**kcrs(myzone))
                gi_blk_ka(n,nblk) = max(gi_blk_ka(n,nblk),1)
                gi_blk_kb(n,nblk) = min(gi_blk_kb(n,nblk),gi_kb(n))

             end do

          end if

       end do

       allocate (li_blk_ia(ng,nblk), &
                 li_blk_ib(ng,nblk), &
                 li_blk_ja(ng,nblk), &
                 li_blk_jb(ng,nblk), &
                 li_blk_ka(ng,nblk), &
                 li_blk_kb(ng,nblk) )

       ! local processor, interior == exclude ghost points
       ! 
       do nb = 1, nblk
       do n = 1, ng
          li_blk_ia(n,nb) = gi_blk_ia(n,nb) - gi_ia(n) + 1
          li_blk_ib(n,nb) = gi_blk_ib(n,nb) - gi_ia(n) + 1
          li_blk_ja(n,nb) = gi_blk_ja(n,nb) - gi_ja(n) + 1
          li_blk_jb(n,nb) = gi_blk_jb(n,nb) - gi_ja(n) + 1
          li_blk_ka(n,nb) = gi_blk_ka(n,nb) - gi_ka(n) + 1
          li_blk_kb(n,nb) = gi_blk_kb(n,nb) - gi_ka(n) + 1
!#ifdef DEBUG
!print*, 'blanking info: myid, ia, ib, ...  '
!print*, myid, n, li_blk_ia(n,nb),li_blk_ib(n,nb),li_blk_ja(n,nb),li_blk_jb(n,nb),li_blk_ka(n,nb),li_blk_kb(n,nb)
!#endif

#ifdef DEBUG
    myunit = myid + 50
    write (filename, fmt = '(a,i2.2)') 'blk', myunit
    open (unit = myunit, file = trim(filename), form = 'formatted')
    write(unit = myunit, fmt = '(3g)')gi_ia(1),gi_blk_ia(1,1),li_blk_ia(1,1)
    write(unit = myunit, fmt = '(3g)')gi_ib(1),gi_blk_ib(1,1),li_blk_ib(1,1)
    write(unit = myunit, fmt = '(3g)')gi_ja(1),gi_blk_ja(1,1),li_blk_ja(1,1)
    write(unit = myunit, fmt = '(3g)')gi_jb(1),gi_blk_jb(1,1),li_blk_jb(1,1)
    write(unit = myunit, fmt = '(3g)')gi_ka(1),gi_blk_ka(1,1),li_blk_ka(1,1)
    write(unit = myunit, fmt = '(3g)')gi_kb(1),gi_blk_kb(1,1),li_blk_kb(1,1)
    close(unit = myunit)
#endif

       end do
       end do

    end if

    deallocate (gi_blk_ia, gi_blk_ib, &
                gi_blk_ja, gi_blk_jb, &
                gi_blk_ka, gi_blk_kb )

  end subroutine init_blanking

  ! --

  subroutine init_osg_interface ()

    ! init_interface = read and distribute interface information in overset grid
    ! 
    ! n_int(nz)         ; total no. of interface grid node in nz-th block
    ! nodes(1:3,i,nz)   ; n-th interface grid node's i,j,k
    ! hosts(1:4,i,nz)   ; 1:3 = i,j,k of host
    !                     4 = host block number
    ! coef(1:3,i,nz)    ; interpolation coefficient along three-dir 1,2,3

    ! read interface information on root
    ! bcast entire information to all processes
    ! store local piece of informationgrid (incl. ghost layers)
    !

    integer, dimension(:,:), allocatable :: nodest
    integer, dimension(:,:), allocatable :: hostst
    real (kind = rdf), dimension(:,:), allocatable :: coeft

    integer :: ii

    integer :: np
    integer :: n_int_total, n_int_max, n_int_max_check
    integer :: no1, no2, no3, no4
    integer :: ho1, ho2, ho3, ho4

    character (len = 255) :: filename
    integer :: myunit

    allocate (n_int(1:nzone))

    ! root reads grid from disk
    ! 
    if (myid == root+20000) then

       open (unit = 31, file = 'interface.dat', form = 'unformatted')

       ! read total no. of interface nodes in all zones
       read (unit = 31) n_int_total, n_int_max

       close (unit = 31)

    end if 

    call mpi_bcast (n_int_total, 1, MPI_INTEGER, root, mpi_comm_world, ierr)

    allocate (nodest(4,n_int_total), &
              hostst(4,n_int_total), &
               coeft(3,n_int_total) )
    nodest = 0
    hostst = 0
    coeft  = 0

    if (myid == root+20000) then

       open (unit = 31, file = 'interface.dat', form = 'unformatted')

       read (unit = 31) n_int_total, n_int_max

       ! get grid from file
       !
       ii = 0
       do nz = 1, nzone
          read (unit = 31) n_int(nz)
       do i = 1, n_int(nz)
          ii = ii + 1
         ! read (unit = 31) (nodest(j,ii), j = 1, 4)
         ! read (unit = 31) (hostst(j,ii), j = 1, 4)
         ! read (unit = 31) ( coeft(j,ii), j = 1, 3)
          read (unit = 31) nodest(1,ii),nodest(2,ii),nodest(3,ii),nodest(4,ii)
          read (unit = 31) hostst(1,ii),hostst(2,ii),hostst(3,ii),hostst(4,ii)
          read (unit = 31)  coeft(1,ii), coeft(2,ii), coeft(3,ii)
       end do
       end do

       close (unit = 31)

      !n_int_max =  maxval(n_int)

    end if

    !
    call mpi_bcast (n_int_max, 1, MPI_INTEGER, root, mpi_comm_world, ierr)
 
    ! distribute info
    ! 
    call mpi_bcast (nodest, 4*n_int_total, MPI_INTEGER, root, &
                    mpi_comm_world, ierr)
    call mpi_bcast (hostst, 4*n_int_total, MPI_INTEGER, root, &
                    mpi_comm_world, ierr)
    call mpi_bcast (coeft,  3*n_int_total, MPI_REAL,    root, &
                    mpi_comm_world, ierr)

    ! mpi-osg only
    ! read ijk range in each processor
    do np = 0, nproc - 1
    myunit = np !+ 10 !50
    write (filename, fmt = '(a,i3.3)') 'decofiles/gi.ijkz', myunit
    open (unit = myunit, file = trim(filename), form = 'formatted')
    read (unit = myunit, fmt = *) id_ijkz(1,np), id_ijkz(2,np), &
                                  id_ijkz(3,np), id_ijkz(4,np), &
                                  id_ijkz(5,np), id_ijkz(6,np), id_ijkz(7,np)
    close(unit = myunit)
    end do

    ! specify donor and donee precessors
    allocate ( donee(0:nproc-1), &
               donor(0:nproc-1) )

    donee = 0
    donor = 0

    allocate (hosts(1:4,1:n_int_max,0:nproc-1), &
              coef (1:3,1:n_int_max,0:nproc-1), &
              nodes(1:4,1:n_int_max,0:nproc-1) )
    hosts = 0
    coef  = 0
    nodes = 0

    do ii = 1, n_int_total

       ho1 = hostst(1,ii)
       ho2 = hostst(2,ii)
       ho3 = hostst(3,ii)
       ho4 = hostst(4,ii)
       no1 = nodest(1,ii)
       no2 = nodest(2,ii)
       no3 = nodest(3,ii)
       no4 = nodest(4,ii)

       if (myzone == ho4 .and. &
          ho3 >= gi_ka(1) .and. ho3 <= gi_kb(1) .and. &
          ho2 >= gi_ja(1) .and. ho2 <= gi_jb(1) .and. &
          ho1 >= gi_ia(1) .and. ho1 <= gi_ib(1) ) then

          do np = 0, nproc - 1
             if (id_ijkz(7,np) == no4 .and. & 
                 id_ijkz(5,np) <= no3 .and. no3 <= id_ijkz(6,np) .and. &
                 id_ijkz(3,np) <= no2 .and. no2 <= id_ijkz(4,np) .and. &
                 id_ijkz(1,np) <= no1 .and. no1 <= id_ijkz(2,np) ) then

                donor(np) = donor(np) + 1

!    li_ia(:) = 1 ; li_ib(:) = gi_ib(:) - gi_ia(:) + 1

                hosts(1,donor(np),np) = hostst(1,ii) - gi_ia(1) + 1
                hosts(2,donor(np),np) = hostst(2,ii) - gi_ja(1) + 1
                hosts(3,donor(np),np) = hostst(3,ii) - gi_ka(1) + 1
                hosts(4,donor(np),np) = hostst(4,ii)

               ! hosts(1:4,donor(np),np) = hostst(1:4,ii)
                coef (1:3,donor(np),np) = coeft (1:3,ii)
             end if
          end do
       end if

       if (myzone == no4 .and. &
          no3 >= gi_ka(1) .and. no3 <= gi_kb(1) .and. &
          no2 >= gi_ja(1) .and. no2 <= gi_jb(1) .and. &
          no1 >= gi_ia(1) .and. no1 <= gi_ib(1) ) then

          do np = 0, nproc - 1
             if (id_ijkz(7,np) == ho4 .and. &
                 id_ijkz(5,np) <= ho3 .and. ho3 <= id_ijkz(6,np) .and. &
                 id_ijkz(3,np) <= ho2 .and. ho2 <= id_ijkz(4,np) .and. &
                 id_ijkz(1,np) <= ho1 .and. ho1 <= id_ijkz(2,np) ) then

                donee(np) = donee(np) + 1

               ! nodes(1:4,donee(np),np) = nodest(1:4,ii)
                nodes(1,donee(np),np) = nodest(1,ii) - gi_ia(1) + 1
                nodes(2,donee(np),np) = nodest(2,ii) - gi_ja(1) + 1
                nodes(3,donee(np),np) = nodest(3,ii) - gi_ka(1) + 1
                nodes(4,donee(np),np) = nodest(4,ii)
             end if
          end do
       end if

    end do

#ifdef DEBUG
    myunit = myid + 50
    write (filename, fmt = '(a,i2.2)') 'donor_donee', myunit
    open (unit = myunit, file = trim(filename), form = 'formatted')
    write(unit = myunit, fmt = '(a15,3x,g)') 'myid', myid
    write(unit = myunit, fmt = '(a15,3x,g)') 'myzone', myzone
    write(unit = myunit, fmt = '(a15)') 'donor'
    write(unit = myunit, fmt = *) donor
    write(unit = myunit, fmt = '(a15)') 'donee'
    write(unit = myunit, fmt = *) donee
    close(unit = myunit)
#endif

    deallocate(nodest, hostst, coeft)

  end subroutine init_osg_interface

  ! --

  subroutine init_grid ()

    ! init_grid = read and distribute grid to all processes
    ! 
    ! read grid on root
    ! bcast entire grid to all processes
    ! store local piece of grid (incl. ghost layers) in xyz
    !
    integer :: sbuf_size
    integer, dimension(:,:,:,:,:), allocatable :: sbuf_idx

    real (kind = rdf), dimension (:), allocatable :: sbuf

    real (kind = rdf), dimension(:,:,:), allocatable :: xtmp
    real (kind = rdf), dimension(:,:,:), allocatable :: ztmp
    real (kind = rdf), dimension(:,:,:), allocatable :: ytmp
    real (kind = rdf), dimension(:,:,:), allocatable :: dtmp

    integer :: ia, ib
    integer :: ja, jb
    integer :: ka, kb

    integer :: ii
    integer :: jj
    integer :: kk

    integer :: nv

    integer :: imax
    integer :: jmax
    integer :: kmax

    integer :: myunit
    character (len = 256) :: filename

    integer :: num_vars
    integer :: nodes_total


    ! create send buffer index on all processes
    ! 
    
    if ( turbulence ) then
       num_vars = 4
    else
       num_vars = 3
    end if

    nodes_total = 0
    do nz = 1, nzone
       nodes_total = nodes_total + img(1,nz) * jmg(1,nz) * kmg(1,nz)
    end do

    sbuf_size = num_vars * nodes_total

    imax = maxval(img)
    jmax = maxval(jmg)
    kmax = maxval(kmg)

    allocate (sbuf(sbuf_size), &
              sbuf_idx(num_vars,imax,jmax,kmax,nzone))

    l = 0
    do nz = 1, nzone
    do k = 1, kmg(1,nz)
       do j = 1, jmg(1,nz)
          do i = 1, img(1,nz)
             do nv = 1, num_vars
                l = l + 1; sbuf_idx(nv,i,j,k,nz) = l
             end do
          end do
       end do
    end do
    end do
    
    ! root reads grid from disk
    ! 
     if (myid == root) then

       ! get grid from file
       !
       open  (unit = 41, file = 'grid', form='unformatted')
       do nz = 1, nzone

          allocate (xtmp(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz)), &
                    ytmp(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz)), &
                    ztmp(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz)), &
                    dtmp(1:img(1,nz),1:jmg(1,nz),1:kmg(1,nz))  )

          read  (unit = 41) (((xtmp(i,j,k),i=1,img(1,nz)), &
                                           j=1,jmg(1,nz)),k=1,kmg(1,nz))
          read  (unit = 41) (((ytmp(i,j,k),i=1,img(1,nz)), &
                                           j=1,jmg(1,nz)),k=1,kmg(1,nz))
          read  (unit = 41) (((ztmp(i,j,k),i=1,img(1,nz)), &
                                           j=1,jmg(1,nz)),k=1,kmg(1,nz))

          if ( turbulence ) &
          read  (unit = 41) (((dtmp(i,j,k),i=1,img(1,nz)), &
                                           j=1,jmg(1,nz)),k=1,kmg(1,nz))

          do k = 2, kmg(1,nz)
             do j = 2, jmg(1,nz)-1
                do i = 1, img(1,nz)
			if (dtmp(i,j,k).eq.0) print *,i,j,k
		end do
	end do
	end do





             ! pack send buf
          ! pack send buffer
          ! 
          if (nz == 1) l = 0
          do k = 1, kmg(1,nz)
             do j = 1, jmg(1,nz)
                do i = 1, img(1,nz)
                   l = l + 1; sbuf(l) = xtmp(i,j,k)
                   l = l + 1; sbuf(l) = ytmp(i,j,k)
                   l = l + 1; sbuf(l) = ztmp(i,j,k)
                   if ( turbulence ) &
                   l = l + 1; sbuf(l) = dtmp(i,j,k)
                end do
             end do
          end do
       
          ! probe planes part moved to input_file
          !

          deallocate (xtmp, ytmp, ztmp, dtmp)

       end do ! nz = 1, nzone
       close (unit = 41)

     end if
 
    ! distribute info
    ! 
     call mpi_bcast (sbuf, sbuf_size, MPI_REAL, root, mpi_comm_world, ierr)

    ! 3d decomp
    !
    kb = le_kb(1)
    ka = le_ka(1)
    jb = le_jb(1)
    ja = le_ja(1)
    ib = le_ib(1)
    ia = le_ia(1)

    if (mydown  == mpi_proc_null) ka = li_ka(1)
    if (myup    == mpi_proc_null) kb = li_kb(1)
    if (myleft  == mpi_proc_null) ja = li_ja(1)
    if (myright == mpi_proc_null) jb = li_jb(1)
    if (myback  == mpi_proc_null) ia = li_ia(1)
    if (myfront == mpi_proc_null) ib = li_ib(1)

    do k = ka, kb
       do j = ja, jb
          do i = ia, ib
             ii = i + gi_ia(1) - 1
             jj = j + gi_ja(1) - 1
             kk = k + gi_ka(1) - 1
             x(le_idx(i,j,k,1))  = sbuf(sbuf_idx(1,ii,jj,kk,myzone))
             y(le_idx(i,j,k,1))  = sbuf(sbuf_idx(2,ii,jj,kk,myzone))
             z(le_idx(i,j,k,1))  = sbuf(sbuf_idx(3,ii,jj,kk,myzone)) 
             if ( turbulence ) &
             wd(le_idx(i,j,k,1)) = sbuf(sbuf_idx(4,ii,jj,kk,myzone)) 
          end do
       end do
    end do

    deallocate(sbuf, sbuf_idx)

  end subroutine init_grid

  ! --

  subroutine init_solu ()

    integer ::  tmp_size
    integer :: sbuf_size

    integer, dimension(:,:,:,:,:), allocatable :: sbuf_idx

    real (kind = rdf), dimension (:), allocatable :: sbuf

    real (kind = rdf), dimension(:,:), allocatable :: qtmp
    real (kind = rdf), dimension(:,:), allocatable :: qntmp
    real (kind = rdf), dimension(:),   allocatable :: nutmp

    integer :: num_vars, me_start
    integer :: ia, ib
    integer :: ja, jb
    integer :: ka, kb
    integer :: kk, jj, ii

    integer :: nv
    integer :: l1, l2
    integer :: nodes_total

    integer :: myunit
    character (len=256) :: filename

    integer :: imax, jmax, kmax

    num_vars = 2 * me
    if ( turbulence ) num_vars = num_vars + 1

    nodes_total = 0
    do nz = 1, nzone
       nodes_total = nodes_total + img(1,nz) * jmg(1,nz) * kmg(1,nz)
    end do

    sbuf_size = num_vars * nodes_total

    imax = maxval(img)
    jmax = maxval(jmg)
    kmax = maxval(kmg)

    allocate (sbuf(sbuf_size), &
              sbuf_idx(num_vars,imax,jmax,kmax,nzone))

    l = 0
    do nz = 1, nzone
    do k = 1, kmg(1,nz)
       do j = 1, jmg(1,nz)
          do i = 1, img(1,nz)
             do nv = 1, num_vars
                l = l + 1; sbuf_idx(nv,i,j,k,nz) = l
             end do
          end do
       end do
    end do
    end do

    if (myid == root) then


       ! get solution from 'solu' file
       !
       if ( turbulence ) then
		  
		  me_start = me


          open  (unit = 42, file = 'solu', form = 'unformatted')
          do nz = 1, nzone

             ! temporary arrays on root
             ! 
             tmp_size = img(1,nz) * jmg(1,nz) * kmg(1,nz)
       
             allocate(qtmp(1:me,tmp_size), &
                     qntmp(1:me,tmp_size) )

             ! allocate temp array for turbulence viscosity
             ! 
             allocate(nutmp(tmp_size))

             do m = 1, me_start
                read  (unit = 42)  (qtmp(m,l), l = 1, tmp_size)
             end do
                read  (unit = 42)  (nutmp(l),  l = 1, tmp_size)
             do m = 1, me_start
                read  (unit = 42)  (qntmp(m,l),l = 1, tmp_size)
             end do

             ! pack send buf
             ! 
             l1 = 0
             if (nz == 1) l2 = 0
             do k = 1, kmg(1,nz)
                do j = 1, jmg(1,nz)
                   do i = 1, img(1,nz)
                      l1 = l1 + 1
                      do nv = 1, me
                         l2 = l2 + 1; sbuf(l2) =  qtmp(nv,l1)
                      end do
                      do nv = 1, me
                         l2 = l2 + 1; sbuf(l2) = qntmp(nv,l1)
                      end do
                      if ( turbulence ) &
                         l2 = l2 + 1; sbuf(l2) =   nutmp(l1)
                   end do
                end do
             end do

             deallocate (qtmp, qntmp, nutmp)

          end do
          close (unit = 42)

       else ! not turbulence

          open  (unit = 42, file = 'solu', form = 'unformatted')
          do nz = 1, nzone

             ! temporary arrays on root
             ! 
             tmp_size = img(1,nz) * jmg(1,nz) * kmg(1,nz)
       
             allocate(qtmp(1:me,tmp_size), &
                     qntmp(1:me,tmp_size) )

             read  (unit = 42)  ((qtmp(nv,l),l = 1, tmp_size), nv = 1, me)
             read  (unit = 42) ((qntmp(nv,l),l = 1, tmp_size), nv = 1, me)

             ! pack send buf
             ! 
             l1 = 0
             if (nz == 1) l2 = 0
             do k = 1, kmg(1,nz)
                do j = 1, jmg(1,nz)
                   do i = 1, img(1,nz)
                      l1 = l1 + 1
                      do nv = 1, me
                         l2 = l2 + 1; sbuf(l2) =  qtmp(nv,l1)
                      end do
                      do nv = 1, me
                         l2 = l2 + 1; sbuf(l2) = qntmp(nv,l1)
                      end do
                   end do
                end do
             end do

             deallocate (qtmp, qntmp)

          end do
          close (unit = 42)

       end if

    end if

    ! broadcast to all processes
    ! 
    call mpi_bcast (sbuf, sbuf_size, MPI_REAL, root, mpi_comm_world, ierr)

    ! upack subuf including ghost celss into
    ! local copies of q, qn, xnut
    !
    ia = le_ia(1)
    ja = le_ja(1)
    ka = le_ka(1)

    ib = le_ib(1)
    jb = le_jb(1)
    kb = le_kb(1)

    ! 3d decomp
    ! 
    if (myup    == mpi_proc_null) kb = li_kb(1)
    if (mydown  == mpi_proc_null) ka = li_ka(1)
    if (myright == mpi_proc_null) jb = li_jb(1)
    if (myleft  == mpi_proc_null) ja = li_ja(1)
    if (myfront == mpi_proc_null) ib = li_ib(1)
    if (myback  == mpi_proc_null) ia = li_ia(1)

    do k = ka, kb
       do j = ja, jb
          do i = ia, ib
             ii = i + gi_ia(1) - 1
             jj = j + gi_ja(1) - 1
             kk = k + gi_ka(1) - 1

            do nv = 1, me
            q(nv,le_idx(i,j,k,1)) = sbuf(sbuf_idx(nv,ii,jj,kk,myzone))
            end do
            do nv = 1, me
            qn(nv,le_idx(i,j,k,1))= sbuf(sbuf_idx(nv+me,ii,jj,kk,myzone))
            end do
            if ( turbulence ) &
            xnut(le_idx(i,j,k,1)) = sbuf(sbuf_idx(num_vars,ii,jj,kk,myzone))
          end do
       end do
    end do

    deallocate(sbuf, sbuf_idx)
    
#ifdef DEBUG

    allocate (qtmp(me,1:li_idx_b(1)), &
             qntmp(me,1:li_idx_b(1)) )

    ! check q & qn
    ! 
    n = 1
    do k = li_ka(n), li_kb(n)
       do j = li_ja(n), li_jb(n)
          do i = li_ia(n), li_ib(n)
             qtmp(:,li_idx(i,j,k,n)) = q(:,le_idx(i,j,k,n))
            qntmp(:,li_idx(i,j,k,n)) =qn(:,le_idx(i,j,k,n))
          end do
       end do
    end do

    call checksum_2d_par('q  :',  qtmp)
    call checksum_2d_par('qn :', qntmp)

    deallocate (qtmp, qntmp)

#endif

  end subroutine init_solu

  ! --   
  
  subroutine ck_variables ()

    real (kind = rdf), dimension (:,:), allocatable :: ct
    real (kind = rdf), dimension (:,:), allocatable :: et
    real (kind = rdf), dimension (:,:), allocatable :: zett, qt
    real (kind = rdf), dimension (:), allocatable   :: at, xt, yt, zt
    real (kind = rdf), dimension (:), allocatable   :: dtt

    real (kind = rdf), dimension (:,:,:), allocatable :: mait, n1it, n2it, mct
    real (kind = rdf), dimension (:,:), allocatable :: sprt

    integer :: la, lb

    character (len = 255) :: filename

    integer :: myunit
    integer :: myunit1
    integer :: myunit2
    integer :: myunit3
    integer :: myunit4

    allocate (ct(1:3,1:li_idx_b(ng)), & 
              et(1:3,1:li_idx_b(ng)), &
            zett(1:3,1:li_idx_b(ng)), &
                  at(1:li_idx_b(ng)), &
                 dtt(1:li_idx_b(ng)), &
            sprt(1:3,1:li_idx_b(ng)), &
        mait(1:4,1:4,1:li_idx_b(ng)), &
        n1it(1:4,1:4,1:li_idx_b(ng)), &
        n2it(1:4,1:4,1:li_idx_b(ng)), &
         mct(1:4,1:4,1:li_idx_b(ng)), &
                  xt(1:li_idx_b(ng)), &
                  yt(1:li_idx_b(ng)), &
                  zt(1:li_idx_b(ng)), &
              qt(1:me,1:li_idx_b(ng)) )

    do n = 1, ng
       do k = li_ka(n), li_kb(n)
          do j = li_ja(n), li_jb(n)
             do i = li_ia(n), li_ib(n)
                la = li_idx(i,j,k,n)
                lb = le_idx(i,j,k,n)
!                 ct(:,la) = csi(:,lb)
!                 et(:,la) = eta(:,lb)
                 zett(:,la) = zet(:,lb)
                xt(la) = x(lb)
                yt(la) = y(lb)
                zt(la) = z(lb)
                qt(:,la) = q(:,lb)
                at(la) = aj(lb)
                dtt(la) = xnut(lb)
!               sprt(:,la) = spr(:,lb)
!               mait(:,:,la) =  mai(:,:,lb)
!               n1it(:,:,la) =  n1i(:,:,lb)
!               n2it(:,:,la) =  n2i(:,:,lb)
!                mct(:,:,la) =   mc(:,:,lb)
             end do
          end do
       end do
    end do
    
!     call checksum_2d ('csi 1:', ct(1:1,:))
!     call checksum_2d ('csi 2:', ct(2:2,:))
!     call checksum_2d ('csi 3:', ct(3:3,:))
!     call checksum_2d ('eta 1:', et(1:1,:))
!     call checksum_2d ('eta 2:', et(2:2,:))
!     call checksum_2d ('eta 3:', et(3:3,:))
!     call checksum_2d ('zet 1:', zt(1:1,:))
!     call checksum_2d ('zet 2:', zt(2:2,:))
!     call checksum_2d ('zet 3:', zt(3:3,:))
!     call checksum_1d ('aj :', at)

!     call checksum_2d_par ('csi', ct)
!     call checksum_2d_par ('eta', et)
!     call checksum_2d_par ('zet', zt)

    call checksum_1d_par ('init.x', xt)    
    call checksum_1d_par ('init.y', yt)    
    call checksum_1d_par ('init.z', zt)    
    call checksum_2d_par ('init.zet', zett)
    call checksum_1d_par ('init.aj', at)    
    call checksum_2d_par ('init.q', qt)
    !call checksum_1d_par ('init.xnut', dtt)
    !call checksum_1d_par ('init.dtev', dtevt)

!     call checksum_1d_par ('dtau', dtt)    
!     call checksum_2d_par ('spr', sprt)
!     call checksum_3d_par ('mai', mait)
!     call checksum_3d_par ('n1i', n1it)
!     call checksum_3d_par ('n2i', n2it)
!     call checksum_3d_par ('mc', mct)

!     ! check via output file
!     !
!     myunit1 = myid + 60
!     myunit2 = myid + 70
!     myunit3 = myid + 80
!     myunit4 = myid + 90

!     write (filename, fmt = '(a,i2.2)') 'ck_cs', myunit1
!     open (unit = myunit1, file = trim(filename), form = 'formatted')
!     write (filename, fmt = '(a,i2.2)') 'ck_et', myunit2
!     open (unit = myunit2, file = trim(filename), form = 'formatted')
!     write (filename, fmt = '(a,i2.2)') 'ck_zt', myunit3
!     open (unit = myunit3, file = trim(filename), form = 'formatted')
!     write (filename, fmt = '(a,i2.2)') 'ck_aj', myunit4
!     open (unit = myunit4, file = trim(filename), form = 'formatted')


!     do n = 1, 1
!        do k = gi_ka(n), gi_kb(n)
!           do j = gi_ja(n), gi_jb(n)
!              do i = gi_ia(n), gi_ib(n)
!                 lb = le_idx(i,j,k-gi_ka(n)+1,n)
!                 write (unit = myunit1,'(4(i2,1x),3(g15.7,1x))') n,i,j,k,csi(:,lb)
!                 write (unit = myunit2,'(4(i2,1x),3(g15.7,1x))') n,i,j,k,eta(:,lb)
!                 write (unit = myunit3,'(4(i2,1x),3(g15.7,1x))') n,i,j,k,zet(:,lb)
!                 write (unit = myunit4,'(4(i2,1x),3(g15.7,1x))') n,i,j,k, aj(lb)
!              end do
!           end do
!        end do
!     end do

!     close (unit = myunit1)
!     close (unit = myunit2)
!     close (unit = myunit3)
!     close (unit = myunit4)

    deallocate (ct, et, zett, at, dtt, sprt, mait, n1it, n2it, mct, qt, xt, yt, zt)


#ifdef DEBUG

    ! check via output file
    !
    myunit = myid + 50

    write (filename, fmt = '(a,i2.2)') 'ce_input', myunit

    open (unit = myunit, file = trim(filename), form = 'formatted')

    nz = myzone
    write (unit = myunit, fmt = '(a15,3x,g)') 'ns', ns
    write (unit = myunit, fmt = '(a15,3x,g)') 'icycle', icycle
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'img', img(:,nz)
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'jmg', jmg(:,nz)
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'kmg', kmg(:,nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'nt2', nt2
    write (unit = myunit, fmt = '(a15,3x,g)') 'delti', delti
    write (unit = myunit, fmt = '(a15,3x,g)') 'e_source', e_source
    write (unit = myunit, fmt = '(a15,3x,g)') 'it_min', it_min
    write (unit = myunit, fmt = '(a15,3x,g)') 'er_min', er_min
    write (unit = myunit, fmt = '(a15,3x,g)') 'eo_min', eo_min
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'iter', iter
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'itm', itm
    write (unit = myunit, fmt = '(a15,3x,g)') 'ren', ren
    write (unit = myunit, fmt = '(a15,3x,g)') 'cfl1', cfl1(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'vnn1', vnn1(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'beta', beta
    write (unit = myunit, fmt = '(a15,3x,g)') 'icn', icn
    write (unit = myunit, fmt = '(a15,3x,g)') 'icnw', icnw
    write (unit = myunit, fmt = '(a15,3x,g)') 'irk', irk
    write (unit = myunit, fmt = '(a15,3x,4(g,1x))') 'alfa', alfa
    write (unit = myunit, fmt = '(a15,3x,4(g,1x))') 'eps', eps
    write (unit = myunit, fmt = '(a15,3x,g)') 'itk', itk
    write (unit = myunit, fmt = '(a15,3x,g)') 'cfl2', cfl2(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'vnn2', vnn2(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'cdes', cdes
    write (unit = myunit, fmt = '(a15,3x,g)') 'ckin', ckin
    write (unit = myunit, fmt = '(a15,3x,g)') 'cein', cein
    write (unit = myunit, fmt = '(a15,3x,g)') 'cvin', cvin
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'ep', ep(:,:,nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'pdiss_coef', pdiss_coef
    write (unit = myunit, fmt = '(a15,3x,g)') 'icrs', icrs(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'jcrs', jcrs(nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'kcrs', kcrs(nz)
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'btype', btype(:,nz)
    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'bdir', bdir(:,nz)
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'local_ifix', local_ifix
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'local_jfix', local_jfix
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'local_kfix', local_kfix
    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'pfix_proc', pfix_proc

    write (unit = myunit, fmt = '(a15,3x,3(g,1x))') 'dims', dims(:,nz)
    write (unit = myunit, fmt = '(a15,3x,g)') 'nproc_nz', nproc_nz(nz)

    write (unit = myunit, fmt = '(a15,3x,6(g,1x))') 'blanking', blanking(:,:,nz)

    close (unit = myunit)

#endif

  end subroutine ck_variables
  
  ! --

  subroutine init_outlet_bc ()

    integer :: l0, l1, l2
    integer :: b

    integer :: imax
    integer :: jmax
    integer :: kmax
    integer :: ijkmax

    imax = le_ib(1)-le_ia(1)+1
    jmax = le_jb(1)-le_ja(1)+1
    kmax = le_kb(1)-le_ka(1)+1
    ijkmax = max(imax,jmax,kmax)

    allocate (rat(ijkmax,ijkmax,6))

    !if(myid ==0) print *, ijkmax
    ! csi direction
    b = 1
    if (myback == mpi_proc_null) then
    if (btype(b,myzone) == 4 .or. btype(b,myzone) == 5) then
       do k = le_ka(1), le_kb(1)
       do j = le_ja(1), le_jb(1)
          l0 = le_idx(1,j,k,1)
          l1 = le_idx(2,j,k,1)
          l2 = le_idx(3,j,k,1)
!          if (bdir(b,myzone) == 1) &
!              rat(j,k,b) = (x(l1) - x(l0)) / (x(l2) - x(l1))
          if (bdir(b,myzone) == 1) then !modificado para que no se vaya fuera de rango el rat
              if(j .ge. 1 .and. k.ge.1) &
               rat(j,k,b) = (x(l1) - x(l0)) / (x(l2) - x(l1))

              !if(myid == 0) print *,rat(j,k,b), j, k
          end if
          if (bdir(b,myzone) == 2) &
              rat(j,k,b) = (y(l1) - y(l0)) / (y(l2) - y(l1))
          if (bdir(b,myzone) == 3) &
              rat(j,k,b) = (z(l1) - z(l0)) / (z(l2) - z(l1))
       end do
       end do
    end if
    end if

    b = 2
    if (myfront == mpi_proc_null) then
    if (btype(b,myzone) == 4 .or. btype(b,myzone) == 5) then
       do k = le_ka(1), le_kb(1)
       do j = le_ja(1), le_jb(1)
          l0 = le_idx(img(1,myzone)  ,j,k,1)
          l1 = le_idx(img(1,myzone)-1,j,k,1)
          l2 = le_idx(img(1,myzone)-2,j,k,1)
          if (bdir(b,myzone) == 1) &
              rat(j,k,b) = (x(l0) - x(l1)) / (x(l1) - x(l2))
          if (bdir(b,myzone) == 2) &
              rat(j,k,b) = (y(l0) - y(l1)) / (y(l1) - y(l2))
          if (bdir(b,myzone) == 3) &
              rat(j,k,b) = (z(l0) - z(l1)) / (z(l1) - z(l2))
       end do
       end do
    end if
    end if

    ! eta direction
    b = 3
    if (myleft == mpi_proc_null) then
    if (btype(b,myzone) == 4 .or. btype(b,myzone) == 5) then
       do k = le_ka(1), le_kb(1)
       do i = le_ia(1), le_ib(1)
          l0 = le_idx(i,1,k,1)
          l1 = le_idx(i,2,k,1)
          l2 = le_idx(i,3,k,1)
          if (bdir(b,myzone) == 1) &
              rat(i,k,b) = (x(l1) - x(l0)) / (x(l2) - x(l1))
          if (bdir(b,myzone) == 2) &
              rat(i,k,b) = (y(l1) - y(l0)) / (y(l2) - y(l1))
          if (bdir(b,myzone) == 3) &
              rat(i,k,b) = (z(l1) - z(l0)) / (z(l2) - z(l1))
       end do
       end do
    end if
    end if

    b = 4
    if (myright == mpi_proc_null) then
    if (btype(b,myzone) == 4 .or. btype(b,myzone) == 5) then
       do k = le_ka(1), le_kb(1)
       do i = le_ia(1), le_ib(1)
          l0 = le_idx(i,jmg(1,myzone)  ,k,1)
          l1 = le_idx(i,jmg(1,myzone)-1,k,1)
          l2 = le_idx(i,jmg(1,myzone)-2,k,1)
          if (bdir(b,myzone) == 1) &
              rat(i,k,b) = (x(l0) - x(l1)) / (x(l1) - x(l2))
          if (bdir(b,myzone) == 2) &
              rat(i,k,b) = (y(l0) - y(l1)) / (y(l1) - y(l2))
          if (bdir(b,myzone) == 3) &
              rat(i,k,b) = (z(l0) - z(l1)) / (z(l1) - z(l2))
       end do
       end do
    end if
    end if

    ! zet direction
    b = 5
    if (mydown == mpi_proc_null) then
    if (btype(b,myzone) == 4 .and. btype(b,myzone) == 5) then
       do j = le_ja(1), le_jb(1)
       do i = le_ia(1), le_ib(1)
          l0 = le_idx(i,j,1,1)
          l1 = le_idx(i,j,2,1)
          l2 = le_idx(i,j,3,1)
          if (bdir(b,myzone) == 1) &
              rat(i,j,b) = (x(l1) - x(l0)) / (x(l2) - x(l1))
          if (bdir(b,myzone) == 2) &
              rat(i,j,b) = (y(l1) - y(l0)) / (y(l2) - y(l1))
          if (bdir(b,myzone) == 3) &
              rat(i,j,b) = (z(l1) - z(l0)) / (z(l2) - z(l1))
       end do
       end do
    end if
    end if

    b = 6
    if (myup == mpi_proc_null) then
    if (btype(b,myzone) == 4 .and. btype(b,myzone) == 5) then
       do j = le_ja(1), le_jb(1)
       do i = le_ia(1), le_ib(1)
          l0 = le_idx(i,j,kmg(1,myzone)  ,1)
          l1 = le_idx(i,j,kmg(1,myzone)-1,1)
          l2 = le_idx(i,j,kmg(1,myzone)-2,1)
          if (bdir(b,myzone) == 1) &
              rat(i,j,b) = (x(l0) - x(l1)) / (x(l1) - x(l2))
          if (bdir(b,myzone) == 2) &
              rat(i,j,b) = (y(l0) - y(l1)) / (y(l1) - y(l2))
          if (bdir(b,myzone) == 3) &
              rat(i,j,b) = (z(l0) - z(l1)) / (z(l1) - z(l2))
       end do
       end do
    end if
    end if

  end subroutine init_outlet_bc

  ! --

!  subroutine init_mg_old ()

!    integer :: imax, jmax, kmax

    ! keep so that everything else will compile
    !
!    imax = maxval(img)
!    jmax = maxval(jmg)
!    kmax = maxval(kmg)

!    allocate ( ls(ng,nzone), &
!               le(ng,nzone), &
!              icg(ng,nzone), &
!              ieg(ng,nzone), &
!              izg(ng,nzone), &
!              ln(imax,jmax,kmax,ng,nzone) )

!    do nz = 1, nzone
!       if(nz == 1) then
!         ls(1,nz) = 1
!         le(1,nz) = img(1,nz)*jmg(1,nz)*kmg(1,nz)
!       else
!         ls(1,nz) = le(ng,nz-1) + 1
!         le(1,nz) = le(ng,nz-1) + img(1,nz)*jmg(1,nz)*kmg(1,nz)
!       end if
!       do n = 2, ng
!          ls(n,nz) = le(n-1,nz) + 1
!          le(n,nz) = ls(n,nz) + img(n,nz)*jmg(n,nz)*kmg(n,nz) - 1
!       end do
!    end do
!
!    do nz = 1, nzone
!    do n = 1, ng
!       icg(n,nz) =  1
!       ieg(n,nz) = img(n,nz)
!       izg(n,nz) = img(n,nz) * jmg(n,nz)
!    end do
!    end do

!    do nz = 1, nzone
!    do n = 1, ng
!       do k = 1, kmg(n,nz)
!          do j = 1, jmg(n,nz)
!             do i = 1, img(n,nz)
!                ln(i,j,k,n,nz) = (k - 1) * img(n,nz) * jmg(n,nz) + &
!                                 (j - 1) * img(n,nz) + &
!                                  i + &
!                                  ls(n,nz) - 1
!             end do
!          end do
!       end do
!    end do
!    end do

!  end subroutine init_mg_old

  ! --

  subroutine exchng3_mg_1d (n, var)
    
    ! current grid level
    ! 
    integer :: n

    ! limited to a 1d decomposition
    ! 
    integer :: buf_kmx           ! total size of buf
    integer :: buf_jmx           ! total size of buf
    integer :: buf_imx           ! total size of buf
    
    ! variable to exchange
    ! 
    real (kind = rdf), dimension(le_idx_a(1):le_idx_b(ng)), intent(inout) :: var

    ! send and receive buffers
    ! 
    real (kind = rdf), dimension (:), allocatable :: bufu
    real (kind = rdf), dimension (:), allocatable :: bufd
    real (kind = rdf), dimension (:), allocatable :: ubuf
    real (kind = rdf), dimension (:), allocatable :: dbuf

    real (kind = rdf), dimension (:), allocatable :: buff
    real (kind = rdf), dimension (:), allocatable :: bufb
    real (kind = rdf), dimension (:), allocatable :: fbuf
    real (kind = rdf), dimension (:), allocatable :: bbuf

    real (kind = rdf), dimension (:), allocatable :: bufl
    real (kind = rdf), dimension (:), allocatable :: bufr
    real (kind = rdf), dimension (:), allocatable :: lbuf
    real (kind = rdf), dimension (:), allocatable :: rbuf

    integer :: ka
    integer :: kb

    integer :: ja
    integer :: jb

    integer :: ia
    integer :: ib

    ! set size (based on ghost layers)
    ! 
    buf_kmx = (le_ix(n) * le_jx(n) * kgp(n))
    buf_jmx = (le_ix(n) * le_kx(n) * jgp(n))
    buf_imx = (le_jx(n) * le_kx(n) * igp(n))
    
    allocate(bufu(buf_kmx), bufd(buf_kmx), &
             ubuf(buf_kmx), dbuf(buf_kmx) )

    allocate(bufl(buf_jmx), bufr(buf_jmx), &
             lbuf(buf_jmx), rbuf(buf_jmx) )

    allocate(bufb(buf_imx), buff(buf_imx), &
             bbuf(buf_imx), fbuf(buf_imx) )

    ! pack data going to process below
    !
    if (mydown /= mpi_proc_null) then
       
       ka = li_ka(n)
       kb = li_ka(n) + kgp(n) - 1
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; bufd(l) = var(le_idx(i,j,k,n))
             end do
          end do
       end do
    end if

    ! pack data going to process above
    ! 
    if (myup /= mpi_proc_null) then
       
       ka = li_kb(n) - kgp(n) + 1
       kb = li_kb(n)
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; bufu(l) = var(le_idx(i,j,k,n))
             end do
          end do
       end do
    end if

    ! pack data going to process left
    !
    if (myleft /= mpi_proc_null) then
       
       ja = li_ja(n)
       jb = li_ja(n) + jgp(n) - 1
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; bufl(l) = var(le_idx(i,j,k,n))
             end do
          end do
       end do
    end if

    ! pack data going to process right
    ! 
    if (myright /= mpi_proc_null) then
       
       ja = li_jb(n) - jgp(n) + 1
       jb = li_jb(n)
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; bufr(l) = var(le_idx(i,j,k,n))
             end do
          end do
       end do
    end if

    ! pack data going to process back
    !
    if (myback /= mpi_proc_null) then
       
       ia = li_ia(n)
       ib = li_ia(n) + igp(n) - 1
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                l = l + 1; bufb(l) = var(le_idx(i,j,k,n))
             end do
          end do
       end do
    end if

    ! pack data going to process front
    ! 
    if (myfront /= mpi_proc_null) then
       
       ia = li_ib(n) - igp(n) + 1
       ib = li_ib(n)
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                l = l + 1; buff(l) = var(le_idx(i,j,k,n))
             end do
          end do
       end do
    end if

    
    ! sync before communzcation
    ! 
    call mpi_barrier(mpi_comm_world, ierr)

    ! down & up
    ! send to mydown & receive from myup
    ! 
    call mpi_sendrecv (bufd, buf_kmx, mpi_real, mydown, 99,& 
                     & ubuf, buf_kmx, mpi_real, myup, 99,&
                     & mpi_comm_world, status, ierr)

    ! send to myup and receive from mydown
    ! 
    call mpi_sendrecv (bufu, buf_kmx, mpi_real, myup, 100,&
                     & dbuf, buf_kmx, mpi_real, mydown, 100,&
                     & mpi_comm_world, status, ierr)

    ! left & right
    ! send to myleft & receive from myright
    ! 
    call mpi_sendrecv (bufl, buf_jmx, mpi_real, myleft, 199,&
                     & rbuf, buf_jmx, mpi_real, myright, 199,&
                     & mpi_comm_world, status, ierr)

    ! send to myright and receive from myleft
    ! 
    call mpi_sendrecv (bufr, buf_jmx, mpi_real, myright, 200,&
                     & lbuf, buf_jmx, mpi_real, myleft, 200,&
                     & mpi_comm_world, status, ierr)

    ! back & front
    ! send to myback & receive from myfront
    ! 
    call mpi_sendrecv (bufb, buf_imx, mpi_real, myback, 299,&
                     & fbuf, buf_imx, mpi_real, myfront, 299,&
                     & mpi_comm_world, status, ierr)

    ! send to myfront and receive from myback
    ! 
    call mpi_sendrecv (buff, buf_imx, mpi_real, myfront, 300,&
                     & bbuf, buf_imx, mpi_real, myback, 300,&
                     & mpi_comm_world, status, ierr)


    ! unpack arrray from mydown
    ! 
     if (mydown /= mpi_proc_null) then 
       ka = le_ka(n) 
       kb = li_ka(n) - 1 
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; var(le_idx(i,j,k,n)) = dbuf(l)
             end do
          end do
       end do
     end if
 
    ! unpack array from myup
    ! 
     if (myup /= mpi_proc_null) then
       ka = li_kb(n) + 1
       kb = le_kb(n)
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; var(le_idx(i,j,k,n)) = ubuf(l)
             end do
          end do
       end do
     end if

      deallocate (bufu, bufd, ubuf, dbuf)

    ! unpack arrray from myleft
    ! 
     if (myleft /= mpi_proc_null) then 
       ja = le_ja(n) 
       jb = li_ja(n) - 1 
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; var(le_idx(i,j,k,n)) = lbuf(l)
             end do
          end do
       end do
     end if
 
    ! unpack array from myright
    ! 
     if (myright /= mpi_proc_null) then
       ja = li_jb(n) + 1
       jb = le_jb(n)
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                l = l + 1; var(le_idx(i,j,k,n)) = rbuf(l)
             end do
          end do
       end do
     end if

      deallocate (bufl, bufr, lbuf, rbuf)

    ! unpack arrray from myback
    ! 
     if (myback /= mpi_proc_null) then 
       ia = le_ia(n) 
       ib = li_ia(n) - 1 
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                l = l + 1; var(le_idx(i,j,k,n)) = bbuf(l)
             end do
          end do
       end do
     end if
 
    ! unpack array from myfront
    ! 
     if (myfront /= mpi_proc_null) then
       ia = li_ib(n) + 1
       ib = le_ib(n)
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                l = l + 1; var(le_idx(i,j,k,n)) = fbuf(l)
             end do
          end do
       end do
     end if

      deallocate (bufb, buff, bbuf, fbuf)
       
  end subroutine exchng3_mg_1d

  !--

  subroutine exchng3_mg_2d (n, var)
    
    ! current grid level
    ! 
    integer :: n

    ! limited to a 1d decomposition
    ! 
    integer :: buf_kmx           ! total size of buf
    integer :: buf_jmx           ! total size of buf
    integer :: buf_imx           ! total size of buf
    
    ! variable to exchange
    ! 
    real (kind = rdf), dimension (:,:), intent(inout) :: var
   !real (kind = rdf), dimension(le_idx_a(1):le_idx_b(ng)), intent(inout) :: var

    ! send and receive buffers
    ! 
    real (kind = rdf), dimension (:), allocatable :: bufu
    real (kind = rdf), dimension (:), allocatable :: bufd
    real (kind = rdf), dimension (:), allocatable :: ubuf
    real (kind = rdf), dimension (:), allocatable :: dbuf

    real (kind = rdf), dimension (:), allocatable :: buff
    real (kind = rdf), dimension (:), allocatable :: bufb
    real (kind = rdf), dimension (:), allocatable :: fbuf
    real (kind = rdf), dimension (:), allocatable :: bbuf

    real (kind = rdf), dimension (:), allocatable :: bufl
    real (kind = rdf), dimension (:), allocatable :: bufr
    real (kind = rdf), dimension (:), allocatable :: lbuf
    real (kind = rdf), dimension (:), allocatable :: rbuf

    integer :: ka
    integer :: kb

    integer :: ja
    integer :: jb

    integer :: ia
    integer :: ib

    integer :: ma
    integer :: mb
    integer :: mx

    integer :: la
    integer :: lb

    ma = lbound(var,1)
    mb = ubound(var,1)
    mx = mb - ma + 1

    ! check bounds of var,2
    !
    if (lbound(var,2) /= le_idx_a(1)) print *, 'lbound problem in exchange_m_2d'
    if (ubound(var,2) /= le_idx_b(ng))print *, 'ubound problem in exchange_m_2d'

    ! set size (based on ghost layers)
    ! 
    buf_kmx = (mx * le_ix(n) * le_jx(n) * kgp(n))
    buf_jmx = (mx * le_ix(n) * le_kx(n) * jgp(n))
    buf_imx = (mx * le_jx(n) * le_kx(n) * igp(n))
    
    allocate(bufu(buf_kmx), bufd(buf_kmx), &
             ubuf(buf_kmx), dbuf(buf_kmx) )

    allocate(bufl(buf_jmx), bufr(buf_jmx), &
             lbuf(buf_jmx), rbuf(buf_jmx) )

    allocate(bufb(buf_imx), buff(buf_imx), &
             bbuf(buf_imx), fbuf(buf_imx) )


    ! pack data going to process below
    !
    if (mydown /= mpi_proc_null) then
       
       ka = li_ka(n)
       kb = li_ka(n) + kgp(n) - 1
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; bufd(l) = var(m, le_idx(i,j,k,n))
                end do
             end do
          end do
       end do
    end if

    ! pack data going to process above
    ! 
    if (myup /= mpi_proc_null) then
       
       ka = li_kb(n) - kgp(n) + 1
       kb = li_kb(n)
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; bufu(l) = var(m, le_idx(i,j,k,n))
                end do
             end do
          end do
       end do
    end if

    ! pack data going to process left
    !
    if (myleft /= mpi_proc_null) then
       
       ja = li_ja(n)
       jb = li_ja(n) + jgp(n) - 1
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; bufl(l) = var(m, le_idx(i,j,k,n))
                end do
             end do
          end do
       end do
    end if

    ! pack data going to process right
    ! 
    if (myright /= mpi_proc_null) then
       
       ja = li_jb(n) - jgp(n) + 1
       jb = li_jb(n)
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; bufr(l) = var(m, le_idx(i,j,k,n))
                end do
             end do
          end do
       end do
    end if

    ! pack data going to process back
    !
    if (myback /= mpi_proc_null) then
       
       ia = li_ia(n)
       ib = li_ia(n) + igp(n) - 1
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                do m = ma, mb
                   l = l + 1; bufb(l) = var(m, le_idx(i,j,k,n))
                end do
             end do
          end do
       end do
    end if

    ! pack data going to process front
    ! 
    if (myfront /= mpi_proc_null) then
       
       ia = li_ib(n) - igp(n) + 1
       ib = li_ib(n)
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                do m = ma, mb
                   l = l + 1; buff(l) = var(m, le_idx(i,j,k,n))
                end do
             end do
          end do
       end do
    end if

    ! sync before communzcation
    ! 
    call mpi_barrier(mpi_comm_world, ierr) !???

    ! send to mydown & receive from myup
    ! 
    call mpi_sendrecv (bufd, buf_kmx, mpi_real, mydown, 99,&
                     & ubuf, buf_kmx, mpi_real, myup, 99,&
                     & mpi_comm_world, status, ierr)

    ! send to myup and receive from mydown
    ! 
    call mpi_sendrecv (bufu, buf_kmx, mpi_real, myup, 100,&
                     & dbuf, buf_kmx, mpi_real, mydown, 100,&
                     & mpi_comm_world, status, ierr)

    ! send to myleft & receive from myright
    ! 
    call mpi_sendrecv (bufl, buf_jmx, mpi_real, myleft, 199,&
                     & rbuf, buf_jmx, mpi_real, myright, 199,&
                     & mpi_comm_world, status, ierr)

    ! send to myup and receive from mydown
    ! 
    call mpi_sendrecv (bufr, buf_jmx, mpi_real, myright, 200,&
                     & lbuf, buf_jmx, mpi_real, myleft, 200,&
                     & mpi_comm_world, status, ierr)

    ! send to myback & receive from myfront
    ! 
    call mpi_sendrecv (bufb, buf_imx, mpi_real, myback, 299,&
                     & fbuf, buf_imx, mpi_real, myfront, 299,&
                     & mpi_comm_world, status, ierr)

    ! send to myfront and receive from myback
    ! 
    call mpi_sendrecv (buff, buf_imx, mpi_real, myfront, 300,&
                     & bbuf, buf_imx, mpi_real, myback, 300,&
                     & mpi_comm_world, status, ierr)


    ! unpack arrray from mydown
    ! 
!     if (mydown /= mpi_proc_null) then 
       ka = le_ka(n) 
       kb = li_ka(n) - 1 
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; var(m, le_idx(i,j,k,n)) = dbuf(l)
                end do
             end do
          end do
       end do
!     end if
    
    ! unpack array from myup
    ! 
!     if (myup /= mpi_proc_null) then
       ka = li_kb(n) + 1
       kb = le_kb(n)
       l = 0
       do k = ka, kb
          do j = le_ja(n), le_jb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; var(m, le_idx(i,j,k,n)) = ubuf(l)
                end do
             end do
          end do
       end do
!     end if

      deallocate (ubuf, dbuf, bufu, bufd)

    ! unpack arrray from myleft
    ! 
!     if (myleft /= mpi_proc_null) then 
       ja = le_ja(n) 
       jb = li_ja(n) - 1 
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; var(m, le_idx(i,j,k,n)) = lbuf(l)
                end do
             end do
          end do
       end do
!     end if

    ! unpack array from myright
    ! 
!     if (myright /= mpi_proc_null) then
       ja = li_jb(n) + 1
       jb = le_jb(n)
       l = 0
       do j = ja, jb
          do k = le_ka(n), le_kb(n)
             do i = le_ia(n), le_ib(n)
                do m = ma, mb
                   l = l + 1; var(m, le_idx(i,j,k,n)) = rbuf(l)
                end do
             end do
          end do
       end do
!     end if

      deallocate (rbuf, lbuf, bufr, bufl)

    ! unpack arrray from myback
    ! 
!     if (myback /= mpi_proc_null) then 
       ia = le_ia(n) 
       ib = li_ia(n) - 1 
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                do m = ma, mb
                   l = l + 1; var(m, le_idx(i,j,k,n)) = bbuf(l)
                end do
             end do
          end do
       end do
!     end if
    
    ! unpack array from myfront
    ! 
!     if (myfront /= mpi_proc_null) then
       ia = li_ib(n) + 1
       ib = le_ib(n)
       l = 0
       do i = ia, ib
          do k = le_ka(n), le_kb(n)
             do j = le_ja(n), le_jb(n)
                do m = ma, mb
                   l = l + 1; var(m, le_idx(i,j,k,n)) = fbuf(l)
                end do
             end do
          end do
       end do
!     end if

      deallocate (fbuf, bbuf, buff, bufb)

 
  end subroutine exchng3_mg_2d

  ! --
  !
  subroutine decomp1d_cg (crs, a, b)

    implicit none

    integer, intent(in) :: crs
    integer, dimension(ng), intent(inout) :: a
    integer, dimension(ng), intent(inout) :: b
    
    integer :: n

    do n = 2, ng
       if (crs == 0) then
          a(n) = a(1)
          b(n) = b(1)
       else
          if (mod(a(n-1),2) == 0) then
             a(n) = (a(n-1) + 2) / 2
          else
             a(n) = (a(n-1) + 1) / 2
          end if
          if (mod(b(n-1),2) == 0) then
             b(n) = b(n-1) / 2
          else
             b(n) = (b(n-1) + 1) / 2
          end if
       end if
    end do

  end subroutine decomp1d_cg


!================================================================================
! subroutines search and vector were extracted from donor_cell_general.f90 for  !
! blade elements location in the grid.                                          !
! J. Sandoval, 2018                                                             !
!================================================================================

  subroutine search(point,cell,inside,dis)
    implicit none
    real::point(3),cell(8,3),dis(6)
    integer::inside,i,p1(6),p2(6),p3(6),p4(6)
    inside=0
  
    p1(5)=1
    p2(5)=4
    p3(5)=3
    p4(5)=2
  
    p1(3)=1
    p2(3)=2
    p3(3)=6
    p4(3)=5
  
    p1(1)=1
    p2(1)=5
    p3(1)=8
    p4(1)=4
  
    p1(2)=2
    p2(2)=3
    p3(2)=7
    p4(2)=6
  
    p1(6)=5
    p2(6)=6
    p3(6)=7
    p4(6)=8
  
    p1(4)=4
    p2(4)=8
    p3(4)=7
    p4(4)=3
  
    do i=1,6
       call vector(cell(p1(i),:),cell(p2(i),:),cell(p3(i),:),cell(p4(i),:),  &
            & point,dis(i))
       !if(abs(dis(i))<1.0E-06) dis(i)=0.
  
       if(dis(i).lt.0.0)then
          inside=0
          exit
       else
          inside=1
       end if
  
    end do
  
  end subroutine search
  
  subroutine vector(p1,p2,p3,p4,point,dis)
    implicit none
    real::dis,p1(3),p2(3),p3(3),p4(3),point(3)
    real::a1v,a2v,a3v,b1v,b2v,b3v,c1v,c2v,c3v,ab1,ab2,ab3
  
    a1v=p4(1)-p2(1)
    a2v=p4(2)-p2(2)
    a3v=p4(3)-p2(3)
    
    b1v=p3(1)-p1(1)
    b2v=p3(2)-p1(2)
    b3v=p3(3)-p1(3)
  
    c1v=point(1)-(p1(1)+p2(1)+p3(1)+p4(1))/4.
    c2v=point(2)-(p1(2)+p2(2)+p3(2)+p4(2))/4.
    c3v=point(3)-(p1(3)+p2(3)+p3(3)+p4(3))/4.
  
    ab1=a2v*b3v-a3v*b2v
    ab2=a3v*b1v-a1v*b3v
    ab3=a1v*b2v-a2v*b1v
  
    dis=c1v*ab1+c2v*ab2+c3v*ab3
  end subroutine vector

! --
  ! --
end subroutine init


