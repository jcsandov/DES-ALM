
!============================================================
!
!
! collect all of solution to root process and write to file
!
!
!============================================================

  subroutine collect_solu ()

    use global
    use global_param
    use global_app
    use global_mpi

    integer ::  tmp_mx
    integer :: sbuf_mx

    integer, dimension(:,:,:,:), allocatable :: sbuf_idx

    real (kind = rdf), dimension (:), allocatable :: sbuf

    real (kind = rdf), dimension(:,:), allocatable :: qtmp
    real (kind = rdf), dimension(:,:), allocatable :: qntmp
    real (kind = rdf), dimension (:), allocatable :: nutmp
    real (kind = rdf), dimension(:,:), allocatable :: uijtmp

    integer :: nvar
    integer :: ia, ib
    integer :: ja, jb
    integer :: ka, kb
    integer :: kk, jj, ii

    integer :: nv
    integer :: l1, l2

    nvar = 2 * me
    if ( turbulence ) nvar = nvar + 1
    if ( nlinc ) nvar = nvar + 6

    sbuf_mx = nvar * img(1) * jmg(1) * kmg(1)

    allocate(sbuf(sbuf_mx), &
         sbuf_idx(nvar,img(1),jmg(1),kmg(1)))

    l = 0
    do k = 1, kmg(1)
       do j = 1, jmg(1)
          do i = 1, img(1)
             do nv = 1, nvar
                l = l + 1; sbuf_idx(nv,i,j,k) = l
             end do
          end do
       end do
    end do

    if (myid == root) then

       ! temporary arrays on root
       ! 
       tmp_mx = img(1) * jmg(1) * kmg(1)
       
       allocate(qtmp(1:me,me*tmp_mx), &
               qntmp(1:me,me*tmp_mx) )

       ! get solution from 'solu' file
       !
       if ( turbulence ) then

          ! allocate temp array for turbulence viscosity
          ! 
          allocate(nutmp(tmp_mx))

          open  (unit = 42, file = 'solu', form = 'unformatted')
          read  (unit = 42)  ((qtmp(nv,l),l = 1, tmp_mx), nv = 1, me)
          read  (unit = 42)     (nutmp(l),l = 1, tmp_mx)
          read  (unit = 42) ((qntmp(nv,l),l = 1, tmp_mx), nv = 1, me)
          close (unit = 42)
       else
          open  (unit = 42, file = 'solu', form = 'unformatted')
          read  (unit = 42)  ((qtmp(nv,l),l = 1, tmp_mx), nv = 1, me)
          read  (unit = 42) ((qntmp(nv,l),l = 1, tmp_mx), nv = 1, me)
          close (unit = 42)
       end if

       ! pack send buf
       ! 
       l1 = 0
       l2 = 0
       do k = 1, kmg(1)
          do j = 1, jmg(1)
             do i = 1, img(1)
                l1 = l1 + 1
                do nv = 1, me
                   l2 = l2 + 1; sbuf(l2) =  qtmp(nv,l1)
                end do
                do nv = 1, me
                   l2 = l2 + 1; sbuf(l2) = qntmp(nv,l1)
                end do
                if ( turbulence ) then
                   l2 = l2 + 1; sbuf(l2) =   nutmp(l1)
                end if
             end do
          end do
       end do

       if ( turbulence ) then
          deallocate (qtmp, qntmp, nutmp)
       else
          deallocate (qtmp, qntmp)
       end if

    end if

    ! broadcast to all processes
    ! 
    call mpi_bcast (sbuf, sbuf_mx, MPI_REAL, root,&
         & mpi_comm_world, ierr)

    ! upack subuf including ghost celss into
    ! local copies of q, qn, xnut
    !
    ia = le_ia(1)
    ja = le_ja(1)
    ka = le_ka(1)

    ib = le_ib(1)
    jb = le_jb(1)
    kb = le_kb(1)

    ! 1d decomp only; need more for 2d or 3d
    ! 
    if (myup == mpi_proc_null) kb = li_kb(1)
    if (mydown == mpi_proc_null) ka = li_ka(1)

    do k = ka, kb
       do j = ja, jb
          do i = ia, ib
             ii = i + gi_ia(1) - 1
             jj = j + gi_ja(1) - 1
             kk = k + gi_ka(1) - 1

             do nv = 1, me
                q(nv,le_idx(i,j,k,1)) = sbuf(sbuf_idx(nv,ii,jj,kk))                
             end do
             do nv = 1, me
                qn(nv,le_idx(i,j,k,1)) = sbuf(sbuf_idx(nv+me,ii,jj,kk))                
             end do
            if ( turbulence ) xnut(le_idx(i,j,k,1)) = sbuf(sbuf_idx(nvar,ii,jj,kk))
          end do
       end do
    end do

    deallocate(sbuf, sbuf_idx)
    
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

!     call checksum_2d_par('q  :',  qtmp)
!     call checksum_2d_par('qn :', qntmp)

    deallocate (qtmp, qntmp)

  end subroutine collect_solu
