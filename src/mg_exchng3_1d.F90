
  subroutine mg_exchng3_1d (n, var)

    use global
    use global_param
    use global_app
    use global_mpi

    implicit none

    ! current grid level
    ! 
    integer :: n

    ! limited to a 1d decomposition
    ! 
    integer :: buf_kmx           ! total size of buf
    integer :: buf_jmx           ! total size of buf
    integer :: buf_imx           ! total size of buf
    
    ! loop counters, etc
    ! 
    integer :: l
    integer :: i
    integer :: j
    integer :: k
    integer :: m    

    ! variable to exchange
    ! 
    real (kind = rdf), dimension (:), intent(inout) :: var

    ! send and receive buffers
    ! 
    real (kind = rdf), dimension (:), allocatable :: bufu
    real (kind = rdf), dimension (:), allocatable :: bufd
    real (kind = rdf), dimension (:), allocatable :: ubuf
    real (kind = rdf), dimension (:), allocatable :: dbuf

    real (kind = rdf), dimension (:), allocatable :: bufr
    real (kind = rdf), dimension (:), allocatable :: bufl
    real (kind = rdf), dimension (:), allocatable :: rbuf
    real (kind = rdf), dimension (:), allocatable :: lbuf

    real (kind = rdf), dimension (:), allocatable :: buff
    real (kind = rdf), dimension (:), allocatable :: bufb
    real (kind = rdf), dimension (:), allocatable :: fbuf
    real (kind = rdf), dimension (:), allocatable :: bbuf

    integer :: ka
    integer :: kb
    integer :: ja
    integer :: jb
    integer :: ia
    integer :: ib

    integer :: la
    integer :: lb


    ! check bounds of var,2
    !
    if (lbound(var,1) /= le_idx_a(1))  print *, 'lbound problem in exchange_m_1d'
    if (ubound(var,1) /= le_idx_b(ng)) print *, 'ubound problem in exchange_m_1d'

    ! set size (based on ghost layers)
    ! 
    buf_kmx = (le_ix(n) * le_jx(n) * kgp(n))
    buf_jmx = (le_ix(n) * le_kx(n) * jgp(n))
    buf_imx = (le_jx(n) * le_kx(n) * igp(n))
    
    allocate(bufu(1:buf_kmx), bufd(1:buf_kmx), &
             ubuf(1:buf_kmx), dbuf(1:buf_kmx) )
    allocate(bufr(1:buf_jmx), bufl(1:buf_jmx), &
             rbuf(1:buf_jmx), lbuf(1:buf_jmx) )
    allocate(buff(1:buf_imx), bufb(1:buf_imx), &
             fbuf(1:buf_imx), bbuf(1:buf_imx) )

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
    
    ! sync before communication
    ! 
    call mpi_barrier(mpi_comm_world, ierr)

    ! send to mydown & receive from myup
    ! 
    call mpi_sendrecv (bufd, buf_kmx, mpi_real, mydown, 99, &
                       ubuf, buf_kmx, mpi_real, myup, 99, &
                       mpi_comm_world, status, ierr)

    ! send to myup and receive from mydown
    ! 
    call mpi_sendrecv (bufu, buf_kmx, mpi_real, myup, 100, &
                       dbuf, buf_kmx, mpi_real, mydown, 100, &
                       mpi_comm_world, status, ierr)

    ! send to myleft & receive from myright
    ! 
    call mpi_sendrecv (bufl, buf_jmx, mpi_real, myleft, 199, &
                       rbuf, buf_jmx, mpi_real, myright, 199, &
                       mpi_comm_world, status, ierr)

    ! send to myright and receive from myleft
    ! 
    call mpi_sendrecv (bufr, buf_jmx, mpi_real, myright, 200, &
                       lbuf, buf_jmx, mpi_real, myleft, 200, &
                       mpi_comm_world, status, ierr)

    ! send to myback & receive from myfront
    ! 
    call mpi_sendrecv (bufb, buf_imx, mpi_real, myback, 299, &
                       fbuf, buf_imx, mpi_real, myfront, 299, &
                       mpi_comm_world, status, ierr)

    ! send to myfront and receive from myback
    ! 
    call mpi_sendrecv (buff, buf_imx, mpi_real, myfront, 300, &
                       bbuf, buf_imx, mpi_real, myback, 300, &
                       mpi_comm_world, status, ierr)

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

      deallocate(ubuf, dbuf, bufu, bufd)

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

      deallocate(rbuf, lbuf, bufr, bufl)

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

      deallocate(fbuf, bbuf, buff, bufb)


     end subroutine mg_exchng3_1d
