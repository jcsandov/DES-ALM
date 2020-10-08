
subroutine rhs_diag_solver ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! implicit, diagonal approximate factorization solver
  !
  ! input
  !     rh(4,ijk)
  !     dtau(ijk)
  !     csi(ijk), eta(ijk), zet(ijk)
  !     aj(ijk)
  !     dc, de, dz
  !     dc2, de2, dz2
  !
  ! output
  !     rh(4,ijk)
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  implicit none
  
  integer :: info

  ! dimension lhs and temporary rhs
  real (kind = rdf), dimension(:,:), allocatable :: aw, ap, ae, aq

  ! boundary condition array
  real (kind = rdf), dimension(:,:,:,:), allocatable :: brh

  ! local dummy variables
  real (kind = rdf), dimension(4) :: lambda, r, rp
  real (kind = rdf), dimension(4,4) :: a, s

  real (kind = rdf) :: dc2, de2, dz2
  real (kind = rdf) :: dcsq, desq, dzsq
  real (kind = rdf) :: dtc, dds
  real (kind = rdf) :: tmp
  real (kind = rdf) :: alpha, alpha_mult

  integer :: ilength, jlength, klength
  integer :: va, vb

  ! eigenvalue multipliers
  !
  lambda(1)=zero; lambda(2)=zero; lambda(3)=one; lambda(4)=-one

  ! grid spacings
  !
  dc2=pt5*dc; de2=pt5*de; dz2=pt5*dz
  dcsq=dc*dc; desq=de*de; dzsq=dz*dz

  ! compute rh = (-dt*S^(-1)*R)
  !
  alpha_mult = e_source * one_pt_five / delti
  s = zero

  do k=k_mysta,k_myend
  do j=j_mysta,j_myend
  do i=i_mysta,idend !i_myend

     alpha =one+alpha_mult*dtau(i,j,k)
     s(1,1)=beta
     s(2,2)=one/alpha
     s(3,3)=one/alpha
     s(4,4)=one/alpha

     r(:) =-dtau(i,j,k)*rh(:,i,j,k)*aj(i,j,k)
     rp(:)=matmul(s,r)
     rh(:,i,j,k)= rp(:)

  end do
  end do
  end do

  !==============================
  !
  ! CSI SWEEP
  !
  !==============================
  va=i_mysta
  vb=idend
  ilength = vb-va+1
  if (allocated(aw)) deallocate (aw,ap,ae,aq)
  allocate (aw(va:vb,1:4),ap(va:vb,1:4),ae(va:vb,1:4),aq(va:vb,1:4) )
  aw=zero; ap=one; ae=zero; aq=zero

  ! save rh for (semi-explicit) bc treatment
  ! 
   if (n == 1) then

     if (allocated(brh)) deallocate (brh)
     allocate (brh(1:4,1:4,j_mysta:j_myend,k_mysta:k_myend)) 
     brh=zero

     if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
        do k=k_mysta,k_myend
        do j=j_mysta,j_myend
           brh(:,1,j,k)=sa(:,1)*rh(:,i_mysta,j,k)+sb(:,1)*rh(:,i_mysta+1,j,k)
        end do
        end do
     end if

     if (myfront /= mpi_proc_null .and. btype(2,myzone) /= 5) then
        do k=k_mysta,k_myend
        do j=j_mysta,j_myend
           brh(:,2,j,k)=sa(:,2)*rh(:,i_myend,j,k)+sb(:,2)*rh(:,i_myend-1,j,k)
        end do
        end do
     end if
  end if

  ! compute Ma^(-1) * (R)
  !
  do k=k_mysta,k_myend
  do j=j_mysta,j_myend
  do i=i_mysta,idend !i_myend
     a(:,:)=mai(:,:,i,j,k)
     r(:)  =rh(:,i,j,k)
     rp(:) =matmul(a, r)
     rh(:,i,j,k)= rp(:)
  end do
  end do
  end do

  if (n == 1) then
     if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           a(:,:)= mai(:,:,i_mysta-1,j,k)
           r(:)  = brh(:,1,j,k)
           rp(:) = matmul(a, r)
           brh(:,1,j,k) = rp(:)
        end do
        end do
     end if
     if (myfront /= mpi_proc_null .and. btype(2,myzone) /= 5) then
        do k = k_mysta, k_myend
        do j = j_mysta, j_myend
           a(:,:)= mai(:,:,i_myend+1,j,k)
           r(:)  = brh(:,2,j,k)
           rp(:) = matmul(a, r)
           brh(:,2,j,k) = rp(:)
        end do
        end do
     end if
  end if
  
  ! csi-sweep each jk line
  ! 
  do k=k_mysta,k_myend
  do j=j_mysta,j_myend

     do i=i_mysta,i_myend ! don't include boundary (idend)
                          ! need special treatment for boundary

        ! left-hand side

        ! center pt (i)            
        !
        dtc=one
        if (dtau(i,j,k) == zero) dtc=zero
        tmp = dtau(i,j,k) * ep(1,1,myzone) * spr(1,i,j,k) * dcsq
        ap(i,1) = one + two * tmp
        ap(i,2) = ap(i,1)
        ap(i,3) = ap(i,1)
        ap(i,4) = ap(i,1)

        ! left pt (i-1)
        !
        dds     = dtc * dtau(i-1,j,k) * dc2 * spr(1,i-1,j,k)
        aw(i,1) = - tmp
        aw(i,2) = - tmp
        aw(i,3) = - dds * lambda(3) - tmp
        aw(i,4) = - dds * lambda(4) - tmp

        ! right pt (i+1)
        !
        dds     = dtc * dtau(i+1,j,k) * dc2 * spr(1,i+1,j,k)
        ae(i,1) = - tmp
        ae(i,2) = - tmp
        ae(i,3) =   dds * lambda(3) - tmp
        ae(i,4) =   dds * lambda(4) - tmp

        ! right-hand side
        !
        aq(i,1:4) = rh(1:4,i,j,k)

     end do

     ! semi-explicit bc treatment
     !
     if (n == 1) then

        ! back-side bc
        ! 
        if (myback == mpi_proc_null .and. btype(1,myzone) /= 5) then
           i = i_mysta
           aq(i,1:4)= aq(i,1:4) - aw(i,1:4) * brh(1:4,1,j,k)
        end if

        ! front-side bc (adjust equations so we can solve this plane)
        ! 
        if (myfront == mpi_proc_null .and. btype(2,myzone) /= 5) then
            i = i_myend
            aq(i,1:4)= aq(i,1:4) - ae(i,1:4) * brh(1:4,2,j,k)
        end if
          
        if ( myfront == mpi_proc_null .and. btype(2,myzone) == 5) then
           dds=dtau(idend,j,k)*dc*spr(1,idend,j,k)
           ap(idend,1)=one
           ap(idend,2)=one
           ap(idend,3)=one+dds*lambda(3)
           ap(idend,4)=one

           dds=dtau(idend-1,j,k)*dc*spr(1,idend-1,j,k)
           aw(idend,1)=zero
           aw(idend,2)=zero
           aw(idend,3)=-dds * lambda(3)
           aw(idend,4)=zero
           
           aq(idend,1:4) = rh(1:4,idend,j,k)
        end if

     end if

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
            li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) ) then
            aw(li_blk_ia(n,nb):li_blk_ib(n,nb),1:4) = zero
            ap(li_blk_ia(n,nb):li_blk_ib(n,nb),1:4) = zero
            ae(li_blk_ia(n,nb):li_blk_ib(n,nb),1:4) = zero
           if (li_blk_ia(n,nb) > i_mysta) ae(li_blk_ia(n,nb)-1,1:4)=zero
           if (li_blk_ib(n,nb) < i_myend) aw(li_blk_ib(n,nb)+1,1:4)=zero
        end if
     end do
     end if

     ! solve linear tridiagonal eqns
     !
     call sgtsv(ilength,1,aw(va+1:vb,1),ap(va:vb,1),ae(va:vb-1,1),aq(va:vb,1), &
                ilength, info)
     call sgtsv(ilength,1,aw(va+1:vb,2),ap(va:vb,2),ae(va:vb-1,2),aq(va:vb,2), &
                ilength, info)
     call sgtsv(ilength,1,aw(va+1:vb,3),ap(va:vb,3),ae(va:vb-1,3),aq(va:vb,3), &
                ilength, info)
     call sgtsv(ilength,1,aw(va+1:vb,4),ap(va:vb,4),ae(va:vb-1,4),aq(va:vb,4), &
                ilength, info)

     ! put right-hand side back
     !
     do i = i_mysta, idend
        rh(1:4,i,j,k)=aq(i,1:4)
     end do

  end do
  end do

  !==============================
  !
  ! ETA SWEEP
  !
  !==============================
  va = j_mysta
  vb = j_myend
  jlength = vb-va+1
  if (allocated(aw)) deallocate (aw,ap,ae,aq)
  allocate (aw(va:vb,1:4),ap(va:vb,1:4),ae(va:vb,1:4),aq(va:vb,1:4) )
  aw=zero; ap=one; ae=zero; aq=zero

  !
  ! save rh for semi-explicit bc treatment
  !
  ! apply bc on fine grid only
  !
  if (n == 1) then

     if (allocated(brh)) deallocate (brh)
     allocate (brh(1:4,1:4,i_mysta:idend,k_mysta:k_myend)) 
     brh=zero

     if (myleft == mpi_proc_null .and. btype(3,myzone) /= 5) then
        do k = k_mysta, k_myend
        do i = i_mysta, idend !i_myend
           brh(:,1,i,k)=sa(:,3)*rh(:,i,j_mysta,k)+sb(:,3)*rh(:,i,j_mysta+1,k)
        end do
        end do
     end if
     
     if (myright == mpi_proc_null .and. btype(4,myzone) /= 5) then
        do k = k_mysta, k_myend
        do i = i_mysta, idend !i_myend
           brh(:,2,i,k)=sa(:,4)*rh(:,i,j_myend,k)+sb(:,4)*rh(:,i,j_myend-1,k)
        end do
        end do
     end if
  end if

  ! compute Mb^-1 * ( Ma * rh) => (N1)^-1*rh
  !
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, idend !i_myend
     r(:) = rh(:,i,j,k)
     a(:,:) = n1i(:,:,i,j,k)
     rp(:) = matmul(a, r)
     rh(:,i,j,k) = rp(:)
  end do
  end do
  end do

  ! semi-explicit bc treatment
  ! 
  if (n == 1) then   

     if (myleft == mpi_proc_null .and. btype(3,myzone) /= 5) then
        do k = k_mysta, k_myend
        do i = i_mysta, idend !i_myend
           r(:)   = brh(:,1,i,k)
           a(:,:) = n1i(:,:,i,1,k)
           rp(:)  = matmul(a, r)
           brh(:,1,i,k) = rp(:)
        end do
        end do
     end if
     
     if (myright == mpi_proc_null .and. btype(4,myzone) /= 5) then
        do k = k_mysta, k_myend
        do i = i_mysta, idend !i_myend
           r(:)   = brh(:,2,i,k)
           a(:,:) = n1i(:,:,i,j_myend+1,k)
           rp(:)  = matmul(a, r)
           brh(:,2,i,k) = rp(:)
        end do
        end do
     end if
  end if

  ! eta-sweep each ik line
  !
  do k = k_mysta, k_myend
  do i = i_mysta, idend !i_myend

     do j = j_mysta, j_myend

        ! left-hand side

        ! center pt (j)            
        !
        dtc = one
        if (dtau(i,j,k) == zero) dtc=zero
        tmp = dtau(i,j,k) * ep(1,2,myzone) * spr(2,i,j,k) * desq
        ap(j,1) = one + two * tmp
        ap(j,2) = ap(j,1)
        ap(j,3) = ap(j,1)
        ap(j,4) = ap(j,1)

        ! left pt (j-1)
        !
        dds     = dtc * dtau(i,j-1,k) * de2 * spr(2,i,j-1,k)
        aw(j,1) = - tmp
        aw(j,2) = - tmp
        aw(j,3) = - dds * lambda(3) - tmp
        aw(j,4) = - dds * lambda(4) - tmp

        ! right pt (j+1)
        !
        dds     = dtc * dtau(i,j+1,k) * de2 * spr(2,i,j+1,k)
        ae(j,1) = - tmp
        ae(j,2) = - tmp
        ae(j,3) = dds * lambda(3) - tmp
        ae(j,4) = dds * lambda(4) - tmp

        ! right-hand side
        !
        aq(j,1:4) = rh(1:4,i,j,k)

     end do

     ! semi-explicit bc treatment
     !
     if (n == 1) then

        ! left-side bc
        ! 
        if (myleft == mpi_proc_null .and. btype(3,myzone) /= 5) then
           j = j_mysta
           aq(j,1:4)= aq(j,1:4) - aw(j,1:4) * brh(1:4,1,i,k)
        end if

        ! right-side bc
        ! 
        if (myright == mpi_proc_null .and. btype(4,myzone) /= 5) then
           j = j_myend
           aq(j,1:4)= aq(j,1:4) - ae(j,1:4) * brh(1:4,2,i,k)
        end if

     end if

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
            li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
            aw(li_blk_ja(n,nb):li_blk_jb(n,nb),1:4) = zero
            ap(li_blk_ja(n,nb):li_blk_jb(n,nb),1:4) = zero
            ae(li_blk_ja(n,nb):li_blk_jb(n,nb),1:4) = zero
           if (li_blk_ja(n,nb) > j_mysta) ae(li_blk_ja(n,nb)-1,1:4) = zero
           if (li_blk_jb(n,nb) < j_myend) aw(li_blk_jb(n,nb)+1,1:4) = zero
        end if
     end do
     end if

     ! solve linear tridiagonal eqns
     !
     call sgtsv(jlength,1,aw(va+1:vb,1),ap(va:vb,1),ae(va:vb-1,1),aq(va:vb,1), &
                jlength, info)
     call sgtsv(jlength,1,aw(va+1:vb,2),ap(va:vb,2),ae(va:vb-1,2),aq(va:vb,2), &
                jlength, info)
     call sgtsv(jlength,1,aw(va+1:vb,3),ap(va:vb,3),ae(va:vb-1,3),aq(va:vb,3), &
                jlength, info)
     call sgtsv(jlength,1,aw(va+1:vb,4),ap(va:vb,4),ae(va:vb-1,4),aq(va:vb,4), &
                jlength, info)

     ! put right-hand side back
     !
     do j = j_mysta, j_myend
        rh(1:4,i,j,k) = aq(j,1:4)
     end do

  end do
  end do

  !==============================
  !
  ! ZET SWEEP
  !
  !==============================
  va = k_mysta
  vb = k_myend
  klength = vb-va+1
  if (allocated(aw)) deallocate (aw,ap,ae,aq)
  allocate (aw(va:vb,1:4),ap(va:vb,1:4),ae(va:vb,1:4),aq(va:vb,1:4) )
  aw=zero; ap=one; ae=zero; aq=zero

  !
  ! save rh for semi-explicit bc treatment
  !
  ! apply bc on fine grid only
  !
  if (n == 1) then

     if (allocated(brh)) deallocate (brh)
      allocate (brh(1:4,1:4,i_mysta:idend,j_mysta:j_myend)) 
      brh=zero

     if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
        do j = j_mysta, j_myend
        do i = i_mysta, idend !i_myend
           brh(:,1,i,j)=sa(:,5)*rh(:,i,j,k_mysta)+sb(:,5)*rh(:,i,j,k_mysta+1)
        end do
        end do
     end if
     
     if (myup == mpi_proc_null .and. btype(6,myzone) /= 6) then
        do j = j_mysta, j_myend
        do i = i_mysta, idend !i_myend
           brh(:,2,i,j)=sa(:,6)*rh(:,i,j,k_myend)+sb(:,6)*rh(:,i,j,k_myend-1)
        end do
        end do
     end if
  end if

  ! compute Mb^-1 * ( Ma * rh) => (N1)^-1*rh
  !
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, idend !i_myend
     r(:)  = rh(:,i,j,k)
     a(:,:)= n2i(:,:,i,j,k)
     rp(:) = matmul(a, r)
     rh(:,i,j,k) = rp(:)
  end do
  end do
  end do

  ! semi-explicit bc treatment
  ! 
  if (n == 1) then   

     if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
        do j = j_mysta, j_myend
        do i = i_mysta, idend !i_myend
           r(:)   = brh(:,1,i,j)
           a(:,:) = n2i(:,:,i,j,1)
           rp(:)  = matmul(a, r)
           brh(:,1,i,j) = rp(:)
        end do
        end do
     end if
     
     if (myup == mpi_proc_null .and. btype(6,myzone) /= 5) then
        do j = j_mysta, j_myend
        do i = i_mysta, idend !i_myend
           r(:)   = brh(:,2,i,j)
           a(:,:) = n2i(:,:,i,j,k_myend+1)
           rp(:)  = matmul(a, r)
           brh(:,2,i,j) = rp(:)
        end do
        end do
     end if

  end if

  ! zet-sweep each ik line
  !
  do j = j_mysta, j_myend
  do i = i_mysta, idend !i_myend

     do k = k_mysta, k_myend

        ! left-hand side

        ! center pt (k)            
        !
        dtc = one
        if (dtau(i,j,k) == zero) dtc=zero
        tmp = dtau(i,j,k) * ep(1,3,myzone) * spr(3,i,j,k) * dzsq
        ap(k,1) = one + two * tmp
        ap(k,2) = ap(k,1)
        ap(k,3) = ap(k,1)
        ap(k,4) = ap(k,1)

        ! left pt (k-1)
        !
        dds     = dtc * dtau(i,j,k-1) * dz2 * spr(3,i,j,k-1)
        aw(k,1) = - tmp
        aw(k,2) = - tmp
        aw(k,3) = - dds * lambda(3) - tmp
        aw(k,4) = - dds * lambda(4) - tmp

        ! right pt (k+1)
        !
        dds     = dtc * dtau(i,j,k+1) * dz2 * spr(3,i,j,k+1)
        ae(k,1) = - tmp
        ae(k,2) = - tmp
        ae(k,3) = dds * lambda(3) - tmp
        ae(k,4) = dds * lambda(4) - tmp

        ! right-hand side
        !
        aq(k,1:4) = rh(1:4,i,j,k)

     end do

     ! semi-explicit bc treatment
     !
     if (n == 1) then

        ! left-side bc
        ! 
        if (mydown == mpi_proc_null .and. btype(5,myzone) /= 5) then
           k = k_mysta
           aq(k,1:4)= aq(k,1:4) - aw(k,1:4) * brh(1:4,1,i,j)
        end if

        ! right-side bc
        ! 
        if (myup == mpi_proc_null .and. btype(6,myzone) /= 5) then
           k = k_myend
           aq(k,1:4)= aq(k,1:4) - ae(k,1:4) * brh(1:4,2,i,j)
        end if

     end if

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) .and. &
            li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
            aw(li_blk_ka(n,nb):li_blk_kb(n,nb),1:4) = zero
            ap(li_blk_ka(n,nb):li_blk_kb(n,nb),1:4) = zero
            ae(li_blk_ka(n,nb):li_blk_kb(n,nb),1:4) = zero
            if (li_blk_ka(n,nb) > k_mysta) ae(li_blk_ka(n,nb)-1,1:4) = zero
            if (li_blk_kb(n,nb) < k_myend) aw(li_blk_kb(n,nb)+1,1:4) = zero
        end if
     end do
     end if

     ! solve linear tridiagonal eqns
     !
     klength = k_myend - k_mysta + 1
     call sgtsv(klength,1,aw(va+1:vb,1),ap(va:vb,1),ae(va:vb-1,1),aq(va:vb,1), &
                klength, info)
     call sgtsv(klength,1,aw(va+1:vb,2),ap(va:vb,2),ae(va:vb-1,2),aq(va:vb,2), &
                klength, info)
     call sgtsv(klength,1,aw(va+1:vb,3),ap(va:vb,3),ae(va:vb-1,3),aq(va:vb,3), &
                klength, info)
     call sgtsv(klength,1,aw(va+1:vb,4),ap(va:vb,4),ae(va:vb-1,4),aq(va:vb,4), &
                klength, info)

     ! put right-hand side back
     !
     do k = k_mysta, k_myend
        rh(1:4,i,j,k)=aq(k,1:4)
     end do

  end do
  end do

  if (allocated(aw))  deallocate (aw,ap,ae,aq)
  if (allocated(brh)) deallocate (brh)

  !==============================
  !
  ! Final Update
  !
  !==============================
  ! dQ = Mc * dQ^***

  do k=k_mysta, k_myend
  do j=j_mysta, j_myend
  do i=i_mysta, idend !i_myend
     a(:,:)= mc(:,:,i,j,k)
     r(:)  = rh(:,i,j,k)
     rp(:) = matmul(a,r)
     rh(:,i,j,k) = rp(:)
  end do
  end do
  end do


end subroutine rhs_diag_solver


