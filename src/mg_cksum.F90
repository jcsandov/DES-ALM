subroutine mg_cksum (comment, n, var)

  use checksum

  integer , intent(in) :: n
  character (len = *), intent(in) :: comment
  real (kind = rdf), dimension (:,:), intent(in) :: var

  real (kind = rdf), dimension (:,:,:,:), allocatable :: var_4d

  allocate(var_4d(1:ubound(var,1),li_ia(n):li_ib(n),li_ja(n):li_jb(n)&
       &,li_ka(n):li_kb(n)))

  do k = li_ka(n), li_kb(n)
  do j = li_ja(n), li_jb(n)
  do i = li_ia(n), li_ib(n)
     l = le_idx(i,j,k,n)
     var_4d(:,i,j,k) = var(:,l)
  end do
  end do
  end do

  call checksum_4d_par (comment, var_4d)

end subroutine mg_cksum
  
