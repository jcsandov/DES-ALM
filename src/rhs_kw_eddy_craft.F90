
subroutine rhs_kw_eddy_craft
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Eddy viscosity based on Craft, Launder, and Suga (1995)    
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! local variable
  real (kind = rdf) :: sso, wwo, xke, sinv, winv
  real (kind = rdf) :: tum0, tum1, tum2, tum3, Cm, Rtur, fmu

  integer :: ista, &
             iend

  ista = i_mysta
  iend = i_myend

  if ( myback  == mpi_proc_null ) ista = i_mysta - 1
  if ( myfront == mpi_proc_null ) iend = i_myend + 1 

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = ista, iend

     if (q(5,i,j,k) == zero) then

        xnut(i,j,k) = zero

     else

        ! Calculate Cm
        sso = four*so(i,j,k)
        wwo = four*wo(i,j,k)

        xke  = one / ck / q(6,i,j,k)
        sinv = xke * sqrt(pt5*sso)
        winv = xke * sqrt(pt5*wwo)
        tum0 = max(sinv,winv)

        tum1 = exp(-0.75_rdf * tum0)
        tum2 = one - exp(-0.36_rdf/tum1)
        tum3 = one + 0.35_rdf * tum0**onept5

        Cm  = 0.3_rdf * tum2 / tum3 / ck
        Rtur= ren * q(5,i,j,k) / ck / q(6,i,j,k)
        fmu = min(one, 0.2_rdf+0.8_rdf*Rtur/50.0_rdf)

        xnut(i,j,k) = Cm * fmu * q(5,i,j,k) / q(6,i,j,k)

     end if

  end do
  end do
  end do

end subroutine rhs_kw_eddy_craft

