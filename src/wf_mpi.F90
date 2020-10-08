!
! Module : wfunction
!
!
! Joongcheol Paik <joongcheol.paik@ce.gatech.edu>
!
! written : 18 June 2005
! ----------------------
!
module wf_mpi
  use global
  use global_param
  use global_mpi
  implicit none

  private

  public :: wf_csi, wf_eta, wf_zet

  public :: newt, fdjac
  public :: fdjac1, funcv
  public :: lnsrch, ludcmp, lubksb

  !real, public :: fmin

  real, parameter :: kapa = 0.41
  real, parameter :: Blog = 5.45
  real, parameter :: yplm =11.30
  real, parameter :: qmin = 1.0E-08
  real, parameter ::   ks = 0.00 ! roughness should read from input data

  integer :: i, j, k
  integer :: i_mystart, j_mystart, k_mystart
  integer :: i_myend,   j_myend,   k_myend
  integer :: l, l1, l2
  integer :: lu, ld

  real :: yp0, u0, re0
  real :: qc1, qc2, q2

  real :: yplus
  real :: ypks
  real :: vmod
  real :: qstar
  real :: cgama
  real :: fbeta

contains

  ! The idea to implement the wall function is as follows:
  ! Provided layer #2 & #3 are all in the range of wall function area.
  ! Assume the direction on layer #2 is same as layer #3 and
  !    velocity normal to wall are zero on layer #2.
  ! First calculate wall share stress based on the magnitude
  !    of the velocity component parallel to the wall on layer #3.
  ! Use this value to calculate the velocity magnitude on layer #2.
  ! Based on the direction information from layer #3, calculate
  !    velocity components along x, y & z directions on layer #2.

  ! Layer number is defined in the following way:
  !    the wall is defined as layer #1 and the others are defined as
  !    layer#2(first layer away from the wall), layer#3 respectively.

  ! Variables
  !  wd(L,im,jm)   : distance to the wall of the first two close-to-wall layers
  !     where the first dimension means the (layer number-1)
  !  csi(3,im,jm,3): the metrics of first three close-to-wall layers
  !  eta           : 1st dimension stores dx,dy,dz respectively
  !  zet           : the last dimension represents the layer number
  !  aj(:,:  )     : determinant of Jacobian matrix on wall
  !  xc,yc,zc....  : geometry information on the wall
  !  ut0           : friction velocity

  ! =======================
  subroutine wf_zet (kb)
  ! =======================
  integer, intent(in) :: kb
  integer :: k1, k2
  integer, parameter :: MAXITS=200
  integer :: its
  integer :: itss
  logical :: check

  real, dimension(4) :: ut
  real, dimension(:,:), allocatable :: ut0

  real :: qnor, qpar
  real :: deltb
  real :: sg11, sg22
  real :: xc, yc, zc
  real :: xe, ye, ze

  i_mystart=gi_ia(1)
  j_mystart=gi_ja(1)
  k_mystart=gi_ka(1)

  i_myend=gi_ib(1)
  j_myend=gi_jb(1)
  k_myend=gi_kb(1)

  allocate (ut0(i_mystart:i_myend,j_mystart:j_myend))

  if (kb == 1) then
     k =k_mystart
     k1=k_mystart+1
     k2=k_mystart+2
  else
     k =k_myend
     k1=k_myend-1
     k2=k_myend-2
  end if

  re0=ren
  ut(1)=0.03

  do j=j_mystart+1,j_myend-1
  do i=i_mystart+1,i_myend-1
     l =gi_2_le_idx(i,j,k ,1)
     l1=gi_2_le_idx(i,j,k1,1)
     l2=gi_2_le_idx(i,j,k2,1)

     vmod=sqrt(zet(1,l)**2+zet(2,l)**2+zet(3,l)**2)

     ! qnor = velocity normal to the wall
     ! qpar = velocity parallel to the wall

     qnor=(q(2,l2)*zet(1,l)+q(3,l2)*zet(2,l)+q(4,l2)*zet(3,l))/vmod
     qpar=(q(2,l2)**2+q(3,l2)**2+q(4,l2)**2)
     qpar=sqrt(qpar-qnor**2.)

     u0=qpar
     yp0=wd(l2)

  !   if(abs(qpar) > qmin) then
        ! call the nonlinear root find solver to calculate u_tao
        call newt(ut,1,check,MAXITS,its,itss)
        ut0(i,j)=ut(1)
  !   else
  !      ut0(i,j)=0.
  !   end if
  end do
  end do

  !do i=i_mystart,i_myend
  !   ut0(i,j_mystart)=ut0(i,j_mystart+1)
  !   ut0(i,j_myend)  =ut0(i,j_myend-1)
  !end do

  do j=j_mystart+1,j_myend-1
  do i=i_mystart+1,i_myend-1
     l =gi_2_le_idx(i,j,k ,1)
     l1=gi_2_le_idx(i,j,k1,1)
     l2=gi_2_le_idx(i,j,k2,1)

     qpar=(q(2,l2)**2+q(3,l2)**2+q(4,l2)**2)

     if (qpar > qmin) then
        yplus=wd(l1)*ren*ut0(i,j)

        if (ks > zero)then
           ypks = ut0(i,j)*ren*ks
           if (ypks <= 90)then
              deltb=(Blog-8.5+1/kapa*log(ypks))*sin(0.4258*(log(ypks)-0.811))
           else if(ks > 0.01) then
              deltb=(Blog-8.5+1/kapa*log(ypks))
           end if
        else
            deltb=zero
        end if

        !Apply wall function on layer 2
        if (yplus > yplm) then
           q2=(log(yplus)/kapa+Blog-deltb)*ut0(i,j)
        else
           q2=yplus*ut0(i,j)
        endif

        lu=gi_2_le_idx(i+1,j,k,1)
        ld=gi_2_le_idx(i-1,j,k,1)
        xc=pt5*(x(lu)-x(ld))
        yc=pt5*(y(lu)-y(ld))
        zc=pt5*(z(lu)-z(ld))

        lu=gi_2_le_idx(i,j+1,k,1)
        ld=gi_2_le_idx(i,j-1,k,1)
        xe=pt5*(x(lu)-x(ld))
        ye=pt5*(y(lu)-y(ld))
        ze=pt5*(z(lu)-z(ld))

        !Square root of g11,g22
        sg11=sqrt(xc**2+yc**2+zc**2)
        sg22=sqrt(xe**2+ye**2+ze**2)

        !  qc1,qc2 are the velocity magnitude along the csi+,eta+ directions
        !  on layer #3. Here the velocity direction on layer #2 is assumed
        !  to be same as that on layer #3.
        qc1=(csi(1,l2)*q(2,l2)+csi(2,l2)*q(3,l2)+csi(3,l2)*q(4,l2))*sg11
        qc2=(eta(1,l2)*q(2,l2)+eta(2,l2)*q(3,l2)+eta(3,l2)*q(4,l2))*sg22

        ! cgama : angle between csi+ and eta+ directions
        ! fbeta  : angle between u_3 and csi+ direction
        cgama=(xc*xe+yc*ye+zc*ze)/(sg11*sg22)

        fbeta=sqrt(qc1**2+qc2**2+2.*qc1*qc2*cgama)
        fbeta=((qc1+qc2*cgama)/fbeta)
        fbeta=min(fbeta, 1.0)
        fbeta=max(fbeta,-1.0)

        if (qc2 > zero) then
           fbeta= acos(fbeta)
        else
           fbeta=-acos(fbeta)
        endif

        ! qstar : velocity magnitude on layer #2
        qstar=q2

        ! qc1,qc2 are the velocity magnitude along csi+, eta+ directions
        qc2=qstar*sin(fbeta)/sin(acos(cgama))
        qc1=qstar*cos(fbeta)-qc2*cgama

        q(2,l1)=qc1*xc/sg11+qc2*xe/sg22
        q(3,l1)=qc1*yc/sg11+qc2*ye/sg22
        q(4,l1)=qc1*zc/sg11+qc2*ze/sg22


        !if (k_e .or. k_w)
        !   qout(4,i,j)=max((ut0(i,j)**2)/(cmu**0.5),ckin)
        !   qout(5,i,j)=max((ut0(i,j)**3)/(kapa*wd(i,j)),cein)
        !else
           q(5,l1)=kapa*yplus/re0
        !end if
     else
        q(2:4,l1)=zero
        q(5,l1)=0.1/re0
     end if
  end do
  end do

  deallocate (ut0)

  end subroutine wf_zet


  ! =======================
  subroutine wf_eta (jb)
  ! =======================
  integer, intent(in) :: jb
  integer :: j1, j2
  integer, parameter :: MAXITS=200
  integer :: its
  integer :: itss
  logical :: check

  real, dimension(4) :: ut
  real, dimension(:,:), allocatable :: ut0

  real :: qnor, qpar
  real :: deltb
  real :: sg11, sg33
  real :: xc, yc, zc
  real :: xz, yz, zz

  i_mystart=gi_ia(1)
  j_mystart=gi_ja(1)
  k_mystart=gi_ka(1)

  i_myend=gi_ib(1)
  j_myend=gi_jb(1)
  k_myend=gi_kb(1)

  allocate (ut0(i_mystart:i_myend,k_mystart:k_myend))

  if (jb == 1) then
     j =j_mystart
     j1=j_mystart+1
     j2=j_mystart+2
  else
     j =j_myend
     j1=j_myend-1
     j2=j_myend-2
  end if

  re0=ren
  ut(1)=0.03

  do k=k_mystart+1,k_myend-1
  do i=i_mystart+1,i_myend-1
     l =gi_2_le_idx(i,j ,k,1)
     l1=gi_2_le_idx(i,j1,k,1)
     l2=gi_2_le_idx(i,j2,k,1)

     vmod=sqrt(eta(1,l)**2+eta(2,l)**2+eta(3,l)**2)

     ! qnor = velocity normal to the wall
     ! qpar = velocity parallel to the wall

     qnor=(q(2,l2)*eta(1,l)+q(3,l2)*eta(2,l)+q(4,l2)*eta(3,l))/vmod
     qpar=(q(2,l2)**2+q(3,l2)**2+q(4,l2)**2)
     qpar=sqrt(qpar-qnor**2)

     u0=qpar
     yp0=wd(l2)

!     if(abs(qpar) > qmin) then
        ! call the nonlinear root find solver to calculate u_tao
        call newt(ut,1,check,MAXITS,its,itss)
        ut0(i,k)=ut(1)
!     else
!        ut0(i,k)=0.
!     end if
  end do
  end do

  !do i=i_mystart,i_myend
  !   ut0(i,k_mystart)=ut0(i,k_mystart+1)
  !   ut0(i,k_myend)  =ut0(i,k_myend-1)
  !end do

  do k=k_mystart+1,k_myend-1
  do i=i_mystart+1,i_myend-1
     l =gi_2_le_idx(i,j ,k,1)
     l1=gi_2_le_idx(i,j1,k,1)
     l2=gi_2_le_idx(i,j2,k,1)

     qpar=(q(2,l2)**2+q(3,l2)**2+q(4,l2)**2)

     if (qpar > qmin) then
        yplus=wd(l1)*ren*ut0(i,k)

        if (ks > zero)then
           ypks = ut0(i,k)*ren*ks
           if (ypks <= 90)then
              deltb=(Blog-8.5+1/kapa*log(ypks))*sin(0.4258*(log(ypks)-0.811))
           else if(ks > 0.01) then
              deltb=(Blog-8.5+1/kapa*log(ypks))
           end if
        else
            deltb=zero
        end if

        !Apply wall function on layer 2
        if (yplus > yplm) then
           q2=(log(yplus)/kapa+Blog-deltb)*ut0(i,k)
        else
           q2=yplus*ut0(i,k)
        endif

        lu=gi_2_le_idx(i+1,j,k,1)
        ld=gi_2_le_idx(i-1,j,k,1)
        xc=pt5*(x(lu)-x(ld))
        yc=pt5*(y(lu)-y(ld))
        zc=pt5*(z(lu)-z(ld))

        lu=gi_2_le_idx(i,j,k+1,1)
        ld=gi_2_le_idx(i,j,k-1,1)
        xz=pt5*(x(lu)-x(ld))
        yz=pt5*(y(lu)-y(ld))
        zz=pt5*(z(lu)-z(ld))

        !Square root of g11,g22
        sg11=sqrt(xc**2+yc**2+zc**2)
        sg33=sqrt(xz**2+yz**2+zz**2)

        !  qc1,qc2 are the velocity magnitude along the csi+,eta+ directions
        !  on layer #3. Here the velocity direction on layer #2 is assumed
        !  to be same as that on layer #3.
        qc1=(csi(1,l2)*q(2,l2)+csi(2,l2)*q(3,l2)+csi(3,l2)*q(4,l2))*sg11
        qc2=(zet(1,l2)*q(2,l2)+zet(2,l2)*q(3,l2)+zet(3,l2)*q(4,l2))*sg33

        ! cgama : angle between csi+ and eta+ directions
        ! fbeta  : angle between u_3 and csi+ direction
        cgama=(xc*xz+yc*yz+zc*zz)/(sg11*sg33)

        fbeta=sqrt(qc1**2+qc2**2+2.*qc1*qc2*cgama)
        fbeta=((qc1+qc2*cgama)/fbeta)
        fbeta=min(fbeta, 1.0)
        fbeta=max(fbeta,-1.0)

        if (qc2 > zero) then
           fbeta= acos(fbeta)
        else
           fbeta=-acos(fbeta)
        endif

        ! qstar : velocity magnitude on layer #2
        qstar=q2

        ! qc1,qc2 are the velocity magnitude along csi+, eta+ directions
        qc2=qstar*sin(fbeta)/sin(acos(cgama))
        qc1=qstar*cos(fbeta)-qc2*cgama

        q(2,l1)=qc1*xc/sg11+qc2*xz/sg33
        q(3,l1)=qc1*yc/sg11+qc2*yz/sg33
        q(4,l1)=qc1*zc/sg11+qc2*zz/sg33


        !if (k_e .or. k_w)
        !   qout(4,i,j)=max((ut0(i,k)**2)/(cmu**0.5),ckin)
        !   qout(5,i,j)=max((ut0(i,k)**3)/(kapa*wd(i,j)),cein)
        !else
           q(5,l1)=kapa*yplus/re0
        !end if
     else
        q(2:4,l1)=zero
        q(5,l1)=0.1/re0
     end if
  end do
  end do

  deallocate (ut0)

  end subroutine wf_eta

  ! =======================
  subroutine wf_csi (ib)
  ! =======================
  integer, intent(in) :: ib
  integer :: i1, i2
  integer, parameter :: MAXITS=200
  integer :: its
  integer :: itss
  logical :: check

  real, dimension(4) :: ut
  real, dimension(:,:), allocatable :: ut0

  real :: qnor, qpar
  real :: deltb
  real :: sg22, sg33
  real :: xe, ye, ze
  real :: xz, yz, zz

  i_mystart=gi_ia(1)
  j_mystart=gi_ja(1)
  k_mystart=gi_ka(1)

  i_myend=gi_ib(1)
  j_myend=gi_jb(1)
  k_myend=gi_kb(1)

  allocate (ut0(j_mystart:j_myend,k_mystart:k_myend))

  if (ib == 1) then
     i =i_mystart
     i1=i_mystart+1
     i2=i_mystart+2
  else
     i =i_myend
     i1=i_myend-1
     i2=i_myend-2
  end if

  re0=ren
  ut(1)=0.03

  do k=k_mystart+1,k_myend-1
  do j=j_mystart+1,j_myend-1
     l =gi_2_le_idx(i ,j,k,1)
     l1=gi_2_le_idx(i1,j,k,1)
     l2=gi_2_le_idx(i2,j,k,1)

     vmod=sqrt(csi(1,l)**2+csi(2,l)**2+csi(3,l)**2)

     ! qnor = velocity normal to the wall
     ! qpar = velocity parallel to the wall

     qnor=(q(2,l2)*csi(1,l)+q(3,l2)*csi(2,l)+q(4,l2)*csi(3,l))/vmod
     qpar=(q(2,l2)**2+q(3,l2)**2+q(4,l2)**2)
     qpar=sqrt(qpar-qnor**2.)

     u0=qpar
     yp0=wd(l2)

!     if(abs(qpar) > qmin) then
        ! call the nonlinear root find solver to calculate u_tao
        call newt(ut,1,check,MAXITS,its,itss)
        ut0(j,k)=ut(1)
!     else
!        ut0(j,k)=0.
!     end if
  end do
  end do

  !do j=j_mystart,j_myend
  !   ut0(j,k_mystart) =ut0(j,k_mystart+1)
  !   ut0(j,k_myend)   =ut0(j,k_myend-1)
  !end do

  do k=k_mystart+1,k_myend-1
  do j=j_mystart+1,j_myend-1
     l =gi_2_le_idx(i ,j,k,1)
     l1=gi_2_le_idx(i1,j,k,1)
     l2=gi_2_le_idx(i2,j,k,1)

     qpar=(q(2,l2)**2+q(3,l2)**2+q(4,l2)**2)

     if (qpar > qmin) then
        yplus=wd(l1)*ren*ut0(j,k)

        if (ks > zero)then
           ypks = ut0(j,k)*ren*ks
           if (ypks <= 90)then
              deltb=(Blog-8.5+1/kapa*log(ypks))*sin(0.4258*(log(ypks)-0.811))
           else if(ks > 0.01) then
              deltb=(Blog-8.5+1/kapa*log(ypks))
           end if
        else
            deltb=zero
        end if

        !Apply wall function on layer 2
        if (yplus > yplm) then
           q2=(log(yplus)/kapa+Blog-deltb)*ut0(j,k)
        else
           q2=yplus*ut0(j,k)
        endif

        lu=gi_2_le_idx(i,j+1,k,1)
        ld=gi_2_le_idx(i,j-1,k,1)
        xe=pt5*(x(lu)-x(ld))
        ye=pt5*(y(lu)-y(ld))
        ze=pt5*(z(lu)-z(ld))

        lu=gi_2_le_idx(i,j,k+1,1)
        ld=gi_2_le_idx(i,j,k-1,1)
        xz=pt5*(x(lu)-x(ld))
        yz=pt5*(y(lu)-y(ld))
        zz=pt5*(z(lu)-z(ld))

        !Square root of g22,g33
        sg22=sqrt(xe**2+ye**2+ze**2)
        sg33=sqrt(xz**2+yz**2+zz**2)

        !  qc1,qc2 are the velocity magnitude along the csi+,eta+ directions
        !  on layer #3. Here the velocity direction on layer #2 is assumed
        !  to be same as that on layer #3.
        qc1=(eta(1,l2)*q(2,l2)+eta(2,l2)*q(3,l2)+eta(3,l2)*q(4,l2))*sg22
        qc2=(zet(1,l2)*q(2,l2)+zet(2,l2)*q(3,l2)+zet(3,l2)*q(4,l2))*sg33

        ! cgama : angle between csi+ and eta+ directions
        ! fbeta  : angle between u_3 and csi+ direction
        cgama=(xe*xz+ye*yz+ze*zz)/(sg22*sg33)

        fbeta=sqrt(qc1**2+qc2**2+2.*qc1*qc2*cgama)
        fbeta=((qc1+qc2*cgama)/fbeta)
        fbeta=min(fbeta, 1.0)
        fbeta=max(fbeta,-1.0)

        if (qc2 > zero) then
           fbeta= acos(fbeta)
        else
           fbeta=-acos(fbeta)
        endif

        ! qstar : velocity magnitude on layer #2
        qstar=q2

        ! qc1,qc2 are the velocity magnitude along csi+, eta+ directions
        qc2=qstar*sin(fbeta)/sin(acos(cgama))
        qc1=qstar*cos(fbeta)-qc2*cgama

        q(2,l1)=qc1*xe/sg22+qc2*xz/sg33
        q(3,l1)=qc1*ye/sg22+qc2*yz/sg33
        q(4,l1)=qc1*ze/sg22+qc2*zz/sg33


        !if (k_e .or. k_w)
        !   qout(4,i,j)=max((ut0(j,k)**2)/(cmu**0.5),ckin)
        !   qout(5,i,j)=max((ut0(j,k)**3)/(kapa*wd(i,j)),cein)
        !else
           q(5,l1)=kapa*yplus/re0
        !end if
     else
        q(2:4,l1)=zero
        q(5,l1)=0.1/re0
     end if
  end do
  end do

  deallocate (ut0)

  end subroutine wf_csi

!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fdjac1(n,x,fjac)
  implicit none
  integer::n
  real::x(n),fvec(n),fjac(n,n)
  if(re0*yp0*x(1)>11.3)then
     fjac(1,1)=log(re0*yp0*x(1))/kapa+Blog+1./kapa
  else
     fjac(1,1)=re0*2*x(1)*yp0
  end if

END SUBROUTINE fdjac1

SUBROUTINE funcv(n,x,fvec)
! Wall function subroutine
! x(1) : u_tao
! x(2) : |u|
! x(3) : distance to wall
! x(4) : Reynolds #
!  USE wf_variables
  implicit none
  integer::n
  real::x(n),fvec(n),ksp,deltb

!  fvec(2)=0.
!  fvec(3)=0.
!  fvec(4)=0.
!  if(x(1).lt.1.e-4) x(1)=1.e-4
  ksp=x(1)*ks*re0
  deltb=0.
  if(ksp>=2.25)then
  if(ksp.le.90)then
     deltb=(Blog-8.5+1/kapa*log(ksp))*sin(0.4258*(log(ksp)-0.811))
  else
     deltb=Blog-8.5+1/kapa*log(ksp)
  end if
  end if
  if(re0*yp0*x(1).gt.11.3)then
     fvec(1)=(1.*log(re0*yp0*x(1))/kapa+Blog-deltb)*x(1)-u0
  else
     fvec(1)=re0*x(1)**2*yp0-u0
  end if
  if(abs(x(1)*re0*yp0-11.3).lt.1.e-4) fvec(1)=0
END SUBROUTINE funcv


SUBROUTINE newt(x,n,check,MAXITS,its,itss)
  
  INTEGER n,nn,NP
  LOGICAL check
  REAL x(n),fvec,TOLF,TOLMIN,TOLX,STPMX

  PARAMETER (NP=40,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7, &
       &  STPMX=100.)
  COMMON /newtv/ fvec(NP),nn
  SAVE /newtv/
!     USES fdjac,fmin,lnsrch,lubksb,ludcmp

  integer :: its, itss
  INTEGER i,j,indx(NP),MAXITS
  REAL :: d,den,f,fold,stpmax,sum,temp,test,fjac(NP,NP),g(NP),p(NP)
  REAL :: xold(NP) !,fmin
!  EXTERNAL fmin
  nn=n
  f=fmin(x)
  test=0.
  do i=1,n
     if(abs(fvec(i)).gt.test)test=abs(fvec(i))
  end do
  if(test.lt..01*TOLF) then
     check=.false.
     return
  endif
  sum=0.
  do i=1,n
     sum=sum+x(i)**2
  end do
  stpmax=STPMX*max(sqrt(sum),float(n))
  do its=1,MAXITS
!     write(*,*) its
!     call fdjac(n,x,fvec,NP,fjac)
     call funcv(n,x,fvec)
     call fdjac1(n,x,fjac(1:n,1:n))
!     print *,fjac(1:n,1:n)
     do i=1,n
        sum=0.
        do j=1,n
           sum=sum+fjac(j,i)*fvec(j)
        end do
        g(i)=sum
     end do
     do i=1,n
        xold(i)=x(i)
     end do
     fold=f
     do i=1,n
        p(i)=-fvec(i)
     end do
!     write(*,*) '********'
     call ludcmp(fjac,n,NP,indx,d,itss)
!     write(*,*) '*******1'
     call lubksb(fjac,n,NP,indx,p)
!     write(*,*) '*******2'
     call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
!     write(*,*) '*******3'
     test=0.
     do i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
     end do
     if(test.lt.TOLF)then
        check=.false.
        return
     endif
     if(check)then
        test=0.
        den=max(f,.5*n)
        do i=1,n
           temp=abs(g(i))*max(abs(x(i)),1.)/den
           if(temp.gt.test)test=temp
        end do
        if(test.lt.TOLMIN)then
           check=.true.
        else
           check=.false.
        endif
        return
     endif
     test=0.
     do i=1,n
        temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.)
        if(temp.gt.test)test=temp
     end do
     if(test.lt.TOLX)return
  end do
!c      pause 'MAXITS exceeded in newt'
END SUBROUTINE newt

SUBROUTINE fdjac(n,x,fvec,np,df)
  INTEGER n,np,NMAX
  REAL df(np,np),fvec(n),x(n),EPS
  PARAMETER (NMAX=40,EPS=1.e-4)
!    USES funcv
  INTEGER i,j
  REAL h,temp,f(NMAX)
  do j=1,n
     temp=x(j)
     h=EPS*abs(temp)
     if(h.eq.0.)h=EPS
     x(j)=temp+h
     h=x(j)-temp
     call funcv(n,x,f)
     x(j)=temp
     do i=1,n
        df(i,j)=(f(i)-fvec(i))/h
     end do
  end do
  return
END SUBROUTINE fdjac

  FUNCTION fmin(x)
  INTEGER n,NP
  REAL fmin, x(*),fvec
  PARAMETER (NP=40)
  COMMON /newtv/ fvec(NP),n
  SAVE /newtv/
!CU    USES funcv
  INTEGER i
  REAL sum
  call funcv(n,x,fvec)
  sum=0.
  do i=1,n
     sum=sum+fvec(i)**2
  end do
  fmin=0.5*sum
  return
END FUNCTION fmin

SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
  INTEGER n
  LOGICAL check
  REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
  PARAMETER (ALF=1.e-3,TOLX=1.e-7)
  EXTERNAL func
!CU    USES func
  INTEGER i
  REAL a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp, &
       &  test,tmplam
  check=.false.
  sum=0.
  do i=1,n
     sum=sum+p(i)*p(i)
  end do

  sum=sqrt(sum)
  if(sum.gt.stpmax)then
     do i=1,n
        p(i)=p(i)*stpmax/sum
     end do
  endif
  slope=0.
  do i=1,n
     slope=slope+g(i)*p(i)
  end do
  if (slope.ge.0.) then

     stop 'roundoff problem in lnsrch'
  end if
  test=0.
  do i=1,n
     temp=abs(p(i))/max(abs(xold(i)),1.)
     if(temp.gt.test)test=temp
  end do
  alamin=TOLX/test
  alam=1.

  alam2=2.
  f2=2.
  fold2=fold

1     continue
  do i=1,n
     x(i)=xold(i)+alam*p(i)
  end do
  f=func(x)

  if(alam.lt.alamin)then
     do i=1,n
        x(i)=xold(i)
     end do
     check=.true.
     return
  else if(f.le.fold+ALF*alam*slope)then
     return
  else
     if(alam.eq.1.)then
        tmplam=-slope/(2.*(f-fold-slope))
     else
        rhs1=f-fold-alam*slope
        rhs2=f2-fold-alam2*slope
        a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if(a.eq.0.)then
           tmplam=-slope/(2.*b)
        else
           disc=b*b-3.*a*slope
           if(disc.lt.0.) then
              tmplam=0.5*alam
           else if (b.le.0.) then
              tmplam=(-b+sqrt(disc))/(3.*a)
           else
              tmplam=-slope/(b+sqrt(disc))
           endif
!           C              tmplam=(-b+sqrt(disc))/(3.*a)
        endif
        if(tmplam.gt..5*alam)tmplam=.5*alam
     endif
  endif

  alam2=alam
  f2=f
!  fold2=fold
  alam=max(tmplam,.1*alam)
  goto 1
END SUBROUTINE lnsrch

SUBROUTINE ludcmp(a,n,np,indx,d,itss)
  INTEGER n,np,indx(n),NMAX
  REAL d,a(np,np),TINY
  PARAMETER (NMAX=500,TINY=1.0e-20)
  INTEGER i,imax,j,k
  integer :: its, itss
  REAL aamax,dum,sum,vv(NMAX)
  itss=0
  d=1.
  do i=1,n
     aamax=0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
     end do
     if (aamax.eq.0.) stop 'singular matrix in ludcmp'
     if (aamax.eq.0.) itss=1
     if (aamax.eq.0.) return
     vv(i)=1./aamax
  end do
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
     end do
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     end do
     if (j.ne.imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        end do
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if(a(j,j).eq.0.)a(j,j)=TINY
     if(j.ne.n)then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        end do
     endif
  end do
  return
END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER n,np,indx(n)
  REAL a(np,np),b(n)
  INTEGER i,ii,j,ll
  REAL sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0)then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii=i
     endif
     b(i)=sum
  end do
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do

END SUBROUTINE lubksb


end module wf_mpi
