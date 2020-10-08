
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

  real :: cy, cz, radius, dA, LX, Xcent, SIGMA
  integer :: posx, ref_vel, counter, counter2, i
  real :: rad1, rad2, aux1, aux2, m

  real (kind = rdf), dimension(:,:,:), allocatable :: xd
  real (kind = rdf), dimension(:,:,:), allocatable :: zd
  real (kind = rdf), dimension(:,:,:), allocatable :: yd
  real (kind = rdf), dimension(:,:), allocatable :: blade_geom

  character airfoil_name(30)*10

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

    !----------------------------------------------------------
    ! Geometry of disk center
    !----------------------------------------------------------
    Xcent = 0.5*(xd(itur1(n_t),jtur1(n_t),ktur1(n_t))+xd(itur2(n_t),jtur1(n_t),ktur1(n_t)))
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

    	kernel(i-itur1(n_t)+1,n_t) = exp(-0.5*((xd(i,jtur1(n_t),ktur1(n_t))-Xcent)/SIGMA)**2)/sqrt(2*pi*(SIGMA*length_scale)**2)

      ! print*,'i:',i,'kernel',kernel(i-itur1(n_t),n_t)
    end do

    !----------------------------------------------------------
    !Since every x-plane has the same j,k distribution we can use itur1 for defining
    !every other coordinate
    !----------------------------------------------------------
    posx = itur1(n_t)
    cy   = (yd(posx,jtur1(n_t),ktur1(n_t))+yd(posx,jtur2(n_t),ktur1(n_t)))/2
    cz   = (zd(posx,jtur1(n_t),ktur1(n_t))+zd(posx,jtur1(n_t),ktur2(n_t)))/2

    do j=1,jmg(1,nz)
      do k=1,kmg(1,nz)
        radius = sqrt((yd(posx,j,k)-cy)**2+(zd(posx,j,k)-cz)**2)

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

          disk_dat(j,k,1,n_t) = radius

          !frontal area element
          disk_dat(j,k,2,n_t) = ABS(yd(posx,j-1,k)-yd(posx,j+1,k))*ABS(zd(posx,j,k-1)-zd(posx,j,k+1))/4

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

          !disk_dat(j,k,4,n_t)= dA*ABS(xd(posx+1,j,k)-xd(posx-1,j,k))/2
        	disk_dat(j,k,4,n_t) = 0.0

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

              !Interpolation for chord length
              aux1                    = blade_geom(counter,3)
              aux2                    = blade_geom(counter+1,3)
              m                       = (aux2-aux1)/(rad2-rad1)
              disk_dat(j,k,5,n_t)     = m*(radius-rad1)+aux1

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
	!		print*,'DEBUG INSIDE R1<R<�R2'
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
!		print*,'Air aft:�',airfoil_after(j,k,n_t)
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

  deallocate (xd, yd, zd, blade_geom)

	if (myid == root) then
  		print*,"Finished initializing actuator disks."
 		do n_t=1,n_tur
  			print*,'Disk N�',n_t, ' has ', disk_node_coordinates(1,1,n_t), ' points.'
  		end do
  		print*,'--------------------------------------------------------'
	end if

end subroutine init_lines


!===================================================================================================
!init_blades:
!----------
! Initializates blades variables.
! ..................................................................................................
! n_elem                 = number of elements in each blade
! theta_angles           = array that contains the angle of the position of each blade
! blade_elemts           = array that contains the i,j,k coordinate of each of the n_elem blade elements per blade.
!
!
! ..................................................................................................
!===================================================================================================

subroutine init_blades
  integer, parameter :: n_elem = 10;    
  real :: Pi = 4.0*ATAN(1.0);
  real, dimension(1:3) :: theta_angles
  real, dimension(0:n_elem,3,1:3):: blade_elemts;
  !real :: tur_rad = 10.0, hub_rad = 2.0;
  !real :: omega = 2.0, dt = 0.5
  theta_angles =  [0.0, -(2.0/3.0)*Pi, -(4.0/3.0)*Pi];

  !!!Init blades!!!

  do i = 0, n_elem, 1
    blade_elemts(i, 1, 1) = 0.0;                                              !Componente i
    blade_elemts(i, 2, 1) = 0.0;                                              !Componente j
    blade_elemts(i, 3, 1) = 2.0 + i*(tur_rad - hub_rad)/n_elem                !Componente k
  end do

  do i = 2, 3, 1
    do j = 0, n_elem, 1
      blade_elemts(j, 1, i) = 0.0;
      blade_elemts(j, 2, i) = cos(theta_angles(i))*blade_elemts(j, 2, 1) - sin(theta_angles(i))*blade_elemts(j, 3, 1)
      blade_elemts(j, 3, i) = sin(theta_angles(i))*blade_elemts(j, 2, 1) + cos(theta_angles(i))*blade_elemts(j, 3, 1)
    end do
  end do

  !!!Init blades!!!
  do i = 1, 3, 1
    theta_angles(i) = theta_angles(i) + omega*dt
  end do

  do i = 0, n_elem, 1
    blade_elemts(i, 1, 1) = 0.0;                                                        !Componente i
    blade_elemts(i, 2, 1) = blade_elemts(i,3,1)*sin(theta_angles(1));                   !Componente j
    blade_elemts(i, 3, 1) = blade_elemts(i,3,1)*cos(theta_angles(1))                    !Componente k
  end do


end subroutine init_blades
