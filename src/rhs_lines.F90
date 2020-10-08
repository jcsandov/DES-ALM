!=================================================================================================================
! rhs_fd:
! -------
! Calculates the source term associated with the actuator disks, and estimates C_p and C_T.
! It's output is rhs(:,i,j,k) -> p, u, v, w
! .................................................................................................................
! AUX1, AUX2                         = Dummy variables for calculations.
! rho (kg/m3)                        = Water's specific weight. We're using the value of 1178.62 kg/m3
! n_t                                = Number of the actuator disks being evaluated.
! n_tur                              = Number of devices being studied.
! counter                            = Dummy variable. Counts blade point in blade.txt being evaluated.
! counter2                           = Dummy variable. Indicates which point inside the disk is being studied.
! number_disk_nodes                  = Number of nodes inside given actuator disk.
! radius (m)                         = Distance between node and disk center.
! theta (rad)                        = Angle between node and disk center.
! accm_P (W)                         = Accumulated power produced by the farm. Used for estimating C_p.
! accm_T (N)                         = Accumulated thrust produced by the farm. Used for estimating C_T.
! accum_A (m2)                       = Accumulated frontal area of the farm.
! disk_node_coordinates(n+1,1:3,n_t) = (i,j,k) coordinates of the n-th point inside the n_t-th disk.
! i                                  = i-coordinate of the node in study.
! j                                  = j-coordinate of the node in study.
! k                                  = k-coordinate of the node in study.
! radius (m)                         = Radius between disk center and node.
! elm_area (m2)                      = Frontal element of area of node.
! theta    (rad)                     = Azimuthal angle between disk center and node.
! chord    (m)                       = Chord lenght at node.
! twist    (deg)                     = Twist angle at node.
! inter (-)                          = Interpolation coefficient at node.
! relthick (-)                       = Relative thickness of nth blade point in blade.txt.
! omega (rad/s)                      = Angular velocity of the disk.
! phi (rad)                          = Angle between relative velocity and (-theta) direction in theta-x plane.
! AOA (deg)                          = Angle of attack.
! Cl_b                               = Lift coefficient for airfoil before node. Interpolated using node's AOA.
! Cl_a                               = Lift coefficient for airfoil after node. Interpolated using node's AOA.
! Cd_b                               = Drag coefficient for airfoil before node. Interpolated using node's AOA.
! Cd_a                               = Drag coefficient for airfoil after node. Interpolated using node's AOA.
! angle1,angle2,Cl1,Cl2,Cd1,Cd2,m    = Dummy variables for reading airfoil data and interpolating.
! F_prandtl (-)                      = Prandtl's tip-loss correction factor.
! Cl, Cd (-)                         = Node's lift and drag coefficients.
! Cn_2D, Ct_2D (-)                   = Projections of lift and drag to theta-x plane. (n: normal, t:tangential).
! Cn_3D, Ct_3D (-)                   = Lift and drag coefficients in x,y,z plane, they include tip-loss correction.
! .................................................................................................................
! Note: disk_node_coordinates(1,1,n_t) = number of nodes inside n_t disk.
! Note2: J. Sandoval modified init.F90, specifically the routine init_blades, to read the relative thickness
!        and incorporate it as a work variable for lift and drag computation for ALM.
!==================================================================================================================
subroutine rhs_lines ()

  use global_param
  use global_app
  use global_mpi
  use checksum
  use global_turb

  real :: rho, angle_aux
  integer ::   n_t, counter, counter2, number_disk_nodes,number_disk_nodes_ker, i_t, bld, elem, nang !bld, elem & nang from ALM
  real :: radius, theta, elm_area, twist, chord, phi, AOA, AUX1, AUX2, dV, relthick
  real ::  Cl, Cd, vel_scale, omega, W_rel2, F_prandtl, F_hub
  real :: m, inter, angle1, angle2, Cl1, Cl2, Cd1, Cd2, Cl_b, Cd_b, Cl_a, Cd_a
  real :: u_theta, u_x, u_y, u_z, Cn_2D, Ct_2D, Cn_3D, Ct_3D, dF_theta, dF_x
  real :: Urel, fl, fD, fx, fy, fz, epsALM !JS, 2018
  real :: x000,x100,x010,x001,x101,x011,x110,x111,y000,y100,y010,y001,y101,y011,y110,y111
  real :: z000,z100,z010,z001,z101,z011,z110,z111 !variables for trilinear interpolation
  real :: udown,uup,vdown,vup,wdown,wup,uinterp,vinterp,winterp
  real :: xker, yker, zker, fxaccum, fyaccum, fzaccum,xf,yf,zf,xff,yff,zff,fxx,fyy,fzz,rdis,deltar, kerALM!JS, 2018
  integer :: counter_nodes, iker, jker, kker !JS, 2018
  integer :: icyl, jcyl, kcyl !JS, 2018
  real :: dF_cel_theta, dF_cel_x, ker
  real :: accm_T, accm_P, accm_A, sumU, sumV, sumW
  real :: dFX, dFY, dFZ, dT, dP, Vcell, time_aux
  real :: dFXprint, dFYprint, dFZprint ! Auxiliar variables for print the rhs changes
  real, dimension(2, 3) :: dist_search
  real, dimension(:), allocatable :: I_left
  real, dimension(:), allocatable :: I_right
  real, dimension(:), allocatable :: J_left
  real, dimension(:), allocatable :: J_right
  real, dimension(:), allocatable :: K_left
  real, dimension(:), allocatable :: K_right

  rho    = 1178.62
  vel_scale = ren*10E-7/length_scale
  !vel_scale =0.45 !m/s ! ren*10E-6/length_scale
  accm_A = 0.0
  angle_aux = 0.0

  forces_lines=0.0
  xf=0.0

  do n_t=1,n_tur ! do over each declared turbine in the model

    !angle_aux = angle(n_t)
    accm_P = 0.0
    accm_T = 0.0

    !------------------------------------------------------
    ! Accumulated frontal area: Is the circular area
    ! accumulated at specific radius (it is for BEM)
    !------------------------------------------------------
    accm_A = pi*(tur_rad(n_t))**2
    !------------------------------------------------------
    ! Obtain number of nodes for evaluation.
    !------------------------------------------------------
    number_disk_nodes = disk_node_coordinates(1,1,n_t)

    !------------------------------------------------------
    ! Local angular velocity is obtained supposing constant
    ! TSR and using free stream velocity as a reference.
    !------------------------------------------------------
    omega = (vel_infty*TSR(n_t))/tur_rad(n_t)

    do bld = 1,nblades(n_t)

      angle_aux=angle(n_t)+(bld-1)*2*pi/nblades(n_t)

      if(angle_aux.gt.2*pi) then
      	angle_aux=angle_aux-2*pi
      end if

      nang=nint(angle_aux/dtheta_interp(n_t))+1
      !if(myid==root) then
      !  print *, 'bld: ', bld
      !  print *, 'angle: ', angle_aux
      !  print *, '_________________________________'
      !end if
      
      do elem = 1,nelems(n_t)

        i  = node_indexes(n_t,elem,nang,1)
        j  = node_indexes(n_t,elem,nang,2)
        k  = node_indexes(n_t,elem,nang,3)

        xf = node_indexes_pos(n_t,elem,nang,1)
        yf = node_indexes_pos(n_t,elem,nang,2)
        zf = node_indexes_pos(n_t,elem,nang,3)

        x000 = node_indexes_pos(n_t,elem,nang,4) 
        y000 = node_indexes_pos(n_t,elem,nang,5) 
        z000 = node_indexes_pos(n_t,elem,nang,6) 
        x100 = node_indexes_pos(n_t,elem,nang,7) 
        y100 = node_indexes_pos(n_t,elem,nang,8) 
        z100 = node_indexes_pos(n_t,elem,nang,9) 
        x110 = node_indexes_pos(n_t,elem,nang,10)
        y110 = node_indexes_pos(n_t,elem,nang,11)
        z110 = node_indexes_pos(n_t,elem,nang,12)
        x010 = node_indexes_pos(n_t,elem,nang,13)
        y010 = node_indexes_pos(n_t,elem,nang,14)
        z010 = node_indexes_pos(n_t,elem,nang,15)
        x001 = node_indexes_pos(n_t,elem,nang,16)
        y001 = node_indexes_pos(n_t,elem,nang,17)
        z001 = node_indexes_pos(n_t,elem,nang,18)
        x101 = node_indexes_pos(n_t,elem,nang,19)
        y101 = node_indexes_pos(n_t,elem,nang,20)
        z101 = node_indexes_pos(n_t,elem,nang,21)
        x111 = node_indexes_pos(n_t,elem,nang,22)
        y111 = node_indexes_pos(n_t,elem,nang,23)
        z111 = node_indexes_pos(n_t,elem,nang,24)
        x011 = node_indexes_pos(n_t,elem,nang,25)
        y011 = node_indexes_pos(n_t,elem,nang,26)
        z011 = node_indexes_pos(n_t,elem,nang,27)

          !-------------------------------------------------------------------------------------
          ! Check if node is inside local domain
          ! WARNING:
          !   - Code only works if the whole cylinder length is contained in the same processor
          !-------------------------------------------------------------------------------------

        if ( ((i-gi_ia+1).ge.i_mysta.and.(i-gi_ia+1).le.i_myend).and.((j-gi_ja+1).ge.j_mysta.and.(j-gi_ja+1).le.j_myend) &
        .and.((k-gi_ka+1).ge.k_mysta.and.(k-gi_ka+1).le.k_myend)) then

          radius   = disk_dat(j,k,1,n_t)*length_scale
          elm_area = disk_dat(j,k,2,n_t)*(length_scale**2)
          !theta    = disk_dat(j,k,3,n_t)
          dV       = disk_dat(j,k,4,n_t)*(length_scale**3)
          chord    = disk_dat(j,k,5,n_t)
          twist    = disk_dat(j,k,6,n_t)
          inter    = disk_dat(j,k,7,n_t)
          relthick = disk_dat(j,k,8,n_t)

          ! Change to local coordinates
          i = i - gi_ia + 1
          j = j - gi_ja + 1
          k = k - gi_ka + 1
          
          !! Local velocities in [m/s]
          !u_x = q(2,i,j,k)*vel_scale ![m/s]
          !u_y = q(3,i,j,k)*vel_scale ![m/s]
          !u_z = q(4,i,j,k)*vel_scale ![m/s]

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Local interpolation of velocities for actuator points
          ! (JS, 2019. ALM_v2_trilinear_interp)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          udown = 1/((x100-x000)*(y110-y100))*((x100-xf)*((y110-yf)*(q(2,i,j,k)+q(2,i,j+1,k)))+ &
                  (xf-x000)*((yf-y100)*(q(2,i+1,j,k)+q(2,i+1,j+1,k))))

          uup = 1/((x101-x001)*(y111-y101))*((x101-xf)*((y111-yf)*(q(2,i,j,k+1)+q(2,i,j+1,k+1)))+ &
                  (xf-x001)*((yf-y101)*(q(2,i+1,j,k+1)+q(2,i+1,j+1,k+1))))

          uinterp = udown+ (uup-udown)/(z001-z000)*(zf-z000)


          vdown = 1/((x100-x000)*(y110-y100))*((x100-xf)*((y110-yf)*(q(3,i,j,k)+q(3,i,j+1,k)))+ &
                  (xf-x000)*((yf-y100)*(q(3,i+1,j,k)+q(3,i+1,j+1,k))))

          vup = 1/((x101-x001)*(y111-y101))*((x101-xf)*((y111-yf)*(q(3,i,j,k+1)+q(3,i,j+1,k+1)))+ &
                  (xf-x001)*((yf-y101)*(q(3,i+1,j,k+1)+q(3,i+1,j+1,k+1))))

          vinterp = vdown+ (vup-vdown)/(z001-z000)*(zf-z000)


          wdown = 1/((x100-x000)*(y110-y100))*((x100-xf)*((y110-yf)*(q(4,i,j,k)+q(4,i,j+1,k)))+ &
                  (xf-x000)*((yf-y100)*(q(4,i+1,j,k)+q(4,i+1,j+1,k))))

          wup = 1/((x101-x001)*(y111-y101))*((x101-xf)*((y111-yf)*(q(4,i,j,k+1)+q(4,i,j+1,k+1)))+ &
                  (xf-x001)*((yf-y101)*(q(4,i+1,j,k+1)+q(4,i+1,j+1,k+1))))

          winterp = wdown+ (wup-wdown)/(z001-z000)*(zf-z000)


          ! Local velocities in [m/s]
          u_x = uinterp*vel_scale
          u_y = vinterp*vel_scale
          u_z = winterp*vel_scale

          print *,'Velocidades interpoladas:'
          print *,'ux,uy,uz',u_x,u_y,u_z
          print *,'!__________________________________!'
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Geometrical features of the blade at r and theta reading
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Computation of fL and fD
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! Cl and Cd reading from the blades file

            u_theta = omega*radius+u_y*sin(angle_aux)-u_z*cos(angle_aux) !m/s
      
            phi = atan2(u_x,u_theta) !rad
      
            AOA = phi*180.0/pi - twist !deg
      
            !------------------------------------------------------
            ! Interpolation with airfoil data for lift and drag
            ! coefficients.
            !------------------------------------------------------
            open  (unit = 61, FILE = adjustl(adjustr(airfoil_before(j+gi_ja-1,k+gi_ka-1,n_t))//'.txt'), STATUS='OLD')
            angle1 = 1000
            angle2 = 1000
            counter = 1
      
            do while(angle2.ne.-999)
              if (counter.eq.1) then
                read(61,*) angle1, Cl1, Cd1
              else
                angle1 = angle2
                Cl1 = Cl2
                Cd1 = Cd2
              end if
      
              read(61,*) angle2, Cl2, Cd2
      
              if (AOA.le.angle1) then
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                !print*,'Error calculating angle of attack ir airfoil before (AOA < minimum)'
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                angle2 = -999
      
              else if (angle2.eq.-999) then
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                !print*,'Error calculating angle of attack in airfoil before (AOA > maximum)'
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      
              else if ( (AOA <= angle2) .and. (AOA >= angle1) ) then
                m = (Cl2-Cl1)/(angle2-angle1)
                Cl_b = m*(AOA-angle1)+Cl1
                m = (Cd2-Cd1)/(angle2-angle1)
                Cd_b = m*(AOA-angle1)+Cd1
                angle2 = -999
              else
                counter = counter + 1
              end if
            end do
      
            close(unit=61)
      
            open  (unit = 61, FILE = adjustl(adjustr(airfoil_after(j+gi_ja-1,k+gi_ka-1,n_t))//'.txt'), STATUS='OLD')
            angle1 = 1000
            angle2 = 1000
            counter = 1
      
            do while(angle2.ne.-999)
              if (counter.eq.1) then
                read(61,*) angle1, Cl1, Cd1
              else
                angle1 = angle2
                Cl1 = Cl2
                Cd1 = Cd2
              end if
      
              read(61,*) angle2, Cl2, Cd2
              if (AOA.le.angle1) then
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                !print*,'Error calculating angle of attack in airfoil after (AOA < minimum)'
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                angle2 = -999
              else if (angle2.eq.-999) then
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                !print*,'Error calculating angle of attack in airfoil after (AOA > maximum)'
                !print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              else if ( (AOA <= angle2) .and. (AOA >= angle1) ) then
                m = (Cl2-Cl1)/(angle2-angle1)
                Cl_a = m*(AOA-angle1)+Cl1
                m = (Cd2-Cd1)/(angle2-angle1)
                Cd_a = m*(AOA-angle1)+Cd1
                angle2 = -999
              else
               counter = counter + 1
              end if
            end do
      
            close(unit=61)
      
            !------------------------------------------------------
            ! Node's lift and drag coefficients.
            !------------------------------------------------------
            Cl = (Cl_a-Cl_b)*inter+Cl_b
            Cd = (Cd_a-Cd_b)*inter+Cd_b

            ! Computation of fl and fD

            Urel = sqrt(u_theta**2+u_x**2) ![m/s]

            fL=0.5*Cl*rho*Urel**2*chord!*(chord*relthick) ![N/m]
            fD=0.5*Cd*rho*Urel**2*chord!*(chord*relthick) ![N/m]

            fx=fL*cos(phi)+fD*sin(phi) ![N]
            fy=fL*sin(angle_aux)*sin(phi)-fD*sin(angle_aux)*sin(phi) ![N]
            fz=fL*cos(angle_aux)*sin(phi)-fD*cos(angle_aux)*cos(phi) ![N]

            forces_lines(n_t,bld,elem,1) = xf
            forces_lines(n_t,bld,elem,2) = yf
            forces_lines(n_t,bld,elem,3) = zf
            forces_lines(n_t,bld,elem,4) = fx
            forces_lines(n_t,bld,elem,5) = fy
            forces_lines(n_t,bld,elem,6) = fz

            !print *, 'xf = ', xf
            if (xf.gt.0.0) then
            !  print *, '__________________________________'
            !  print *, 'angle_aux: ', angle_aux
            !  print *, 'phi: ', phi
            !  print *, 'blade #: ', bld, 'of ', nblades(n_t)
            !  print *, 'elem #: ', elem, 'of ', nelems(n_t)
            !  print *, 'xf = ', xf, '[-]'
            !  print *, 'u_x: ', u_x, '[m/s]'
            !  print *, 'u_y: ', u_y, '[m/s]'
            !  print *, 'u_z: ', u_z, '[m/s]'
            !  print *, 'u_theta', u_theta, '[m/s]'
            !  print *, 'Urel: ', Urel, '[m/s]'
            !  print *, 'chord length: ', chord, '[m]' 
            !  print *, 'Cl: ', Cl, '[-]' 
            !  print *, 'Cd: ', Cd, '[-]' 
            !  print *, 'fL: ', fL, '[N]'
            !  print *, 'fD: ', fD, '[N]'
            !print *, 'fx = ', forces_lines(n_t,bld,elem,4), '[N]'
            !print *, 'fy = ', forces_lines(n_t,bld,elem,5), '[N]'
            !print *, 'fz = ', forces_lines(n_t,bld,elem,6), '[N]'
            !  print *, 'fx = ', fx, '[N]'
            !  print *, 'fy = ', fy, '[N]'
            !  print *, 'fz = ', fz, '[N]'
            end if


        end if !if inside myid domain
      end do !loop for blade elements
    end do !loop for blades

    number_disk_nodes_ker=disk_node_coordinates_ker(1,1,n_t)

    do counter2=1,number_disk_nodes_ker

      iker = disk_node_coordinates_ker(counter2+1,1,n_t)
      jker = disk_node_coordinates_ker(counter2+1,2,n_t)
      kker = disk_node_coordinates_ker(counter2+1,3,n_t)

      xker = disk_node_coordinates_pos_ker(counter2+1,1,n_t)
      yker = disk_node_coordinates_pos_ker(counter2+1,2,n_t)
      zker = disk_node_coordinates_pos_ker(counter2+1,3,n_t)

      radius   = disk_dat_ker(iker,jker,kker,1,n_t)*length_scale
      elm_area = disk_dat_ker(iker,jker,kker,2,n_t)*(length_scale**2)
      dV       = disk_dat_ker(iker,jker,kker,4,n_t)*(length_scale**3)

      if ( ((iker-gi_ia+1).ge.i_mysta.and.(iker-gi_ia+1).le.i_myend).and.((jker-gi_ja+1).ge.j_mysta &
            .and.(jker-gi_ja+1).le.j_myend).and.((kker-gi_ka+1).ge.k_mysta.and.(kker-gi_ka+1).le.k_myend)) then

      if (radius .gt. hub_rad(n_t)) then
      
        !--------------------------------
        ! Change to  local coordinates
        !--------------------------------
        iker = iker - gi_ia + 1
        jker = jker - gi_ja + 1
        kker = kker - gi_ka + 1
      
        fxaccum = 0.0
        fyaccum = 0.0
        fzaccum = 0.0

        do bld=1,nblades(n_t)
          do elem=1,nelems(n_t)
            xff=forces_lines(n_t,bld,elem,1)
            yff=forces_lines(n_t,bld,elem,2)
            zff=forces_lines(n_t,bld,elem,3)
            fxx=forces_lines(n_t,bld,elem,4)
            fyy=forces_lines(n_t,bld,elem,5)
            fzz=forces_lines(n_t,bld,elem,6)

            !print *, 'xff = ', xff, '[N]'
            !print *, 'yff = ', yff, '[N]'
            !print *, 'zff = ', zff, '[N]'
            !print *, 'fxx = ', fxx, '[N]'
            !print *, 'fyy = ', fyy, '[N]'
            !print *, 'fzz = ', fzz, '[N]'
            !print *, '**********************************************************************'


            rdis = SQRT((xker-xff)**2+(yker-yff)**2+(zker-zff)**2)*length_scale ![m]
            epsALM = 2*(dV)**(1.0/3.0) ! [m]
            kerALM = 1/((epsALM**3)*(pi**(3.0/2.0)))*EXP(-rdis**2/epsALM**2) ! [1/m3]
            deltar = ABS(tur_rad(n_t)-hub_rad(n_t))/nelems(n_t) ! [m]
            fxaccum = fxaccum+kerALM*fxx*deltar
            fyaccum = fyaccum+kerALM*fyy*deltar
            fzaccum = fzaccum+kerALM*fzz*deltar
          end do
        end do
        
        !if(ABS(fxaccum).gt.0.00001) then
        !print *, 'fxaccum = ', fxaccum, '[N]'
        !print *, 'fyaccum = ', fyaccum, '[N]'
        !print *, 'fzaccum = ', fzaccum, '[N]'
        !print *, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'      
        !end if

        dFXprint = fxaccum*(length_scale/(rho*vel_scale**2))*(1.0/aj(iker,jker,kker)) 
        dFYprint = fyaccum*(length_scale/(rho*vel_scale**2))*(1.0/aj(iker,jker,kker)) 
        dFZprint = fzaccum*(length_scale/(rho*vel_scale**2))*(1.0/aj(iker,jker,kker)) 
      
        rh(2,iker,jker,kker) = rh(2,iker,jker,kker)+dFXprint
        rh(3,iker,jker,kker) = rh(3,iker,jker,kker)+dFYprint
        rh(4,iker,jker,kker) = rh(4,iker,jker,kker)+dFZprint

        C_P(n_t) = 1.0! accm_P/(0.5*rho*accm_A*vel_infty**3)
        C_T(n_t) = 1.0! accm_T/(0.5*rho*accm_A*vel_infty**2)

        !if (dFXprint.ne.0.0.or.dFYprint.ne.0.0 .or.dFZprint.ne.0.0) then
        !  print *, 'xff, xker = ', xff, xker, '[-]'
        !  print *, 'dFXprint', dFXprint, '[-]'
        !  print *, 'dFYprint', dFYprint, '[-]'
        !  print *, 'dFZprint', dFZprint, '[-]'
        !  print *, '////////////////////////////////////////////'
        !end if
            !if (myid==root) then
            !   print *, '__________________________________'
            !   print *, 'xff, xker, yff, yker, zff, zker, length_scale', xff, xker, yff, yker, zff, zker, length_scale
            !   print *, 'rdis, epsALM, kerALM, number_disk_nodes: ',rdis, epsALM, kerALM, number_disk_nodes
            !   print *,'The rhs modification aplied on the node: ',iker,jker,kker
            !   print *, 'is', dFXprint, dFYprint, dFZprint
            !   print *, 'and the previous values of rhs were: ', rh(2,i,j,k), rh(3,i,j,k), rh(4,i,j,k)
            !end if
      end if ! if I'm outside the hub
      end if ! inside the processor

    end do ! loop over every cylinder element

  end do !loop for every device

  end subroutine rhs_lines

