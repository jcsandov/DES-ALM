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
! inter                              = Interpolation coefficient at node.
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
!==================================================================================================================
subroutine rhs_fd ()

  use global_param
  use global_app
  use global_mpi
  use checksum
  use global_turb

  real :: rho
  integer ::   n_t, counter, counter2, number_disk_nodes, i_t
  real :: radius, theta, elm_area, twist, chord, phi, AOA, AUX1, AUX2, dV
  real ::  Cl, Cd, vel_scale, omega, W_rel2, F_prandtl, F_hub
  real :: m, inter, angle1, angle2, Cl1, Cl2, Cd1, Cd2, Cl_b, Cd_b, Cl_a, Cd_a
  real :: u_theta, u_x, u_y, u_z, Cn_2D, Ct_2D, Cn_3D, Ct_3D, dF_theta, dF_x
  real :: dF_cel_theta, dF_cel_x, ker
  real :: accm_T, accm_P, accm_A, sumU, sumV, sumW
  real :: dFX, dFY, dFZ, dT, dP, Vcell
  rho    = 1178.62 
  vel_scale = ren*10E-7/length_scale
  accm_A = 0.0

  do n_t=1,n_tur

    accm_P = 0.0
    accm_T = 0.0

    !------------------------------------------------------
    ! Accumulated frontal area.
    !------------------------------------------------------
    accm_A = pi*(tur_rad(n_t))**2
    !------------------------------------------------------
    ! Obtain number of nodes for evaluation.
    !------------------------------------------------------
    number_disk_nodes = disk_node_coordinates(1,1,n_t)

    !------------------------------------------------------
    ! Angular velocity is obtained supposing constant TSR
    ! and using free stream velocity as a reference.
    !------------------------------------------------------
    omega = (vel_infty*TSR(n_t))/tur_rad(n_t)

    do counter2 = 1,number_disk_nodes
      !------------------------------------------------------
      ! Obtain node information in global coordinates.
      ! Note: variables are dimensionalized.
      !------------------------------------------------------
      i        = disk_node_coordinates(counter2+1,1,n_t)
      j        = disk_node_coordinates(counter2+1,2,n_t)
      k        = disk_node_coordinates(counter2+1,3,n_t)

      !------------------------------------------------------
      ! Check if node is inside local domain
      ! WARNING:
      !   - Code only works if the whole cylinder length is contained in the same processor
      !------------------------------------------------------
      if ( ((i-gi_ia+1).ge.i_mysta.and.(i-gi_ia+1).le.i_myend).and.((j-gi_ja+1).ge.j_mysta.and.(j-gi_ja+1).le.j_myend) &
        .and.((k-gi_ka+1).ge.k_mysta.and.(k-gi_ka+1).le.k_myend)) then

        radius   = disk_dat(j,k,1,n_t)*length_scale  
        elm_area = disk_dat(j,k,2,n_t)*(length_scale**2) 
        theta    = disk_dat(j,k,3,n_t)
        dV       = disk_dat(j,k,4,n_t)*(length_scale**3)
        chord    = disk_dat(j,k,5,n_t)        
        twist    = disk_dat(j,k,6,n_t)
        inter    = disk_dat(j,k,7,n_t)        

        if (radius .gt. hub_rad(n_t)) then

        !------------------------------------------------------
        ! Change to  local coordinates
        !------------------------------------------------------
          !i = i - gi_ia + 1
          j = j - gi_ja + 1
          k = k - gi_ka + 1



      !------------------------------------------------------ 
      ! Volume average for velocity field in every x-plane
      ! at this j,k position
      ! WARNING:
      !   - Code only works if the whole cylinder length is contained in the same processor
      !------------------------------------------------------
        sumU = 0.0
        sumV = 0.0
        sumW = 0.0
        dV   = 0.0
       
        do i_t = itur1(n_t),itur2(n_t)

          i = i_t - gi_ia + 1

          sumU = sumU + q(2,i,j,k)/aj(i,j,k)
          sumV = sumV + q(3,i,j,k)/aj(i,j,k)
          sumW = sumW + q(4,i,j,k)/aj(i,j,k)
          dV = dV + 1.0/aj(i,j,k)

        end do

        u_x = sumU/dV*vel_scale
        u_y = sumV/dV*vel_scale
        u_z = sumW/dV*vel_scale



        !------------------------------------------------------
        ! Calculations for phi and angle of attack.
        ! Check velocity direction and use vel at cell's face
        ! for force calculations.
        !------------------------------------------------------
        ! if ((q(2,i-1,j,k)+2*q(2,i,j,k)+q(2,i+1,j,k))/4.ge.0) then
        !   u_x     = 0.5*(q(2,i-1,j,k)+q(2,i,j,k))*vel_scale
        !   u_y     = 0.5*(q(3,i-1,j,k)+q(3,i,j,k))*vel_scale
        !   u_z     = 0.5*(q(4,i-1,j,k)+q(4,i,j,k))*vel_scale
        ! else
        !   u_x     = 0.5*(q(2,i,j,k)+q(2,i+1,j,k))*vel_scale
        !   u_y     = 0.5*(q(3,i,j,k)+q(3,i+1,j,k))*vel_scale
        !   u_z     = 0.5*(q(4,i,j,k)+q(4,i+1,j,k))*vel_scale
        ! end if

        u_theta = omega*radius+u_y*sin(theta)-u_z*cos(theta) !m/s

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
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'Error calculating angle of attack ir airfoil before (AOA < minimum)'
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            angle2 = -999

          else if (angle2.eq.-999) then
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'Error calculating angle of attack in airfoil before (AOA > maximum)'
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

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
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'Error calculating angle of attack in airfoil after (AOA < minimum)'
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            angle2 = -999
          else if (angle2.eq.-999) then
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'Error calculating angle of attack in airfoil after (AOA > maximum)'
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
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

        !------------------------------------------------------
        ! Projection of lift and drag to theta-x plane.
        !------------------------------------------------------
        Cn_2D = Cl*cos(phi)+Cd*sin(phi)
        Ct_2D = Cl*sin(phi)-Cd*cos(phi)

        !------------------------------------------------------
        ! Prandtl Tip Loss factor: F_prandtl
        !------------------------------------------------------
        AUX1 = num_blades*(1-tur_rad(n_t)/radius)
        AUX2 = 2.0*sin(phi)

        if (AUX2 > tiny(1.0)) then
          F_prandtl = 2.0/pi*acos(exp(AUX1/AUX2))
          if (F_prandtl.lt.0 .or. F_prandtl.gt.1.0) then
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'ERROR in F_prandtl'
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            F_prandtl = 0.0
          end if
        else
          F_prandtl = 1.0
        end if

        ! !------------------------------------------------------
        ! ! Hub Loss factor: F_hub
        ! !------------------------------------------------------
        ! AUX1 = num_blades*(hub_rad(n_t)/radius-1)
        ! AUX2 = 2.0*sin(phi)

        ! if (radius.le.hub_rad(n_t)) then
        !   F_hub = 0.0
        ! else if (AUX2 > tiny(1.0)) then
        !   F_hub = 2.0/pi*acos(exp(AUX1/AUX2))
        !   if (F_hub.lt.0 .or. F_hub.gt.1.0) then
        !     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        !     print*,'ERROR in F_hub'
        !     print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        !     F_hub = 0.0
        !   end if
        ! else
        !   F_hub = 1.0
        ! end if


        !------------------------------------------------------
        ! A.C. says that correction factors are not necessary
        !------------------------------------------------------
        F_prandtl = 1.0
        F_hub = 1.0

        Cn_3D = Cn_2D*F_prandtl*F_hub
        Ct_3D = Ct_2D*F_prandtl*F_hub

        !------------------------------------------------------
        ! Forces in normal and tangential coordinates
        !------------------------------------------------------
        W_rel2       = u_theta**2+u_x**2
        dF_x         = 0.5*chord*rho*Cn_3D*W_rel2 !N/m
        dF_theta     = 0.5*chord*rho*Ct_3D*W_rel2 !N/m

        dF_cel_x     = num_blades*dF_x/(2.0*pi*radius)     !N/m2
        dF_cel_theta = num_blades*dF_theta/(2.0*pi*radius) !N/m2


        !------------------------------------------------------
        ! Forces in fixed coordinates
        ! Note: David checks for 'iflow', which can change direction of normal force.
        ! The sign indicates the direction of the force felt by the blade
        !------------------------------------------------------
        dFX =  dF_cel_x
        dFY = -1.0*dF_cel_theta*sin(theta)
        dFZ =  dF_cel_theta*cos(theta)

        !---------------------------------------------------------------------------
        !SOURCE TERM CALCULATION:
        ! N-S eqs:
        ! 	rho * Du/Dt = ... + F 	[kg/m3 * m/s2 = N / m3 = Force / Vol,cartesian] 
        !		Du/Dt = ... + F/rho 	[m/s2 = Force /(rho * Vol,cartesian)]
        !	
        !***************************************************************************	
        !	In our DES application F/rho must be divided by the jacobian aj(i,k,k)
        !	to correctly compute de N-S eqs, and also non-dimensional
        !***************************************************************************
        ! 
        !   S = (F/rho)*(L/U^2)*(1/aj)
        !
        !		Vol,cartesian [m3] = L^3/aj(i,j,k)
        !		Vol,cartesian [-]  = 1/aj(i,j,k)
        !
        !---------------------------------------------------------------------------

        ! print*,'---------------------------------'
        ! print*,'RHS PRE FUERZA:'
        ! print*,'Disk N:',n_t,' Proc:',myid
        ! print*,'i:',i+gi_ia-1,'j:',j+gi_ja-1,'k:',k+gi_ka-1
        ! print*,'rhs_i: ',rh(2,i,j,k)
        ! print*,'rhs_j: ',rh(3,i,j,k)
        ! print*,'rhs_k: ',rh(4,i,j,k)
        ! print*,'*****************'

        !------------------------------------------------------
        ! Accumulate Torque and Power
        !-----------------------------------------------------

        do i_t = itur1(n_t),itur2(n_t)

          ker = kernel(i_t - itur1(n_t) + 1,n_t)

          i = i_t - gi_ia + 1

          Vcell = length_scale**3/aj(i,j,k)

          dT     = ker*dF_cel_x*Vcell
          dP     = radius*omega*ker*dF_cel_theta*Vcell

          accm_T = accm_T + dT
          accm_P = accm_P + dP

          rh(2,i,j,k) = rh(2,i,j,k) + (ker*dFX)*(length_scale/(rho*vel_scale**2))*(1.0/aj(i,j,k))
          rh(3,i,j,k) = rh(3,i,j,k) + (ker*dFY)*(length_scale/(rho*vel_scale**2))*(1.0/aj(i,j,k))
          rh(4,i,j,k) = rh(4,i,j,k) + (ker*dFZ)*(length_scale/(rho*vel_scale**2))*(1.0/aj(i,j,k))

        end do

        ! print*,'---------------------------------------------------------------'
        ! print*,'MODEL:'
        ! print*,'rho(kg/m3) :',rho
        ! print*,'pi         :',pi
        ! print*,'Nb         :',num_blades
        ! print*,'Lscale(m)  :',length_scale
        ! print*,'Vscale(m/s):',vel_scale,'Uinf(m/s):',vel_infty
        print*,'TSR(-)     :',TSR(n_t)
        ! print*,'Re         :',ren
        ! print*,'Tur_rad',tur_rad(n_t), 'hub_rad',hub_rad(n_t)
        ! print*,'*******************************************************'
        ! print*,'COORDINATES:'
        ! print*,'Disk N:',n_t,' Proc:',myid
        ! print*,'i:',i+gi_ia-1,'j:',j+gi_ja-1,'k:',k+gi_ka-1
        ! print*,'radius(m):',radius,'theta(deg):', theta*180.0/pi
        ! print*,'*******************************************************'
        ! print*,'GEOMETRY:'
        ! print*,'chord(m):',chord,'twist(deg):',twist
        ! print*,'dA[-]:',elm_area/length_scale**2,'dV[-]:',1/aj(i,j,k)
        ! print*,'*******************************************************'
        ! print*,'VELOCITY:'
        ! print*,'ux(m/s):',u_x,'uy(m/s):',u_y,'uz(m/s):',u_z
        ! print*,'W2(m/s):',sqrt(W_rel2)
        ! print*,'Omega(rad/s):',omega,'u_th(m/s):',u_theta
        ! print*,'*******************************************************'
        ! print*,'BLADE PHYSICS:'
        ! print*,'phi(deg)',phi*180.0/pi,'AOA(deg)', AOA
        ! print*,'F_prandtl:',F_prandtl,'F_hub:',F_hub
        ! print*,'Cl:', Cl, 'Cd:', Cd
        ! print*,'Cn2D:',Cn_2D,'Ct2D:',Ct_2D
        ! print*,'Cn3D:',Cn_3D,'Ct3D:',Ct_3D
        ! print*,'Cx:', Cn_3D,'Cy:', -1.0*Ct_3D*sin(theta),'Cz:',Ct_3D*cos(theta)
        ! print*,'*******************************************************'
        ! print*,'FORCES:'
        ! print*,'dF_x        :' ,dF_x         
        ! print*,'dF_theta    : ',dF_theta     
        ! print*,'dF_cel_x    :' ,dF_cel_x     
        ! print*,'dF_cel_theta:' ,dF_cel_theta 
        ! print*,'dFX(kg/m)   :' ,dFX
        ! print*,'dFY(kg/m)   :' ,dFY
        ! print*,'dFZ(kg/m)   :' ,dFZ
        ! print*,'*******************************************************'
        ! print*,'RHS POST FUERZA:'
        ! print*,'Disk N:',n_t,' Proc:',myid
        ! print*,'i:',i+gi_ia-1,'j:',j+gi_ja-1,'k:',k+gi_ka-1
        ! print*,'drhx:  ',(dFX*W_rel2)/(rho*length_scale**2*vel_scale**2)
        ! print*,'drhy:  ',(dFY*W_rel2)/(rho*length_scale**2*vel_scale**2)
        ! print*,'drhz:  ',(dFZ*W_rel2)/(rho*length_scale**2*vel_scale**2)
        ! print*,'rhs_i: ',rh(2,i,j,k)
        ! print*,'rhs_j: ',rh(3,i,j,k)
        ! print*,'rhs_k: ',rh(4,i,j,k)
        ! print*,'---------------------------------------------------------------'          

      end if !if outside hub

      end if !if inside myid domain
    end do  !loop for every point inside the n_t device

    C_P(n_t) = accm_P/(0.5*rho*accm_A*vel_infty**3)
    C_T(n_t) = accm_T/(0.5*rho*accm_A*vel_infty**2)

  end do !loop for every device

  end subroutine rhs_fd







