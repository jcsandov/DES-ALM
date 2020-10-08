!===================================================================================================
! global_turb:
! ------------
! Defines all the global parameters for the actuator disks.
!***************************************************************************************************
! Also defines some features of turbine blades for ALM (J. Sandoval, 2018)
!***************************************************************************************************
! ..................................................................................................
! length_scale (m)                     = Length scale used in the whole code. Defined as water depth.
! TSR (-)                              = Tip speed ratio, assumed constant for every device.
! vel_infty (m/s)                      = Free stream velocity.
!
!****************************************************************************************************
! n_tur (-)                            = Number of devices beign modelled.
! itur1, itur2 (-)                     = Nodes that define disk center in x-direction.
! jtur1, jtur2 (-)                     = Nodes that define disk center in y-direction.
! ktur1, ktur2 (-)                     = Nodes that define disk center in z-direction.
! tur_rad (m)                          = Array that stores each rotor's radius.
! disk_dat(j,k,1,n_t) (-)              = Distance to disk center.
! disk_dat(j,k,2,n_t) (-)              = Frontal area element.
! disk_dat(j,k,3,n_t) (rad)            = Azimuthal angle between node and disk center.
! disk_dat(j,k,4,n_t) (-)              = Velocity at reference plane i=42.
! disk_dat(j,k,5,n_t) (m)              = Chord length at node.
! disk_dat(j,k,6,n_t) (deg)            = Twist angle at node.
! disk_dat(j,k,7,n_t) (-)              = Interpolation coefficient based on radial distance.
! airfoil_before(j,k,n_t) (-)          = Name of airfoil before node.
! airfoil_after(j,k,n_t) (-)           = Name of airfoil after node.
! disk_node_coordinates(n+1,1,n_t) (-) = i-coordinate for nth node inside n_t disk area.
! disk_node_coordinates(n+1,2,n_t) (-) = j-coordinate for nth node inside n_t disk area.
! disk_node_coordinates(n+1,3,n_t) (-) = k-coordinate for nth node inside n_t disk area.
! **************************************************************************************************
! ALM variables (J. Sandoval, 2018) **
! ************************************
! angle(n_t) (rad)                           = currently angle of reference blade of n_t turbine
! nblades(n_t) (-)                           = Number of blades of each device
! nelems(n_t) (-)                            = Number of elements to discretize the blades
! blade_angles(n_t,bld) (-)                  = angle with respect to the horizontal of each blade in each 
!                                              turbine 
! blade_elemts(n_t,bld,elem,1:3)			       = [xe, ye,ze]: position of elem element of bld blade of
!											                         n_t turbine
! dtheta(n_t) (rad)								           = angle movement of n_t turbine in each time step.
!											                         dtheta = omega*dt
! nagles(n_t)	(-)						                 = number of discretization angles. It asumes that between
!											                         each dtheta, we will discretize 100 elements
! dtheta_interp(n_t) (rad)						       = 2*pi/nangles(n_t). Angle of circle discretization.
!											                         It is used to calculate elems positions in init.F90
! node_indexes(n_t,elem,nang,1:6)            = icell, jcell, kcell: indexes of a blade element at
!                                              nang (angle counter) position and geometrical positions
! 
! forces_lines(n_t,bld,elem,1:6)             = instantaneous position and magnitude of forces applied
!                                              by the actuator lines
!									   
! ..................................................................................................
! Note: disk_node_coordinates(1,1,n_t) =  of nodes inside n_t disk.
! ****************************************************************************************************
! Note: This variable was deleted and replaced by nblades (modified by J. Sandoval, 18/05/2018)
! num_blades (-)                       = Number of blades for every device.
!=====================================================================================================
module global_turb

  use precision
  implicit none

  real (kind = rdf):: length_scale    
  real (kind = rdf):: vel_infty !, num_blades  
  integer :: n_tur
  integer, dimension (:),allocatable ::  itur1, itur2, jtur1, jtur2, ktur1, ktur2, nelems, nblades, nangles
  real (kind = rdf), dimension(:), allocatable :: tur_rad, hub_rad, TSR
  real (kind = rdf), dimension(:), allocatable :: angle, dtheta, dtheta_interp
  real (kind = rdf), dimension(:,:), allocatable :: kernel, blade_angles
  integer, dimension(:,:,:,:), allocatable :: node_indexes
  real (kind = rdf), dimension(:,:,:,:), allocatable :: node_indexes_pos! JS, 2018
  real (kind = rdf), dimension(:,:,:,:), allocatable :: forces_lines ! JS, 2018
  real (kind = rdf), dimension(:,:,:,:), allocatable :: disk_dat, blade_elemts
  real (kind = rdf), dimension(:,:,:,:,:), allocatable :: disk_dat_ker
  character , dimension(:,:,:), allocatable :: airfoil_before*10, airfoil_after*10
  integer (kind = rdf), dimension(:,:,:), allocatable :: disk_node_coordinates
  real (kind = rdf), dimension(:,:,:), allocatable :: disk_node_coordinates_pos

  integer (kind = rdf), dimension(:,:,:), allocatable :: disk_node_coordinates_ker !JS, 2019. ALM_src_v3
  real (kind = rdf), dimension(:,:,:), allocatable :: disk_node_coordinates_pos_ker !JS, 2019. ALM_src_v3

  real, dimension (:), allocatable :: C_T, C_P
  logical :: turbines, node_indexes_load

end module global_turb



