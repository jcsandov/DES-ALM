rm *.o *.mod
rm makefiles.*

#. /etc/profile.d/modules.sh

module load openmpi/1.6.5/gcc

mpif90 -O3 -c precision.F90 

mpif90 -O3 -c global.F90  

mpif90 -O3 -c global_app.F90  

mpif90 -O3 -c global_mpi.F90  

mpif90 -O3 -c global_param.F90  

mpif90 -O3 -c  global_osg.F90

mpif90 -O3 -c global_turb.F90  

mpif90 -O3 -c checksum.F90  

mpif90 -O3 -c wf_mpi.F90

mpif90 -O3 -c rfg.F90

mpif90 -O3 -c rfg_inlet.F90  

mpif90 -O3 -c  main.F90 

mpif90 -O3 -c init.F90

mpif90 -O3 -c conservation.F90 

mpif90 -O3 -c   mg_driver.F90 

mpif90 -O3 -c mg_nlevel.F90 

mpif90 -O3 -c  	mg_inject.F90  

mpif90 -O3 -c  	mg_metrics.F90  

mpif90 -O3 -c   solver_daf.F90 

mpif90 -O3 -c  	bcond_fm.F90 

mpif90 -O3 -c   mg_brhs.F90 

mpif90 -O3 -c   nonreflect_ibc.F90 

mpif90 -O3 -c   des_eddy.F90 
		
mpif90 -O3 -c  	rsmooth.F90 

mpif90 -O3 -c   sgtsv.f

mpif90 -O3 -c   xerbla.f

mpif90 -O3 -o channel precision.o global.o global_param.o global_app.o \
                global_mpi.o global_osg.o checksum.o wf_mpi.o  global_turb.o \
    rfg.o rfg_inlet.o \
	main.o \
	init.o \
	mg_nlevel.o \
	conservation.o \
	mg_driver.o \
	mg_inject.o \
	mg_metrics.o \
	mg_brhs.o \
	nonreflect_ibc.o \
	solver_daf.o \
	bcond_fm.o \
	des_eddy.o \
        sgtsv.o \
        xerbla.o \
	rsmooth.o \
	-llapack -lblas

mv channel ../

echo Fin!
