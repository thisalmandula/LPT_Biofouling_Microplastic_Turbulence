program main_code

implicit none

integer :: n, p, i, j, proc, ierr, master, tlim,tot,l,nproc,numpar,ierror,tx,m,ts_print,timestep_day,par_proc
double precision :: t_init, w_sink_max, t, itmax, deltat, r_pl, d_pl, rho_pl, sd, lat, timestep, t0,turb_dis_rate_dim,tdr,beta_alg_shear_mean,beta_alg_shear_turb, z0,r_pl1,rho_pl1,xfinal,x_proc_sum,beta_rng_stoc,gauss_rng_stoc,r0,error,collision,growth,mortality,respiration,collision_avg,growth_avg,mortality_avg,respiration_avg,nalg_avg,nalgfinal
double precision :: K,pdt,kappa,u_star_w,phi_sf,theta_lg,MLD,u_star_a,u_10_wind,g,beta_wa,X3,X2,dK,z_0,gturb,ran,stoc,y,var,tt,xfinal_sum,xfinal_avg,r_tot,d_pl2,turbshrmean_proc_sum,turbshrmean_sum, lastday,max_stoc,gturb2,collision_proc_sum,growth_proc_sum,mortality_proc_sum,respiration_proc_sum,nalg_proc_sum,nalg_sum,turbshr_proc_sum,turbshr_sum,turbshr_avg,turbshrfinal,turbshrmean_avg,turbshrmeanfinal, dia_inc,dia_bins(101),pdf_dia(100),depth_bins(21),pdf_depth(20),dia_bins_mid(100),depth_bins_mid(20),min_dia,max_dia,pdf_par_sum ,z_end,min_depth,max_depth,pdf_depth_sum,depth_inc,pdt2,collision_sum,growth_sum,mortality_sum,respiration_sum,collisionfinal,growthfinal,mortalityfinal,respirationfinal
double precision, dimension(365) :: sha
double precision, allocatable, dimension(:) :: beta_rng,gauss_rng,z_end_proc,d_pl_proc
double precision, allocatable, dimension(:,:) :: xfinaltot,beta_rng3,gauss_rng3, x_proc,d_pl_tot,z_end_tot,d_pl_proc2,z_end_proc2
double precision, dimension(2) :: x, xdot
include "mpif.h"  ! define MIP constants
integer status(MPI_STATUS_SIZE)  ! size defined in mpif.h
character(len=20) :: filename,filename2,filename3
real(8) :: start_time, end_time, total_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  DEFINE MAIN LOOP FOR ITERATION   !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, proc, ierr)		! define current process id
call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)		! define number of processes

if(proc .eq. master) then      ! do following only on master ...

write(*,*) "The simulation has started"

end if

K = 0.                  ! Diffusivity constant
pdt = 0.1               ! Random mixing layer limit top
kappa = 0.4             ! Von Karman constant
u_star_w = 0.0029       ! Friction velocity near sea surface
phi_sf = 0.9            ! Stability function
theta_lg = 1.           ! Langmuir circulation enchancement factor
MLD = 100.              ! Mixed layer depth
u_star_a = 0.4          ! Friction velocity of air
u_10_wind = 10.         ! Wind speed at 10 m
g = 9.81                ! Gravitational acceleration
beta_wa = 35.           ! Assume constant wave age (cp/u_star_a)
X3 = 1.                 ! Order of magnitude of dissipation
var=1.			! Variance

call cpu_time(start_time)

open(unit=10,file='param.in',status='old')
       	read(10,*) d_pl
       	read(10,*) rho_pl1
       	read(10,*) sd
       	read(10,*) lat
       	read(10,*) timestep
	read(10,*) tlim
	read(10,*) i
	read(10,*) numpar
	read(10,*) ts_print
	close(10)

r_pl1 = d_pl/2.0				!calculate particle radius
r0=r_pl1                     			! characteristic scale of all quantities related to plastic size
t_init=0.
w_sink_max=0.
par_proc = numpar/nproc				!divide number of particles equally among processors. 
tlim=tlim*24.*3600. 				!conversion from days to seconds
timestep_day = 24.*3600./timestep

tx = tlim/timestep

! NOTE: Change lastday to the starting timestep of your lastday			
lastday = 85536		
			
allocate(x_proc(par_proc,2))
allocate(beta_rng(tx))
allocate(gauss_rng(tx))
allocate(beta_rng3(tx,par_proc))
allocate(gauss_rng3(tx,par_proc))
allocate(z_end_proc(par_proc))
allocate(d_pl_proc(par_proc))
allocate(z_end_tot(par_proc,nproc))
allocate(d_pl_tot(par_proc,nproc))
allocate(d_pl_proc2(par_proc,timestep_day))
allocate(z_end_proc2(par_proc,timestep_day))

t0=1.0*24.0*3600.0 				!characteristic time scale 
z0=10.0 					!characteristic space scale
x_proc(:,:) = 0.0				!Initialize a zero matrix for x_proc
d_pl_proc(:) = 0.0				!Initialize a zero matrix for d_pl_proc
z_end_proc(:) = 0.0				!Initialize a zero matrix for d_pl_proc
z_end_tot(:,:) = 0.0				!Initialize a zero matrix for z_end_tot
d_pl_tot(:,:) = 0.0				!Initialize a zero matrix for d_pl_tot
z_end_proc2(:,:) = 0.0				!Initialize a zero matrix for z_end_tot
d_pl_proc2(:,:) = 0.0				!Initialize a zero matrix for d_pl_tot

deltat=timestep/t0 				!Define non-dimensional timestep
tlim=tlim/t0 					!dimensionless limit time  						!LHS Initializing
t=0.0

itmax=NINT(tlim/deltat)      	 		!calculation of total timesteps required

do m = 1,par_proc

call random_seed()				!create a random seed for every timestep
call random_number(beta_rng)			!generate a random number between 0 and 1

call gauss_rng2(gauss_rng,itmax)

gauss_rng3(:,m) = gauss_rng
beta_rng3(:,m) = beta_rng

end do

max_stoc = minval(gauss_rng3)
y = 0
error = 0.02
z_0 = (3.5153e-5)*((beta_wa*u_star_a/u_10_wind)**(-0.42))*(u_10_wind**2/g) 	! Roughness scale

do while (error >= 0.001)		! Determine random mixed layer depth at the ocean surface for the boundary condition
dK =  ((kappa*u_star_w*theta_lg/phi_sf)*((1-((y*z0)/MLD))**2)+(-2/MLD)*(kappa*u_star_w*theta_lg/phi_sf)*((y*z0)+z_0)*((1-((y*z0)/MLD))))*(t0/z0) ! Diffusivity gradient
K =  (kappa*u_star_w*theta_lg/phi_sf)*(((y*z0))+z_0)*((1-(((y*z0))/MLD))**2)	! Diffusivity
gturb2 = sqrt(2*(K*(t0/z0**2))/(deltat*var))*max_stoc            ! Stochastic term calculation 
pdt2 = - dK*deltat + gturb2*deltat
error = abs(y-abs(pdt2))
y = abs(pdt2)
end do

call solar_tracker(sha,lat)     		!solar hour angles
									
dK =  ((kappa*u_star_w*theta_lg/phi_sf)*((1-((0*z0)/MLD))**2)+(-2/MLD)*(kappa*u_star_w*theta_lg/phi_sf)*((0*z0)+z_0)*((1-((0*z0)/MLD))))*(t0/z0)	! Diffusivity gradient (starting value at sea surface)
K =  (kappa*u_star_w*theta_lg/phi_sf)*(((0*z0)+0.5*dK*(z0/t0)*deltat*t0)+z_0)*((1-(((0*z0)+0.5*dK*(z0/t0)*deltat*t0)/MLD))**2)				! Diffusivity (starting value at sea surface)

if(proc .eq. master) then 

filename='results/par.out'

open(unit=110,file=filename,status='replace', action='write', iostat=l)	! Open the file for writing
write(110,*) 'Variables="time","z","collision","growth","mortality","respiration","nalg_turbshear","nalg_turbshear_mean","nalg"'
end if

do j=1,tx

x_proc_sum = 0.
collision_proc_sum = 0.
growth_proc_sum = 0.
mortality_proc_sum = 0.
respiration_proc_sum = 0.
nalg_proc_sum = 0.
turbshr_proc_sum = 0.
turbshrmean_proc_sum = 0.	

	do m = 1,par_proc	! loop over number of processes	

	gauss_rng_stoc = gauss_rng3(j,m)
	beta_rng_stoc = beta_rng3(j,m) 

	x(1) = x_proc(m,1)					! Allocate value for nalg
	x(2) = x_proc(m,2)					! Allocate value for z

	stoc=gauss_rng_stoc;                       	! Call the gaussian rng value 

	call biofouling2(xdot,t,x,r_pl1,rho_pl1,i,sha,sd,K,deltat,stoc,X3,X2,r_tot,collision,growth,mortality,respiration,turb_dis_rate_dim,tdr,beta_alg_shear_mean,beta_alg_shear_turb)
	
	X3 = X2

	if (x(2)>-MLD/z0) then                ! Setup diffusivity values for mixed layer depth location of particle
    		y = abs(x(2))

		dK =  ((kappa*u_star_w*theta_lg/phi_sf)*((1-((y*z0)/MLD))**2)+(-2/MLD)*(kappa*u_star_w*theta_lg/phi_sf)*((y*z0)+z_0)*((1-((y*z0)/MLD))))*(t0/z0) ! Diffusivity gradient
		K =  (kappa*u_star_w*theta_lg/phi_sf)*(((y*z0))+z_0)*((1-(((y*z0))/MLD))**2)	! Diffusivity
	else
		dK = 0
		K = 0
	end if

	gturb=sqrt(2*(K*(t0/z0**2))/(deltat*var))*stoc               ! Stochastic term calculation 

	x(1) = x(1) + xdot(1)*deltat                 ! Make sure turbulence does not affect number of algae
        x(2) = x(2) - dK*deltat + gturb*deltat + xdot(2)*deltat   ! Find new z depth with turbulence      

	ran = beta_rng_stoc              ! Random value from 0 to 1 for boundary conditions

	if(x(1)<0.) then                 ! nalg<0 is not physical
    		x(1)=0.
	end if

	if (x(2)>pdt2) then                ! plastic can not go above the sea surface
    		x(2)=ran*pdt2
	end if
	
	turbshrfinal = beta_alg_shear_turb/(t0*r0**2)
	turbshrmeanfinal = beta_alg_shear_mean/(t0*r0**2)
	xfinal = x(2)*z0	
	z_end_proc(m) = x(2)*z0			! Calculate final particle location
	nalgfinal = x(1)/(t0*r0**2)
	collisionfinal = collision/(t0*r0**2)
	growthfinal = growth/(t0*r0**2)
	mortalityfinal = mortality/(t0*r0**2)
	respirationfinal = respiration/(t0*r0**2)
	x_proc(m,1)=x(1)		! carry on nalg to next timestep using global variable
	x_proc(m,2)=x(2)		! carry on z to next timestep using global variable
	x_proc_sum=x_proc_sum+xfinal 	! sum all particle locations in a processor
	collision_proc_sum = collision_proc_sum+collisionfinal
	growth_proc_sum = growth_proc_sum+growthfinal
	mortality_proc_sum = mortality_proc_sum+mortalityfinal
	respiration_proc_sum = respiration_proc_sum+respirationfinal
	nalg_proc_sum=nalg_proc_sum+nalgfinal
	turbshr_proc_sum=turbshr_proc_sum+turbshrfinal
	turbshrmean_proc_sum=turbshrmean_proc_sum+turbshrmeanfinal
	
	if(j > lastday) then 
	d_pl_proc2(m,j-lastday) = r_tot*r0*2		! Calculate final particle diameter	
	!Activate z_end_proc2 only if you want average particle location for the last day (In the paper its for particles of size 1e-3 m)
	!z_end_proc2(m,j-lastday) = x(2)*z0			! Calculate final particle location
	end if
	
	end do
	!write(*,*) "d_pl_proc is", d_pl_proc2
	t=t+deltat
	!write(*,*) "collision_proc_sum is", respiration_sum
	call MPI_REDUCE(x_proc_sum, xfinal_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(collision_proc_sum, collision_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(growth_proc_sum, growth_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(mortality_proc_sum, mortality_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(respiration_proc_sum, respiration_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(nalg_proc_sum, nalg_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(turbshr_proc_sum, turbshr_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	call MPI_REDUCE(turbshrmean_proc_sum, turbshrmean_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
	!write(*,*) "x_proc_sum is", x_proc_sum
	if(proc .eq. master) then 
	xfinal_avg = xfinal_sum/numpar
	collision_avg = collision_sum/numpar
	growth_avg = growth_sum/numpar
	mortality_avg = mortality_sum/numpar
	respiration_avg = respiration_sum/numpar
	nalg_avg = nalg_sum/numpar
	turbshr_avg = turbshr_sum/numpar
	turbshrmean_avg = turbshrmean_sum/numpar
	if (mod(j, ts_print) == 0) then
	write(110, 79) t,xfinal_avg,collision_avg,growth_avg,mortality_avg,respiration_avg,turbshr_avg,turbshrmean_avg,nalg_avg 	! The * means that data is written in free format
	end if
	end if
end do

do m = 1,par_proc
	do j = 1,timestep_day
	d_pl_proc(m) = d_pl_proc(m) + d_pl_proc2(m,j)
	!Activate z_end_proc only if you want average particle location for the last day (In the paper its for particles of size 1e-3 m)
	!z_end_proc(m) = z_end_proc(m) + z_end_proc2(m,j)
	end do
	d_pl_proc(m) = d_pl_proc(m)/timestep_day
	!Activate z_end_proc only if you want average particle location for the last day (In the paper its for particles of size 1e-3 m)
	!z_end_proc(m) = z_end_proc(m)/timestep_day
end do	

if(proc .eq. master) then 
 close(110)
79 format(3x,9(f25.15,3x))
 end if

call MPI_Send(d_pl_proc, par_proc, MPI_DOUBLE_PRECISION,master,2, MPI_COMM_WORLD, ierr)	! Send a message to all processors to start collecting particle motion data
call MPI_Send(z_end_proc, par_proc, MPI_DOUBLE_PRECISION,master,3, MPI_COMM_WORLD, ierr)	! Send a message to all processors to start collecting particle motion data

if(proc .eq. master) then      ! do following only on master ...

write(*,*) "Collecting results from all processors"
d_pl_tot(:,1) = d_pl_proc		! Collect final particle diameters from master processor
z_end_tot(:,1) = abs(z_end_proc)	! Collect final particle diameters from master processor

do proc=1,nproc-1   ! loop on processors to collect local sum
call MPI_RECV(d_pl_proc, par_proc, MPI_DOUBLE_PRECISION, proc, 2, MPI_COMM_WORLD, status, ierr)	! Send data from all processors to master processor for data collection
call MPI_RECV(z_end_proc, par_proc, MPI_DOUBLE_PRECISION, proc, 3, MPI_COMM_WORLD, status, ierr)	! Send data from all processors to master processor for data collection
d_pl_tot(:,proc+1) = d_pl_proc		! Allocate particle diameters from each processor to one matrix
z_end_tot(:,proc+1) = abs(z_end_proc)	! Allocate particle final location from each processor to one matrix
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!  PDF CALCULATION   !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

min_dia = minval(d_pl_tot)		! Calculate minimum value of particle diameters
max_dia = maxval(d_pl_tot)		! Calculate maximum value of particle diameters
dia_inc = (max_dia-min_dia)/100		! Calculate increments to particle diameters
dia_bins(1) = min_dia			! Allocate a bin to collect particle numbers for each diameter range, set first value
dia_bins(101) = max_dia			! Allocate a bin to collect particle numbers for each diameter range, set final value

do j = 1,99
dia_bins(j+1) = dia_bins(j) + dia_inc
end do

do m = 1,par_proc			! Allocate a particle diameters to specific bins
 	do j = 1, nproc
		do p = 2, 101
			if ((d_pl_tot(m,j).gt.dia_bins(p-1)).and.(d_pl_tot(m,j).lt.dia_bins(p))) then
			pdf_dia(p-1) = pdf_dia(p-1) + 1
			end if
		end do
	end do
end do

dia_bins_mid(1) = min_dia + dia_inc/2		! Allocate a new bin with mid point values of data ranges

do j = 2,100
dia_bins_mid(j) = dia_bins_mid(j-1) + dia_inc
end do

do j = 1,100
pdf_dia(j) = pdf_dia(j)/(numpar*dia_inc)	! Calculate PDF
end do

pdf_par_sum = sum(pdf_dia)			! Sum total number of particles to see if all particles accounted for

min_depth = 0					! Calculate minimum value of particle depth
max_depth = MLD					! Calculate maximum value of particle depth
depth_inc = (max_depth-min_depth)/20		! Calculate increments to particle depth
depth_bins(1) = min_depth			! Allocate a bin to collect particle numbers for each depth range, set first value
depth_bins(21) = max_depth			! Allocate a bin to collect particle numbers for each depth range, set final value

do j = 1,19
depth_bins(j+1) = depth_bins(j) + depth_inc	! Allocate bins for particle depths
end do

do m = 1,par_proc				! Allocate particles to bins of specific depths
 	do j = 1, nproc
		do p = 2, 21
			if ((z_end_tot(m,j).gt.depth_bins(p-1)).and.(z_end_tot(m,j).lt.depth_bins(p))) then
			pdf_depth(p-1) = pdf_depth(p-1) + 1
			end if
		end do
	end do
end do

pdf_depth_sum = sum(pdf_depth)

depth_bins_mid(1) = min_depth + depth_inc/2

do j = 2,20
depth_bins_mid(j) = depth_bins_mid(j-1) + depth_inc	! Allocate a new bin with mid point of all depths
end do

do j = 1,20
pdf_depth(j) = pdf_depth(j)/(numpar*depth_inc)		! Calculate PDF
end do

pdf_depth_sum = sum(pdf_depth)

filename2='results/par_dia.out'			! Dedicate a file name to write data 
filename3='results/par_depth.out'		! Dedicate a file name to write data 

open(unit=110,file=filename2,status='replace', action='write', iostat=l)	! Open the file for writing
write(10, *) '"diameter","pdf"'					! Write data to the file
 do j = 1, 100
      	write(110, *) dia_bins_mid(j),pdf_dia(j)		! The * means that data is written in free format
 end do
 close(110)

open(unit=110,file=filename3,status='replace', action='write', iostat=l)	! Open the file for writing
write(10, *) '"pdf","depth"'					! Write data to the file
 do j = 1, 20
      	write(110, *) pdf_depth(j),depth_bins_mid(j)		! The * means that data is written in free format
 end do
 close(110)

end if

call MPI_Finalize(ierr)
call cpu_time(end_time)

total_time = end_time - start_time
write(*,*) 'Total simulation time:', total_time, 'seconds'

end program main_code
     
