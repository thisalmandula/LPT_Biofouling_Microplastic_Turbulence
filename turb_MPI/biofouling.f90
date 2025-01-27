subroutine biofouling2(xdot,t,x,r_pl1,rho_pl1,i,sha,sd,K,deltat,stoc,X3,X2,r_tot,collision,growth,mortality,respiration,turb_dis_rate_dim,tdr,density,rho_tot)

implicit none
integer i
double precision, dimension(365) ::sha
double precision, dimension(2) ::x, xdot
double precision :: r0, rho0, t0, mu0, nu0, na_amb0, Gama, g, comh, X2, X3, z0, kappa, Ca_cell, vol_alg, rho_bf, shear, temp_inf, v_tot, beta_alg_shear_turb, tdr
double precision :: ext_water, light_Im, ext_algae, mu_max, alpha, Q_10, temp_min, temp_opt, temp_max, dD, density, chla_C, turb_dis_rate_dim, mu_tot, rho_tot
double precision :: sec_per_day, mort_alg, resp_alg, light_Iopt, nalg, z, temp, dT, salt, dS, chla, a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, light_Iz, r_alg, beta_alg_shear_t
double precision :: viscosity, a11, b11, viscosity_dyn, t_s, t_f, sha_v, theta, hour, pi, light_I0, mu_opt, na_amb, respiration, dnalg, r_tot,sd,t, r_pl, rho_pl
double precision :: teta_pl, v_pl, v_bf, t_bf, viscosity_kin, d_ast, om_ast, log_d, w_sink, dz, T00, beta_alg_brownian, beta_alg_settling, beta_alg_shear_l,stoc
double precision :: beta_alg_shear_lam, N_2, turb_dis_rate, K, m, t_k, tau, tke, vel, tayscale, Re_ta, T_inf, T_x, sigma_x, Xbar, deltat_s, CA1, CA2, dX,deltat
double precision :: shear_turb_dim, shear_turb, shear_turb_m_dim, shear_turb2, beta_alg_shear_mean, beta_alg_shear_m, collision, growth, mortality, r_pl1, rho_pl1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  DEFINE CHARACTERISTIC SCALES   !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

r0=r_pl1                     	! characteristic scale of all quantities related to plastic size
rho0=1000.                   	! characteristic scale of all densities (water and plastic)
z0=10.                       	! characteristic scale of depth
t0=1.*24.*3600.                	! characteristic scale of time
mu0=1.0e-3                    	! characteristic scale of dynamic viscosity
nu0=1.0e-6                    	! characteristic scale of kinematic viscosity
T00=20.+273.                  	! characteristic scale of kelvin temperature
na_amb0=1.0e-7                	! characteristic scale of the inverse of ambient algae concentration
Gama=0.2                    	! turbulent shear coefficient
g=9.81                      	! gravitational constant
comh = 2.0**(1.0/2.0)           ! mixing coefficient
X2 = X3                      	! order of magnitude of dissipation

!!!!!!!!!!!!!!  INITIALIZE DIMENSIONLESS PARAMETERS   !!!!!!!!!!!!!!!!	

r_pl=r_pl1/r0                	! dimensionless plastic radius
rho_pl=rho_pl1/rho0          	! dimensionless plastic density
g=g*(t0**2./z0)               	! dimensionless gravitational constant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  DEFINE PARAMETERS   !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

kappa = 1.38064852e-23      	! Boltzman constant

!!!!!!!!!!!!!!!!!!!!!!  LIGHT PARAMETERS   !!!!!!!!!!!!!!!!!!!!!!!!!!!

light_Im = 1.2e+8           	! Surface light intensity at noon
ext_water = 0.2             	! Extinction coefficient water
ext_algae = 0.02            	! Extinction coefficient algae

!!!!!!!!!!!!!!!!!!!  ALGAE GROWTH PARAMETERS   !!!!!!!!!!!!!!!!!!!!!!!

mu_max = 1.85               	! Maximum growth rate
alpha = 0.12                	! Initial slope in growth equation
mort_alg = 0.39             	! Mortality rate
resp_alg = 0.1              	! Respiration rate
Q_10 = 2.0                    	! Temperature coefficient for respiration             
temp_min = 0.2              	! Minimum temperature for algae growth
temp_opt = 26.7             	! Optimal temperature for algae growth
temp_max = 33.3             	! Maximum temperature for algae growth
light_Iopt = 1.75392e+13    	! Optimal light intensity for algae growth

!!!!!!!!!!!!!!!!!!!!!!  ALGAE PROPERTIES   !!!!!!!!!!!!!!!!!!!!!!!!!!!

Ca_cell = 2726e-9           	! mg carbon for a single cell
vol_alg = 2.0e-16/r0**3      	! dimensionless volume individual algae
rho_bf = 1388.0/rho0        	! dimensionless biofilm density 
shear = 1.7e+5              	! shear used in encounter rate

!!!!!!!!!!!!!!!  CONVERSION FROM DAY-1 to SEC-1   !!!!!!!!!!!!!!!!!!!!

sec_per_day = 86400.               ! seconds per day
light_Im=light_Im/sec_per_day      ! light intensity per second
mu_max=mu_max/sec_per_day          ! max growth rate per second
alpha=alpha/sec_per_day            ! initial slope in growth eq. per second
mort_alg=mort_alg/sec_per_day      ! mortality rate per second
resp_alg=resp_alg/sec_per_day      ! respiration rate per second
light_Iopt=light_Iopt/sec_per_day  ! opt. lig. int. al. grw. per second
shear=shear/sec_per_day            ! shear in encounter rate per second

!!!!!!!!!!!  CONVERSION TO DIMENSIONLESS PARAMETERS   !!!!!!!!!!!!!!!!

light_Im=light_Im*t0        ! light intensity dimensionalized
mu_max=mu_max*t0            ! max growth rate dimensionalized
alpha=alpha*t0              ! initial slope in growth eq. dimensionalized
mort_alg=mort_alg*t0        ! mortality rate dimensionalized
resp_alg=resp_alg*t0        ! respiration rate dimensionalized
light_Iopt=light_Iopt*t0    ! opt. lig. int. al. grw. dimensionalized
shear=shear*t0              ! shear in encounter rate dimensionalized

!!!!!!!!!!!!!!!!  DEFINE DIFFERENTIAL EQUATIONS   !!!!!!!!!!!!!!!!!!!!

nalg=x(1)                   ! x is the vector of vectorial 
z=x(2)                      ! differential equation x'=f(x)

!!!!!!!!!!!!!!!!!!!!  BOUNDARY CONDITIONS  !!!!!!!!!!!!!!!!!!!!!!!!!!!

if(nalg<0.) then                 ! nalg<0 is not physical
    nalg=0.
end if

if (z>0.) then                  ! plastic can not go above the sea surface
    z=-z
end if

!!!!!!!!!!  READ TEMPERATURE SALINITY CHLOROPHYLL DATA   !!!!!!!!!!!!!

call temperature(z0,z,i,temp,dT)
call salinity(z0,z,i,salt,dS)
call chlorophyll(z0,z,i,chla)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  CALCULATE DENSITY AND VISCOSITY   !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  FROM EQ. OF STATE (SEA WATER)   !!!!!!!!!!!!!!!!!!!
!!!!!  ATTENTION: convert salinity from g/kg (PSU) to kg/kg    !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

salt=salt*1.0e-3
dS = dS*1.0e-3
a1=9.999e+2
a2=2.034e-2
a3=-6.162e-3
a4=2.261e-5
a5=-4.657e-8
b1=8.020e+2
b2=-2.001
b3=1.677e-2
b4=2.261e-5
b5=-4.657e-5

!!!!!!!!!!!!!!!!!!!!  DENSITY CALCULATION  !!!!!!!!!!!!!!!!!!!!!!!!!!!

density=a1+a2*temp+a3*(temp**2.)+a4*(temp**3.)+a5*(temp**4.)+b1*salt+b2*salt*temp+b3*salt*(temp**2.)+b4*salt*(temp**3.)+b5*(salt**2.)*(temp**2.)

density=density/rho0        	! dimensionless sea water density

dD=a2*dT+2.*a3*dT*temp+3.*a4*dT*(temp**2.)+ 4.*a5*dT*(temp**3.)+b1*dS+b2*dS*temp+b2*salt*dT+b3*dS*(temp**2)+ 2.*b3*salt*dT*temp+b4*dS*(temp**3.)+ 3.*b4*salt*dT*(temp**2.)+ 2.*b5*dS*(salt)*(temp**2.)+ 2.*b5*dT*(salt**2.)*(temp)

dD=dD*(z0/rho0)   		! dimensionless sea water density gradient

!!!!!!!!!!!!!!!!!!  VISCOSITY CALCULATION  !!!!!!!!!!!!!!!!!!!!!!!!!!!

viscosity=0.156*(temp+64.993)**2.	
viscosity=4.2844e-5+1.0/(viscosity-91.296)	

!!!!!!!!!!!!!!  DYNAMIC VISCOSITY CALCULATION  !!!!!!!!!!!!!!!!!!!!!!!

a11=1.541+1.998e-2*temp-9.52e-5*(temp**2)
b11=7.974-7.561e-2*temp+4.724e-4*(temp**2)

viscosity_dyn=viscosity*(1.0+a11*salt+b11*(salt**2))/mu0  !  dimensionless

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!  SUN TRACKER   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pi = acos(-1.0)   !  = 3.14159...
t_s = ceiling(t)
t_s = t_s+sd-1.0
if (t_s > 364.0) then
 t_f=modulo(t_s,364.0)
else
 t_f=t_s
end if

sha_v = sha(t_f+1.)

hour=modulo((t*t0)/3600.,24.0)
theta = (hour-(12.-sha_v))*pi/(sha_v*2.)

if (theta < -pi) then
    theta = theta + pi             ! Need for day lengths smaller than 8 hours
end if

if (theta > 2.*pi) then
    theta = theta - pi             ! Need for day lengths smaller than 8 hours
end if

light_I0=light_Im*sin(theta) 

!!!!!!!!!!!!!!!!!  LIGHT BOUNDARY CONDITIONS  !!!!!!!!!!!!!!!!!!!!!!!!

if(light_I0<0.) then
  light_I0 = 0.
end if

!!!!!!!!!!!!!!  LIGHT INTENSITY AT OCEAN DEPTH   !!!!!!!!!!!!!!!!!!!!!

light_Iz=light_I0*exp((ext_water+ext_algae*chla)*z*z0)

!!!!!!!!!!!!! GROWTH RATE UNDER OPT. TEMPERATURE   !!!!!!!!!!!!!!!!!!!

mu_opt=mu_max*light_Iz/(light_Iz+t0*(mu_max/alpha)*(light_Iz/light_Iopt-1.0)**2.)

!!!!!!!!!!!!!!! TEMP. INFLUENCE ON GROWTH RATE   !!!!!!!!!!!!!!!!!!!!!

temp_inf=((temp-temp_max)*(temp-temp_min)**2.)/((temp_opt-temp_min)*((temp_opt-temp_min)*(temp-temp_opt)-(temp_opt-temp_max)*(temp_opt+temp_min-2.*temp)));

!!!!!!!!!!!!!!!!!!!!!!  ALGAE GROWTH RATE  !!!!!!!!!!!!!!!!!!!!!!!!!!!

mu_tot=mu_opt*temp_inf

!!!!!!!!!!!!!! CONVERSION FACTOR: CHLR TO CARBON  !!!!!!!!!!!!!!!!!!!!

chla_C=0.003+0.0154*exp(0.05*temp)*exp(-0.059*light_Iz/t0*1.0e-6*sec_per_day)*mu_tot/mu_max

!!!!!!!!!!!!!!!   AMBIENT ALGAE CONCENTRATION  !!!!!!!!!!!!!!!!!!!!!!!

na_amb=na_amb0*(chla/chla_C)/Ca_cell

!!!!!!!!!!!!!!!!!     BOUNDARY CONDITIONS    !!!!!!!!!!!!!!!!!!!!!!!!!

if (light_Im*exp((ext_water+ext_algae*chla)*z*z0)<=0.01*light_Im) then
     na_amb=0.
end if

!!!!!!!!!!!!!!!!!    PARTICLE TOTAL RADIUS   !!!!!!!!!!!!!!!!!!!!!!!!!

teta_pl=4.*pi*(r_pl**2.)
v_pl=pi*(r_pl**3)*4./3.
v_bf=vol_alg*nalg*teta_pl  
v_tot=v_bf+v_pl
          
t_bf=(v_tot*3./(4.*pi))**(1./3.)-r_pl

if (t_bf<0) then
   t_bf = 0.
end if

r_tot=(v_tot*3./(4.*pi))**(1./3.)
r_alg=(vol_alg*3./(4.*pi))**(1./3.)
         
rho_tot=(v_pl*rho_pl+v_bf*rho_bf)/v_tot

!!!!!!!!!!!!!  DIMENSIONLESS PARTICLE DIAMETER   !!!!!!!!!!!!!!!!!!!!!

viscosity_kin=viscosity_dyn/density  ! dimensionless kinematic viscosity        

d_ast=(rho_tot-density)*9.81*r0**3./nu0**2.*(2.*r_tot)**3.
d_ast=abs(d_ast/(density*viscosity_kin**2.))

!!!!!!!!!!!!!!!!!!!    SINKING VELOCITY   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (d_ast<5.0e-2) then
    om_ast=(d_ast**2.)/5832.
else
    log_d=log10(d_ast)
    om_ast=10.**(-3.76715+1.92944*log_d-0.09815*(log_d**2.)-0.00575*(log_d**3.)+0.00056*(log_d**4.))
end if

w_sink=9.81*nu0*om_ast*viscosity_kin*(rho_tot-density)/density
w_sink=-(w_sink)**(1./3.)*t0/z0      !    dimensionless velocity

!!!!!!!!!!!!!!!!!   BOUNDARY CONDITIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (z*z0>=0. .AND. w_sink>0.) then                        ! plastic cannot fly
    w_sink=0.
end if
              
if (z*z0<=-4000. .AND. w_sink<0.) then                    ! ocean floor set at -4000m 
    w_sink=0.
end if
          
dz=w_sink  

temp=(temp+273.16)/T00                          	! dimensionless kelvin temperature

!!!!!!!!!!!!!!!!!  ENCOUNTER KERNEL RATE   !!!!!!!!!!!!!!!!!!!!!!!!!!

beta_alg_brownian=2./3.*(temp/viscosity_dyn)*(1./r_tot+1./r_alg)*(r_tot+r_alg)      	! dimensionless beta_alg_brownian
beta_alg_settling=0.5*pi*r_tot**2.*abs(w_sink)   					! dimensionless beta_alg_settling
beta_alg_shear_l=4./3.*shear*(r_tot+r_alg)**3.      					! dimensionless beta_shear (Laminar)
beta_alg_shear_lam =beta_alg_shear_l*r0**3.*na_amb/(na_amb0*teta_pl)

!!!!!!!!!!!!!!!!  TURBULENT SHEAR FACTOR   !!!!!!!!!!!!!!!!!!!!!!!!!!

if(K <= 0.) then                 			! Stability requirement
    	K=1e-6
end if

N_2 = (g/density)*abs(dD)                             	! Bouyancy frequency
turb_dis_rate = (K*(t0/z0**2))/Gama*N_2               	! Turbulence dissipation rate dimensionless
m = turb_dis_rate*(z0**2/t0**3)            		! Turbulence dissipation rate dimensional

if(m<1e-10) then                 			! for stability
	m=1e-10
end if

tdr = m						      	! for plotting

t_k = (nu0/m)**(1./2.)                                  ! Kolmogorov time scale   
tau = comh/(N_2**(1./2.)/t0)                            ! Characteristic time scale
tke = K/tau                                           	! Turbulent kinetic energy
vel = (2.*tke)**(1./2.)

tayscale = (10.*(nu0*tke)/m)**(1./2.)                   ! Taylor microscale
Re_ta = abs(vel)*tayscale/nu0                         	! Taylor Scale Reynolds Number
T_inf = (3./20.)**(1./2.)*Re_ta*t_k                     ! Intergral time scale
T_x = T_inf*(0.055+3.55/Re_ta**0.7)                   	! Decorrelation time scale

sigma_x = -0.863 + (3.*0.25/2.)*log(Re_ta)              ! Variance of X
Xbar = -1./(2.*sigma_x**2.);                            ! Mean order of magnitude of dissipation
deltat_s = deltat*t0;

CA1 = -(X2-Xbar)*deltat_s/T_x
CA2 = ((2.*deltat_s*sigma_x**2.)/T_x)**0.5*stoc
dX = CA1 + CA2
X2 = X2 + dX

turb_dis_rate_dim = m*exp(X2)

if(turb_dis_rate_dim<1e-12) then                 			! for stability
	turb_dis_rate_dim=1e-12
end if

!!!!!!!!!!!!!!!!  FLUCTUATING TURBULENCE   !!!!!!!!!!!!!!!!!!!!!!!!!!!

shear_turb_dim = (turb_dis_rate_dim/(viscosity_kin*nu0))**(1./2.)    		! Calculate dimensional turbulence shear rate
shear_turb = t0*shear_turb_dim                                     		! Non dimensionalize turbulence shear rate
beta_alg_shear_t=1.3*shear_turb*(r_tot+r_alg)**3.                     		! dimensionless beta_shear (Turbulent)
beta_alg_shear_turb = (beta_alg_shear_t*r0**3.)*na_amb/(na_amb0*teta_pl)	
!!!!!!!!!!!!!!!!!!!!  MEAN TURBULENCE   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

shear_turb_m_dim = (m/(viscosity_kin*nu0))**(1./2.)     			! Calculate dimensional turbulence shear rate
shear_turb2 = t0*shear_turb_m_dim                                     		! Non dimensionalize turbulence shear rate
beta_alg_shear_m=1.3*shear_turb2*(r_tot+r_alg)**3.                      	! dimensionless beta_shear (Turbulent)
beta_alg_shear_mean = (beta_alg_shear_m*r0**3.)*na_amb/(na_amb0*teta_pl)

!!!!!!!!!!!  R.H.S VALUES. FOR ALGAE GROWTH EQU.   !!!!!!!!!!!!!!!!!!!

! Important!!!! dimensionless collision term (change the shear term here based on mean [beta_alg_shear_m] or fluctuating [beta_alg_shear_t])
! If you are using fluctuating turbulence term, make sure to use a time step of 1 second or lower.

collision=(beta_alg_brownian*(kappa*T00*t0)/mu0+beta_alg_settling*z0*r0**2.+beta_alg_shear_m*r0**3.)*na_amb/(na_amb0*teta_pl) 	

growth=mu_tot*nalg                             				! dimensionless growth term

mortality=mort_alg*nalg                        				! dimensionless mortality term

respiration=resp_alg*nalg*Q_10**((temp*T00-20.-273.16)*0.1) 		! dimensionless respiration term

dnalg=collision+growth-mortality-respiration   				! dimensionless dA/dt

xdot=[dnalg,dz]                                				! new value for x'=f(x)

end subroutine biofouling2

