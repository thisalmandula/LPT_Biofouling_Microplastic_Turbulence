subroutine temperature(z0,z,i,temp,dT)

implicit none
integer i
double precision,dimension(3) :: zc, Tsurf, Tbot, p 
double precision temp, dT, z, z0

z = z*z0

p=[2., 1., 2.]
zc=[-300., -400., -1200.]
Tsurf=[25., 16., 8.]
Tbot=[1.5, 1., 2.]

temp=Tsurf(i)+(Tbot(i)-Tsurf(i))*z**p(i)/(z**p(i)+zc(i)**p(i));
dT = (Tbot(i)-Tsurf(i))*(p(i)*z**(p(i)-1)*zc(i)**p(i))/((z**p(i)+zc(i)**p(i))**2);
dT = -dT;
z = z/z0
end subroutine temperature

