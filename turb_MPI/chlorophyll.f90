subroutine chlorophyll(z0,z,i,chla)

implicit none
integer i
real(8),dimension(3) :: Cab, chlazbase, s, Camax, zmax, Dz
real(8) chla, z, z0

z = z*z0

chlazbase=[0.151,0.185,2.95]
Cab=[0.533,0.428,0.188]
s=[1.72e-3,1.67e-5,0.0];
Camax=[1.194,1.015,0.885];
zmax=[92.01,82.360,9.870];
Dz=[43.46,57.330,28.210];

chla=Cab(i)-s(i)*-z+Camax(i)*exp(-((-z-zmax(i))/Dz(i))**2);
chla=chla*chlazbase(i);
z = z/z0

end subroutine chlorophyll

