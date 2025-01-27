subroutine salinity(z0,z,i,salt,dS)

implicit none
integer i
double precision,dimension(3) :: Sfix, b1, b2, b3, b4, b5, b6 
double precision salt, dS, z, zfix, z0

z = z*z0

zfix=-1000.
Sfix=[34.6, 34.6, 34.5]
b1=[9.9979979767e-17, -1.2695423698e-17, -2.4306209411e-16]
b2=[1.0536246487e-12, -6.3765898788e-14, -2.3825929024e-12]
b3=[3.9968286066e-09, 1.2655205499e-10, -8.0374560926e-09]
b4=[6.5411526250e-06, 1.0818886978e-06, -1.0613797323e-05]
b5=[4.1954014008e-03, 1.5454960921e-03, -4.0153966208e-03]
b6=[3.5172984035e+01, 3.4995207081e+01, 3.4908138893e+01]

if (z > zfix) then
    salt=b1(i)*z**5+b2(i)*z**4+b3(i)*z**3+b4(i)*z**2+b5(i)*z+b6(i)
    dS = 5*b1(i)*z**4+4*b2(i)*z**3+3*b3(i)*z**2+2*b4(i)*z+b5(i)
else 
    salt=Sfix(i)
    dS = 0
end if

dS = -dS
z = z/z0

end subroutine salinity

