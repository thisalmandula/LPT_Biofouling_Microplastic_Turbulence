subroutine solar_tracker(sha,lat)

implicit none
integer i
double precision,dimension(365) :: decl_ang, sha
double precision :: pi, lat, B
real, parameter :: deg2rad = 0.01745329251994329576923690768489 ! Conversion factor
pi = acos(-1.0)   			!  = 3.14159...


do i=1,365
B = 2.*pi*(i-1.)*(1./365.)
decl_ang(i) = (0.006918-0.399912*cos(B)+0.070257*sin(B)-0.006758*cos(2.*B)+0.000907*sin(2.*B)-0.002697*cos(3.*B)+0.00148*sin(3.*B))/deg2rad
end do
!write(*,*) "The value of B is ",decl_ang

sha = acos(-tan(deg2rad*lat)*tan(deg2rad*decl_ang))
sha = sha/deg2rad
sha = sha*4./60.
!write(*,*) "The value of B is ",sha

end subroutine solar_tracker

