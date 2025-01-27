subroutine gauss_rng2(gauss_rng,itmax)

implicit none

integer :: i
double precision :: pi,itmax
double precision :: u1(itmax), u2(itmax)
double precision :: gauss_rng(itmax)

call random_seed()
call random_number(u1)
call random_number(u2)

! Generate a gaussian distributed random number generator based on the box muller method. 
pi = acos(-1.0)   !  = 3.14159...
do i = 1, itmax   
gauss_rng(i) = sqrt(-2*log(u1(i)))*cos(2*pi*u2(i))
end do

end subroutine gauss_rng2
     

