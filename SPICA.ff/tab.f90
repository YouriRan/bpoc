program gentable
implicit none

real, parameter :: delr=0.002, rcut=2.5
real :: r
integer :: nbins, j

nbins=int( (rcut+1)/delr ) + 1
do j=0,nbins
    r=delr*j
    write(6,*) r, 1/r, 1/(r*r), -1/(r**6), -6/(r**7), 1/(r**9), 9/(r**10)
end do

end program