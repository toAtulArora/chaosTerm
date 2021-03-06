program lorenz                  
use gnuplot_fortran
implicit none
!dt=4.0E-2,
real, parameter :: sigma=10,rho=29,rho2=35, beta=8.0/3, gamma=20.0 , delta=5.0 ,dt=4.0E-2,epsilon=2.85,maxT=500, maxEps=4,deltaEps=0.2, edg=2.0
integer, parameter :: maxI = maxT/dt, snakeSize=10, maxSpheres=200, maxEpsI=maxEps/deltaEps
integer :: i,j,k,l
real :: relativeDistance,t,eps
logical :: plotGraphs = .false.

real, dimension(3) :: x0 = (/ 0.8, 0.3068, 7.0 /), y0 = (/ -2.8, -0.3068, -7.0 /), yf= (/ 20,0,0/)
real, dimension(maxI,3) :: xp,yp !yf = reshape( (/ (40,i=1,3*maxI) /), (/ maxI, 3 /) )
real, dimension(2*maxI,3) :: xy
real, dimension (3) :: x,xc,k1,k2,k3,k4
real, dimension (3) :: y,yc,m1,m2,m3,m4
real, dimension (maxI,2) :: dis
real, dimension (maxEpsI,2) :: epscount

type sphereLike
   integer :: count
   real, dimension(3):: origin
end type sphereLike

type(sphereLike), dimension(maxSpheres)::sphere
   
call startPlot()

call setXrange(0.0,100.0)
call setYrange(0.0,100.0)
call setZrange(0.0,100.0)

call srand(100)

write (*,*) "Iterating through..(will do ", maxI, "iterations)"

i=1
t=0
do while (t<maxT) !i<3000)
   t=t+dt

   k1=xDot(x0)
   k2=xDot(x0 + 0.5*dt*k1)
   k3=xDot(x0 + 0.5*dt*k2)
   k4=xDot(x0 + dt*k3)
   xc = x0 + dt*(1.0/6)*(k1+2*k2+2*k3+k4)
   x0=xc

   k1=yDot(xc(1),y0)
   k2=yDot(xc(1),y0 + 0.5*dt*k1)
   k3=yDot(xc(1),y0 + 0.5*dt*k2)
   k4=yDot(xc(1),y0 + dt*k3)
   yc = y0 + dt*(1.0/6)*(k1+2*k2+2*k3+k4)
   y0=yc



   xp(i,:)=xc
   yp(i,:)=yc + yf
   dis(i,2)=sqrt(sum((xc-yc)*(xc-yc)))
   dis(i,1)=t

   xy(1:i,:)=xp(1:i,:)
   xy(i:2*i,:)=yp(1:i,:)
   if(plotGraphs) then
      call nextPlot3d(xy(1:2*i,1),xy(1:2*i,2),xy(1:2*i,3))
      call nextPlot2d(dis(:i,1),dis(:i,2))
   end if
   i=i+1

end do


write (*,*) "Finding the box dimension.."

j=0
eps=0.0
call initOrigin(sphere(:))
!Initialization is correct
! write (*,*) "Spheres:", sphere(1:maxSpheres)%origin(1)

!isInsideSphere is working too
! write(*,*) "Should be 0: ", isInsideSphere( (/ 0.0,0.0,0.0 /) , 1.0, (/ 2.0,2.0,2.0 /))
! write(*,*) "Should be 1: ", isInsideSphere( (/ 1.0,1.0,1.0 /) , 2.0, (/ 2.0,2.0,2.0 /))

do while (eps<maxEps)
   j=j+1
   eps=eps + deltaEps
   do k=1,maxSpheres
      sphere(k)%count=0
      do i=1,maxI
         sphere(k)%count=sphere(k)%count + isInsideSphere(sphere(k)%origin,eps,xp(i,:))
      end do
   end do
   !write (*,*) "For eps:", eps, " counts:", sphere(1:maxSpheres)%count
   epscount(j,1)=eps
   epscount(j,2)=sum(sphere(1:maxSpheres)%count)/real(maxSpheres)
end do

call plot2dSave(log(epscount(5:j,1)),log(epscount(5:j,2)+1.0E-13),filename="dimensionLogLog.pdf",picFormat=1)
call plot2dSave(epscount(1:j,1),epscount(1:j,2),filename="dimension.pdf",picFormat=1)

open(file="dimension.dat",unit=4)
do i=1,j
   write (4,*) epscount(i,1),epscount(i,2)
end do

call endPlot()
! call system("xdg-open result3d.avi")


contains
  function xDot(x)
    real, dimension(3)::xDot,x
    xDot(1)=sigma*(x(2)-x(1))
    xDot(2)=rho*x(1) - x(2) - gamma*x(1)*x(3)
    !xDot(2)=(x(1)*(rho-x(3))) - x(2)
    xDot(3)=(delta*x(1)*x(2)) - (beta*x(3))
  end function xDot

  function yDot(u,y)
    real, dimension(3)::yDot,y
    real :: u
    yDot(1)=sigma*(y(2)-y(1))
    yDot(2)=rho*u - y(2) - gamma*u*y(3)
    yDot(3)=(delta*u*y(2)) - (beta*y(3))
  end function yDot

  function isInsideSphere(r0,eps,r)
    real, dimension(3) :: r0,r
    real :: eps
    integer:: isInsideSphere
    isInsideSphere=0
    if(eps*eps > sum( (r0-r)*(r0-r)) )  then
       isInsideSphere=1
    end if
  end function isInsideSphere

  subroutine initOrigin(spheres)
    type(sphereLike), dimension(:) :: spheres
    do l=1,maxSpheres
       spheres(l)%origin(1)=(0.5-rand())*2.0*edg
       spheres(l)%origin(2)=(0.5-rand())*2.0*edg
       spheres(l)%origin(3)=(0.5-rand())*2.0*edg
    end do
  end subroutine
  
end program lorenz
