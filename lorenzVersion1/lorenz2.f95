program lorenz                  
use gnuplot_fortran
implicit none

!real, parameter :: sigma=10,rho=29,rho2=35, beta=8.0/3, dt=4.0E-2,epsilon=2.85,maxT=30
real, parameter :: sigma=10,rho=29,rho2=35, beta=8.0/3, gamma=20.0 , delta=5.0 ,dt=4.0E-2,epsilon=2.85,maxT=20
integer, parameter :: maxI = maxT/dt, snakeSize=10
integer :: i,j
real :: relativeDistance,t

real, dimension(3) :: x0 = (/ 0.8, 0.3068, 7.0 /), y0 = (/ -2.8, -0.3068, -7.0 /), yf= (/ 20,0,0/)
real, dimension(maxI,3) :: xp,yp !yf = reshape( (/ (40,i=1,3*maxI) /), (/ maxI, 3 /) )
real, dimension(2*maxI,3) :: xy
real, dimension (3) :: x,xc,k1,k2,k3,k4
real, dimension (3) :: y,yc,m1,m2,m3,m4
real, dimension (maxI,2) :: dis

call startPlot()

call setXrange(0.0,100.0)
call setYrange(0.0,100.0)
call setZrange(0.0,100.0)

i=1
j=0
t=0
do while (t<maxT) !i<3000)

   t=t+dt
   k1=xDot(x0)
   k2=xDot(x0 + 0.5*dt*k1)
   k3=xDot(x0 + 0.5*dt*k2)
   k4=xDot(x0 + dt*k3)
   xc = x0 + dt*(1.0/6)*(k1+2*k2+2*k3+k4)
   !xc=x0 + xDot(x0)*dt
   x0=xc

   k1=yDot(xc(1),y0)
   k2=yDot(xc(1),y0 + 0.5*dt*k1)
   k3=yDot(xc(1),y0 + 0.5*dt*k2)
   k4=yDot(xc(1),y0 + dt*k3)
   yc = y0 + dt*(1.0/6)*(k1+2*k2+2*k3+k4)
   !yc=y0 + yDot(y0)*dt
   y0=yc



   ! j=j+1
   ! if(mod(j,1)==0) then
   !    j=0
      xp(i,:)=xc
      yp(i,:)=yc + yf
      dis(i,2)=sqrt(sum((xc-yc)*(xc-yc)))
      dis(i,1)=t

      xy(1:i,:)=xp(1:i,:)
      xy(i:2*i,:)=yp(1:i,:)
      
      !yp(i,:)=yc
      ! if(i>snakeSize) then
      !      call nextPlot3d(xp(i-snakeSize:i,1),xp(i-snakeSize:i,2),xp(i-snakeSize:i,3))
      ! end if
      !call nextPlot3d(xp(:2*i + 1,1),xp(:2*i + 1,2),xp(:2*i + 1,3))      
      call nextPlot3d(xy(1:2*i,1),xy(1:2*i,2),xy(1:2*i,3))
      call nextPlot2d(dis(:i,1),dis(:i,2))

      !call nextPlot3d(yp(:i,1),yp(:i,2),yp(:i,3))
      ! call nextPlot3d(yp(:i,1)+yf(:i,1),yp(:i,2)+yf(:i,1),yp(:i,3)+yf(:i,1))

      !call nextPlot3d((/ xp(:1), yp(:i,1)+yf(:i,1) /),(/ xp(:1), yp(:i,2)+yf(:i,1) /),(/ xp(:1), yp(:i,3)+yf(:i,1) /))

      i=i+1
      !write(*,*) "Gla"

   ! end if

end do


call endPlot()
call system("xdg-open result3d.avi")


! type vector
!    real, dimension(3) :: x
! end type vector

!type vector q
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


end program lorenz
