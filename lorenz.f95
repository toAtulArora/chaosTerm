program lorenz
!chaoticatmospheres.deviantart.com/art/Strange-Attractors-The-Coupled-Lorenz-Attractor-376065911

  real, parameter :: sigma=10, rho1=29, rho2=35, beta=8.0/3,dt=1.0E-3,epsilon=2.85
  real :: relativeDistance
  ! real :: x=0.008,y=0.103068,z=7,t
  real :: xO=0.008,yO=0.103068,zO=7.0,t
  real :: x,y,z

  ! real :: x1,y1,z1
  ! real :: x1o,y1o,z1o

  real :: x2o=0.008,y2o=100.1030681,z2o=7.0
  real :: x2,y2,z2
  !real :: x2o,y2o,z2o

  logical :: yoyoSucks=.true.
  !write(*,*) "Yoyo sucks :("
  do while (t<1000)
     t=t+dt
     x=xO + xDot(xO,yO,zO)*dt
     y=yO + yDot(xO,yO,zO,rho1)*dt
     z=zO + zDot(xO,yO,zO)*dt
     
     x2=x2o + xDotCoupled(x2o,y2o,z2o,xO)*dt
     y2=y2o + yDot(x2o,y2o,z2o,rho2)*dt
     z2=z2o + zDot(x2o,y2o,z2o)*dt
     
     xO=x
     yO=y
     zO=z
     
     x2o=x2
     y2o=y2
     z2o=z2

     if(yoyoSucks) then
        yoyoSucks=.false.
        relativeDistance= sqrt((x-x2)**2 + (y-y2)**2 + (z-z2)**2)
        write (*,*) x,y,z,x2,y2,z2,relativeDistance,t        
     else
        yoyoSucks=.true.
     end if
  end do

contains
  function xDotCoupled(x,y,z,x2)
    real:: x,y,z,xDotCoupled,x2
    xDotCoupled=sigma*(y-x) + (epsilon*(x2-x))
  end function xDotCoupled

  function xDot(x,y,z)
    real:: x,y,z,xDot
    xDot=sigma*(y-x)
  end function xDot
  
  function yDot(x,y,z,rho)
    real:: x,y,z,yDot,rho
    yDot=(x*(rho-z)) - y
  end function yDot
  
  function zDot(x,y,z)
    real:: x,y,z,zDot
    zDot=(x*y)-(beta*z)
  end function zDot
end program lorenz
