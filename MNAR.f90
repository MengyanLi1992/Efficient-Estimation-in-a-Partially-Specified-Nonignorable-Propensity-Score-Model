
module commondata 
  implicit none
  integer, parameter ::  sim=100, n=1000, L=12, lt=2, LL=L+lt, len1 = 20, len2=13,len3=20, len4=20
  integer, parameter :: iseed=1, lw=lt*(lt*3+13)/2, lw1=L*(L*3+13)/2
  double precision, parameter :: beta0= -1.d0, alpha0=2.d0
  double precision, dimension(L) :: gamma0, gamma, gamma_int
  double precision :: alpha, beta
  double precision :: tol, tol1
  double precision,dimension(lw) :: wv
  double precision,dimension(lw1) :: wv1
  double precision, dimension(L) :: h
  double precision, dimension(len1) :: yy1, ww1
  double precision, dimension(len2) :: yy2, ww2
  double precision, dimension(len3) :: uu3, ww3
  double precision, dimension(len4) :: zz4, ww4
  double precision, dimension(n) :: u, z, y, u0, z0, y0, u1, z1, y1
  integer, dimension(n) :: r
  double precision, dimension(L) :: d
  integer :: n1

end module commondata

!start of program

program mnar
  use commondata 
  double precision, dimension(lt) :: theta, theta0
  double precision :: m0
  double precision, dimension(2) :: endpts
  double precision, dimension(len1) :: bb1
  double precision, dimension(len2) :: bb2
  double precision, dimension(len3) :: bb3
  double precision, dimension(len4) :: bb4
  integer :: eflag1, eflag
  integer :: i,iter
  double precision, dimension(lt) :: esteq_out
  double precision :: temp1, temp
  double precision, dimension(LL,LL) :: V
  double precision, dimension(lt,lt) :: V1, RRinv

  !double precision :: out_ij, out0_i, out1_i
  !double precision, dimension(n) :: expcta0_out,expcta1_out, expcta0_out0,expcta1_out0

  external esteqtheta 
  external esteqgamma

  open(1, file='out11.dat')
  open(2, file='Vout11.dat')
  open(3, file='VVout11.dat')
  !open(4, file='TGout4.dat')
  !open(5, file='TGout5.dat')
  !open(6, file='out6.dat')
  !open(7, file='out7.dat')
  !open(8, file='out8.dat')


  call gaussq(4,len1,0.d0,0.d0,0,endpts,bb1,yy1,ww1)
  call gaussq(4,len2,0.d0,0.d0,0,endpts,bb2,yy2,ww2)
  ! numeric integral on [0,1]
  !JACOBI QUADRATURE, W(X) = (1-X)**ALPHA * (1+X)**BETA ON (-1, 1), ALPHA, BETA .GT. -1.
  ! need to change variable uu = 2u-1
  call gaussq(5,len3,0.d0,0.d0,0,endpts,bb3,uu3,ww3)
  call gaussq(4,len4,0.d0,0.d0,0,endpts,bb4,zz4,ww4)

  d(1) = -2.5d0
  d(2) = -1.8d0
  d(3) = -1.2d0
  d(4) = -0.7d0
  d(5) = -0.3d0
  d(6) = 0.d0
  d(7) = 0.3d0
  d(8) = 0.7d0
  d(9) = 1.2d0
  d(10) = 1.8d0
  d(11) = 2.5d0
  d(12) = 3.5d0

  h(1)=1.5*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(2)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(3)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(4)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(5)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(6)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(7)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(8)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(9)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(10)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(11)=2*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  h(12)=1.5*sqrt(7.0/3.0)*((n*0.79)**(-0.25))
  !h(13)=1.5*sqrt(7.0/3.0)*((n*0.79)**(-0.25))

  call gFCNVEC(L, d, gamma0)

  ! tolerance for hybrd10
  tol = 0.00001
  tol1 = 0.01
  eflag = 1
  eflag1 = 1
  iter=1

  theta0(1) = alpha0
  theta0(2) = beta0
  
  temp=rand(iseed)

10 if (iter.le.sim) then 
     print*, 'iter=', iter
     m0 = 3.6d0
     call genordata(m0,u,z,y,r,u1,z1,y1,u0,z0,y0)
     n1 = sum(r) 
     !call genordata1(nn,uu,zz)
     !print*, 'n1=', n1
     !n0 = n - n1
     !add random perturbations to the true values of thata
     do i=1,lt
       temp1=(rand(0)/5)-0.1
       theta(i)=theta0(i)*(1.0+temp1) 
     end do

     do i=1,L
       temp1=(rand(0)/5)-0.1
       gamma_int(i)=gamma0(i)*(1.0+temp1)
     end do
    
     call hybrd10(esteqtheta,lt,theta,esteq_out,tol,eflag,wv,lw)

     if ((eflag.ne.1).or.(sum(abs(esteq_out)).gt.1000.)) then
       print*,'fail estimation, eflag=',eflag,'iter=',iter
       goto 10
     end if

     write(1,*) theta, gamma

     call matrixV(theta, gamma, V)
     call matrixV1(theta, gamma, V1, RRinv)

     write(2,*)  V(1,1)/n, V(2,2)/n, V(1,3)/n, V(2,3)/n, V(1,9)/n, V(2,9)/n, V(1,14)/n, V(2,14)/n
     write(3,*)  V1(1,1)/n, V1(2,2)/n, RRinv(1,1)/n, RRinv(2,2)/n

     iter=iter+1
     goto 10
    end if
  !write(1,*), y1
  close(1)
  close(2)
  close(3)
  !close(4)
  !close(5)
  !close(6)
  !close(7)
  !close(8)

return
end program mnar




! data generation

subroutine fyz1(alpha,y,z,prob)
  double precision, intent(in) :: alpha
  double precision, intent(in) :: y,z
  double precision, intent(out) :: prob

  prob = 1/sqrt(2*3.1415926)*exp(-(y-alpha*z)**2/2)

return
end subroutine fyz1

subroutine gFCN(y, gy)
  double precision, intent(in) :: y
  double precision, intent(out) :: gy

  gy = 3.d0*exp(y)/(1+exp(y))

return
end subroutine gFCN


subroutine piFCN(beta,y,u,prob)
  double precision, intent(in) :: beta
  double precision, intent(in) :: y,u
  double precision, intent(out) :: prob
  double precision :: gy

  call gFCN(y,gy)
  prob=exp(beta*u + gy)/(1 + exp(beta*u + gy))

return
end subroutine piFCN


subroutine fygivenx(alpha,beta,y,z,u,prob)
  use commondata, only : len1, yy1, ww1
  double precision, intent(in) :: alpha, beta
  double precision, intent(in) :: y,z,u
  double precision, intent(out) :: prob
  double precision :: pi, prob1, temp,denom, yy, gy
  integer :: l1

  denom=0.d0
  do l1=1,len1
    yy = sqrt(2.0)*yy1(l1)+alpha*z
    call gFCN(yy, gy)
    temp=1.0/sqrt(3.1415926)*(1.0+exp(-beta*u)*exp(-gy))*ww1(l1)
    denom=denom+temp
  end do
  call fyz1(alpha,y,z,prob1)
  call piFCN(beta,y,u,pi)
  prob=(prob1/pi)/denom

  return
end subroutine fygivenx



subroutine genordata(m0,u,z,y,r,u1,z1,y1,u0,z0,y0)
  use commondata,only : n,beta0,alpha0
  double precision, intent(in) :: m0
  double precision,dimension(n),intent(out) :: u,z,y,u1,z1,y1,u0,z0,y0
  integer,dimension(n),intent(out) :: r
  real :: zz0,zz1
  double precision ::yy0,temp,prob1,prob2,pi
  integer :: i,itr,j,k

  j=1
  k=1
  do i=1,n
    !print*, 'i=', i
    ! generate uniform(0,1) u
    u(i)=rand(0)
    ! generate normal(u,1) z
    call rnorm(zz0,zz1)
    z(i)=u(i)+zz0*0.5
    ! generate y using rejection sampling
    yy0=zz1+alpha0*z(i)
    temp=rand(0)
    call fygivenx(alpha0,beta0,yy0,z(i),u(i),prob1)
    call fyz1(alpha0,yy0,z(i),prob2)
    itr=1
    do while (temp.gt.(prob1/(m0*prob2)))
      itr=itr+1
      !print*, 'itr=', itr 
      call rnorm(zz0,zz1)
      yy0=zz1+alpha0*z(i)
      temp=rand(0)
      call fygivenx(alpha0,beta0,yy0,z(i),u(i),prob1)
      call fyz1(alpha0,yy0,z(i),prob2)
    end do
    y(i)=yy0
    ! generate r
    temp = rand(0)
    call piFCN(beta0,y(i),u(i),pi)
    if (temp.lt.pi) then 
      r(i)=1
      u1(j)=u(i)
      z1(j)=z(i)
      y1(j)=y(i)
      j=j+1
    else
      r(i)=0
      u0(k)=u(i)
      z0(k)=z(i)
      y0(k)=y(i)
      k=k+1
    end if
  end do 
  return
end subroutine genordata


! estimate function g using local linear and cubic interpolation

! kernel epanechnikov

subroutine kh(nn,x,h,y)
  integer,intent(in) :: nn
  double precision,dimension(nn),intent(in) :: x
  double precision,intent(in) :: h
  double precision,dimension(nn),intent(out) :: y
  double precision,dimension(nn) :: tmp
  integer::i

  tmp=x/h
  y=0.0
  do i=1,nn
    if (abs(tmp(i)).lt.1) then
      y(i)=(1.-tmp(i)**2.)/h*3./4.
    end if
  end do
  return
end subroutine kh


subroutine kh_Gaussian(nn,x,h,y)
  integer,intent(in) :: nn
  double precision,dimension(nn),intent(in) :: x
  double precision,intent(in) :: h
  double precision,dimension(nn),intent(out) :: y
  double precision,dimension(nn) :: tmp
  integer :: i

  tmp=x/h
  do i=1,nn
    y(i)=1/sqrt(2*3.1415926)*exp(-(tmp(i)**2)/2)/h
  end do
  return
end subroutine kh_Gaussian


subroutine fyz1VEC(alpha,n,y,z,prob)
  double precision, intent(in) :: alpha
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: y,z
  double precision, dimension(n), intent(out) :: prob

  prob = 1/sqrt(2*3.1415926)*exp(-(y-alpha*z)**2/2)

return
end subroutine fyz1VEC


subroutine gFCNVEC(p,y, gy)
  integer, intent(in) :: p
  double precision, dimension(p), intent(in) :: y
  double precision, dimension(p), intent(out) :: gy

  gy = 3*exp(y)/(1+exp(y))

return
end subroutine gFCNVEC



subroutine cubicinterp(xx,yy,xx0,yy0)
  double precision, dimension(4), intent(in) :: xx, yy
  double precision, intent(in) :: xx0
  double precision, intent(out) :: yy0
  double precision, dimension(4,4) :: A 
  double precision :: cond
  integer,dimension(4) :: ipvt
  double precision,dimension(4) :: wv, yytemp

  A(:,1) = 1.d0
  A(:,2) = xx
  A(:,3) = xx**2.0
  A(:,4) = xx**3.0

  yytemp = yy ! this step is necessary. Without it, when we call this step second time, 
  !yy is the root of AX=b of first call. And if we keep using the same yy for the two calls,
  !even if we redefine yy, there can be problems. 

  call decomp(4,4,A,cond,ipvt,wv)
  call solvels(4,4,A,yytemp,ipvt)
  yy0 = yytemp(1) + yytemp(2)*xx0 + yytemp(3)*xx0**2.0 + yytemp(4)*xx0**3.0

  return
end subroutine cubicinterp


subroutine quadrainterp(xx,yy,xx0,yy0)
  double precision, dimension(3), intent(in) :: xx, yy
  double precision, intent(in) :: xx0
  double precision, intent(out) :: yy0
  double precision, dimension(3,3) :: A 
  double precision :: cond
  integer,dimension(3) :: ipvt
  double precision,dimension(3) :: wv, yytemp

  A(:,1) = 1.d0
  A(:,2) = xx
  A(:,3) = xx**2.0

  yytemp = yy ! this step is necessary. Without it, when we call this step second time, 
  !yy is the root of AX=b of first call. And if we keep using the same yy for the two calls,
  !even if we redefine yy, there can be problems. 

  call decomp(3,3,A,cond,ipvt,wv)
  call solvels(3,3,A,yytemp,ipvt)
  yy0 = yytemp(1) + yytemp(2)*xx0 + yytemp(3)*xx0**2.0 

return 
end subroutine quadrainterp



subroutine cubicinterpVEC(l1,xx,yy,l2,cubicinterp_in,cubicinterp_out)
  !l1: L, xx:d, yy:gamma, l2:len2, cubicinterp_in:yy2
  integer, intent(in) :: l1, l2
  double precision, dimension(l1), intent(in) :: xx, yy
  double precision, dimension(l2), intent(in) :: cubicinterp_in
  double precision, dimension(l2), intent(out) :: cubicinterp_out
  double precision, dimension(4) :: tempx,tempy
  double precision, dimension(3) :: tempx1,tempy1
  integer :: i, j

  call dsort(xx, yy, l1, 2)
  j=1
  do i=1,l2
  !print*, 'i=', i
    ! Find the interval that cubicinterp_in(i) resides 
    do while ((cubicinterp_in(i).gt.xx(j)).and.(j.lt.(l1+1)))
      j=j+1
    end do
    !print*, 'j=', j
    ! cubicinterp_in(l2) is in interval (xx(j-1), xx(j)]
    if ((j.gt.2).and.(j.lt.l1)) then 
      tempx = xx((j-2):(j+1))
      tempy = yy((j-2):(j+1))
      call cubicinterp(tempx, tempy, cubicinterp_in(i), cubicinterp_out(i))
    else if (j.eq.2) then 
      !tempx = xx(1:4)
      !tempy = yy(1:4)
      !call cubicinterp(tempx, tempy, cubicinterp_in(i), cubicinterp_out(i))
      tempx1 = xx(1:3)
      tempy1 = yy(1:3)
      call quadrainterp(tempx1, tempy1,cubicinterp_in(i), cubicinterp_out(i))
    else if (j.eq.l1) then
      !k=l1-1
      !cubicinterp_out(i) = (yy(l1)-yy(k))/(xx(l1)-xx(k))*cubicinterp_in(i) + (yy(k)*xx(l1)-yy(l1)*xx(k))/(xx(l1) - xx(k))
      !tempx = xx((l1-3):l1)
      !tempy = yy((l1-3):l1)
      !call cubicinterp(tempx, tempy, cubicinterp_in(i), cubicinterp_out(i))
      tempx1 = xx((l1-2):l1)
      tempy1 = yy((l1-2):l1)
      call quadrainterp(tempx1, tempy1, cubicinterp_in(i), cubicinterp_out(i))
    ! boudary (may use different methods)
    else if (j.eq.1) then
      cubicinterp_out(i) = yy(1)
      !cubicinterp_out(i) = ((yy(2) - yy(1))*cubicinterp_in(i) + yy(1)*xx(2) - yy(2)*xx(1))/(xx(2) - xx(1))
    else if (j.eq.(l1+1)) then
      cubicinterp_out(i) = yy(l1)
    end if 
    !print*, 'cubicinterp_out(i)=', cubicinterp_out(i)         
  end do
return 
end subroutine cubicinterpVEC



subroutine appexpct(alpha,gamma,n0,z,appexpct_out)
  use commondata, only: L, d, len2, yy2, ww2
  double precision, intent(in) :: alpha
  double precision, dimension(L), intent(in) :: gamma
  integer, intent(in) :: n0
  double precision, dimension(n0), intent(in) :: z
  double precision, dimension(n0), intent(out) :: appexpct_out
  integer ::  k
  double precision, dimension(len2) :: temp
  double precision, dimension(len2) :: yy, gy0, gy

  do k=1,n0
    yy = sqrt(2.0)*yy2+alpha*z(k)
    !call gFCNVEC(len2, yy, gy0)
    call cubicinterpVEC(L,d,gamma,len2,yy,gy)
    !print*, 'dif=', gy0-gy
    temp = exp(-gy)/sqrt(3.1415926)*ww2
    appexpct_out(k) = sum(temp)
  end do 

  !temp1 = 0.d0
  !do l2=1,len2
    !temp = exp(-gamma(l2))/sqrt(3.1415926)*exp(-(1/2)*(sqrt(2.d0)*yy2(l2)*(alpha - 2*alpha*z) + (alpha/2 - alpha*z)**2))*ww2(l2)
    !temp1 = temp1 + temp
  !end do 
  !appexpct_out = temp1
  
return
end subroutine appexpct


subroutine appexpctkh(alpha,gamma,h,yy0,n0,z,appexpctkh_out)
  use commondata, only: L, d, len2, yy2, ww2
  double precision, intent(in) :: alpha
  double precision, dimension(L), intent(in) :: gamma
  double precision, intent(in) :: h, yy0
  integer, intent(in) :: n0
  double precision, dimension(n0), intent(in) :: z
  double precision, dimension(n0), intent(out) :: appexpctkh_out
  integer ::  k
  double precision, dimension(len2) :: temp
  double precision, dimension(len2) :: yy, gy0, gy, kh_out


  do k=1,n0
    yy = sqrt(2.0)*yy2+alpha*z(k)
    !call gFCNVEC(len2, yy, gy0)
    call cubicinterpVEC(L,d,gamma,len2,yy,gy)
    !print*, 'dif=', gy0-gy
    !call kh_Gaussian(len2,yy-yy0,h,kh_out)
    call kh(len2,yy-yy0,h,kh_out)
    temp = kh_out*exp(-gy)/sqrt(3.1415926)*ww2
    appexpctkh_out(k) = sum(temp)
  end do 

  !temp1 = 0.d0
  !do l2=1,len2
    !temp = exp(-gamma(l2))/sqrt(3.1415926)*exp(-(1/2)*(sqrt(2.d0)*yy2(l2)*(alpha - 2*alpha*z) + (alpha/2 - alpha*z)**2))*ww2(l2)
    !temp1 = temp1 + temp
  !end do 
  !appexpct_out = temp1
  
return
end subroutine appexpctkh


subroutine appexpctkh_G(alpha,gamma,h,yy0,n0,z,appexpctkh_out)
  use commondata, only: L, d, len2, yy2, ww2
  double precision, intent(in) :: alpha
  double precision, dimension(L), intent(in) :: gamma
  double precision, intent(in) :: h, yy0
  integer, intent(in) :: n0
  double precision, dimension(n0), intent(in) :: z
  double precision, dimension(n0), intent(out) :: appexpctkh_out
  integer ::  k
  double precision, dimension(len2) :: temp
  double precision, dimension(len2) :: yy, gy0, gy, kh_out


  do k=1,n0
    yy = sqrt(2.0)*yy2+alpha*z(k)
    !call gFCNVEC(len2, yy, gy0)
    call cubicinterpVEC(L,d,gamma,len2,yy,gy)
    !print*, 'dif=', gy0-gy
    call kh_Gaussian(len2,yy-yy0,h,kh_out)
    !call kh(len2,yy-yy0,h,kh_out)
    temp = kh_out*exp(-gy)/sqrt(3.1415926)*ww2
    appexpctkh_out(k) = sum(temp)
  end do

!temp1 = 0.d0
!do l2=1,len2
!temp = exp(-gamma(l2))/sqrt(3.1415926)*exp(-(1/2)*(sqrt(2.d0)*yy2(l2)*(alpha - 2*alpha*z) + (alpha/2 - alpha*z)**2))*ww2(l2)
!temp1 = temp1 + temp
!end do
!appexpct_out = temp1

  return
end subroutine appexpctkh_G


subroutine gammaFNC(alpha, beta, L, gamma, d, h, n, n1, y1, u1, z0, u0,gammaFNC_out)
  use commondata, only: len2, yy2, ww2
  double precision, intent(in) :: alpha,  beta
  integer, intent(in) :: L
  double precision, dimension(L), intent(in) :: gamma, d, h 
  integer, intent(in) :: n, n1
  double precision, dimension(n), intent(in) :: y1, u1
  double precision, dimension(n), intent(in) :: z0, u0
  double precision, dimension(L), intent(out) :: gammaFNC_out
  integer :: j, n0
  double precision, dimension(n) :: temp 
  double precision, dimension(n1) :: temp1, kh_out
  double precision :: left, right
  double precision, dimension(n-n1) ::  temp2, temp3, prob0, trueexpct_out, appexpct_out, appexpctkh_out
  double precision, dimension(n-n1) :: pi0
  double precision, dimension(n1) :: pi1
  
  n0 = n-n1

  !call trueexpct(alpha,n0,z0(1:n0),trueexpct_out)
  call appexpct(alpha,gamma,n0,z0(1:n0),appexpct_out)
  !print*, 'diff=', sum((appexpct_out/trueexpct_out) - 1.0)/n0
  
  do j=1,L
    temp = y1-d(j)
    temp1 = temp(1:n1)
    call kh(n1,temp1,h(j),kh_out)
    if ((j.lt.1.5).or.(j.gt.12.5)) then
      call kh_Gaussian(n1,temp1,h(j),kh_out)
    end if
    pi1 = 1.0/(1.0 + exp(-beta*u1(1:n1)-gamma(j)))
    !temp1 = kh_out/(1+exp(gamma(j) + beta*u1(1:n1)))
    temp1 = kh_out*(1.0-pi1)
    left = sum(temp1)/n

    !temp2 = d(j)
    !call fyz1VEC(alpha,n0,temp2,z0(1:n0),prob0)
    !pi0 = 1/(1 + exp(-beta*u0(1:n0)-gamma(j)))
    !temp3 = exp(-gamma(j))*prob0/appexpct_out*pi0 !or change trueexpct_out to appexpct_out
    !right = sum(temp3)/n

    call appexpctkh(alpha,gamma,h(j),d(j),n0,z0(1:n0),appexpctkh_out)
    if ((j.lt.1.5).or.(j.gt.12.5)) then
      call appexpctkh_G(alpha,gamma,h(j),d(j),n0,z0(1:n0),appexpctkh_out)
    end if
    pi0 = 1.0/(1.0 + exp(-beta*u0(1:n0)-gamma(j)))
    temp3 = appexpctkh_out/appexpct_out*pi0 !or change trueexpct_out to appexpct_out
    right = sum(temp3)/n

    gammaFNC_out(j) = left-right
  end do

return 
end subroutine gammaFNC



subroutine gammaFNC1(alpha, beta, L, gamma, d, h, n, n1, y1, u1, z0, u0,gammaFNC_out)
  use commondata, only: len2, yy2, ww2
  double precision, intent(in) :: alpha,  beta
  integer, intent(in) :: L
  double precision, dimension(L), intent(in) :: gamma, d, h 
  integer, intent(in) :: n, n1
  double precision, dimension(n), intent(in) :: y1, u1
  double precision, dimension(n), intent(in) :: z0, u0
  double precision, dimension(L,n), intent(out) :: gammaFNC_out
  integer :: j, n0
  double precision, dimension(n) :: temp 
  double precision, dimension(n1) :: temp1, kh_out
  double precision, dimension(n) :: left, right
  double precision, dimension(n-n1) ::  temp2, temp3, prob0, trueexpct_out, appexpct_out, appexpctkh_out
  double precision, dimension(n-n1) :: pi0
  double precision, dimension(n1) :: pi1
  
  n0 = n-n1

  !call trueexpct(alpha,n0,z0(1:n0),trueexpct_out)
  call appexpct(alpha,gamma,n0,z0(1:n0),appexpct_out)
  !print*, 'diff=', sum((appexpct_out/trueexpct_out) - 1.0)/n0
  
  do j=1,L
    temp = y1-d(j)
    temp1 = temp(1:n1)
    call kh(n1,temp1,h(j),kh_out)
    if ((j.lt.1.5).or.(j.gt.12.5)) then
      call kh_Gaussian(n1,temp1,h(j),kh_out)
    end if
    pi1 = 1.0/(1.0 + exp(-beta*u1(1:n1)-gamma(j)))
    !temp1 = kh_out/(1+exp(gamma(j) + beta*u1(1:n1)))
    temp1 = kh_out*(1.0-pi1)
    left = temp1

    !temp2 = d(j)
    !call fyz1VEC(alpha,n0,temp2,z0(1:n0),prob0)
    !pi0 = 1/(1 + exp(-beta*u0(1:n0)-gamma(j)))
    !temp3 = exp(-gamma(j))*prob0/appexpct_out*pi0 !or change trueexpct_out to appexpct_out
    !right = sum(temp3)/n

    call appexpctkh(alpha,gamma,h(j),d(j),n0,z0(1:n0),appexpctkh_out)
    if ((j.lt.1.5).or.(j.gt.12.5)) then
      call appexpctkh_G(alpha,gamma,h(j),d(j),n0,z0(1:n0),appexpctkh_out)
    end if
    pi0 = 1/(1 + exp(-beta*u0(1:n0)-gamma(j)))
    temp3 = appexpctkh_out/appexpct_out*pi0 !or change trueexpct_out to appexpct_out
    right = temp3

    gammaFNC_out(j,:) = left-right
  end do

return 
end subroutine gammaFNC1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integral equation for a1 and a2



subroutine appexpctsg(alpha,gamma,z,appexpctsg_out)
  use commondata, only: L, d, len2, yy2, ww2
  double precision, intent(in) :: alpha
  double precision, dimension(L), intent(in) :: gamma
  double precision, intent(in) :: z
  double precision, intent(out) :: appexpctsg_out
  double precision, dimension(len2) :: yy
  double precision, dimension(len2) :: temp1
  double precision, dimension(len2) :: gy

  yy = yy2*sqrt(2.0) + z*alpha
  call cubicinterpVEC(L,d,gamma,len2,yy,gy)
  !l1: L, xx:d, yy:gamma, l2:len2, cubicinterp_in:yy2

  temp1 = exp(-gy)*sqrt(2.0)*yy2*z/sqrt(3.1415926)*ww2
 
  appexpctsg_out = sum(temp1)

return
end subroutine appexpctsg



subroutine integralopt_hat(alpha, beta, gamma, i, j, n, z, u, out_ij, out0_i, out1_i)
  ! here i and j should be <= len2
  ! out_ij = k(yi, tj); out0_i = phi_0(yi); out1_i = phi_1(yi)
  use commondata, only: L, yy2, d, len2, ww2
  double precision, intent(in) :: alpha,  beta
  double precision, dimension(L), intent(in) :: gamma
  integer, intent(in) :: i, j, n
  double precision, dimension(n), intent(in) :: z, u
  double precision, intent(out) :: out_ij
  double precision, intent(out) :: out0_i, out1_i
  double precision, dimension(n) :: w_inv
  double precision, dimension(n) :: yi, tj
  double precision, dimension(n) :: trueexpct_out, appexpct_out, trueexpctsg_out, appexpctsg_out
  double precision, dimension(n) :: probi, probj, temp1, temp2, temp3
  integer :: k

  !num, dem, temp1, temp2, temp3, temp4, temp5, probi, probj
  !double precision :: tempu_ij, tempz_ij, tempz0_i, tempu1_i,tempz1_i
  !double precision :: trueexpctsg_out, appexpctsg_out
  !double precision, dimension(len4) :: trueexpct_out, appexpct_out 
  !integer :: k,l4

  yi = sqrt(2.0)*yy2(i) + alpha/2
  tj = sqrt(2.0)*yy2(j) + alpha/2

  !call trueexpct(alpha,n,z,trueexpct_out)
  call appexpct(alpha,gamma,n,z,appexpct_out)
  !print*, 'diffexpctg=', trueexpct_out-appexpct_out
  w_inv  = 1.0 + exp(-beta*u)*appexpct_out

  call fyz1VEC(alpha,n,yi,z,probi)
  call fyz1VEC(alpha,n,tj,z,probj)
  temp1 = probi*probj/appexpct_out*exp(-beta*u)/(w_inv**2.0)
  out_ij = sum(temp1)/n
  
  do k=1,n
    !call trueexpctsg(alpha,z(i),trueexpctsg_out(i))   
    call appexpctsg(alpha,gamma,z(k),appexpctsg_out(k))
  end do

  temp2 = probi*appexpctsg_out/appexpct_out*exp(-beta*u)/(w_inv**2.0)
  out0_i = -sum(temp2)/n

  temp3 = u*probi*exp(-beta*u)/(w_inv**2.0)
  out1_i = sum(temp3)/n

return
end subroutine integralopt_hat



subroutine MandPhi(alpha, beta, gamma, n, z, u, M, Phi0, Phi1)
  !!! treat a1(y)exp(-g(y)) as the unknown function !!!
  use commondata, only: L, yy2, len2, d, ww2
  double precision, intent(in) :: alpha,  beta
  double precision, dimension(L), intent(in) :: gamma
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: z, u
  double precision, dimension(len2, len2), intent(out) :: M
  double precision, dimension(len2), intent(out) :: Phi0, Phi1
  double precision :: out_ij
  integer :: i, j

  do i=1,len2
    do j=1,len2
      !call integralopt(alpha, beta, gamma, i, j, out_ij, Phi0(i), Phi1(i))
      call integralopt_hat(alpha, beta, gamma, i, j, n, z, u, out_ij, Phi0(i), Phi1(i))
      ! Use Hermite-Gauss quadrature
      !M(i,j) = out_ij*exp(yy2(j)**2.0)*ww2(j)
      ! without Hermite-Gauss quadrature
      M(i,j) = out_ij
    end do
  end do
  return
end subroutine MandPhi



subroutine aFCN(alpha, beta, gamma, n, z, u, a0VEC, a1VEC, cond)
  use commondata, only: L, yy2, len2, d, ww2
  double precision, intent(in) :: alpha,  beta
  double precision, dimension(L), intent(in) :: gamma
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: z, u
  double precision, dimension(len2), intent(out) :: a0VEC, a1VEC
  double precision,intent(out) :: cond
  double precision, dimension(len2, len2) :: M, MM
  double precision, dimension(len2) :: Phi0, Phi1
  integer,dimension(len2) :: ipvt
  double precision,dimension(len2) :: wv
  double precision :: delta

  !delta= 0.d0
  delta = 1e-6
  call MandPhi(alpha, beta, gamma, n, z, u, M, Phi0, Phi1)
  do i=1,len2
    M(i,i)=M(i,i)+delta
  end do
 
  MM = M 

  ! cond is condition number, the algorithm breaks down if it's too large
  call decomp(len2,len2,M,cond,ipvt,wv)
  !print*, 'cond=', cond
  !print*, 'Phi0=', Phi0
  !print*, 'M=', M
  call solvels(len2,len2,M,Phi0,ipvt)
  a0VEC = Phi0

  call decomp(len2,len2,MM,cond,ipvt,wv)
  call solvels(len2,len2,MM,Phi1,ipvt)
  a1VEC = Phi1

  return
end subroutine aFCN


! construct estimating equationa for alpha and beta

subroutine expcta1a2(alpha,beta,gamma,n,z,u,expcta0_out,expcta1_out)
  use commondata, only: L, yy2, len2, d, ww2
  double precision, intent(in) :: alpha, beta
  double precision, dimension(L), intent(in) :: gamma
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: z, u
  double precision, dimension(n), intent(out) :: expcta0_out, expcta1_out
  integer :: l2
  double precision, dimension(n) :: temp0, temp1, temp
  double precision, dimension(len2) :: a0VEC, a1VEC
  double precision :: cond
  double precision, dimension(n) :: temp2

  call aFCN(alpha, beta, gamma, n, z, u, a0VEC, a1VEC, cond)

  temp0 = 0.d0
  temp1 = 0.d0
  do l2=1, len2
    temp2 = exp(-((alpha/2.0 - z*alpha)**2.0)/2.0)*exp(-sqrt(2.0)*(alpha/2.0 - z*alpha)*yy2(l2))
    temp = a0VEC(l2)/sqrt(3.1415926)*temp2*ww2(l2)
    temp0 = temp0 + temp
    temp = a1VEC(l2)/sqrt(3.1415926)*temp2*ww2(l2)
    temp1 = temp1 + temp
  end do
  expcta0_out = temp0
  expcta1_out = temp1

  return
end subroutine expcta1a2



subroutine effsc(alpha, beta, gamma, n, y, z, r, u, effscalpha_out, effscbeta_out)
  use commondata, only: L, d, len2, ww2, yy2
  double precision, intent(in) :: alpha,  beta
  double precision, dimension(L), intent(in) :: gamma
  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: z, u, y
  integer, dimension(n), intent(in) ::r
  double precision, dimension(n), intent(out) :: effscbeta_out, effscalpha_out
  double precision, dimension(n) :: expcta0_out, expcta1_out, trueexpctsg_out, appexpctsg_out, trueexpct_out, appexpct_out
  integer :: i
  double precision, dimension(n) :: w

  do i=1,n
    !call trueexpctsg(alpha,z(i),trueexpctsg_out(i))
    call appexpctsg(alpha,gamma,z(i),appexpctsg_out(i))
  end do 

  !call trueexpct(alpha,n,z,trueexpct_out)
  call appexpct(alpha,gamma,n,z,appexpct_out)
  !print*, 'appexpct_out=', sum(appexpct_out)

  call expcta1a2(alpha,beta,gamma,n,z,u,expcta0_out,expcta1_out)
  !print*, 'expcta0_out=', sum(expcta0_out)
  !print*, 'expcta1_out=', sum(expcta1_out)
  w = 1.0/(1.0+exp(-beta*u)*appexpct_out)
  effscalpha_out = r*(y-alpha*z)*z - (r - w)/appexpct_out*(appexpctsg_out + expcta0_out)
  effscbeta_out = (r-w)*u - (r-w)/appexpct_out*expcta1_out
  return
end subroutine effsc



subroutine esteqgamma(L, gamma, esteqgamma_out, eflag)
  use commondata, only : h, d, n, n1, y1, u1, z0, u0, alpha, beta,len2, yy2, ww2
  integer,intent(in) :: L
  double precision,dimension(L),intent(in) :: gamma
  double precision,dimension(L),intent(out) :: esteqgamma_out
  integer, intent(out) :: eflag
  double precision :: temp, temp0
  integer :: j

  !print*, 'alpha=', alpha
  !print*, 'beta=', beta

  esteqgamma_out = 10000000.d0

  temp0 = 1.d0
  do j=1,L
    temp=0.d0
    if ((-6.d0.lt.gamma(j)).and.(gamma(j).lt.6.d0)) then
      temp=1.d0
    end if
    temp0 = temp * temp0
  end do

  !do j=1,1
    !temp=0.d0
    !if (abs(gamma(j)).lt.2.d0) then
      !temp=1.d0
    !end if
    !temp0 = temp * temp0
  !end do

  if (temp0.gt.0.5d0) then
    call gammaFNC(alpha, beta, L, gamma, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out)
  end if

  !print*, 'gamma=', gamma
  !print*, 'esteqgamma_out=', esteqgamma_out

return 
end subroutine esteqgamma



subroutine esteqtheta(lt, theta, esteq_out, eflag)
  use commondata, only : alpha,beta,gamma_int,gamma,L,d,tol1,wv1,lw1,len2,ww2,yy2,n,y,z,r,u,h,n1,y1,u1,z0,u0
  integer,intent(in) :: lt
  double precision,dimension(lt),intent(in) :: theta
  double precision,dimension(lt),intent(out) :: esteq_out
  integer, intent(out) :: eflag
  double precision, dimension(L) :: esteqgamma_out
  integer :: i, eflag1
  double precision, dimension(n) :: effscbeta_out, effscalpha_out
  double precision :: temp, temp0
  integer :: j
  external esteqgamma

  alpha = theta(1)
  beta = theta(2)
  
  !print*, 'theta=', theta

  esteq_out = 10000000.d0

  temp0 = 1.d0
  do j=1,2
    temp=0.d0
    if ((-10.d0.lt.theta(j)).and.(theta(j).lt.10.d0)) then
      temp=1.d0
    end if
    temp0 = temp * temp0
  end do

  if (temp0.gt.0.5d0) then
    gamma = gamma_int
    call hybrd10(esteqgamma,L,gamma,esteqgamma_out,tol1,eflag1,wv1,lw1)
    !print*, 'esteqgamma_out=', esteqgamma_out
    !print*, 'eflag1 =', eflag1
    !print*, 'gamma=', gamma

    call effsc(alpha, beta, gamma, n, y, z, r, u, effscalpha_out, effscbeta_out)
    esteq_out(1) = sum(effscalpha_out)/n
    esteq_out(2) = sum(effscbeta_out)/n
  end if

  !print*, 'theta=', theta
  !print*, 'esteq_out=', esteq_out

return 
end subroutine esteqtheta




subroutine matrixV(theta, gamma, V)
  use commondata, only : lt, L, LL, n, z, r, y, u, n1, d, h,  y1, u1, z0, u0, len2, yy2, ww2
  double precision,dimension(lt),intent(in) :: theta
  double precision,dimension(L), intent(in) :: gamma
  double precision,dimension(LL,LL),intent(out) :: V
  double precision :: dd
  double precision :: delta
  double precision, dimension(L) :: gammal, gammar
  double precision, dimension(n) :: effscbeta_out_l, effscalpha_out_l, effscbeta_out_r, effscalpha_out_r
  double precision, dimension(L) :: esteqgamma_out_l, esteqgamma_out_r
  double precision,dimension(LL,LL) :: Q, Qinv, RR, Q1
  double precision, dimension(n) :: effscbeta_out, effscalpha_out
  double precision, dimension(L,n) :: esteqgamma_out
  double precision,dimension(LL,n) :: effsc_out 
  integer :: i, j
  double precision :: cond
  integer,dimension(LL) :: ipvt
  double precision,dimension(LL) :: wv

  dd=1e-10
  delta=1e-8
  !delta=0.d0
  gammal = gamma
  gammar = gamma
  
  call effsc(theta(1)-dd, theta(2), gamma, n, y, z, r, u, effscalpha_out_l, effscbeta_out_l)
  call effsc(theta(1)+dd, theta(2), gamma, n, y, z, r, u, effscalpha_out_r, effscbeta_out_r)
  call gammaFNC(theta(1)-dd, theta(2), L, gamma, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out_l)
  call gammaFNC(theta(1)+dd, theta(2), L, gamma, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out_r)

  Q(1,1) = sum(effscalpha_out_r-effscalpha_out_l)/2/dd/n
  Q(2,1) = sum(effscbeta_out_r - effscbeta_out_l)/2/dd/n
  Q(3:LL,1) = (esteqgamma_out_r - esteqgamma_out_l)/2/dd

  call effsc(theta(1), theta(2)-dd, gamma, n, y, z, r, u, effscalpha_out_l, effscbeta_out_l)
  call effsc(theta(1), theta(2)+dd, gamma, n, y, z, r, u, effscalpha_out_r, effscbeta_out_r)
  call gammaFNC(theta(1), theta(2)-dd, L, gamma, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out_l)
  call gammaFNC(theta(1), theta(2)+dd, L, gamma, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out_r)
  Q(1,2) = sum(effscalpha_out_r-effscalpha_out_l)/2/dd/n
  Q(2,2) = sum(effscbeta_out_r - effscbeta_out_l)/2/dd/n
  Q(3:LL,2) = (esteqgamma_out_r - esteqgamma_out_l)/2/dd

  do i=1,L
    gammal(i) = gamma(i)-dd
    gammar(i) = gamma(i)+dd
    call effsc(theta(1), theta(2), gammal, n, y, z, r, u, effscalpha_out_l, effscbeta_out_l)
    call effsc(theta(1), theta(2), gammar, n, y, z, r, u, effscalpha_out_r, effscbeta_out_r)
    call gammaFNC(theta(1), theta(2), L, gammal, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out_l)
    call gammaFNC(theta(1), theta(2), L, gammar, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out_r)
    Q(1,(i+2)) = sum(effscalpha_out_r-effscalpha_out_l)/2/dd/n
    Q(2,(i+2)) = sum(effscbeta_out_r - effscbeta_out_l)/2/dd/n
    Q(3:LL,(i+2)) = (esteqgamma_out_r - esteqgamma_out_l)/2/dd
    gammal(i) = gamma(i)
    gammar(i) = gamma(i)
  end do 

  do j=1,LL
    Q(j,j)=Q(j,j)+delta
  end do

  !print*, 'Q=', Q
  
  !Q1 = Q
  !call decomp(LL,LL,Q1,cond,ipvt,wv)
  !print*, 'condQ1=', cond

  call inv(LL, Q, Qinv)
  !print*, 'Q=', Q
  !print*, 'Qinv=', Qinv

  call effsc(theta(1), theta(2), gamma, n, y, z, r, u, effscalpha_out, effscbeta_out)
  call gammaFNC1(theta(1), theta(2), L, gamma, d, h, n, n1, y1, u1, z0, u0, esteqgamma_out)
  effsc_out(1,:) = effscalpha_out
  effsc_out(2,:) = effscbeta_out
  effsc_out(3:LL,:) = esteqgamma_out
  RR=matmul(effsc_out,transpose(effsc_out))/n
  V=matmul(matmul(Qinv,RR),transpose(Qinv))

return
end subroutine matrixV




subroutine matrixV1(theta, gamma, V, RRinv)
  use commondata, only : lt, n, L, z, r, y, u,len2,yy2,ww2, d
  double precision,dimension(lt),intent(in) :: theta
  double precision,dimension(L), intent(in) :: gamma
  double precision,dimension(lt,lt),intent(out) :: V
  double precision,dimension(lt,lt),intent(out) :: RRinv
  double precision :: dd
  double precision :: delta
  double precision,dimension(lt,lt) :: Q, Qinv, RR, Q1
  double precision, dimension(n) :: effscbeta_out_l, effscalpha_out_l, effscbeta_out_r, effscalpha_out_r
  double precision, dimension(n) :: effscbeta_out, effscalpha_out
  double precision,dimension(lt,n) :: effsc_out
  integer :: j
  
  double precision :: cond
  integer,dimension(2) :: ipvt
  double precision,dimension(2) :: wv
  
  dd=1e-10
  delta=0.d0
  !delta = 1e-10

  call effsc(theta(1)-dd, theta(2), gamma, n, y, z, r, u, effscalpha_out_l, effscbeta_out_l)
  call effsc(theta(1)+dd, theta(2), gamma, n, y, z, r, u, effscalpha_out_r, effscbeta_out_r)
  Q(1,1) = sum(effscalpha_out_r-effscalpha_out_l)/2/dd/n
  Q(2,1) = sum(effscbeta_out_r - effscbeta_out_l)/2/dd/n

  call effsc(theta(1), theta(2)-dd, gamma, n, y, z, r, u, effscalpha_out_l, effscbeta_out_l)
  call effsc(theta(1), theta(2)+dd, gamma, n, y, z, r, u, effscalpha_out_r, effscbeta_out_r)
  Q(1,2) = sum(effscalpha_out_r-effscalpha_out_l)/2/dd/n
  Q(2,2) = sum(effscbeta_out_r - effscbeta_out_l)/2/dd/n

  do j=1,lt 
    Q(j,j)=Q(j,j)+delta
  end do

  !Q1 = Q
  !call decomp(lt,lt,Q1,cond,ipvt,wv)
  !print*, 'condQ2=', cond

  call inv(lt, Q, Qinv)
  !print*, 'Q=', Q
  !print*, 'Qinv=', Qinv

  call effsc(theta(1), theta(2), gamma, n, y, z, r, u, effscalpha_out, effscbeta_out)
  effsc_out(1,:) = effscalpha_out
  effsc_out(2,:) = effscbeta_out
  RR=matmul(effsc_out,transpose(effsc_out))/n
  V=matmul(matmul(Qinv,RR),transpose(Qinv))

  call inv(lt, RR, RRinv)

return
end subroutine matrixV1






