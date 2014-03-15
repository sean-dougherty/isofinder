  ! Pedro A. Bernaola-Galvan, Dept. of Applied Physics II
  ! University of Malaga, Malaga (SPAIN)
  ! E-MAIL: rick@ctima.uma.es
  !
  ! Last updated: Mar-13-99
  ! This module is used to compute the significance
  ! level using the function gammp (chi-square)
  !(Numerical Recipes in Fortran)
  ! includes also BETA FUNCTION
  !
  !Last updated: Mar-13-1999


  module chi_c

  implicit none
  
  contains

 !*************************************************
  function normal(x)
    real*8 x,normal

	if (x>=0) then 
	  normal=0.5d0*(1.d0+gammp(0.5d0,abs(x*x)))
	else
	  normal=0.5d0*(1.d0-gammp(0.5d0,abs(x*x)))
	endif

  end function normal
 !*************************************************

  function gammp_approx(a,x)
    real*8  x2,x,gammp_approx,a

      x2=dsqrt(2.d0*x)-dsqrt(2.d0*a-1.d0)
	  gammp_approx=normal(x2)
  end function gammp_approx
  
 !**************************************************
  
  function gammln(xx)
  real*8 gammln,xx
  integer*4 j
  double precision ser,stp,tmp,x,y,cof(6)
  
  save cof,stp
  DATA cof,stp/76.18009172947146d0, -86.50532032941677d0, &
               24.01409824083091d0,-1.231739572450155d0, &
               .1208650973866179d-2,-.5395239384953d-5,   &
               2.5066282746310005d0/
               
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  enddo               
  gammln=tmp+log(stp*ser/x)
  return
  end function gammln
                  
  
  !******************************************************
  function gammp(a,x)
  real*8 a,gammp,x
  
  ! uses GCF and GSER
  real*8 gammcf,gamser,gln
!  if (x.lt.0.or.a.le.0.) pause 'bad argument in GAMMP'
  if (x.lt.a+1.) then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=1.-gammcf
  endif
  return
  end function gammp    
  
  !******************************************************
  
  subroutine gser(gamser,a,x,gln)
  integer*4 itmax
  real*8 a,gamser,gln,x,eps
  parameter (itmax=100,eps=3.e-7)
  ! uses GAMMLN
  integer*4 n
  real*8 ap,del,sum !,gammln
  gln=gammln(a)
  if (x.le.0) then
!    if (x.lt.0) pause 'x<0 in GSER'
    gamser=0.
    return
  endif
  ap=a
  sum=1./a
  del=sum  
  do n=1,itmax
    ap=ap+1
    del=del*x/ap
    sum=sum+del
    if (abs(del).lt.abs(sum)*eps) goto 1  
  enddo
!  pause 'a too large, ITMAX too small in gser'
  1 gamser=sum*exp(-x+a*log(x)-gln)
  return
  end subroutine gser
  
  !*********************************************************
  
  
  subroutine gcf(gammcf,a,x,gln)
  real*8 a,gammcf,gln,x,eps,fpmin
  parameter (eps=3.e-7,fpmin=1.e-30)
  integer*4, parameter :: itmax=100
  !uses GAMLN
  integer*4 i
  real*8 an,b,c,d,del,h  !,gammln
  gln=gammln(a)
  b=x+1.-a
  c=1./fpmin
  d=1./b
  h=d
  do i=1,itmax
    an=-i*(i-a)
    b=b+2
    d=an*d+b
    if (abs(d).lt.fpmin) d=fpmin
    c=b+an/c
    if (abs(c).lt.fpmin) c=fpmin
    d=1./d
    del=d*c
    h=h*del
    if (abs(del-1.).lt.eps) goto 1  
  enddo
!  pause 'a too large ITMAX too small in gcf'
  1 gammcf=exp(-x+a*log(x)-gln)*h
  return
  end subroutine gcf
  
  !**********************************************************
  
  
  FUNCTION betai2(a,b,x)
    real*8  a,b,x,betai2
    real*8  bt
    
!    if (x<0.0 .or. x>1.0) pause 'bad argument x in beta i'
    if (x==0.0 .or. x==1.0) then
      bt=0.0
    else
      bt=dexp(gammln(a+b)-gammln(a)-gammln(b) &
         +a*dlog(x) + b*dlog(1.0-x))
    
    endif
    if (x<(a+1.0)/(a+b+2.0)) then
      betai2=bt*betacf(a,b,x)/a
      return
    else
      betai2=1.0-bt*betacf(b,a,1.0-x)/b
      return
    endif     
  END FUNCTION betai2  
  
  
 !************************************************************
 FUNCTION betacf(a,b,x)
   integer*4   maxit
   real*8      betacf,a,b,x,eps,fpmin
   parameter   (maxit=100,eps=3.d-7,fpmin=1.0d-30)
   integer*4   m,m2
   real*8      aa,c,d,del,h,qab,qam,qap
   
   qab=a+b
   qap=a+1.d0
   qam=a-1.d0
   c=1.0
   d=1.0-qab*x/qap
   if (dabs(d)<fpmin) d=fpmin
   d=1.0/d
   h=d
  
   do m=1,maxit
     m2=2*m
     aa=m*(b-m)*x/((qam+m2)*(a+m2))
     d=1.0+aa*d
     if (dabs(d)<fpmin) d=fpmin
     c=1.0+aa/c
     if (dabs(c)<fpmin) c=fpmin
     d=1.0/d
     h=h*d*c
     aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
     d=1.0+aa*d
     if (dabs(d)<fpmin) d=fpmin
     c=1.0+aa/c
     if (dabs(c)<fpmin) c=fpmin
     d=1.0/d
     del=d*c
     h=h*del
     if (dabs(del-1.0)<eps) goto 1
   enddo
!   pause 'a b too big, or MAXIT too small in betacf'
 1 betacf=h
   return
  
 END FUNCTION betacf
 
 
 !************************************************************ 
  
 
  
 end module chi_c 
 
