module special

implicit none

contains

!*****************************************************

!***** CALCULATION OF SPECIAL FUNCTIONS (Numerical Recipes) *****
!****************************************************************


!**************************************************
FUNCTION  ttest0(t,df)
  real*8  df,t,ttest0

   ttest0=betai(0.5d0*df,0.5d0,df/(df+t**2))

END function ttest0
!**************************************************

!**************************************************
!*** Test de Student para varianzas diferentes ****
!**************************************************
FUNCTION ttest_difvar(ave1,ave2,var1,var2,n1,n2)
   REAL*8     n1,n2
   REAL*8     prob,t,ave1,ave2,var1,var2,ttest_difvar
   REAL*8     df
   
   t=(ave1-ave2)/sqrt(var1/n1+var2/n2)
   df=(var1/n1+var2/n2)**2/((var1/n1)**2/(n1-1)+(var2/n2)**2/(n2-1))
   prob=betai(0.5d0*df,0.5d0,df/(df+t**2))
   ttest_difvar=prob
END function ttest_difvar
!****************************************************



!*****************************************************
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

  !***********************************************************

  FUNCTION betai(a,b,x)
    real*8  a,b,x,betai
    real*8  bt

!   if (x<0.0 .or. x>1.0) pause 'bad argument x in beta i'
    if (x==0.0 .or. x==1.0) then
      bt=0.0
    else
      bt=dexp(gammln(a+b)-gammln(a)-gammln(b) &
         +a*dlog(x) + b*dlog(1.0-x))

    endif
    if (x<(a+1.0)/(a+b+2.0)) then
      betai=bt*betacf(a,b,x)/a
      return
    else
      betai=1.0-bt*betacf(b,a,1.0-x)/b
      return
    endif
  END FUNCTION betai


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
!  pause 'a b too big, or MAXIT too small in betacf'
 1 betacf=h
   return

 END FUNCTION betacf


 !************************************************************

END MODULE special
