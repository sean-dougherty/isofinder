! Pedro A. Bernaola-Galvan, Dept. of Applied Physics II    
! University of Malaga, Malaga (SPAIN)
! E-MAIL: rick@ctima.uma.es
!
!
! Shannon Entropy, Jensen-Shannon divergence
! and related functions
!
! Last updated: Jun-10-99


module diver

use chi_c

implicit none

real*8, parameter :: log2=0.6931471805599453d0
real*8, parameter :: log10=2.3025851d0
real*8, parameter:: An4=4.165,Bn4=0.596
real*8, parameter:: Ak4=4.837,Bk4=0.168  

real*8, parameter:: An2=8.082,Bn2=0.155
real*8, parameter:: Ak2=1.952,Bk2=0.157  

real*8, parameter:: Anew4 =-3.701,Bnew4 =2.831,factor4 =0.824
real*8, parameter:: Anew2 =-1.988,Bnew2 =3.086,factor2 =0.824 !it is the same

!real*8, parameter:: Anew12=-9.8657,Bnew12=2.4483,factor12=0.835 !DNA

real*8, parameter:: Anew12=-7.490,Bnew12=2.090,factor12=0.833 !this is for chi2 OTHERS

real*8, parameter:: Anew64=-1.607,Bnew64=0.492,factor64=0.350 !this is for chi2
real*8, parameter:: Anew7 =-9.208,Bnew7 =2.486,factor7 =0.833 !this is for chi2

!real*8, parameter:: Anew444=-19.60,Bnew444=4.993,factor444=0.842 !this is for chi2
real*8, parameter:: Anew444=-15.28,Bnew444=4.386,factor444=0.842 !this is for chi2


INTERFACE entro  
  module procedure entroi,entror
end interface

INTERFACE jsdiver
  module procedure diveri,diverr
end interface


INTERFACE sig_level
  module procedure sl_chi_d,sl_chi_i,sl_chi_r
end interface  


!******* WARNING ************
!This routine only works with ALPHABET=4

INTERFACE sig_levelm
  module procedure sl_m_d,sl_m_i,sl_m_r
end interface

!***********************************************

CONTAINS

!*************  ENTROPY ************************

!Shannon entropy for a relative frequency vector
!Note that the parameter is a real array

function entror(x)
   real*8,dimension(:) ::x
   real*8  entror

   where (x==0) x=1.d0  !avoids log(0)

   entror=-sum(x*dlog(x))/log2
end function entror

!**********************************************

!Shannon entropy for an absolute frequency vector
!The parameter is an integer array

function entroi(x)
   integer*4, dimension(:) ::x
   real*8 entroi
   if (sum(x)==0) then
      print*,'ERROR: Null frequency vector'
      entroi=0
   else   
     entroi=entror(1.d0*x/sum(x))
    endif 
end function entroi

!*********************************************


!**********  JS DIVERGENCE *******************

!JS divergence for a relative frequency vector
!The parameter is a real array
 
function diverr(p1,p2,n1,n2)
  integer*4 n1,n2
  real*8, dimension(:) :: p1,p2
  real*8  diverr,n

  n=n1+n2
  diverr=entro((n1*p1+n2*p2)/n)- &
  & (n1*entro(p1)+n2*entro(p2))/n
end function diverr

!**********************************************

! JS divergence for a relative frequency vector
! The paramter is an integer array

function diveri(f1,f2)
  real*8 n1,n2,n
  integer*4, dimension(:) :: f1,f2
  real*8 diveri 
  n1=sum(f1)
  n2=sum(f2)
  n=n1+n2
  diveri=entro(f1+f2)-(n1*entro(f1)+n2*entro(f2))/n
end function diveri

!**********************************************


!*************  SIGNIFICANCE LEVEL ************

function sl_chi_d(dim,d,n)

  !Given the alphabet size (dim), the divergence 
  !value (d) and the size of the sequence to be
  !splited (n) this function returns the signficance
  !level.

  integer*4 dim,n
  real*8 sl_chi_d,d
  if (dim<100) then
    sl_chi_d=gammp((dim-1.0)/2.d0,d*log2*n)
  else
    sl_chi_d=gammp_approx((dim-1.0)/2.d0,d*log2*n)
  endif
end function sl_chi_d

!**********************************************

function sl_chi_i(dim,f1,f2)

  !The same as the previous one but the divergence
  !is computed inside the function

  integer*4 dim,n
  integer*4, dimension(:) ::f1,f2
  real*8 d,sl_chi_i
  n=sum(f1)+sum(f2)
  d=jsdiver(f1,f2)
  sl_chi_i=sl_chi_d(dim,d,n)
end function sl_chi_i


!***********************************************

function sl_chi_r(dim,p1,p2,n1,n2)

  ! The same as the previous one but using
  ! relative frequency vectors instead of absolute

  integer*4 dim,n1,n2
  real*8, dimension(:) :: p1,p2
  real*8 d,sl_chi_r
  d=jsdiver(p1,p2,n1,n2)
  sl_chi_r=sl_chi_d(dim,d,n1+n2)
end function sl_chi_r

!***********************************************


!*******  SIGNIFICANCE LEVEL Maximum************

!***********************************************
function Neff2(n,dim)
  real*8    Neff2,n
  integer*4 dim

  if (dim==4) then
    Neff2=An4+Bn4*dlog(n)
  else
    Neff2=An2+Bn2*dlog(n)
  endif    

end function Neff2


!*******************************
function Keff(n,dim)
  real*8    Keff,n
  integer*4 dim

  if (dim==4) then
    Keff=Ak4+Bk4*dlog(n)
  else
    keff=Ak2+Bk2*dlog(n)
  endif  

end function Keff

!*******************************

function Nnew(n,dim)
  real*8     Nnew,n
  integer*4  dim
 
  if (dim==4) then
    Nnew=Anew4+Bnew4*dlog(n)
  else if (dim==2) then
    Nnew=Anew2+Bnew2*dlog(n)
  else if (dim==64) then
    Nnew=Anew64+Bnew64*dlog(n)  
  else if (dim==12) then
    Nnew=Anew12+Bnew12*dlog(n) 
  else if (dim==7) then
    Nnew=Anew7+Bnew7*dlog(n)    
  else if (dim==444) then
    Nnew=Anew444+Bnew444*dlog(n)   
  endif

end function Nnew

!*******************************


 !computes chi-sqare dist. with
 !k-1 degrees of freedom, to the 
 !power of n

function gammp_n(x,n,k)
  real*8    k,n
  real*8    x,gammp_n,inter

   inter=gammp((k-1.d0)/2.d0,x/2.d0)
   if (inter<1.0d-10) then
     gammp_n=0.0
   else
     gammp_n=inter**n
   endif
end function gammp_n

!*******************************



 !computes chi-sqare dist. with
 !k-1 degrees of freedom, to the 
 !power of n, incluiding a factor
 

function gammp_n_f(x,n,k,factor)
  real*8    k,n,factor
  real*8    x,gammp_n_f,inter

  ! print*, factor, n
   inter=gammp((k-1.d0)/2.d0,factor*x/2.d0)
   if (inter<1.0d-10) then
     gammp_n_f=0.0
   else
     gammp_n_f=inter**n
   endif
end function gammp_n_f

!*******************************


function sl_m_d(dim,d,n)

  !Given the alphabet size (dim), the divergence 
  !value (d) and the size of the sequence to be
  !splited (n) this function returns the signficance
  !level considered as maximum value. It uses the
  !fitting approximation of Neff2 and Keff

  integer*4   dim,n
  real*8      sl_m_d,d
  real*8      ne

  ne=Nnew(1.d0*n,dim)
  !ne=Neff2(1.d0*n,dim)
  !ke=Keff2(1.d0*n,dim)

 
 if (dim==2) then 
   sl_m_d=gammp_n_f(2*d*log2*n,ne,1.d0*dim,factor2)
 else if  (dim==4) then
   sl_m_d=gammp_n_f(2*d*log2*n,ne,1.d0*dim,factor4)
 else if  (dim==64) then
   sl_m_d=gammp_n_f(2*d*log2*n,ne,16.d0,factor64)
 else if  (dim==12) then
   !sl_m_d=gammp_n_f(2*d*log2*n,ne,10.d0,factor12)  !DNA
   sl_m_d=gammp_n_f(2*d*log2*n,ne,12.d0,factor12)  !OTHERS
 else if  (dim==7) then
   sl_m_d=gammp_n_f(2*d*log2*n,ne,7.d0,factor7)   
 else if  (dim==444) then 
   sl_m_d=gammp_n_f(2*d*log2*n,ne,4.d0,factor444) 
 endif

end function sl_m_d

!**********************************************

function sl_m_i(dim,f1,f2)

  !The same as the previous one but the divergence
  !is computed inside the function

  integer*4 dim,n
  integer*4, dimension(:) ::f1,f2
  real*8 d,sl_m_i
  n=sum(f1)+sum(f2)
  d=jsdiver(f1,f2)
  sl_m_i=sl_m_d(dim,d,n)
end function sl_m_i


!***********************************************

function sl_m_r(dim,p1,p2,n1,n2)

  ! The same as the previous one but using
  ! relative frequency vectors instead of absolute

  integer*4 dim,n1,n2
  real*8, dimension(:) :: p1,p2
  real*8 d,sl_m_r
  d=jsdiver(p1,p2,n1,n2)
  sl_m_r=sl_m_d(dim,d,n1+n2)
end function sl_m_r

!***********************************************





END MODULE diver

