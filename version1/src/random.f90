module random

CONTAINS

!*****************************
FUNCTION ran3(idum)
  integer*4  idum
  integer*4 mbig,mseed,mz
  real*8    ran3,fac

  parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)

  integer*4 i,iff,ii,inext,inextp,k
  integer*4 mj,mk,ma(55)

  save iff,inext,inextp,ma
  data iff /0/
  if (idum<0 .or. iff==0) then
    iff=1
    mj=mseed-iabs(idum)
    mj=mod(mj,mbig)
    ma(55)=mj
    mk=1
    do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if (mk<mz) mk=mk+mbig
      mj=ma(ii)
    enddo
    do k=1,4
      do i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if (ma(i)<mz) ma(i)=ma(i)+mbig
      enddo
    enddo
    inext=0
    inextp=31
    idum=1
  endif
  inext=inext+1
  if (inext==56) inext=1
  inextp=inextp+1
  if (inextp==56) inextp=1
  mj=ma(inext)-ma(inextp)
  if (mj<mz) mj=mj+mbig
  ma(inext)=mj
  ran3=mj*fac
  return

end function ran3
!*****************************



end MODULE random